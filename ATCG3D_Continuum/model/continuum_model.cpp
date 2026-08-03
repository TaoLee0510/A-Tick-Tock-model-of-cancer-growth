#include "model/continuum_model.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <type_traits>

#include "common/density_growth_rule.hpp"
#include "engine/simulation.hpp"
#include "geometry/footprint.hpp"

namespace atcg3d::continuum {
namespace {

constexpr std::array<char, 8> kCheckpointMagic{
    {'A', 'T', 'C', 'G', 'C', 'P', 'D', '1'}};
constexpr std::uint32_t kCheckpointVersion = 1;
constexpr double kFixed26DiffusionFactor = 9.0 / 26.0;

bool same_time(double lhs, double rhs) noexcept {
    return std::abs(lhs - rhs) <=
        1.0e-10 * std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

std::uint64_t hash_mix(std::uint64_t state, std::uint64_t value) noexcept {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    value ^= value >> 31U;
    return state ^ (value + (state << 6U) + (state >> 2U));
}

template <class T>
void write_pod(std::ostream& stream, const T& value) {
    static_assert(std::is_trivially_copyable_v<T>);
    stream.write(reinterpret_cast<const char*>(&value), sizeof(value));
    if (!stream) throw std::runtime_error("unable to write continuum checkpoint");
}

template <class T>
T read_pod(std::istream& stream) {
    static_assert(std::is_trivially_copyable_v<T>);
    T result{};
    stream.read(reinterpret_cast<char*>(&result), sizeof(result));
    if (!stream) throw std::runtime_error("truncated continuum checkpoint");
    return result;
}

double beta_mean(const BetaRateConfig& config) noexcept {
    return config.scale * config.alpha / (config.alpha + config.beta);
}

double positive_part(double value) noexcept {
    return std::max(0.0, value);
}

}  // namespace

ContinuumModel3D::ContinuumModel3D(ContinuumModelConfig3D config)
    : config_(std::move(config)) {
    config_.validate();
    voxel_count_ = static_cast<std::size_t>(config_.grid.shape[0]) *
        static_cast<std::size_t>(config_.grid.shape[1]) *
        static_cast<std::size_t>(config_.grid.shape[2]);
    const int dimensions = config_.base.thin_layer ? 2 : 3;
    voxel_measure_ = std::pow(config_.grid.spacing_voxels, dimensions);
    large_cell_volume_ = std::pow(config_.base.large_footprint_edge, dimensions);
    for (auto& field : populations_) field.assign(voxel_count_, 0.0);
    for (auto& field : work_) field.assign(voxel_count_, 0.0);
    nutrient_.assign(voxel_count_, 0.0);
    nutrient_next_.assign(voxel_count_, 0.0);
    vessel_.assign(voxel_count_, 0.0);

    double maximum_diffusion = 0.0;
    for (std::size_t field = 0; field < kPopulationFieldCount3D; ++field) {
        maximum_diffusion = std::max(
            maximum_diffusion,
            base_diffusion(static_cast<PopulationField3D>(field), 1.0));
    }
    const double cfl = config_.time_step_hours * 2.0 * dimensions *
        maximum_diffusion /
        (config_.grid.spacing_voxels * config_.grid.spacing_voxels);
    if (cfl > 0.45 + 1.0e-12) {
        throw std::invalid_argument(
            "continuum time step violates the explicit migration CFL bound");
    }
}

std::size_t ContinuumModel3D::index(int x, int y, int z) const noexcept {
    return (static_cast<std::size_t>(z) * config_.grid.shape[1] +
            static_cast<std::size_t>(y)) * config_.grid.shape[0] +
        static_cast<std::size_t>(x);
}

bool ContinuumModel3D::grid_coordinate(Vec3i site,
                                       int& x,
                                       int& y,
                                       int& z) const noexcept {
    const std::array<double, 3> coordinate{
        static_cast<double>(site.x), static_cast<double>(site.y),
        static_cast<double>(site.z)};
    std::array<int, 3> result{};
    for (std::size_t axis = 0; axis < 3; ++axis) {
        result[axis] = static_cast<int>(std::floor(
            (coordinate[axis] - config_.grid.origin[axis]) /
            config_.grid.spacing_voxels));
        if (result[axis] < 0 || result[axis] >= config_.grid.shape[axis]) {
            return false;
        }
    }
    x = result[0];
    y = result[1];
    z = result[2];
    return true;
}

void ContinuumModel3D::initialize_from_abm(const Simulation3D& simulation) {
    if (initialized_) throw std::logic_error("continuum model is already initialized");
    for (auto& field : populations_) std::fill(field.begin(), field.end(), 0.0);
    std::fill(vessel_.begin(), vessel_.end(), 0.0);

    const double small_weight = 1.0 / voxel_measure_;
    const double large_site_weight = 1.0 / (large_cell_volume_ * voxel_measure_);
    for (const Slot slot : simulation.cells().alive_slots()) {
        const CellType type = simulation.cells().type(slot);
        const CellStage stage = simulation.cells().stage(slot);
        const PopulationField3D field = type == CellType::r
            ? (stage == CellStage::large ? PopulationField3D::r_large
                                         : PopulationField3D::r_small)
            : (stage == CellStage::large ? PopulationField3D::K_large
                                         : PopulationField3D::K_small);
        auto& target = populations_[static_cast<std::size_t>(field)];
        if (stage == CellStage::large) {
            for (const Vec3i site : large_footprint(simulation.cells().anchor(slot))) {
                int x{}, y{}, z{};
                if (!grid_coordinate(site, x, y, z)) {
                    throw std::runtime_error(
                        "continuum grid does not contain an ABM large-cell footprint");
                }
                target[index(x, y, z)] += large_site_weight;
            }
        } else {
            int x{}, y{}, z{};
            if (!grid_coordinate(simulation.cells().anchor(slot), x, y, z)) {
                throw std::runtime_error(
                    "continuum grid does not contain an ABM cell anchor");
            }
            // The four-field v1 reduction maps optional ultrasmall agents to
            // the small density class. Supplied reference profiles disable
            // ultrasmall cells, so this is primarily an import compatibility path.
            target[index(x, y, z)] += small_weight;
        }
    }

    if (config_.vascular.source_mode == "abm_perfusion" ||
        config_.vascular.source_mode == "abm_plus_synthetic_line") {
        for (const Vec3i site : simulation.vessel_grid().occupied_sites()) {
            if (!simulation.vessel_grid().perfused(site)) continue;
            int x{}, y{}, z{};
            if (!grid_coordinate(site, x, y, z)) continue;
            const std::size_t location = index(x, y, z);
            vessel_[location] = std::min(1.0, vessel_[location] + 1.0 / voxel_measure_);
        }
    }
    if (config_.vascular.source_mode == "synthetic_central_line" ||
        config_.vascular.source_mode == "abm_plus_synthetic_line") {
        add_synthetic_vessel();
    }
    time_hours_ = config_.initialization_mode == "abm_checkpoint"
        ? simulation.clock().time_hours : config_.start_time_hours;
    if (time_hours_ + 1.0e-10 < config_.start_time_hours ||
        time_hours_ >= config_.end_time_hours) {
        throw std::runtime_error("ABM import time is outside continuum run bounds");
    }
    solve_nutrient();
    next_nutrient_refresh_hours_ = time_hours_ + config_.nutrient.refresh_every_hours;
    initialized_ = true;
    validate_state();
}

void ContinuumModel3D::initialize_from_arrays(
    std::array<std::vector<double>, kPopulationFieldCount3D> populations,
    std::vector<double> vessel_fraction,
    double time_hours) {
    if (initialized_) throw std::logic_error("continuum model is already initialized");
    for (const auto& field : populations) {
        if (field.size() != voxel_count_) {
            throw std::invalid_argument("continuum initial population size mismatch");
        }
    }
    if (vessel_fraction.size() != voxel_count_) {
        throw std::invalid_argument("continuum initial vessel size mismatch");
    }
    if (!std::isfinite(time_hours) ||
        time_hours + 1.0e-10 < config_.start_time_hours ||
        time_hours >= config_.end_time_hours) {
        throw std::invalid_argument("continuum initial array time is outside run bounds");
    }
    populations_ = std::move(populations);
    vessel_ = std::move(vessel_fraction);
    if (config_.vascular.source_mode == "synthetic_central_line" ||
        config_.vascular.source_mode == "abm_plus_synthetic_line") {
        add_synthetic_vessel();
    }
    time_hours_ = time_hours;
    solve_nutrient();
    next_nutrient_refresh_hours_ = time_hours_ + config_.nutrient.refresh_every_hours;
    initialized_ = true;
    validate_state();
}

void ContinuumModel3D::add_synthetic_vessel() {
    const int axis = config_.vascular.synthetic_axis == "x" ? 0
        : (config_.vascular.synthetic_axis == "y" ? 1 : 2);
    const double radius_squared =
        config_.vascular.synthetic_radius_voxels *
        config_.vascular.synthetic_radius_voxels;
    for (std::size_t location = 0; location < voxel_count_; ++location) {
        const auto point = coordinate(location);
        double distance_squared = 0.0;
        for (int dimension = 0; dimension < 3; ++dimension) {
            if (dimension == axis || (config_.base.thin_layer && dimension == 2)) {
                continue;
            }
            const double offset = point[dimension] -
                config_.vascular.synthetic_center[dimension];
            distance_squared += offset * offset;
        }
        if (distance_squared <= radius_squared) vessel_[location] = 1.0;
    }
}

double ContinuumModel3D::occupied_fraction(std::size_t location) const noexcept {
    return populations_[0][location] + populations_[2][location] +
        large_cell_volume_ *
            (populations_[1][location] + populations_[3][location]);
}

std::array<double, 3> ContinuumModel3D::coordinate(
    std::size_t location) const noexcept {
    const int nx = config_.grid.shape[0];
    const int ny = config_.grid.shape[1];
    const int x = static_cast<int>(location % static_cast<std::size_t>(nx));
    const std::size_t yz = location / static_cast<std::size_t>(nx);
    const int y = static_cast<int>(yz % static_cast<std::size_t>(ny));
    const int z = static_cast<int>(yz / static_cast<std::size_t>(ny));
    return {
        config_.grid.origin[0] + (x + 0.5) * config_.grid.spacing_voxels,
        config_.grid.origin[1] + (y + 0.5) * config_.grid.spacing_voxels,
        config_.grid.origin[2] + (z + 0.5) * config_.grid.spacing_voxels};
}

double ContinuumModel3D::capacity_multiplier(double value) const noexcept {
    const double local = std::clamp(value, 0.0, config_.nutrient.vessel_value);
    const double raw = local / (config_.nutrient.capacity_half_saturation + local);
    const double at_vessel = config_.nutrient.vessel_value /
        (config_.nutrient.capacity_half_saturation + config_.nutrient.vessel_value);
    const double saturation = at_vessel > 0.0
        ? std::clamp(raw / at_vessel, 0.0, 1.0) : 0.0;
    return 1.0 + (config_.nutrient.maximum_capacity_multiplier - 1.0) *
        saturation;
}

void ContinuumModel3D::solve_nutrient() {
    const int nx = config_.grid.shape[0];
    const int ny = config_.grid.shape[1];
    const int nz = config_.grid.shape[2];
    const int dimensions = config_.base.thin_layer ? 2 : 3;
    const double diffusion = config_.nutrient.diffusion_voxels2_per_hour /
        (config_.grid.spacing_voxels * config_.grid.spacing_voxels);
    const double laplacian_diagonal = 2.0 * dimensions * diffusion;
    for (int iteration = 0; iteration < config_.nutrient.solver_iterations;
         ++iteration) {
        for (int z = 0; z < nz; ++z) {
            for (int y = 0; y < ny; ++y) {
                for (int x = 0; x < nx; ++x) {
                    const std::size_t here = index(x, y, z);
                    double neighbor_sum = 0.0;
                    if (x > 0) neighbor_sum += nutrient_[index(x - 1, y, z)];
                    if (x + 1 < nx) neighbor_sum += nutrient_[index(x + 1, y, z)];
                    if (y > 0) neighbor_sum += nutrient_[index(x, y - 1, z)];
                    if (y + 1 < ny) neighbor_sum += nutrient_[index(x, y + 1, z)];
                    if (!config_.base.thin_layer) {
                        if (z > 0) neighbor_sum += nutrient_[index(x, y, z - 1)];
                        if (z + 1 < nz) neighbor_sum += nutrient_[index(x, y, z + 1)];
                    }
                    const double old = nutrient_[here];
                    const double r_occupied = populations_[0][here] +
                        large_cell_volume_ * populations_[1][here];
                    const double K_occupied = populations_[2][here] +
                        large_cell_volume_ * populations_[3][here];
                    const double r_sink =
                        config_.nutrient.r_consumption_per_occupied_voxel_hour *
                        r_occupied /
                        (config_.nutrient.r_consumption_half_saturation + old);
                    const double K_sink =
                        config_.nutrient.K_consumption_per_occupied_voxel_hour *
                        K_occupied /
                        (config_.nutrient.K_consumption_half_saturation + old);
                    const double exchange = config_.nutrient.vessel_exchange_per_hour *
                        std::clamp(vessel_[here], 0.0, 1.0);
                    const double denominator = laplacian_diagonal +
                        config_.nutrient.decay_per_hour + exchange + r_sink + K_sink;
                    const double candidate = denominator > 0.0
                        ? (diffusion * neighbor_sum +
                           exchange * config_.nutrient.vessel_value) / denominator
                        : 0.0;
                    nutrient_next_[here] = std::clamp(
                        old + config_.nutrient.relaxation * (candidate - old),
                        0.0, config_.nutrient.vessel_value);
                }
            }
        }
        nutrient_.swap(nutrient_next_);
    }
    ++nutrient_solve_count_;
}

double ContinuumModel3D::mean_growth_rate(CellType type) const noexcept {
    if (config_.base.initial_growth_rate_model == "fixed") {
        return type == CellType::r ? config_.base.initial_r_growth_rate
                                   : config_.base.initial_K_growth_rate;
    }
    const auto& distribution = type == CellType::r
        ? config_.base.initial_r_growth_truncated_normal
        : config_.base.initial_K_growth_truncated_normal;
    return std::clamp(distribution.mean, distribution.minimum,
                      distribution.maximum);
}

double ContinuumModel3D::base_diffusion(
    PopulationField3D field,
    double local_occupied_fraction) const noexcept {
    const bool r_type = field == PopulationField3D::r_small ||
                        field == PopulationField3D::r_large;
    const bool large = field == PopulationField3D::r_large ||
                       field == PopulationField3D::K_large;
    double migration_rate = r_type
        ? beta_mean(config_.base.normal_r_migration_beta)
        : (config_.base.initial_K_migration_rate_model == "fixed"
               ? config_.base.initial_K_migration_rate
               : beta_mean(config_.base.initial_K_migration_beta));
    double diffusion = kFixed26DiffusionFactor * migration_rate *
        config_.migration.diffusion_scale;
    if (r_type && config_.base.migration_activation_enabled &&
        local_occupied_fraction >= config_.base.migration_activation_threshold) {
        diffusion *= config_.migration.activated_r_mobility_multiplier;
    }
    if (large) diffusion *= config_.migration.large_mobility_multiplier;
    return diffusion;
}

void ContinuumModel3D::migrate(double dt) {
    for (std::size_t field = 0; field < kPopulationFieldCount3D; ++field) {
        work_[field] = populations_[field];
    }
    const int nx = config_.grid.shape[0];
    const int ny = config_.grid.shape[1];
    const int nz = config_.grid.shape[2];
    const double inverse_h2 = 1.0 /
        (config_.grid.spacing_voxels * config_.grid.spacing_voxels);
    const auto exchange_face = [&](std::size_t lhs, std::size_t rhs) {
        const double lhs_phi = std::clamp(
            occupied_fraction(lhs) / config_.reaction.maximum_occupied_fraction,
            0.0, 1.0);
        const double rhs_phi = std::clamp(
            occupied_fraction(rhs) / config_.reaction.maximum_occupied_fraction,
            0.0, 1.0);
        const double lhs_vacancy = 1.0 - lhs_phi;
        const double rhs_vacancy = 1.0 - rhs_phi;
        const double mobility = std::pow(
            0.5 * (lhs_vacancy + rhs_vacancy),
            config_.migration.crowding_exponent);
        const double face_phi = 0.5 * (lhs_phi + rhs_phi);
        for (std::size_t field = 0; field < kPopulationFieldCount3D; ++field) {
            const double diffusion = base_diffusion(
                static_cast<PopulationField3D>(field), face_phi);
            const bool large = field ==
                    static_cast<std::size_t>(PopulationField3D::r_large) ||
                field == static_cast<std::size_t>(PopulationField3D::K_large);
            const double vacancy_exponent = large ? large_cell_volume_ : 1.0;
            const double lhs_availability = std::pow(lhs_vacancy, vacancy_exponent);
            const double rhs_availability = std::pow(rhs_vacancy, vacancy_exponent);
            const double rate = diffusion * inverse_h2 * mobility *
                (populations_[field][lhs] * rhs_availability -
                 populations_[field][rhs] * lhs_availability);
            work_[field][lhs] -= dt * rate;
            work_[field][rhs] += dt * rate;
        }
    };
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                const std::size_t here = index(x, y, z);
                if (x + 1 < nx) exchange_face(here, index(x + 1, y, z));
                if (y + 1 < ny) exchange_face(here, index(x, y + 1, z));
                if (!config_.base.thin_layer && z + 1 < nz) {
                    exchange_face(here, index(x, y, z + 1));
                }
            }
        }
    }
    for (std::size_t field = 0; field < kPopulationFieldCount3D; ++field) {
        for (double& value : work_[field]) {
            if (value < -1.0e-10) {
                throw std::runtime_error("continuum migration produced negative density");
            }
            value = std::max(0.0, value);
        }
        populations_[field].swap(work_[field]);
    }
}

void ContinuumModel3D::build_local_counts(
    std::vector<double>& r_counts,
    std::vector<double>& K_counts) const {
    const int nx = config_.grid.shape[0];
    const int ny = config_.grid.shape[1];
    const int nz = config_.grid.shape[2];
    const int px = nx + 1;
    const int py = ny + 1;
    const int pz = nz + 1;
    const std::size_t prefix_size = static_cast<std::size_t>(px) * py * pz;
    std::vector<double> r_prefix(prefix_size, 0.0);
    std::vector<double> K_prefix(prefix_size, 0.0);
    const auto prefix_index = [px, py](int x, int y, int z) {
        return (static_cast<std::size_t>(z) * py + y) * px + x;
    };
    for (int z = 1; z <= nz; ++z) {
        for (int y = 1; y <= ny; ++y) {
            for (int x = 1; x <= nx; ++x) {
                const std::size_t source = index(x - 1, y - 1, z - 1);
                const double r_value =
                    (populations_[0][source] + populations_[1][source]) *
                    voxel_measure_;
                const double K_value =
                    (populations_[2][source] + populations_[3][source]) *
                    voxel_measure_;
                const auto update = [&](std::vector<double>& prefix, double value) {
                    prefix[prefix_index(x, y, z)] = value
                        + prefix[prefix_index(x - 1, y, z)]
                        + prefix[prefix_index(x, y - 1, z)]
                        + prefix[prefix_index(x, y, z - 1)]
                        - prefix[prefix_index(x - 1, y - 1, z)]
                        - prefix[prefix_index(x - 1, y, z - 1)]
                        - prefix[prefix_index(x, y - 1, z - 1)]
                        + prefix[prefix_index(x - 1, y - 1, z - 1)];
                };
                update(r_prefix, r_value);
                update(K_prefix, K_value);
            }
        }
    }
    const int radius = std::max(0, static_cast<int>(std::floor(
        0.5 * config_.base.growth_density_window_edge /
        config_.grid.spacing_voxels)));
    r_counts.assign(voxel_count_, 0.0);
    K_counts.assign(voxel_count_, 0.0);
    const auto box_sum = [&](const std::vector<double>& prefix,
                             int x0, int y0, int z0,
                             int x1, int y1, int z1) {
        return prefix[prefix_index(x1, y1, z1)]
            - prefix[prefix_index(x0, y1, z1)]
            - prefix[prefix_index(x1, y0, z1)]
            - prefix[prefix_index(x1, y1, z0)]
            + prefix[prefix_index(x0, y0, z1)]
            + prefix[prefix_index(x0, y1, z0)]
            + prefix[prefix_index(x1, y0, z0)]
            - prefix[prefix_index(x0, y0, z0)];
    };
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                const int x0 = std::max(0, x - radius);
                const int y0 = std::max(0, y - radius);
                const int z0 = config_.base.thin_layer ? 0 : std::max(0, z - radius);
                const int x1 = std::min(nx, x + radius + 1);
                const int y1 = std::min(ny, y + radius + 1);
                const int z1 = config_.base.thin_layer ? 1 : std::min(nz, z + radius + 1);
                const std::size_t here = index(x, y, z);
                r_counts[here] = box_sum(r_prefix, x0, y0, z0, x1, y1, z1);
                K_counts[here] = box_sum(K_prefix, x0, y0, z0, x1, y1, z1);
            }
        }
    }
}

void ContinuumModel3D::react(double dt) {
    std::vector<double> r_counts;
    std::vector<double> K_counts;
    build_local_counts(r_counts, K_counts);
    const double r_inherent = mean_growth_rate(CellType::r);
    const double K_inherent = mean_growth_rate(CellType::K);
    const double r_limit = config_.base.thin_layer
        ? config_.base.legacy_mapping.source_r_limit : config_.base.r_limit;
    const double K_limit = config_.base.thin_layer
        ? config_.base.legacy_mapping.source_K_limit : config_.base.K_limit;
    const double r_capacity = config_.base.thin_layer
        ? config_.base.legacy_mapping.source_carrying_capacity_r
        : config_.base.carrying_capacity_r;
    const double K_capacity = config_.base.thin_layer
        ? config_.base.legacy_mapping.source_carrying_capacity_K
        : config_.base.carrying_capacity_K;
    const std::array<double, 4> volume_weight{1.0, large_cell_volume_,
                                               1.0, large_cell_volume_};

    for (std::size_t location = 0; location < voxel_count_; ++location) {
        const double multiplier = capacity_multiplier(nutrient_[location]);
        const double r_count = r_counts[location] / multiplier;
        const double K_count = K_counts[location] / multiplier;
        const double total_count = r_count + K_count;
        const double r_growth = calculate_density_growth_rate_continuous(
            static_cast<int>(CellType::r), r_inherent,
            r_count, K_count, total_count, r_limit, K_limit,
            config_.base.alpha, config_.base.beta, r_capacity, K_capacity);
        const double K_growth = calculate_density_growth_rate_continuous(
            static_cast<int>(CellType::K), K_inherent,
            r_count, K_count, total_count, r_limit, K_limit,
            config_.base.alpha, config_.base.beta, r_capacity, K_capacity);
        const double r_division_rate = positive_part(r_growth) /
            (config_.base.division_timing.base_cycle_hours *
             std::max(1.0e-12, r_inherent));
        const double K_division_rate = positive_part(K_growth) /
            (config_.base.division_timing.base_cycle_hours *
             std::max(1.0e-12, K_inherent));
        const double r_death_rate = r_growth <= config_.base.death_growth_rate_threshold
            ? 1.0 / config_.base.r_death_delay_hours : 0.0;
        const double K_death_rate = K_growth <= config_.base.death_growth_rate_threshold
            ? 1.0 / config_.base.K_death_delay_hours : 0.0;
        const double occupied = occupied_fraction(location);
        const double vacancy = std::clamp(
            1.0 - occupied / config_.reaction.maximum_occupied_fraction,
            0.0, 1.0);
        const double large_success = std::pow(
            vacancy, config_.reaction.large_daughter_vacancy_exponent);
        const double small_success = std::pow(
            vacancy, config_.reaction.small_daughter_vacancy_exponent);
        const double conversion_probability =
            config_.base.r_to_K_conversion.enabled &&
            occupied >= config_.base.r_to_K_conversion.density_threshold
            ? config_.base.r_to_K_conversion.probability_per_division : 0.0;

        const std::array<double, 4> old{
            populations_[0][location], populations_[1][location],
            populations_[2][location], populations_[3][location]};
        std::array<double, 4> delta{
            -r_death_rate * old[0], -r_death_rate * old[1],
            -K_death_rate * old[2], -K_death_rate * old[3]};

        const double r_large_events = r_division_rate * old[1];
        const double r_large_daughters = large_success * r_large_events;
        const double r_shape_reductions = (1.0 - large_success) * r_large_events;
        delta[1] += (1.0 - conversion_probability) * r_large_daughters
            - r_shape_reductions;
        delta[3] += conversion_probability * r_large_daughters;
        delta[0] += (2.0 - conversion_probability) * r_shape_reductions;
        delta[2] += conversion_probability * r_shape_reductions;

        const double K_large_events = K_division_rate * old[3];
        const double K_large_daughters = large_success * K_large_events;
        const double K_shape_reductions = (1.0 - large_success) * K_large_events;
        delta[3] += K_large_daughters - K_shape_reductions;
        delta[2] += 2.0 * K_shape_reductions;

        const double r_small_events = r_division_rate * old[0];
        const double r_small_births = small_success * r_small_events;
        delta[0] += (1.0 - conversion_probability) * r_small_births
            - (1.0 - small_success) * r_small_events *
                config_.reaction.failed_r_division_death_fraction;
        delta[2] += conversion_probability * r_small_births;
        delta[2] += small_success * K_division_rate * old[2];

        std::array<double, 4> base{};
        std::array<double, 4> positive{};
        double base_occupied = 0.0;
        double positive_occupied = 0.0;
        for (std::size_t field = 0; field < 4; ++field) {
            base[field] = std::max(0.0, old[field] + dt * std::min(0.0, delta[field]));
            positive[field] = dt * std::max(0.0, delta[field]);
            base_occupied += volume_weight[field] * base[field];
            positive_occupied += volume_weight[field] * positive[field];
        }
        const double available = std::max(
            0.0, config_.reaction.maximum_occupied_fraction - base_occupied);
        const double positive_scale = positive_occupied > available &&
                positive_occupied > 0.0
            ? available / positive_occupied : 1.0;
        for (std::size_t field = 0; field < 4; ++field) {
            populations_[field][location] =
                base[field] + positive_scale * positive[field];
        }
    }
}

bool ContinuumModel3D::step() {
    if (!initialized_) throw std::logic_error("continuum model is not initialized");
    if (time_hours_ >= config_.end_time_hours ||
        same_time(time_hours_, config_.end_time_hours)) return false;
    double target = std::min(time_hours_ + config_.time_step_hours,
                             config_.end_time_hours);
    if (next_nutrient_refresh_hours_ < target &&
        !same_time(next_nutrient_refresh_hours_, target)) {
        target = next_nutrient_refresh_hours_;
    }
    const double dt = target - time_hours_;
    migrate(dt);
    react(dt);
    time_hours_ = target;
    ++step_count_;
    if (time_hours_ > next_nutrient_refresh_hours_ ||
        same_time(time_hours_, next_nutrient_refresh_hours_)) {
        solve_nutrient();
        do {
            next_nutrient_refresh_hours_ += config_.nutrient.refresh_every_hours;
        } while (time_hours_ > next_nutrient_refresh_hours_ ||
                 same_time(time_hours_, next_nutrient_refresh_hours_));
    }
    validate_state();
    return true;
}

ContinuumDiagnostics3D ContinuumModel3D::diagnostics() const {
    ContinuumDiagnostics3D result;
    std::array<long double, 2> nutrient_weight{};
    std::array<long double, 2> radius_weight{};
    long double nutrient_sum = 0.0L;
    for (std::size_t location = 0; location < voxel_count_; ++location) {
        const auto point = coordinate(location);
        const double radius = std::sqrt(point[0] * point[0] + point[1] * point[1] +
                                        (config_.base.thin_layer ? 0.0
                                                                 : point[2] * point[2]));
        for (std::size_t field = 0; field < 4; ++field) {
            const double mass = populations_[field][location] * voxel_measure_;
            result.population_mass[field] += mass;
            const std::size_t type = field < 2 ? 0 : 1;
            result.type_mass[type] += mass;
            nutrient_weight[type] += mass * nutrient_[location];
            radius_weight[type] += mass * radius;
        }
        const double occupied = occupied_fraction(location);
        result.occupied_volume += occupied * voxel_measure_;
        result.maximum_occupied_fraction = std::max(
            result.maximum_occupied_fraction, occupied);
        nutrient_sum += nutrient_[location];
        result.maximum_nutrient = std::max(result.maximum_nutrient,
                                           nutrient_[location]);
        result.vessel_volume += vessel_[location] * voxel_measure_;
    }
    result.mean_nutrient = voxel_count_ == 0 ? 0.0
        : static_cast<double>(nutrient_sum / voxel_count_);
    for (std::size_t type = 0; type < 2; ++type) {
        if (result.type_mass[type] > 0.0) {
            result.mean_nutrient_by_type[type] = static_cast<double>(
                nutrient_weight[type] / result.type_mass[type]);
            result.mean_radius_by_type[type] = static_cast<double>(
                radius_weight[type] / result.type_mass[type]);
        }
    }
    return result;
}

void ContinuumModel3D::validate_state() const {
    if (!std::isfinite(time_hours_) || time_hours_ < 0.0 ||
        !std::isfinite(next_nutrient_refresh_hours_) ||
        next_nutrient_refresh_hours_ <= time_hours_) {
        throw std::runtime_error("continuum clock state is invalid");
    }
    for (std::size_t location = 0; location < voxel_count_; ++location) {
        for (const auto& field : populations_) {
            if (!std::isfinite(field[location]) || field[location] < 0.0) {
                throw std::runtime_error("continuum population state is invalid");
            }
        }
        if (occupied_fraction(location) >
                config_.reaction.maximum_occupied_fraction + 1.0e-8 ||
            !std::isfinite(nutrient_[location]) || nutrient_[location] < 0.0 ||
            nutrient_[location] > config_.nutrient.vessel_value + 1.0e-10 ||
            !std::isfinite(vessel_[location]) || vessel_[location] < 0.0 ||
            vessel_[location] > 1.0) {
            throw std::runtime_error("continuum field bounds are invalid");
        }
    }
}

std::uint64_t ContinuumModel3D::state_checksum() const {
    std::uint64_t result = config_.dynamics_fingerprint();
    result = hash_mix(result, std::bit_cast<std::uint64_t>(time_hours_));
    result = hash_mix(result,
                      std::bit_cast<std::uint64_t>(next_nutrient_refresh_hours_));
    result = hash_mix(result, step_count_);
    result = hash_mix(result, nutrient_solve_count_);
    for (const auto& field : populations_) {
        for (const double value : field) {
            result = hash_mix(result, std::bit_cast<std::uint64_t>(value));
        }
    }
    for (const double value : nutrient_) {
        result = hash_mix(result, std::bit_cast<std::uint64_t>(value));
    }
    for (const double value : vessel_) {
        result = hash_mix(result, std::bit_cast<std::uint64_t>(value));
    }
    return result;
}

void ContinuumModel3D::save_checkpoint(const std::filesystem::path& path) const {
    if (!initialized_) throw std::logic_error("cannot checkpoint uninitialized continuum model");
    if (std::filesystem::exists(path)) {
        throw std::runtime_error("refusing to overwrite continuum checkpoint: " +
                                 path.string());
    }
    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    const std::filesystem::path temporary = path.string() + ".tmp";
    std::filesystem::remove(temporary);
    std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
    if (!stream) throw std::runtime_error("unable to create continuum checkpoint");
    stream.write(kCheckpointMagic.data(), kCheckpointMagic.size());
    write_pod(stream, kCheckpointVersion);
    write_pod(stream, config_.dynamics_fingerprint());
    for (const int extent : config_.grid.shape) write_pod(stream, extent);
    write_pod(stream, static_cast<std::uint64_t>(voxel_count_));
    write_pod(stream, std::bit_cast<std::uint64_t>(time_hours_));
    write_pod(stream, std::bit_cast<std::uint64_t>(next_nutrient_refresh_hours_));
    write_pod(stream, step_count_);
    write_pod(stream, nutrient_solve_count_);
    const auto write_field = [&](const std::vector<double>& field) {
        stream.write(reinterpret_cast<const char*>(field.data()),
                     static_cast<std::streamsize>(field.size() * sizeof(double)));
        if (!stream) throw std::runtime_error("unable to write continuum field");
    };
    for (const auto& field : populations_) write_field(field);
    write_field(nutrient_);
    write_field(vessel_);
    stream.flush();
    if (!stream) throw std::runtime_error("unable to finish continuum checkpoint");
    stream.close();
    std::filesystem::rename(temporary, path);
}

void ContinuumModel3D::load_checkpoint(const std::filesystem::path& path) {
    if (initialized_) throw std::logic_error("continuum model is already initialized");
    std::ifstream stream(path, std::ios::binary);
    if (!stream) throw std::runtime_error("unable to open continuum checkpoint: " + path.string());
    std::array<char, 8> magic{};
    stream.read(magic.data(), magic.size());
    if (!stream || magic != kCheckpointMagic ||
        read_pod<std::uint32_t>(stream) != kCheckpointVersion) {
        throw std::runtime_error("unsupported continuum checkpoint format");
    }
    if (read_pod<std::uint64_t>(stream) != config_.dynamics_fingerprint()) {
        throw std::runtime_error("continuum checkpoint configuration mismatch");
    }
    for (const int extent : config_.grid.shape) {
        if (read_pod<int>(stream) != extent) {
            throw std::runtime_error("continuum checkpoint grid mismatch");
        }
    }
    if (read_pod<std::uint64_t>(stream) != voxel_count_) {
        throw std::runtime_error("continuum checkpoint field-size mismatch");
    }
    time_hours_ = std::bit_cast<double>(read_pod<std::uint64_t>(stream));
    next_nutrient_refresh_hours_ =
        std::bit_cast<double>(read_pod<std::uint64_t>(stream));
    step_count_ = read_pod<std::uint64_t>(stream);
    nutrient_solve_count_ = read_pod<std::uint64_t>(stream);
    const auto read_field = [&](std::vector<double>& field) {
        stream.read(reinterpret_cast<char*>(field.data()),
                    static_cast<std::streamsize>(field.size() * sizeof(double)));
        if (!stream) throw std::runtime_error("truncated continuum field");
    };
    for (auto& field : populations_) read_field(field);
    read_field(nutrient_);
    read_field(vessel_);
    if (stream.peek() != std::char_traits<char>::eof()) {
        throw std::runtime_error("continuum checkpoint has trailing data");
    }
    initialized_ = true;
    validate_state();
    if (time_hours_ >= config_.end_time_hours &&
        !same_time(time_hours_, config_.end_time_hours)) {
        throw std::runtime_error("continuum checkpoint is beyond configured end time");
    }
}

}  // namespace atcg3d::continuum
