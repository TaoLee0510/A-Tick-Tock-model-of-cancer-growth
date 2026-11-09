#include "model/structured_pde_model.hpp"
#include "model/operator_math.hpp"
#include "model/beta_duration.hpp"
#include "model/feasible_jump.hpp"
#include "model/shared_resource.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <utility>

#include "common/density_growth_rule.hpp"
#include "engine/parallelism.hpp"
#include "engine/simulation.hpp"
#include "geometry/directions.hpp"
#include "geometry/footprint.hpp"
#include "rules/initial_rates.hpp"

namespace atcg3d::structured_pde {
namespace {

constexpr std::array<char, 8> kCheckpointMagic{
    {'A', 'T', 'C', 'G', 'S', 'P', 'D', '1'}};
constexpr std::uint32_t kRefractoryCheckpointVersion = 3;

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
    if (!stream) throw std::runtime_error("unable to write structured PDE checkpoint");
}

template <class T>
T read_pod(std::istream& stream) {
    static_assert(std::is_trivially_copyable_v<T>);
    T result{};
    stream.read(reinterpret_cast<char*>(&result), sizeof(result));
    if (!stream) throw std::runtime_error("truncated structured PDE checkpoint");
    return result;
}

template <class T>
void write_vector(std::ostream& stream, const T& values) {
    const std::uint64_t size = values.size();
    write_pod(stream, size);
    stream.write(reinterpret_cast<const char*>(values.data()),
                 static_cast<std::streamsize>(size * sizeof(typename T::value_type)));
    if (!stream) throw std::runtime_error("unable to write structured PDE field");
}

template <class T>
void read_vector(std::istream& stream,
                 T& values,
                 std::size_t expected) {
    const std::uint64_t size = read_pod<std::uint64_t>(stream);
    if (size != expected) {
        throw std::runtime_error("structured PDE checkpoint field size mismatch");
    }
    values.resize(expected);
    stream.read(reinterpret_cast<char*>(values.data()),
                static_cast<std::streamsize>(expected * sizeof(typename T::value_type)));
    if (!stream) throw std::runtime_error("truncated structured PDE field");
}

void write_bounds(std::ostream& stream,
                  const StructuredActiveBounds3D& bounds) {
    write_pod(stream, static_cast<std::uint8_t>(bounds.valid));
    for (const int value : {bounds.x0, bounds.y0, bounds.z0,
                            bounds.x1, bounds.y1, bounds.z1}) {
        write_pod(stream, value);
    }
}

StructuredActiveBounds3D read_bounds(std::istream& stream,
                                     int nx,
                                     int ny,
                                     int nz) {
    StructuredActiveBounds3D bounds;
    bounds.valid = read_pod<std::uint8_t>(stream) != 0;
    bounds.x0 = read_pod<int>(stream);
    bounds.y0 = read_pod<int>(stream);
    bounds.z0 = read_pod<int>(stream);
    bounds.x1 = read_pod<int>(stream);
    bounds.y1 = read_pod<int>(stream);
    bounds.z1 = read_pod<int>(stream);
    if (bounds.valid &&
        (bounds.x0 < 0 || bounds.y0 < 0 || bounds.z0 < 0 ||
         bounds.x1 > nx || bounds.y1 > ny || bounds.z1 > nz ||
         bounds.x0 >= bounds.x1 || bounds.y0 >= bounds.y1 ||
         bounds.z0 >= bounds.z1)) {
        throw std::runtime_error("invalid structured checkpoint bounds");
    }
    return bounds;
}

template <class T>
void write_region(std::ostream& stream,
                  const T& values,
                  const StructuredActiveBounds3D& bounds,
                  int nx,
                  int ny) {
    write_pod(stream, bounds_size(bounds));
    if (!bounds.valid) return;
    const std::streamsize row_bytes = static_cast<std::streamsize>(
        (bounds.x1 - bounds.x0) * sizeof(typename T::value_type));
    for (int z = bounds.z0; z < bounds.z1; ++z) {
        for (int y = bounds.y0; y < bounds.y1; ++y) {
            const std::size_t begin =
                (static_cast<std::size_t>(z) * ny + y) * nx + bounds.x0;
            stream.write(reinterpret_cast<const char*>(values.data() + begin),
                         row_bytes);
        }
    }
    if (!stream) throw std::runtime_error("unable to write structured PDE region");
}

template <class T>
void read_region(std::istream& stream,
                 T& values,
                 const StructuredActiveBounds3D& bounds,
                 int nx,
                 int ny) {
    const std::uint64_t stored_size = read_pod<std::uint64_t>(stream);
    if (stored_size != bounds_size(bounds)) {
        throw std::runtime_error("structured PDE checkpoint region size mismatch");
    }
    fill_field(values,typename T::value_type{});
    if (!bounds.valid) return;
    const std::streamsize row_bytes = static_cast<std::streamsize>(
        (bounds.x1 - bounds.x0) * sizeof(typename T::value_type));
    for (int z = bounds.z0; z < bounds.z1; ++z) {
        for (int y = bounds.y0; y < bounds.y1; ++y) {
            const std::size_t begin =
                (static_cast<std::size_t>(z) * ny + y) * nx + bounds.x0;
            stream.read(reinterpret_cast<char*>(values.data() + begin), row_bytes);
        }
    }
    if (!stream) throw std::runtime_error("truncated structured PDE region");
}

}  // namespace

StructuredPdeModel3D::StructuredPdeModel3D(StructuredPdeConfig3D config)
    : config_(std::move(config)), operators_(config_) {
    if (operators_.activation.transported_refractory) {
        config_.continuum.base.static_vasculature =
            config_.continuum.shared_vascular_geometry();
    }
    config_.validate();
    const auto& continuum = config_.continuum;
    voxel_count_ = static_cast<std::size_t>(continuum.grid.shape[0]) *
        static_cast<std::size_t>(continuum.grid.shape[1]) *
        static_cast<std::size_t>(continuum.grid.shape[2]);
    const int dimensions = continuum.base.thin_layer ? 2 : 3;
    double maximum_normal_diffusion = 0.0;
    for (std::size_t stage = 0; stage < 2; ++stage) {
        for (const auto type : {CellType::r, CellType::K}) {
            maximum_normal_diffusion = std::max(maximum_normal_diffusion,
                normal_diffusion(static_cast<StructuredStage3D>(stage), type));
        }
    }
    const double normal_cfl = continuum.time_step_hours * 2.0 * dimensions *
        maximum_normal_diffusion /
        (continuum.grid.spacing_voxels * continuum.grid.spacing_voxels);
    if (!std::isfinite(normal_cfl) || normal_cfl > 0.45 + 1.0e-12) {
        throw std::invalid_argument("structured time step violates ordinary diffusion CFL bound");
    }
    voxel_measure_ = std::pow(continuum.grid.spacing_voxels, dimensions);
    large_cell_volume_ =
        std::pow(continuum.base.large_footprint_edge, dimensions);
    for (DirectionId direction = 1; direction <= 26; ++direction) {
        if (continuum.base.thin_layer && direction_vector(direction).z != 0) continue;
        direction_ids_.push_back(direction);
    }
    const std::size_t buckets = direction_ids_.size() + 1;
    turn_buckets_.resize(buckets);
    for (std::size_t bucket = 1; bucket < buckets; ++bucket) {
        const DirectionId previous = direction_ids_[bucket - 1];
        for (std::size_t candidate = 0; candidate < direction_ids_.size();
             ++candidate) {
            if (candidate + 1 == bucket) continue;
            if (direction_angle_degrees(previous, direction_ids_[candidate]) <=
                continuum.base.turn_half_angle_degrees + 1.0e-10) {
                turn_buckets_[bucket].push_back(candidate + 1);
            }
        }
    }
    if (operators_.transport.directional_sectors && continuum.base.thin_layer) {
        const int edge = operators_.nutrient.transient_resources
            ? config_.migration.direction_nutrient_window_edge
            : config_.migration.direction_density_window_edge;
        const int lower = (edge - 1) / 2;
        const int upper = edge - lower - 1;
        const double minimum_cosine = std::cos(
            continuum.base.direction_density_half_angle_degrees *
            std::acos(-1.0) / 180.0);
        direction_row_spans_.resize(direction_ids_.size());
        direction_sector_site_counts_.assign(direction_ids_.size(), 0U);
        for (std::size_t direction_index = 0;
             direction_index < direction_ids_.size(); ++direction_index) {
            const Vec3i forward = direction_vector(direction_ids_[direction_index]);
            const double forward_length = std::sqrt(
                static_cast<double>(squared_length(forward)));
            for (int dy = -lower; dy <= upper; ++dy) {
                int first = upper + 1;
                int last = -lower - 1;
                std::size_t row_sites = 0;
                for (int dx = -lower; dx <= upper; ++dx) {
                    if (dx == 0 && dy == 0) continue;
                    const Vec3i offset{dx, dy, 0};
                    const double offset_length = std::sqrt(
                        static_cast<double>(squared_length(offset)));
                    const double cosine =
                        static_cast<double>(dot(offset, forward)) /
                        (offset_length * forward_length);
                    if (cosine + 1.0e-12 < minimum_cosine) continue;
                    first = std::min(first, dx);
                    last = std::max(last, dx);
                    ++row_sites;
                    ++direction_sector_site_counts_[direction_index];
                }
                if (first <= last) {
                    if (row_sites != static_cast<std::size_t>(last - first + 1)) {
                        throw std::logic_error(
                            "direction cone row is not a contiguous interval");
                    }
                    direction_row_spans_[direction_index].push_back(
                        {dy, first, last});
                }
            }
        }
    }
    if (config_.sector_mean_model == "prepared_prefix_fft_v2") {
        sector_mean_ = std::make_unique<continuum::SectorMeanField3D>(
            continuum.grid.shape, continuum.base.thin_layer,
            config_.migration.direction_nutrient_window_edge,
            continuum.base.direction_density_half_angle_degrees,
            continuum.base.threads);
    }
    const bool sparse = config_.storage_model != "dense_v1";
    for (std::size_t stage = 0; stage < kStructuredStageCount3D; ++stage) {
        for (auto* field : {&r_normal_[stage], &K_[stage],
                            &r_normal_work_[stage], &K_work_[stage],
                            &r_refractory_[stage], &refractory_clock_[stage],
                            &refractory_work_[stage], &refractory_clock_work_[stage],
                            &activation_density_[stage]}) {
            field->set_sparse(sparse);
        }
        active_total_[stage].set_sparse(sparse);
        activation_cooldown_[stage].set_sparse(sparse);
        activation_armed_[stage].set_sparse(sparse);
        r_normal_[stage].assign(voxel_count_, 0.0);
        if (operators_.activation.transported_refractory) {
            r_refractory_[stage].assign(voxel_count_, 0.0);
            refractory_clock_[stage].assign(voxel_count_, 0.0);
            refractory_work_[stage].assign(voxel_count_, 0.0);
            refractory_clock_work_[stage].assign(voxel_count_, 0.0);
        }
        K_[stage].assign(voxel_count_, 0.0);
        r_normal_work_[stage].assign(voxel_count_, 0.0);
        K_work_[stage].assign(voxel_count_, 0.0);
        activation_density_[stage].assign(voxel_count_, 0.0);
        activation_cooldown_[stage].assign(voxel_count_, 0.0F);
        activation_armed_[stage].assign(voxel_count_, 1U);
        active_total_[stage].assign(voxel_count_, 0.0F);
        active_direction_[stage].resize(buckets);
        active_clock_[stage].resize(buckets);
        for (std::size_t bucket = 0; bucket < buckets; ++bucket) {
            active_direction_[stage][bucket].set_sparse(sparse);
            active_clock_[stage][bucket].set_sparse(sparse);
            active_direction_[stage][bucket].assign(voxel_count_, 0.0F);
            active_clock_[stage][bucket].assign(voxel_count_, 0.0F);
        }
    }
    active_work_.resize(buckets);
    clock_work_.resize(buckets);
    for (std::size_t bucket = 0; bucket < buckets; ++bucket) {
        active_work_[bucket].set_sparse(sparse);
        clock_work_[bucket].set_sparse(sparse);
        active_work_[bucket].assign(voxel_count_, 0.0F);
        clock_work_[bucket].assign(voxel_count_, 0.0F);
    }
    if (operators_.transport.cached_cohort_transport) {
        guidance_weight_cache_.resize(buckets);
        for (auto& field : guidance_weight_cache_) {
            field.set_sparse(sparse);
            field.assign(voxel_count_, 0.0F);
        }
        guidance_weight_cache_stamp_.set_sparse(sparse);
        active_location_stamp_.set_sparse(sparse);
        guidance_weight_cache_stamp_.assign(voxel_count_, 0U);
        active_location_stamp_.assign(voxel_count_, 0U);
    }
    nutrient_.assign(voxel_count_, 0.0);
    nutrient_next_.assign(voxel_count_, 0.0);
    vessel_.assign(voxel_count_, 0.0);
    if(config_.continuum.angiogenesis.model != "disabled") {
        const auto& c=config_.continuum;
        angiogenesis_=std::make_unique<continuum::AngiogenesisField3D>(c.angiogenesis,c.grid.shape,c.grid.spacing_voxels,c.base.thin_layer,c.base.threads);
    }

    tumour_mask_.assign(voxel_count_, 0U);
    if (config_.division_clock_model == "transported_shifted_geometric_v1") {
        renewal_ = std::make_unique<DivisionRenewal3D>(continuum.base.division_timing,
            config_.division_work_bin_width, config_.division_maximum_work,
            std::array<double, 2>{mean_growth_rate(CellType::r), mean_growth_rate(CellType::K)},
            config_.storage_model == "sparse_zero_pages_3d_v2");
    }
    if (config_.migration.activation_clock == "beta_duration_distribution_v2") {
        const auto& base = continuum.base;
        auto kernel = beta_duration_kernel(base.migration_activation_duration_alpha,
            base.migration_activation_duration_beta, base.division_timing.base_cycle_hours / mean_growth_rate(CellType::r),
            config_.migration.activation_time_bin_width_hours, config_.migration.activation_maximum_hours);
        duration_ = std::make_unique<DivisionRenewal3D>(config_.migration.activation_time_bin_width_hours,
            config_.migration.activation_maximum_hours, std::array<std::vector<double>, 2>{kernel, kernel});
    }
    if (config_.migration.activation_rate_model == "beta_rate_distribution_v2") {
        const auto& law = continuum.base.normal_r_migration_beta;
        auto kernel = beta_duration_kernel(law.alpha, law.beta, law.scale,
            law.scale / config_.migration.activation_rate_bins, law.scale);
        velocity_ = std::make_unique<DivisionRenewal3D>(law.scale / config_.migration.activation_rate_bins,
            law.scale, std::array<std::vector<double>, 2>{kernel, kernel});
    }
}

void StructuredPdeModel3D::finish_duration_transport() {
    if (duration_) duration_->finish_transport([&](auto location, auto channel) {
        double mass = 0.0;
        if (channel < 2) for (const auto& bucket : active_direction_[channel]) mass += bucket[location];
        return mass;
    });
}
void StructuredPdeModel3D::finish_velocity_transport() {
    if (velocity_) velocity_->finish_transport([&](auto location, auto channel) {
        double mass = 0.0;
        if (channel < 2) for (const auto& bucket : active_direction_[channel]) mass += bucket[location];
        return mass;
    });
}

double StructuredPdeModel3D::division_channel_mass(std::size_t location,
                                                 std::size_t channel) const {
    const auto stage = channel % 2;
    return channel < 2 ? r_normal_[stage][location] + active_total_[stage][location]
                       : K_[stage][location];
}
void StructuredPdeModel3D::finish_division_transport() {
    if (renewal_) renewal_->finish_transport([&](auto location, auto channel) {
        return division_channel_mass(location, channel);
    });
}

std::size_t StructuredPdeModel3D::index(int x, int y, int z) const noexcept {
    return (static_cast<std::size_t>(z) * config_.continuum.grid.shape[1] +
            static_cast<std::size_t>(y)) * config_.continuum.grid.shape[0] +
        static_cast<std::size_t>(x);
}

bool StructuredPdeModel3D::grid_coordinate(Vec3i site,
                                           int& x,
                                           int& y,
                                           int& z) const noexcept {
    const std::array<double, 3> coordinate{
        static_cast<double>(site.x), static_cast<double>(site.y),
        static_cast<double>(site.z)};
    std::array<int, 3> result{};
    for (std::size_t axis = 0; axis < 3; ++axis) {
        result[axis] = static_cast<int>(std::floor(
            (coordinate[axis] - config_.continuum.grid.origin[axis]) /
            config_.continuum.grid.spacing_voxels));
        if (result[axis] < 0 ||
            result[axis] >= config_.continuum.grid.shape[axis]) return false;
    }
    x = result[0];
    y = result[1];
    z = result[2];
    return true;
}

double StructuredPdeModel3D::r_normal(
    StructuredStage3D stage, std::size_t location) const noexcept {
    return r_normal_[static_cast<std::size_t>(stage)][location];
}

double StructuredPdeModel3D::r_active(
    StructuredStage3D stage, std::size_t location) const noexcept {
    return active_total_[static_cast<std::size_t>(stage)][location];
}

double StructuredPdeModel3D::K(
    StructuredStage3D stage, std::size_t location) const noexcept {
    return K_[static_cast<std::size_t>(stage)][location];
}

double StructuredPdeModel3D::activation_density(
    StructuredStage3D stage, std::size_t location) const noexcept {
    return activation_density_[static_cast<std::size_t>(stage)][location];
}

double StructuredPdeModel3D::occupied_fraction(
    std::size_t location) const noexcept {
    const double small = r_normal_[0][location] + r_active(
        StructuredStage3D::small, location) + K_[0][location];
    const double large = r_normal_[1][location] + r_active(
        StructuredStage3D::large, location) + K_[1][location];
    const double own=small+large_cell_volume_*large;
    return external_occupied_.empty() ? own : own+external_occupied_[location];
}

std::array<double, 3> StructuredPdeModel3D::coordinate(
    std::size_t location) const noexcept {
    const int nx = config_.continuum.grid.shape[0];
    const int ny = config_.continuum.grid.shape[1];
    const int x = static_cast<int>(location % static_cast<std::size_t>(nx));
    const std::size_t yz = location / static_cast<std::size_t>(nx);
    const int y = static_cast<int>(yz % static_cast<std::size_t>(ny));
    const int z = static_cast<int>(yz / static_cast<std::size_t>(ny));
    return {
        config_.continuum.grid.origin[0] +
            (x + 0.5) * config_.continuum.grid.spacing_voxels,
        config_.continuum.grid.origin[1] +
            (y + 0.5) * config_.continuum.grid.spacing_voxels,
        config_.continuum.grid.origin[2] +
            (z + 0.5) * config_.continuum.grid.spacing_voxels};
}

void StructuredPdeModel3D::include_active_location(
    std::size_t stage, int x, int y, int z) noexcept {
    auto& bounds = active_bounds_[stage];
    if (!bounds.valid) {
        bounds = {x, y, z, x + 1, y + 1, z + 1, true};
        return;
    }
    bounds.x0 = std::min(bounds.x0, x);
    bounds.y0 = std::min(bounds.y0, y);
    bounds.z0 = std::min(bounds.z0, z);
    bounds.x1 = std::max(bounds.x1, x + 1);
    bounds.y1 = std::max(bounds.y1, y + 1);
    bounds.z1 = std::max(bounds.z1, z + 1);
}

void StructuredPdeModel3D::include_population_location(
    int x, int y, int z) noexcept {
    if (!population_bounds_.valid) {
        population_bounds_ = {x, y, z, x + 1, y + 1, z + 1, true};
        return;
    }
    population_bounds_.x0 = std::min(population_bounds_.x0, x);
    population_bounds_.y0 = std::min(population_bounds_.y0, y);
    population_bounds_.z0 = std::min(population_bounds_.z0, z);
    population_bounds_.x1 = std::max(population_bounds_.x1, x + 1);
    population_bounds_.y1 = std::max(population_bounds_.y1, y + 1);
    population_bounds_.z1 = std::max(population_bounds_.z1, z + 1);
}

StructuredActiveBounds3D StructuredPdeModel3D::expanded_bounds(
    const StructuredActiveBounds3D& bounds) const noexcept {
    if (!bounds.valid) return {};
    return {
        std::max(0, bounds.x0 - 1),
        std::max(0, bounds.y0 - 1),
        std::max(0, bounds.z0 - (config_.continuum.base.thin_layer ? 0 : 1)),
        std::min(config_.continuum.grid.shape[0], bounds.x1 + 1),
        std::min(config_.continuum.grid.shape[1], bounds.y1 + 1),
        std::min(config_.continuum.grid.shape[2],
                 bounds.z1 + (config_.continuum.base.thin_layer ? 0 : 1)),
        true};
}

void StructuredPdeModel3D::clear_active_work(
    const StructuredActiveBounds3D& bounds) {
    if (!bounds.valid) return;
    const auto clear_row = [&](int y, int z) {
        const std::size_t begin = index(bounds.x0, y, z);
        const std::size_t end = index(bounds.x1 - 1, y, z) + 1;
        for (std::size_t bucket = 0; bucket < active_work_.size(); ++bucket) {
            std::fill(active_work_[bucket].begin() + begin,
                      active_work_[bucket].begin() + end, 0.0F);
            std::fill(clock_work_[bucket].begin() + begin,
                      clock_work_[bucket].begin() + end, 0.0F);
        }
    };
    if (operators_.transport.cached_cohort_transport) {
        const int rows_per_plane = bounds.y1 - bounds.y0;
        const std::size_t row_count = static_cast<std::size_t>(
            rows_per_plane) * static_cast<std::size_t>(bounds.z1 - bounds.z0);
        const int workers = std::max(1, std::min(
            config_.continuum.base.threads, available_worker_threads()));
        deterministic_parallel_for(
            row_count, workers, [&](std::size_t row_offset) {
                const int y = bounds.y0 + static_cast<int>(
                    row_offset % static_cast<std::size_t>(rows_per_plane));
                const int z = bounds.z0 + static_cast<int>(
                    row_offset / static_cast<std::size_t>(rows_per_plane));
                clear_row(y, z);
            });
        return;
    }
    for (int z = bounds.z0; z < bounds.z1; ++z) {
        for (int y = bounds.y0; y < bounds.y1; ++y) {
            clear_row(y, z);
        }
    }
}

void StructuredPdeModel3D::shrink_active_bounds(std::size_t stage) {
    const auto old = active_bounds_[stage];
    StructuredActiveBounds3D next;
    if (!old.valid) return;
    for (int z = old.z0; z < old.z1; ++z) {
        for (int y = old.y0; y < old.y1; ++y) {
            for (int x = old.x0; x < old.x1; ++x) {
                const std::size_t location = index(x, y, z);
                if (active_total_[stage][location] <
                    config_.migration.minimum_density) continue;
                if (!next.valid) {
                    next = {x, y, z, x + 1, y + 1, z + 1, true};
                } else {
                    next.x0 = std::min(next.x0, x);
                    next.y0 = std::min(next.y0, y);
                    next.z0 = std::min(next.z0, z);
                    next.x1 = std::max(next.x1, x + 1);
                    next.y1 = std::max(next.y1, y + 1);
                    next.z1 = std::max(next.z1, z + 1);
                }
            }
        }
    }
    active_bounds_[stage] = next;
}

void StructuredPdeModel3D::shrink_population_bounds() {
    const auto old = population_bounds_;
    StructuredActiveBounds3D next;
    if (!old.valid) return;
    for (int z = old.z0; z < old.z1; ++z) {
        for (int y = old.y0; y < old.y1; ++y) {
            for (int x = old.x0; x < old.x1; ++x) {
                const std::size_t location = index(x, y, z);
                double total = r_normal_[0][location] + r_normal_[1][location] +
                    active_total_[0][location] + active_total_[1][location] +
                    K_[0][location] + K_[1][location];
                if(!external_occupied_.empty()) total+=external_occupied_[location];
                if (total < config_.migration.minimum_density) continue;
                if (!next.valid) {
                    next = {x, y, z, x + 1, y + 1, z + 1, true};
                } else {
                    next.x0 = std::min(next.x0, x);
                    next.y0 = std::min(next.y0, y);
                    next.z0 = std::min(next.z0, z);
                    next.x1 = std::max(next.x1, x + 1);
                    next.y1 = std::max(next.y1, y + 1);
                    next.z1 = std::max(next.z1, z + 1);
                }
            }
        }
    }
    population_bounds_ = next;
}

void StructuredPdeModel3D::initialize_from_abm(
    const Simulation3D& simulation) {
    if (initialized_) throw std::logic_error("structured PDE is already initialized");
    for (auto& field : r_normal_) fill_field(field,0.0);
    for (auto& field : K_) fill_field(field,0.0);
    for (auto& stage : active_direction_) {
        for (auto& bucket : stage) fill_field(bucket,0.0F);
    }
    for (auto& stage : active_clock_) {
        for (auto& bucket : stage) fill_field(bucket,0.0F);
    }
    for (auto& field : active_total_) fill_field(field,0.0F);
    for (auto& field : activation_cooldown_) {
        fill_field(field,0.0F);
    }
    for (auto& field : activation_armed_) {
        fill_field(field,1U);
    }
    active_bounds_ = {};
    population_bounds_ = {};
    fill_field(vessel_,0.0);

    const double small_weight = 1.0 / voxel_measure_;
    const double large_site_weight =
        1.0 / (large_cell_volume_ * voxel_measure_);
    const double now = simulation.clock().time_hours;
    for (const Slot slot : simulation.cells().alive_slots()) {
        const CellType type = simulation.cells().type(slot);
        const CellStage cell_stage = simulation.cells().stage(slot);
        const std::size_t stage = cell_stage == CellStage::large ? 1U : 0U;
        const bool active = type == CellType::r &&
            (simulation.cells().flags(slot) &
             static_cast<std::uint8_t>(kMigrationActive)) != 0;
        std::size_t bucket = 0;
        if (active) {
            const DirectionId last = simulation.cells().last_direction(slot);
            const auto found = std::find(direction_ids_.begin(), direction_ids_.end(), last);
            if (found != direction_ids_.end()) {
                bucket = static_cast<std::size_t>(found - direction_ids_.begin()) + 1;
            }
        }
        const double remaining = active
            ? std::max(0.0,
                  simulation.cells().migration_activation_end_time(slot) - now)
            : 0.0;
        const auto add = [&](Vec3i site, double weight) {
            int x{}, y{}, z{};
            if (!grid_coordinate(site, x, y, z)) {
                throw std::runtime_error(
                    "structured PDE grid does not contain an ABM cell");
            }
            const std::size_t location = index(x, y, z);
            include_population_location(x, y, z);
            if (renewal_) {
                const double completed = std::max(0.0, now - simulation.cells().last_update_time(slot)) *
                    std::max(0.0, double(simulation.cells().density_growth_rate(slot)));
                renewal_->add(location, (type == CellType::r ? 0 : 2) + stage,
                    weight, std::max(0.0, double(simulation.cells().division_work_remaining(slot)) - completed));
            }
            if (type == CellType::K) {
                K_[stage][location] += weight;
            } else if (!active) {
                r_normal_[stage][location] += weight;
            } else {
                if (duration_) duration_->add(location, stage, weight, remaining);
                if (velocity_) velocity_->add(location, stage, weight, simulation.cells().normal_migration_rate(slot));
                active_direction_[stage][bucket][location] +=
                    static_cast<float>(weight);
                active_clock_[stage][bucket][location] +=
                    static_cast<float>(weight * remaining);
                active_total_[stage][location] += static_cast<float>(weight);
                if (operators_.nutrient.transient_resources) {
                    activation_armed_[stage][location] = 0U;
                    activation_cooldown_[stage][location] =
                        static_cast<float>(
                            config_.migration.reactivation_cooldown_hours);
                }
                include_active_location(stage, x, y, z);
            }
        };
        if (cell_stage == CellStage::large) {
            const Vec3i anchor = simulation.cells().anchor(slot);
            for (const Vec3i site : large_footprint(anchor)) {
                if (config_.continuum.base.thin_layer && site.z != anchor.z) continue;
                add(site, large_site_weight);
            }
        } else {
            add(simulation.cells().anchor(slot), small_weight);
        }
    }
    const auto& vascular = config_.continuum.vascular;
    if (vascular.source_mode == "abm_perfusion" ||
        vascular.source_mode == "abm_plus_synthetic_line") {
        for (const Vec3i site : simulation.vessel_grid().occupied_sites()) {
            if (!simulation.vessel_grid().perfused(site)) continue;
            int x{}, y{}, z{};
            if (!grid_coordinate(site, x, y, z)) continue;
            vessel_[index(x, y, z)] = 1.0;
        }
    }
    if (operators_.activation.transported_refractory ||
        vascular.source_mode == "synthetic_central_line" ||
        vascular.source_mode == "abm_plus_synthetic_line") add_synthetic_vessel();
    clear_cells_from_vessels();
    time_hours_ = config_.continuum.initialization_mode == "abm_checkpoint"
        ? now : config_.continuum.start_time_hours;
    if (operators_.nutrient.transient_resources) {
        fill_field(nutrient_,config_.continuum.nutrient.initial_value);
    }
    rebuild_moving_tumour_front();
    solve_nutrient();
    build_activation_density();
    next_nutrient_refresh_hours_ =
        time_hours_ + config_.continuum.nutrient.refresh_every_hours;
    if(angiogenesis_) angiogenesis_->initialize(vessel_);
    initialized_ = true;
    validate_state();
}

void StructuredPdeModel3D::initialize_from_arrays(
    StructuredInitialFields3D fields,
    double time_hours) {
    if (initialized_) throw std::logic_error("structured PDE is already initialized");
    for (std::size_t stage = 0; stage < kStructuredStageCount3D; ++stage) {
        for (const auto* field : {&fields.r_normal[stage], &fields.r_active[stage],
                                  &fields.K[stage],
                                  &fields.active_remaining_hours[stage]}) {
            if (field->size() != voxel_count_) {
                throw std::invalid_argument("structured initial field size mismatch");
            }
        }
    }
    if (fields.vessel_fraction.size() != voxel_count_) {
        throw std::invalid_argument("structured vessel field size mismatch");
    }
    if (!std::isfinite(time_hours) ||
        time_hours + 1.0e-10 < config_.continuum.start_time_hours ||
        time_hours >= config_.continuum.end_time_hours) {
        throw std::invalid_argument("structured initial time is outside run bounds");
    }
    for (std::size_t stage = 0; stage < kStructuredStageCount3D; ++stage) {
        active_bounds_[stage] = {};
        fill_field(active_total_[stage],0.0F);
        fill_field(activation_cooldown_[stage],0.0F);
        fill_field(activation_armed_[stage],1U);
        r_normal_[stage] = std::move(fields.r_normal[stage]);
        if (operators_.activation.transported_refractory) {
            if ((!fields.r_refractory[stage].empty() &&
                 fields.r_refractory[stage].size() != voxel_count_) ||
                (!fields.refractory_remaining_hours[stage].empty() &&
                 fields.refractory_remaining_hours[stage].size() != voxel_count_)) {
                throw std::invalid_argument("structured refractory initial size mismatch");
            }
            for (std::size_t location = 0; location < voxel_count_; ++location) {
                const double mass = fields.r_refractory[stage].empty() ? 0.0 :
                    fields.r_refractory[stage][location];
                const double remaining = fields.refractory_remaining_hours[stage].empty() ? 0.0 :
                    fields.refractory_remaining_hours[stage][location];
                if (!std::isfinite(mass) || !std::isfinite(remaining) || mass < 0.0 ||
                    mass > r_normal_[stage][location] || remaining < 0.0) {
                    throw std::invalid_argument("invalid structured refractory initial state");
                }
                r_refractory_[stage][location] = mass;
                refractory_clock_[stage][location] = mass * remaining;
            }
        }
        K_[stage] = std::move(fields.K[stage]);
        active_direction_[stage][0].assign(voxel_count_, 0.0F);
        active_clock_[stage][0].assign(voxel_count_, 0.0F);
        for (std::size_t location = 0; location < voxel_count_; ++location) {
            const double active = fields.r_active[stage][location];
            const double remaining = fields.active_remaining_hours[stage][location];
            if (active < 0.0 || remaining < 0.0) {
                throw std::invalid_argument("structured active state is negative");
            }
            active_direction_[stage][0][location] = static_cast<float>(active);
            active_clock_[stage][0][location] =
                static_cast<float>(active * remaining);
            active_total_[stage][location] = static_cast<float>(active);
            if (duration_ && active > 0.0) duration_->add(location, stage,
                active_total_[stage][location], remaining);
            if (velocity_ && active > 0.0) velocity_->add_fresh(location, stage, active_total_[stage][location]);
            if (active >= config_.migration.minimum_density) {
                if (operators_.nutrient.transient_resources) {
                    activation_armed_[stage][location] = 0U;
                    activation_cooldown_[stage][location] =
                        static_cast<float>(
                            config_.migration.reactivation_cooldown_hours);
                }
                const int x = static_cast<int>(
                    location % static_cast<std::size_t>(config_.continuum.grid.shape[0]));
                const std::size_t yz = location /
                    static_cast<std::size_t>(config_.continuum.grid.shape[0]);
                const int y = static_cast<int>(yz %
                    static_cast<std::size_t>(config_.continuum.grid.shape[1]));
                const int z = static_cast<int>(yz /
                    static_cast<std::size_t>(config_.continuum.grid.shape[1]));
                include_active_location(stage, x, y, z);
            }
        }
        for (std::size_t bucket = 1; bucket < active_direction_[stage].size();
             ++bucket) {
            fill_field(active_direction_[stage][bucket],0.0F);
            fill_field(active_clock_[stage][bucket],0.0F);
        }
    }
    vessel_ = std::move(fields.vessel_fraction);
    population_bounds_ = {};
    const int nx = config_.continuum.grid.shape[0];
    const int ny = config_.continuum.grid.shape[1];
    for (std::size_t location = 0; location < voxel_count_; ++location) {
        const double total = r_normal_[0][location] + r_normal_[1][location] +
            active_total_[0][location] + active_total_[1][location] +
            K_[0][location] + K_[1][location];
        if (total < config_.migration.minimum_density) continue;
        const int x = static_cast<int>(location % static_cast<std::size_t>(nx));
        const std::size_t yz = location / static_cast<std::size_t>(nx);
        const int y = static_cast<int>(yz % static_cast<std::size_t>(ny));
        const int z = static_cast<int>(yz / static_cast<std::size_t>(ny));
        include_population_location(x, y, z);
    }
    const auto& source = config_.continuum.vascular.source_mode;
    if (operators_.activation.transported_refractory || source == "synthetic_central_line" ||
        source == "abm_plus_synthetic_line") add_synthetic_vessel();
    clear_cells_from_vessels();
    time_hours_ = time_hours;
    if (renewal_) {
        for (std::size_t location = 0; location < voxel_count_; ++location)
            for (std::size_t channel = 0; channel < 4; ++channel) {
                const double mass = division_channel_mass(location, channel);
                if (mass > 0.0) renewal_->add_fresh(location, channel, mass);
            }
    }
    if (operators_.nutrient.transient_resources) {
        fill_field(nutrient_,config_.continuum.nutrient.initial_value);
    }
    rebuild_moving_tumour_front();
    solve_nutrient();
    build_activation_density();
    next_nutrient_refresh_hours_ =
        time_hours_ + config_.continuum.nutrient.refresh_every_hours;
    if(angiogenesis_) angiogenesis_->initialize(vessel_);
    initialized_ = true;
    validate_state();
}

std::array<double, 2> StructuredPdeModel3D::growth_counts_at(Vec3i site) const {
    int x{}, y{}, z{};
    const auto bounds = population_bounds_;
    if (!bounds.valid || !grid_coordinate(site, x, y, z) || x < bounds.x0 || x >= bounds.x1 ||
        y < bounds.y0 || y >= bounds.y1 || z < bounds.z0 || z >= bounds.z1) return {};
    std::vector<double> r_counts, K_counts;
    build_local_counts(r_counts, K_counts);
    const auto offset = (static_cast<std::size_t>(z - bounds.z0) * (bounds.y1 - bounds.y0) +
        (y - bounds.y0)) * (bounds.x1 - bounds.x0) + (x - bounds.x0);
    return {r_counts[offset], K_counts[offset]};
}

double StructuredPdeModel3D::refractory_mass(
    StructuredStage3D stage, std::size_t location) const noexcept {
    return operators_.activation.transported_refractory ? r_refractory_[static_cast<std::size_t>(stage)][location] : 0.0;
}

double StructuredPdeModel3D::refractory_mean_hours(
    StructuredStage3D stage, std::size_t location) const noexcept {
    const double mass = refractory_mass(stage, location);
    return mass > 0.0 ? refractory_clock_[static_cast<std::size_t>(stage)][location] / mass : 0.0;
}

bool StructuredPdeModel3D::step() {
    if (!initialized_) throw std::logic_error("structured PDE is not initialized");
    if (config_.storage_model != "dense_v1" &&
        bounds_size(expanded_bounds(population_bounds_)) > config_.maximum_active_voxels) {
        throw std::runtime_error("sparse PDE active-region memory budget exceeded");
    }
    const auto& continuum = config_.continuum;
    if (time_hours_ >= continuum.end_time_hours ||
        same_time(time_hours_, continuum.end_time_hours)) return false;
    double target = std::min(
        time_hours_ + continuum.time_step_hours, continuum.end_time_hours);
    if (!operators_.nutrient.transient_resources &&
        next_nutrient_refresh_hours_ < target &&
        !same_time(next_nutrient_refresh_hours_, target)) {
        target = next_nutrient_refresh_hours_;
    }
    const double dt = target - time_hours_;
    operators_.vascular.advance(*this, dt);
    if(external_after_vascular_advance_) external_after_vascular_advance_();
    // Density is only the trigger. Active cohorts keep their own remaining
    // clock and therefore do not deactivate when this field later falls.
    if (config_.activation_operator_enabled) operators_.activation.advance(*this, dt);
    if (config_.migration_operator_enabled) {
        operators_.transport.advance(*this, dt);
        if (config_.exchange_operator_enabled) operators_.exchange.advance(*this, dt);
    } else {
        operators_.transport.expire(*this, dt);
    }
    // V4 conversion is evaluated at the post-transport division location.
    // At 200x migration, reusing the pre-transport 70x70 field can otherwise
    // classify mass hundreds of voxels away using its source neighbourhood.
    if (operators_.transport.cached_cohort_transport &&
        continuum.base.r_to_K_conversion.enabled) {
        build_activation_density();
    }
    operators_.reaction.advance(*this, dt);
    // Diffusion expands the conservative work box by one voxel. Remove the
    // numerical halo after all operators have finished so subsequent steps
    // scale with the occupied support rather than the full production grid.
    shrink_population_bounds();
    time_hours_ = target;
    ++step_count_;
    operators_.nutrient.advance(*this, dt);
    validate_state();
    return true;
}

StructuredPdeDiagnostics3D StructuredPdeModel3D::diagnostics() const {
    StructuredPdeDiagnostics3D result;
    if (operators_.vascular.removal_diagnostics)
        result.vascular_removed_mass = vascular_removed_mass_;
    long double nutrient_sum = 0.0L;
    long double r_nutrient = 0.0L;
    long double K_nutrient = 0.0L;
    long double r_radius = 0.0L;
    long double K_radius = 0.0L;
    long double tumour_nutrient = 0.0L;
    long double tumour_front_nutrient = 0.0L;
    double maximum_radius = 0.0;
    for (std::size_t location = 0; location < voxel_count_; ++location) {
        const auto point = coordinate(location);
        maximum_radius = std::max(maximum_radius,
            std::sqrt(point[0] * point[0] + point[1] * point[1] +
                (config_.continuum.base.thin_layer ? 0.0
                                                   : point[2] * point[2])));
    }
    const double width = config_.continuum.grid.spacing_voxels;
    std::vector<long double> radial_mass(
        static_cast<std::size_t>(std::floor(maximum_radius / width)) + 1, 0.0L);
    for (std::size_t location = 0; location < voxel_count_; ++location) {
        if (config_.storage_model != "dense_v1") {
            const int nx = config_.continuum.grid.shape[0];
            const int ny = config_.continuum.grid.shape[1];
            const int x = static_cast<int>(location % nx);
            const int y = static_cast<int>((location / nx) % ny);
            const int z = static_cast<int>(location / (static_cast<std::size_t>(nx) * ny));
            const bool populated = population_bounds_.valid &&
                x >= population_bounds_.x0 && x < population_bounds_.x1 &&
                y >= population_bounds_.y0 && y < population_bounds_.y1 &&
                z >= population_bounds_.z0 && z < population_bounds_.z1;
            if (!populated && tumour_mask_[location] == 0U) {
                // Preserve full-domain resource sums without faulting empty
                // population mappings into resident memory during reporting.
                result.maximum_nutrient = std::max(
                    result.maximum_nutrient, nutrient_[location]);
                nutrient_sum += nutrient_[location];
                result.vessel_volume += vessel_[location] * voxel_measure_;
                continue;
            }
        }
        const auto point = coordinate(location);
        const double radius = std::sqrt(
            point[0] * point[0] + point[1] * point[1] +
            (config_.continuum.base.thin_layer ? 0.0 : point[2] * point[2]));
        double r_here = 0.0;
        double K_here = 0.0;
        for (std::size_t stage = 0; stage < 2; ++stage) {
            const auto stage_value = static_cast<StructuredStage3D>(stage);
            const double normal = r_normal_[stage][location] * voxel_measure_;
            const double active = r_active(stage_value, location) * voxel_measure_;
            const double K_mass = K_[stage][location] * voxel_measure_;
            result.r_normal_mass[stage] += normal;
            result.r_active_mass[stage] += active;
            result.K_mass[stage] += K_mass;
            r_here += normal + active;
            K_here += K_mass;
        }
        result.r_total += r_here;
        result.r_active_total +=
            (r_active(StructuredStage3D::small, location) +
             r_active(StructuredStage3D::large, location)) * voxel_measure_;
        result.K_total += K_here;
        result.occupied_volume += occupied_fraction(location) * voxel_measure_;
        result.maximum_occupied_fraction = std::max(
            result.maximum_occupied_fraction, occupied_fraction(location));
        result.maximum_nutrient =
            std::max(result.maximum_nutrient, nutrient_[location]);
        nutrient_sum += nutrient_[location];
        r_nutrient += r_here * nutrient_[location];
        K_nutrient += K_here * nutrient_[location];
        r_radius += r_here * radius;
        K_radius += K_here * radius;
        radial_mass[static_cast<std::size_t>(std::floor(radius / width))] += r_here;
        result.vessel_volume += vessel_[location] * voxel_measure_;
        if (operators_.nutrient.moving_front && tumour_mask_[location] != 0U) {
            result.tumour_volume += voxel_measure_;
            tumour_nutrient += nutrient_[location];
            const int nx = config_.continuum.grid.shape[0];
            const int ny = config_.continuum.grid.shape[1];
            const int x = static_cast<int>(
                location % static_cast<std::size_t>(nx));
            const int y = static_cast<int>(
                location / static_cast<std::size_t>(nx));
            const bool front =
                (x > 0 && tumour_mask_[index(x - 1, y, 0)] == 0U) ||
                (x + 1 < nx && tumour_mask_[index(x + 1, y, 0)] == 0U) ||
                (y > 0 && tumour_mask_[index(x, y - 1, 0)] == 0U) ||
                (y + 1 < ny && tumour_mask_[index(x, y + 1, 0)] == 0U);
            if (front) {
                result.tumour_front_volume += voxel_measure_;
                tumour_front_nutrient += nutrient_[location];
            }
        }
    }
    result.active_fraction = result.r_total > 0.0
        ? result.r_active_total / result.r_total : 0.0;
    const bool per_cell = config_.continuum.nutrient.consumption_model ==
        "per_cell_ratio_v2";
    const double r_consumers = result.r_normal_mass[0] +
        result.r_active_mass[0] +
        (per_cell ? 1.0 : large_cell_volume_) *
            (result.r_normal_mass[1] + result.r_active_mass[1]);
    const double K_consumers = result.K_mass[0] +
        (per_cell ? 1.0 : large_cell_volume_) * result.K_mass[1];
    result.assembled_r_consumption_per_hour =
        config_.continuum.nutrient.r_consumption_rate_per_hour * r_consumers;
    result.assembled_K_consumption_per_hour =
        config_.continuum.nutrient.K_consumption_rate_per_hour * K_consumers;
    result.mean_nutrient = voxel_count_ > 0
        ? static_cast<double>(nutrient_sum / voxel_count_) : 0.0;
    result.r_mean_nutrient = result.r_total > 0.0
        ? static_cast<double>(r_nutrient / result.r_total) : 0.0;
    result.K_mean_nutrient = result.K_total > 0.0
        ? static_cast<double>(K_nutrient / result.K_total) : 0.0;
    result.tumour_mean_nutrient = result.tumour_volume > 0.0
        ? static_cast<double>(tumour_nutrient /
              (result.tumour_volume / voxel_measure_))
        : 0.0;
    result.tumour_front_mean_nutrient = result.tumour_front_volume > 0.0
        ? static_cast<double>(tumour_front_nutrient /
              (result.tumour_front_volume / voxel_measure_))
        : 0.0;
    result.r_mean_radius = result.r_total > 0.0
        ? static_cast<double>(r_radius / result.r_total) : 0.0;
    result.K_mean_radius = result.K_total > 0.0
        ? static_cast<double>(K_radius / result.K_total) : 0.0;
    const auto quantile = [&](double fraction) {
        const long double target = fraction * result.r_total;
        long double cumulative = 0.0L;
        for (std::size_t bin = 0; bin < radial_mass.size(); ++bin) {
            cumulative += radial_mass[bin];
            if (cumulative >= target) return (bin + 1) * width;
        }
        return maximum_radius;
    };
    if (result.r_total > 0.0) {
        result.r_radius_50 = quantile(0.50);
        result.r_radius_90 = quantile(0.90);
        result.r_radius_99 = quantile(0.99);
    }
    return result;
}

void StructuredPdeModel3D::validate_resources() const {
    const auto& mode = config_.continuum.nutrient.boundary_mode;
    const bool local_validation = initialized_ && operators_.nutrient.moving_front &&
        (mode == "moving_tumor_front_dirichlet_v2" ||
         mode == "moving_tumor_front_and_vessels_dirichlet_v2") &&
        nutrient_update_bounds_.valid;
    const auto bounds = local_validation
        ? nutrient_update_bounds_
        : StructuredActiveBounds3D{
              0, 0, 0,
              config_.continuum.grid.shape[0],
              config_.continuum.grid.shape[1],
              config_.continuum.grid.shape[2], true};
    for (int z = bounds.z0; z < bounds.z1; ++z) {
        for (int y = bounds.y0; y < bounds.y1; ++y) {
            for (int x = bounds.x0; x < bounds.x1; ++x) {
                const std::size_t location = index(x, y, z);
        if (!std::isfinite(nutrient_[location]) || nutrient_[location] < 0.0 ||
            nutrient_[location] >
                config_.continuum.nutrient.vessel_value + 1.0e-10 ||
            !std::isfinite(vessel_[location]) || vessel_[location] < 0.0 ||
            vessel_[location] > 1.0 + 1.0e-10) {
            throw std::runtime_error("invalid structured resource field");
        }
            }
        }
    }
}

void StructuredPdeModel3D::validate_state() const {
    if (!initialized_ && step_count_ != 0) {
        throw std::logic_error("uninitialized structured PDE has steps");
    }
    const double maximum = config_.continuum.reaction.maximum_occupied_fraction;
    if (population_bounds_.valid) {
        const auto bounds = population_bounds_;
        for (int z = bounds.z0; z < bounds.z1; ++z) {
            for (int y = bounds.y0; y < bounds.y1; ++y) {
                for (int x = bounds.x0; x < bounds.x1; ++x) {
                    const std::size_t location = index(x, y, z);
                    for (std::size_t stage = 0; stage < 2; ++stage) {
                        if (operators_.activation.transported_refractory &&
                            (!std::isfinite(r_refractory_[stage][location]) ||
                             !std::isfinite(refractory_clock_[stage][location]) ||
                             r_refractory_[stage][location] < 0.0 ||
                             r_refractory_[stage][location] > r_normal_[stage][location] + 1.0e-10 ||
                             refractory_clock_[stage][location] < 0.0 ||
                             (r_refractory_[stage][location] == 0.0 && refractory_clock_[stage][location] != 0.0))) {
                            throw std::runtime_error("invalid structured cohort refractory state");
                        }
                        if (!std::isfinite(r_normal_[stage][location]) ||
                            r_normal_[stage][location] < -1.0e-10 ||
                            !std::isfinite(K_[stage][location]) ||
                            K_[stage][location] < -1.0e-10 ||
                            !std::isfinite(active_total_[stage][location]) ||
                            active_total_[stage][location] < -1.0e-7F) {
                            throw std::runtime_error(
                                "invalid structured population density");
                        }
                    }
                    if (occupied_fraction(location) > maximum + 5.0e-5) {
                        throw std::runtime_error(
                            "structured PDE exceeded occupied capacity");
                    }
                    if (vessel_blocks_cells(location) &&
                        occupied_fraction(location) > 1.0e-10) {
                        throw std::runtime_error(
                            "structured PDE placed cells inside a vessel");
                    }
                }
            }
        }
    }
    for (std::size_t stage = 0; stage < 2; ++stage) {
        const auto bounds = active_bounds_[stage];
        if (!bounds.valid) continue;
        for (int z = bounds.z0; z < bounds.z1; ++z) {
            for (int y = bounds.y0; y < bounds.y1; ++y) {
                for (int x = bounds.x0; x < bounds.x1; ++x) {
                    const std::size_t location = index(x, y, z);
                    double total = 0.0;
                    for (std::size_t bucket = 0;
                         bucket < active_direction_[stage].size(); ++bucket) {
                        const float mass =
                            active_direction_[stage][bucket][location];
                        const float clock = active_clock_[stage][bucket][location];
                        if (!std::isfinite(mass) || mass < -1.0e-7F ||
                            !std::isfinite(clock) || clock < -1.0e-6F ||
                            (mass == 0.0F && clock != 0.0F)) {
                            throw std::runtime_error("invalid structured active state");
                        }
                        total += mass;
                    }
                    if (std::abs(total - active_total_[stage][location]) >
                        1.0e-5 * std::max(1.0, total)) {
                        throw std::runtime_error("structured active cache mismatch");
                    }
                }
            }
        }
    }
    if (operators_.nutrient.transient_resources) {
        const auto refractory_bounds = initialized_ && step_count_ > 0
            ? population_bounds_
            : StructuredActiveBounds3D{
                  0, 0, 0,
                  config_.continuum.grid.shape[0],
                  config_.continuum.grid.shape[1],
                  config_.continuum.grid.shape[2], true};
        for (std::size_t stage = 0; stage < 2; ++stage) {
            if (!refractory_bounds.valid) continue;
            for (int z = refractory_bounds.z0; z < refractory_bounds.z1; ++z) {
                for (int y = refractory_bounds.y0; y < refractory_bounds.y1; ++y) {
                    for (int x = refractory_bounds.x0;
                         x < refractory_bounds.x1; ++x) {
                const std::size_t location = index(x, y, z);
                if (!std::isfinite(activation_cooldown_[stage][location]) ||
                    activation_cooldown_[stage][location] < 0.0F ||
                    activation_armed_[stage][location] > 1U) {
                    throw std::runtime_error(
                        "invalid structured activation refractory state");
                }
                    }
                }
            }
        }
    }
}

std::uint64_t StructuredPdeModel3D::state_checksum() const {
    std::uint64_t state = config_.dynamics_fingerprint();
    state = hash_mix(state, std::bit_cast<std::uint64_t>(time_hours_));
    state = hash_mix(state,
        std::bit_cast<std::uint64_t>(next_nutrient_refresh_hours_));
    state = hash_mix(state, step_count_);
    state = hash_mix(state, nutrient_solve_count_);
    const auto hash_bounds = [&](const StructuredActiveBounds3D& bounds) {
        state = hash_mix(state, bounds.valid);
        for (const int value : {bounds.x0, bounds.y0, bounds.z0,
                                bounds.x1, bounds.y1, bounds.z1}) {
            state = hash_mix(state, static_cast<std::uint64_t>(value));
        }
    };
    const auto hash_double_region = [&](const auto& values,
                                        const StructuredActiveBounds3D& bounds) {
        if (!bounds.valid) return;
        for (int z = bounds.z0; z < bounds.z1; ++z) {
            for (int y = bounds.y0; y < bounds.y1; ++y) {
                for (int x = bounds.x0; x < bounds.x1; ++x) {
                    state = hash_mix(state, std::bit_cast<std::uint64_t>(
                        static_cast<double>(values[index(x, y, z)])));
                }
            }
        }
    };
    hash_bounds(population_bounds_);
    for (std::size_t stage = 0; stage < 2; ++stage) {
        hash_double_region(r_normal_[stage], population_bounds_);
        hash_double_region(K_[stage], population_bounds_);
        if (operators_.activation.transported_refractory) {
            hash_double_region(r_refractory_[stage], population_bounds_);
            hash_double_region(refractory_clock_[stage], population_bounds_);
        }
        const auto bounds = active_bounds_[stage];
        hash_bounds(bounds);
        if (bounds.valid) {
            for (std::size_t bucket = 0;
                 bucket < active_direction_[stage].size(); ++bucket) {
                for (int z = bounds.z0; z < bounds.z1; ++z) {
                    for (int y = bounds.y0; y < bounds.y1; ++y) {
                        for (int x = bounds.x0; x < bounds.x1; ++x) {
                            const std::size_t location = index(x, y, z);
                            state = hash_mix(state, std::bit_cast<std::uint32_t>(
                                active_direction_[stage][bucket][location]));
                            state = hash_mix(state, std::bit_cast<std::uint32_t>(
                                active_clock_[stage][bucket][location]));
                        }
                    }
                }
            }
        }
    }
    if (operators_.nutrient.transient_resources) {
        for (std::size_t stage = 0; stage < 2; ++stage) {
            for (std::size_t location = 0; location < voxel_count_; ++location) {
                state = hash_mix(state, std::bit_cast<std::uint32_t>(
                    activation_cooldown_[stage][location]));
                state = hash_mix(state,
                    activation_armed_[stage][location]);
            }
        }
    }
    const StructuredActiveBounds3D full_bounds{
        0, 0, 0, config_.continuum.grid.shape[0],
        config_.continuum.grid.shape[1], config_.continuum.grid.shape[2], true};
    hash_double_region(nutrient_, full_bounds);
    hash_double_region(vessel_, full_bounds);
    if(angiogenesis_) state=hash_mix(state,angiogenesis_->checksum());
    if (renewal_) state = hash_mix(state, renewal_->checksum());
    if (duration_) state = hash_mix(state, duration_->checksum());
    if (velocity_) state = hash_mix(state, velocity_->checksum());
    if (operators_.vascular.removal_diagnostics) {
        for (const auto* field : {&vascular_removed_mass_.r_normal,
                                  &vascular_removed_mass_.r_active,
                                  &vascular_removed_mass_.K}) {
            for (const double mass : *field)
                state = hash_mix(state, std::bit_cast<std::uint64_t>(mass));
        }
    }
    if (operators_.nutrient.persistent_workspaces) {
        for (std::size_t stage = 0; stage < 2; ++stage) {
            for (const auto* field : {&r_normal_[stage], &K_[stage],
                                      &r_refractory_[stage], &refractory_clock_[stage]}) {
                hash_double_region(*field, full_bounds);
            }
            for (const auto* fields : {&active_direction_[stage], &active_clock_[stage]}) {
                for (const auto& field : *fields) {
                    hash_double_region(field, full_bounds);
                }
            }
            hash_double_region(active_total_[stage], full_bounds);
        }
        hash_double_region(nutrient_next_, full_bounds);
        hash_bounds(tumour_mask_bounds_);
        hash_bounds(nutrient_update_bounds_);
        for (const auto value : tumour_mask_) {
            state = hash_mix(state, value);
        }
        state = hash_mix(state, tumour_voxel_count_);
        state = hash_mix(state, tumour_front_voxel_count_);
    }
    return state;
}

std::uint32_t StructuredPdeModel3D::checkpoint_version() const noexcept {
    return operators_.checkpoint_version;
}

void StructuredPdeModel3D::save_checkpoint(
    const std::filesystem::path& path) const {
    if (!initialized_) throw std::logic_error("cannot checkpoint uninitialized model");
    if (std::filesystem::exists(path) ||
        std::filesystem::exists(path.string() + ".tmp")) {
        throw std::runtime_error("refusing to overwrite structured checkpoint");
    }
    const std::filesystem::path temporary = path.string() + ".tmp";
    std::filesystem::create_directories(path.parent_path());
    std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
    if (!stream) throw std::runtime_error("unable to create structured checkpoint");
    stream.write(kCheckpointMagic.data(), kCheckpointMagic.size());
    write_pod(stream, checkpoint_version());
    write_pod(stream, config_.dynamics_fingerprint());
    write_pod(stream, static_cast<std::uint64_t>(voxel_count_));
    write_pod(stream, static_cast<std::uint64_t>(direction_ids_.size()));
    write_pod(stream, time_hours_);
    write_pod(stream, next_nutrient_refresh_hours_);
    write_pod(stream, step_count_);
    write_pod(stream, nutrient_solve_count_);
    const int nx = config_.continuum.grid.shape[0];
    const int ny = config_.continuum.grid.shape[1];
    const StructuredActiveBounds3D full_bounds{
        0, 0, 0, nx, ny, config_.continuum.grid.shape[2], true};
    const auto population_fields = operators_.nutrient.persistent_workspaces
        ? full_bounds : population_bounds_;
    write_bounds(stream, population_bounds_);
    for (std::size_t stage = 0; stage < 2; ++stage) {
        write_region(stream, r_normal_[stage], population_fields, nx, ny);
        write_region(stream, K_[stage], population_fields, nx, ny);
        write_bounds(stream, active_bounds_[stage]);
        const auto active_fields = operators_.nutrient.persistent_workspaces
            ? full_bounds : active_bounds_[stage];
        for (const auto& bucket : active_direction_[stage]) {
            write_region(stream, bucket, active_fields, nx, ny);
        }
        for (const auto& bucket : active_clock_[stage]) {
            write_region(stream, bucket, active_fields, nx, ny);
        }
    }
    if (operators_.nutrient.transient_resources) {
        for (std::size_t stage = 0; stage < 2; ++stage) {
            write_vector(stream, activation_cooldown_[stage]);
            write_vector(stream, activation_armed_[stage]);
        }
    }
    if (operators_.activation.transported_refractory) {
        for (std::size_t stage = 0; stage < 2; ++stage) {
            write_region(stream, r_refractory_[stage], population_fields, nx, ny);
            write_region(stream, refractory_clock_[stage], population_fields, nx, ny);
        }
    }
    write_vector(stream, nutrient_);
    write_vector(stream, vessel_);
    if(angiogenesis_) angiogenesis_->save(stream);
    if (renewal_) renewal_->save(stream);
    if (duration_) duration_->save(stream);
    if (velocity_) velocity_->save(stream);
    if (operators_.vascular.removal_diagnostics) {
        for (const auto* field : {&vascular_removed_mass_.r_normal,
                                  &vascular_removed_mass_.r_active,
                                  &vascular_removed_mass_.K}) {
            for (const double mass : *field) write_pod(stream, mass);
        }
    }
    if (operators_.nutrient.persistent_workspaces) {
        for (const auto& field : active_total_) {
            write_vector(stream, field);
        }
        write_vector(stream, nutrient_next_);
        write_bounds(stream, tumour_mask_bounds_);
        write_bounds(stream, nutrient_update_bounds_);
        write_vector(stream, tumour_mask_);
        write_pod(stream, tumour_voxel_count_);
        write_pod(stream, tumour_front_voxel_count_);
    }
    write_pod(stream, state_checksum());
    stream.flush();
    if (!stream) throw std::runtime_error("unable to finish structured checkpoint");
    stream.close();
    std::filesystem::rename(temporary, path);
}

void StructuredPdeModel3D::load_checkpoint(
    const std::filesystem::path& path) {
    if (initialized_) throw std::logic_error("structured PDE is already initialized");
    std::ifstream stream(path, std::ios::binary);
    if (!stream) throw std::runtime_error("unable to open structured checkpoint");
    std::array<char, 8> magic{};
    stream.read(magic.data(), magic.size());
    const std::uint32_t saved_version = read_pod<std::uint32_t>(stream);
    if (magic != kCheckpointMagic || saved_version != checkpoint_version()) {
        throw std::runtime_error("unsupported structured checkpoint format");
    }
    if (read_pod<std::uint64_t>(stream) != config_.dynamics_fingerprint()) {
        throw std::runtime_error("structured checkpoint configuration mismatch");
    }
    if (read_pod<std::uint64_t>(stream) != voxel_count_ ||
        read_pod<std::uint64_t>(stream) != direction_ids_.size()) {
        throw std::runtime_error("structured checkpoint grid mismatch");
    }
    time_hours_ = read_pod<double>(stream);
    next_nutrient_refresh_hours_ = read_pod<double>(stream);
    step_count_ = read_pod<std::uint64_t>(stream);
    nutrient_solve_count_ = read_pod<std::uint64_t>(stream);
    const int nx = config_.continuum.grid.shape[0];
    const int ny = config_.continuum.grid.shape[1];
    const int nz = config_.continuum.grid.shape[2];
    const StructuredActiveBounds3D full_bounds{0, 0, 0, nx, ny, nz, true};
    population_bounds_ = read_bounds(stream, nx, ny, nz);
    const auto population_fields = operators_.nutrient.persistent_workspaces
        ? full_bounds : population_bounds_;
    for (std::size_t stage = 0; stage < 2; ++stage) {
        read_region(stream, r_normal_[stage], population_fields, nx, ny);
        read_region(stream, K_[stage], population_fields, nx, ny);
        active_bounds_[stage] = read_bounds(stream, nx, ny, nz);
        const auto active_fields = operators_.nutrient.persistent_workspaces
            ? full_bounds : active_bounds_[stage];
        for (auto& bucket : active_direction_[stage]) {
            read_region(stream, bucket, active_fields, nx, ny);
        }
        for (auto& bucket : active_clock_[stage]) {
            read_region(stream, bucket, active_fields, nx, ny);
        }
        fill_field(active_total_[stage],0.0F);
        const auto bounds = active_bounds_[stage];
        if (bounds.valid) {
            for (int z = bounds.z0; z < bounds.z1; ++z) {
                for (int y = bounds.y0; y < bounds.y1; ++y) {
                    for (int x = bounds.x0; x < bounds.x1; ++x) {
                        const std::size_t location = index(x, y, z);
                        double total = 0.0;
                        for (const auto& bucket : active_direction_[stage]) {
                            total += bucket[location];
                        }
                        active_total_[stage][location] = static_cast<float>(total);
                    }
                }
            }
        }
    }
    if (saved_version >= kRefractoryCheckpointVersion) {
        for (std::size_t stage = 0; stage < 2; ++stage) {
            read_vector(stream, activation_cooldown_[stage], voxel_count_);
            read_vector(stream, activation_armed_[stage], voxel_count_);
        }
    }
    if (operators_.activation.transported_refractory) {
        for (std::size_t stage = 0; stage < 2; ++stage) {
            read_region(stream, r_refractory_[stage], population_fields, nx, ny);
            read_region(stream, refractory_clock_[stage], population_fields, nx, ny);
        }
    }
    read_vector(stream, nutrient_, voxel_count_);
    read_vector(stream, vessel_, voxel_count_);
    if(angiogenesis_) angiogenesis_->load(stream);
    if (renewal_) renewal_->load(stream, voxel_count_);
    if (duration_) duration_->load(stream, voxel_count_);
    if (velocity_) velocity_->load(stream, voxel_count_);
    if (operators_.vascular.removal_diagnostics) {
        for (auto* field : {&vascular_removed_mass_.r_normal,
                            &vascular_removed_mass_.r_active,
                            &vascular_removed_mass_.K}) {
            for (double& mass : *field) {
                mass = read_pod<double>(stream);
                if (!std::isfinite(mass) || mass < 0.0)
                    throw std::runtime_error("invalid vascular removed mass");
            }
        }
    }
    if (operators_.nutrient.persistent_workspaces) {
        for (auto& field : active_total_) {
            read_vector(stream, field, voxel_count_);
        }
        read_vector(stream, nutrient_next_, voxel_count_);
        tumour_mask_bounds_ = read_bounds(stream, nx, ny, nz);
        nutrient_update_bounds_ = read_bounds(stream, nx, ny, nz);
        read_vector(stream, tumour_mask_, voxel_count_);
        tumour_voxel_count_ = read_pod<decltype(tumour_voxel_count_)>(stream);
        tumour_front_voxel_count_ = read_pod<decltype(tumour_front_voxel_count_)>(stream);
    }
    const std::uint64_t expected_checksum = read_pod<std::uint64_t>(stream);
    if (stream.peek() != std::char_traits<char>::eof()) {
        throw std::runtime_error("structured checkpoint has trailing data");
    }
    if (!operators_.nutrient.persistent_workspaces) {
        rebuild_moving_tumour_front();
        nutrient_next_ = nutrient_;
    }
    validate_resources();
    initialized_ = true;
    validate_state();
    if (state_checksum() != expected_checksum) {
        throw std::runtime_error("structured checkpoint checksum mismatch");
    }
}

}  // namespace atcg3d::structured_pde
