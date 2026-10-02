#include "model/structured_pde_model.hpp"
#include "model/shared_resource.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cmath>
#include <fstream>
#include <limits>
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

namespace atcg3d::structured_pde {
namespace {

constexpr std::array<char, 8> kCheckpointMagic{
    {'A', 'T', 'C', 'G', 'S', 'P', 'D', '1'}};
constexpr std::uint32_t kCohortCheckpointVersion = 4;
constexpr std::uint32_t kLegacyCheckpointVersion = 2;
constexpr std::uint32_t kRefractoryCheckpointVersion = 3;
constexpr double kFixed26DiffusionFactor = 9.0 / 26.0;

bool same_time(double lhs, double rhs) noexcept {
    return std::abs(lhs - rhs) <=
        1.0e-10 * std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

double beta_mean(const BetaRateConfig& config) noexcept {
    return config.scale * config.alpha / (config.alpha + config.beta);
}

double positive_part(double value) noexcept {
    return std::max(0.0, value);
}

int floor_div(int value, int divisor) noexcept {
    int quotient = value / divisor;
    const int remainder = value % divisor;
    if (remainder != 0 && ((remainder < 0) != (divisor < 0))) --quotient;
    return quotient;
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
void write_vector(std::ostream& stream, const std::vector<T>& values) {
    const std::uint64_t size = values.size();
    write_pod(stream, size);
    stream.write(reinterpret_cast<const char*>(values.data()),
                 static_cast<std::streamsize>(size * sizeof(T)));
    if (!stream) throw std::runtime_error("unable to write structured PDE field");
}

template <class T>
void read_vector(std::istream& stream,
                 std::vector<T>& values,
                 std::size_t expected) {
    const std::uint64_t size = read_pod<std::uint64_t>(stream);
    if (size != expected) {
        throw std::runtime_error("structured PDE checkpoint field size mismatch");
    }
    values.resize(expected);
    stream.read(reinterpret_cast<char*>(values.data()),
                static_cast<std::streamsize>(expected * sizeof(T)));
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

std::uint64_t bounds_size(const StructuredActiveBounds3D& bounds) noexcept {
    if (!bounds.valid) return 0;
    return static_cast<std::uint64_t>(bounds.x1 - bounds.x0) *
        static_cast<std::uint64_t>(bounds.y1 - bounds.y0) *
        static_cast<std::uint64_t>(bounds.z1 - bounds.z0);
}

template <class T>
void write_region(std::ostream& stream,
                  const std::vector<T>& values,
                  const StructuredActiveBounds3D& bounds,
                  int nx,
                  int ny) {
    write_pod(stream, bounds_size(bounds));
    if (!bounds.valid) return;
    const std::streamsize row_bytes = static_cast<std::streamsize>(
        (bounds.x1 - bounds.x0) * sizeof(T));
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
                 std::vector<T>& values,
                 const StructuredActiveBounds3D& bounds,
                 int nx,
                 int ny) {
    const std::uint64_t stored_size = read_pod<std::uint64_t>(stream);
    if (stored_size != bounds_size(bounds)) {
        throw std::runtime_error("structured PDE checkpoint region size mismatch");
    }
    std::fill(values.begin(), values.end(), T{});
    if (!bounds.valid) return;
    const std::streamsize row_bytes = static_cast<std::streamsize>(
        (bounds.x1 - bounds.x0) * sizeof(T));
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
    : config_(std::move(config)) {
    if (config_.schema_version >= 7) {
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
    if (config_.schema_version >= 3 && continuum.base.thin_layer) {
        const int edge = config_.schema_version >= 5
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
    for (std::size_t stage = 0; stage < kStructuredStageCount3D; ++stage) {
        r_normal_[stage].assign(voxel_count_, 0.0);
        if (config_.schema_version >= 7) {
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
            active_direction_[stage][bucket].assign(voxel_count_, 0.0F);
            active_clock_[stage][bucket].assign(voxel_count_, 0.0F);
        }
    }
    active_work_.resize(buckets);
    clock_work_.resize(buckets);
    for (std::size_t bucket = 0; bucket < buckets; ++bucket) {
        active_work_[bucket].assign(voxel_count_, 0.0F);
        clock_work_[bucket].assign(voxel_count_, 0.0F);
    }
    if (config_.schema_version >= 4) {
        guidance_weight_cache_.resize(buckets);
        for (auto& field : guidance_weight_cache_) {
            field.assign(voxel_count_, 0.0F);
        }
        guidance_weight_cache_stamp_.assign(voxel_count_, 0U);
        active_location_stamp_.assign(voxel_count_, 0U);
    }
    nutrient_.assign(voxel_count_, 0.0);
    nutrient_next_.assign(voxel_count_, 0.0);
    vessel_.assign(voxel_count_, 0.0);
    tumour_mask_.assign(voxel_count_, 0U);
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
    return small + large_cell_volume_ * large;
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

namespace {

StructuredActiveBounds3D union_bounds(
    const StructuredActiveBounds3D& lhs,
    const StructuredActiveBounds3D& rhs) noexcept {
    if (!lhs.valid) return rhs;
    if (!rhs.valid) return lhs;
    return {
        std::min(lhs.x0, rhs.x0),
        std::min(lhs.y0, rhs.y0),
        std::min(lhs.z0, rhs.z0),
        std::max(lhs.x1, rhs.x1),
        std::max(lhs.y1, rhs.y1),
        std::max(lhs.z1, rhs.z1),
        true};
}

}  // namespace

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
    if (config_.schema_version >= 4) {
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
                const double total = r_normal_[0][location] + r_normal_[1][location] +
                    active_total_[0][location] + active_total_[1][location] +
                    K_[0][location] + K_[1][location];
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
    for (auto& field : r_normal_) std::fill(field.begin(), field.end(), 0.0);
    for (auto& field : K_) std::fill(field.begin(), field.end(), 0.0);
    for (auto& stage : active_direction_) {
        for (auto& bucket : stage) std::fill(bucket.begin(), bucket.end(), 0.0F);
    }
    for (auto& stage : active_clock_) {
        for (auto& bucket : stage) std::fill(bucket.begin(), bucket.end(), 0.0F);
    }
    for (auto& field : active_total_) std::fill(field.begin(), field.end(), 0.0F);
    for (auto& field : activation_cooldown_) {
        std::fill(field.begin(), field.end(), 0.0F);
    }
    for (auto& field : activation_armed_) {
        std::fill(field.begin(), field.end(), 1U);
    }
    active_bounds_ = {};
    population_bounds_ = {};
    std::fill(vessel_.begin(), vessel_.end(), 0.0);

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
            if (type == CellType::K) {
                K_[stage][location] += weight;
            } else if (!active) {
                r_normal_[stage][location] += weight;
            } else {
                active_direction_[stage][bucket][location] +=
                    static_cast<float>(weight);
                active_clock_[stage][bucket][location] +=
                    static_cast<float>(weight * remaining);
                active_total_[stage][location] += static_cast<float>(weight);
                if (config_.schema_version >= 5) {
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
    if (config_.schema_version >= 7 ||
        vascular.source_mode == "synthetic_central_line" ||
        vascular.source_mode == "abm_plus_synthetic_line") add_synthetic_vessel();
    clear_cells_from_vessels();
    time_hours_ = config_.continuum.initialization_mode == "abm_checkpoint"
        ? now : config_.continuum.start_time_hours;
    if (config_.schema_version >= 5) {
        std::fill(nutrient_.begin(), nutrient_.end(),
                  config_.continuum.nutrient.initial_value);
    }
    rebuild_moving_tumour_front();
    solve_nutrient();
    build_activation_density();
    next_nutrient_refresh_hours_ =
        time_hours_ + config_.continuum.nutrient.refresh_every_hours;
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
        std::fill(active_total_[stage].begin(), active_total_[stage].end(), 0.0F);
        std::fill(activation_cooldown_[stage].begin(),
                  activation_cooldown_[stage].end(), 0.0F);
        std::fill(activation_armed_[stage].begin(),
                  activation_armed_[stage].end(), 1U);
        r_normal_[stage] = std::move(fields.r_normal[stage]);
        if (config_.schema_version >= 7) {
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
            if (active >= config_.migration.minimum_density) {
                if (config_.schema_version >= 5) {
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
            std::fill(active_direction_[stage][bucket].begin(),
                      active_direction_[stage][bucket].end(), 0.0F);
            std::fill(active_clock_[stage][bucket].begin(),
                      active_clock_[stage][bucket].end(), 0.0F);
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
    if (config_.schema_version >= 7 || source == "synthetic_central_line" ||
        source == "abm_plus_synthetic_line") add_synthetic_vessel();
    clear_cells_from_vessels();
    time_hours_ = time_hours;
    if (config_.schema_version >= 5) {
        std::fill(nutrient_.begin(), nutrient_.end(),
                  config_.continuum.nutrient.initial_value);
    }
    rebuild_moving_tumour_front();
    solve_nutrient();
    build_activation_density();
    next_nutrient_refresh_hours_ =
        time_hours_ + config_.continuum.nutrient.refresh_every_hours;
    initialized_ = true;
    validate_state();
}

void StructuredPdeModel3D::add_synthetic_vessel() {
    if (config_.schema_version >= 7) {
        const auto geometry = config_.continuum.shared_vascular_geometry();
        const int nx = geometry.shape[0];
        const int ny = geometry.shape[1];
        for (std::size_t location = 0; location < voxel_count_; ++location) {
            const int x = static_cast<int>(location % nx);
            const auto yz = location / nx;
            const int y = static_cast<int>(yz % ny);
            const int z = static_cast<int>(yz / ny);
            if (geometry.source_voxel({x, y, z})) vessel_[location] = 1.0;
        }
        return;
    }
    const auto& vascular = config_.continuum.vascular;
    const int axis = vascular.synthetic_axis == "x" ? 0
        : (vascular.synthetic_axis == "y" ? 1 : 2);
    const double radius_squared = vascular.synthetic_radius_voxels *
        vascular.synthetic_radius_voxels;
    for (std::size_t location = 0; location < voxel_count_; ++location) {
        const auto point = coordinate(location);
        double distance_squared = 0.0;
        for (int dimension = 0; dimension < 3; ++dimension) {
            if (dimension == axis ||
                (config_.continuum.base.thin_layer && dimension == 2)) continue;
            const double offset =
                point[dimension] - vascular.synthetic_center[dimension];
            distance_squared += offset * offset;
        }
        if (distance_squared <= radius_squared) vessel_[location] = 1.0;
    }
}

bool StructuredPdeModel3D::vessel_blocks_cells(
    std::size_t location) const noexcept {
    return config_.migration.vessel_exclusion && vessel_[location] > 0.0;
}

void StructuredPdeModel3D::clear_cells_from_vessels() {
    if (!config_.migration.vessel_exclusion) return;
    for (std::size_t location = 0; location < voxel_count_; ++location) {
        if (!vessel_blocks_cells(location)) continue;
        if (config_.schema_version >= 7 && occupied_fraction(location) > 1.0e-10) {
            throw std::invalid_argument("structured initial population overlaps an excluded vessel");
        }
        for (std::size_t stage = 0; stage < kStructuredStageCount3D; ++stage) {
            r_normal_[stage][location] = 0.0;
            if (config_.schema_version >= 7) {
                r_refractory_[stage][location] = 0.0;
                refractory_clock_[stage][location] = 0.0;
            }
            K_[stage][location] = 0.0;
            active_total_[stage][location] = 0.0F;
            activation_cooldown_[stage][location] = 0.0F;
            activation_armed_[stage][location] = 0U;
            for (std::size_t bucket = 0;
                 bucket < active_direction_[stage].size(); ++bucket) {
                active_direction_[stage][bucket][location] = 0.0F;
                active_clock_[stage][bucket][location] = 0.0F;
            }
        }
    }
    for (std::size_t stage = 0; stage < kStructuredStageCount3D; ++stage) {
        shrink_active_bounds(stage);
    }
    shrink_population_bounds();
}

double StructuredPdeModel3D::capacity_multiplier(double value) const noexcept {
    if (config_.schema_version >= 5) return 1.0;
    const auto& nutrient = config_.continuum.nutrient;
    const double local = std::clamp(value, 0.0, nutrient.vessel_value);
    const double raw = local / (nutrient.capacity_half_saturation + local);
    const double at_vessel = nutrient.vessel_value /
        (nutrient.capacity_half_saturation + nutrient.vessel_value);
    const double saturation = at_vessel > 0.0
        ? std::clamp(raw / at_vessel, 0.0, 1.0) : 0.0;
    return 1.0 + (nutrient.maximum_capacity_multiplier - 1.0) * saturation;
}

void StructuredPdeModel3D::rebuild_moving_tumour_front() {
    if (config_.schema_version < 6) return;
    const auto& nutrient = config_.continuum.nutrient;
    const auto& mode = nutrient.boundary_mode;
    if (mode != "moving_tumor_front_dirichlet_v2" &&
        mode != "moving_tumor_front_and_vessels_dirichlet_v2") return;

    const int nx = config_.continuum.grid.shape[0];
    const int ny = config_.continuum.grid.shape[1];
    const StructuredActiveBounds3D previous_bounds = tumour_mask_bounds_;
    if (tumour_mask_bounds_.valid) {
        for (int y = tumour_mask_bounds_.y0; y < tumour_mask_bounds_.y1; ++y) {
            for (int x = tumour_mask_bounds_.x0; x < tumour_mask_bounds_.x1;
                 ++x) {
                tumour_mask_[index(x, y, 0)] = 0U;
            }
        }
    }
    tumour_mask_bounds_ = {};
    tumour_voxel_count_ = 0U;
    tumour_front_voxel_count_ = 0U;
    if (!population_bounds_.valid) {
        nutrient_update_bounds_ = previous_bounds;
        return;
    }

    const int margin = nutrient.tumor_front_smoothing_radius_voxels + 2;
    const int x0 = std::max(0, population_bounds_.x0 - margin);
    const int y0 = std::max(0, population_bounds_.y0 - margin);
    const int x1 = std::min(nx, population_bounds_.x1 + margin);
    const int y1 = std::min(ny, population_bounds_.y1 + margin);
    const int width = x1 - x0;
    const int height = y1 - y0;
    if (width <= 0 || height <= 0) {
        nutrient_update_bounds_ = previous_bounds;
        return;
    }

    tumour_occupancy_work_.resize(
        static_cast<std::size_t>(width) * height);
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            tumour_occupancy_work_[
                static_cast<std::size_t>(y) * width + x] =
                occupied_fraction(index(x0 + x, y0 + y, 0));
        }
    }
    const auto summary = continuum::build_moving_tumor_front_mask_2d(
        tumour_occupancy_work_, width, height,
        nutrient.tumor_front_smoothing_radius_voxels,
        nutrient.tumor_front_density_threshold,
        tumour_front_workspace_, tumour_local_mask_work_);
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            if (tumour_local_mask_work_[
                    static_cast<std::size_t>(y) * width + x] == 0U) continue;
            tumour_mask_[index(x0 + x, y0 + y, 0)] = 1U;
        }
    }
    tumour_mask_bounds_ = {x0, y0, 0, x1, y1, 1, true};
    nutrient_update_bounds_ = union_bounds(previous_bounds, tumour_mask_bounds_);
    tumour_voxel_count_ = summary.tumour_voxels;
    tumour_front_voxel_count_ = summary.front_voxels;
}

bool StructuredPdeModel3D::nutrient_source(
    int x, int y, int z, std::size_t location) const noexcept {
    if (config_.schema_version < 5) return false;
    const auto& mode = config_.continuum.nutrient.boundary_mode;
    if (config_.schema_version >= 6 &&
        (mode == "moving_tumor_front_dirichlet_v2" ||
         mode == "moving_tumor_front_and_vessels_dirichlet_v2")) {
        // The in-grid host exterior is maintained at Nmax. When tumour reaches
        // the computational box, tumour voxels on the box edge are not sources
        // and the existing ghost-cell rule therefore remains zero-flux.
        const bool use_vessels = mode ==
            "moving_tumor_front_and_vessels_dirichlet_v2";
        return tumour_mask_[location] == 0U ||
            (use_vessels && vessel_[location] > 0.0);
    }
    const auto& grid = config_.continuum.grid;
    const bool planar_edge = x == 0 || x + 1 == grid.shape[0] ||
        y == 0 || y + 1 == grid.shape[1] ||
        (!config_.continuum.base.thin_layer &&
         (z == 0 || z + 1 == grid.shape[2]));
    const bool use_edges = mode != "vessels_dirichlet_v1";
    const bool use_vessels = mode != "planar_edges_dirichlet_v1";
    return (use_edges && planar_edge) ||
        (use_vessels && vessel_[location] > 0.0);
}

void StructuredPdeModel3D::advance_transient_nutrient(double dt) {
    const auto& continuum = config_.continuum;
    const auto& nutrient = continuum.nutrient;
    const int nx = continuum.grid.shape[0];
    const int ny = continuum.grid.shape[1];
    const int nz = continuum.grid.shape[2];
    const int dimensions = continuum.base.thin_layer ? 2 : 3;
    const double inverse_h2 = 1.0 /
        (continuum.grid.spacing_voxels * continuum.grid.spacing_voxels);
    const double mu = nutrient.diffusion_voxels2_per_hour * dt * inverse_h2;
    if (mu > 1.0 / (2.0 * dimensions) + 1.0e-12) {
        throw std::runtime_error(
            "transient nutrient diffusion violates the explicit CFL limit");
    }
    const double source_value = nutrient.vessel_value;
    const double decay = std::exp(-nutrient.decay_per_hour * dt);
    const int workers = std::max(
        1, std::min(continuum.base.threads, available_worker_threads()));
    const bool moving_front = config_.schema_version >= 6 &&
        (nutrient.boundary_mode == "moving_tumor_front_dirichlet_v2" ||
         nutrient.boundary_mode ==
             "moving_tumor_front_and_vessels_dirichlet_v2");
    const StructuredActiveBounds3D update_bounds = moving_front
        ? nutrient_update_bounds_
        : StructuredActiveBounds3D{0, 0, 0, nx, ny, nz, true};
    if (!update_bounds.valid) {
        ++nutrient_solve_count_;
        return;
    }
    const std::size_t bx = static_cast<std::size_t>(
        update_bounds.x1 - update_bounds.x0);
    const std::size_t by = static_cast<std::size_t>(
        update_bounds.y1 - update_bounds.y0);
    const std::size_t bz = static_cast<std::size_t>(
        update_bounds.z1 - update_bounds.z0);
    const std::size_t update_count = bx * by * bz;
    const auto coordinates = [&](std::size_t offset, int& x, int& y, int& z) {
        x = update_bounds.x0 + static_cast<int>(offset % bx);
        const std::size_t local_yz = offset / bx;
        y = update_bounds.y0 + static_cast<int>(local_yz % by);
        z = update_bounds.z0 + static_cast<int>(local_yz / by);
    };
    deterministic_parallel_for(update_count, workers, [&](std::size_t offset) {
        int x{}, y{}, z{};
        coordinates(offset, x, y, z);
        const std::size_t here = index(x, y, z);
        if (nutrient_source(x, y, z, here)) {
            nutrient_next_[here] = source_value;
            return;
        }

        const double old = nutrient_[here];
        // A boundary that is not selected as a fixed nutrient source is
        // reflecting. The centre-valued ghost cell gives zero normal flux and
        // avoids any domain-exterior access in vessel-only controls.
        double neighbor_sum =
            (x > 0 ? nutrient_[index(x - 1, y, z)] : old) +
            (x + 1 < nx ? nutrient_[index(x + 1, y, z)] : old) +
            (y > 0 ? nutrient_[index(x, y - 1, z)] : old) +
            (y + 1 < ny ? nutrient_[index(x, y + 1, z)] : old);
        if (!continuum.base.thin_layer) {
            neighbor_sum +=
                (z > 0 ? nutrient_[index(x, y, z - 1)] : old) +
                (z + 1 < nz ? nutrient_[index(x, y, z + 1)] : old);
        }
        const double diffused = std::clamp(
            old + mu * (neighbor_sum - 2.0 * dimensions * old),
            0.0, source_value) * decay;

        // Schema v3 consumption is per biological cell, independent of
        // phenotype and footprint. Solve the local Michaelis-Menten sink
        // implicitly so consumption cannot drive the field negative.
        const double consumers = r_normal_[0][here] +
            r_active(StructuredStage3D::small, here) +
            r_normal_[1][here] +
            r_active(StructuredStage3D::large, here) +
            K_[0][here] + K_[1][here];
        const double demand =
            nutrient.K_consumption_rate_per_hour * consumers;
        const double half = nutrient.K_consumption_half_saturation;
        const double b = half + dt * demand - diffused;
        const double discriminant = std::max(
            0.0, b * b + 4.0 * half * diffused);
        const double consumed = 0.5 * (-b + std::sqrt(discriminant));
        nutrient_next_[here] = config_.schema_version >= 7
            ? continuum::resource_after_uptake(diffused, consumers, nutrient.K_consumption_rate_per_hour,
                half, dt, source_value)
            : std::clamp(consumed, 0.0, source_value);
    });
    nutrient_.swap(nutrient_next_);
    if (moving_front) {
        // The explicit stencil above must read the previous value of a voxel
        // that has just become exterior. After the swap, synchronize those
        // Dirichlet sites in the inactive buffer so the solve box may shrink
        // again on the next step without leaving stale nutrient behind.
        deterministic_parallel_for(
            update_count, workers, [&](std::size_t offset) {
                int x{}, y{}, z{};
                coordinates(offset, x, y, z);
                const std::size_t here = index(x, y, z);
                if (nutrient_source(x, y, z, here)) {
                    nutrient_next_[here] = source_value;
                }
            });
    }
    ++nutrient_solve_count_;
    validate_resources();
}

void StructuredPdeModel3D::solve_nutrient() {
    const auto& continuum = config_.continuum;
    if (config_.schema_version >= 5) {
        const int nx = continuum.grid.shape[0];
        const int ny = continuum.grid.shape[1];
        for (std::size_t here = 0; here < voxel_count_; ++here) {
            const int x = static_cast<int>(
                here % static_cast<std::size_t>(nx));
            const std::size_t yz = here / static_cast<std::size_t>(nx);
            const int y = static_cast<int>(
                yz % static_cast<std::size_t>(ny));
            const int z = static_cast<int>(
                yz / static_cast<std::size_t>(ny));
            if (nutrient_source(x, y, z, here)) {
                nutrient_[here] = continuum.nutrient.vessel_value;
            }
        }
        ++nutrient_solve_count_;
        validate_resources();
        if (config_.schema_version >= 6) nutrient_next_ = nutrient_;
        return;
    }
    const int nx = continuum.grid.shape[0];
    const int ny = continuum.grid.shape[1];
    const int nz = continuum.grid.shape[2];
    const int dimensions = continuum.base.thin_layer ? 2 : 3;
    const double diffusion = continuum.nutrient.diffusion_voxels2_per_hour /
        (continuum.grid.spacing_voxels * continuum.grid.spacing_voxels);
    const double laplacian_diagonal = 2.0 * dimensions * diffusion;
    const int workers = std::max(
        1, std::min(continuum.base.threads, available_worker_threads()));
    for (int iteration = 0; iteration < continuum.nutrient.solver_iterations;
         ++iteration) {
        deterministic_parallel_for(voxel_count_, workers, [&](std::size_t here) {
            const int x = static_cast<int>(here % static_cast<std::size_t>(nx));
            const std::size_t yz = here / static_cast<std::size_t>(nx);
            const int y = static_cast<int>(yz % static_cast<std::size_t>(ny));
            const int z = static_cast<int>(yz / static_cast<std::size_t>(ny));
            double neighbor_sum = 0.0;
            if (x > 0) neighbor_sum += nutrient_[index(x - 1, y, z)];
            if (x + 1 < nx) neighbor_sum += nutrient_[index(x + 1, y, z)];
            if (y > 0) neighbor_sum += nutrient_[index(x, y - 1, z)];
            if (y + 1 < ny) neighbor_sum += nutrient_[index(x, y + 1, z)];
            if (!continuum.base.thin_layer) {
                if (z > 0) neighbor_sum += nutrient_[index(x, y, z - 1)];
                if (z + 1 < nz) neighbor_sum += nutrient_[index(x, y, z + 1)];
            }
            const double old = nutrient_[here];
            const bool per_cell = continuum.nutrient.consumption_model ==
                "per_cell_ratio_v2";
            const double large_consumption_weight =
                per_cell ? 1.0 : large_cell_volume_;
            const double r_consumers = r_normal_[0][here] +
                r_active(StructuredStage3D::small, here) +
                large_consumption_weight * (r_normal_[1][here] +
                    r_active(StructuredStage3D::large, here));
            const double K_consumers =
                K_[0][here] + large_consumption_weight * K_[1][here];
            const double r_sink =
                continuum.nutrient.r_consumption_rate_per_hour *
                r_consumers /
                (continuum.nutrient.r_consumption_half_saturation + old);
            const double K_sink =
                continuum.nutrient.K_consumption_rate_per_hour *
                K_consumers /
                (continuum.nutrient.K_consumption_half_saturation + old);
            const double exchange = continuum.nutrient.vessel_exchange_per_hour *
                std::clamp(vessel_[here], 0.0, 1.0);
            const double denominator = laplacian_diagonal +
                continuum.nutrient.decay_per_hour + exchange + r_sink + K_sink;
            const double candidate = denominator > 0.0
                ? (diffusion * neighbor_sum +
                   exchange * continuum.nutrient.vessel_value) / denominator
                : 0.0;
            nutrient_next_[here] = std::clamp(
                old + continuum.nutrient.relaxation * (candidate - old),
                0.0, continuum.nutrient.vessel_value);
        });
        nutrient_.swap(nutrient_next_);
    }
    ++nutrient_solve_count_;
    validate_resources();
}

void StructuredPdeModel3D::build_activation_density() {
    const auto& continuum = config_.continuum;
    const int nx = continuum.grid.shape[0];
    const int ny = continuum.grid.shape[1];
    const int nz = continuum.grid.shape[2];
    const int edge = continuum.base.migration_activation_window_edge;
    const int query = continuum.base.migration_activation_block_edge;
    const int lower = (edge - 1) / 2;
    const int upper = edge - lower - 1;
    const int origin_x = static_cast<int>(std::floor(continuum.grid.origin[0]));
    const int origin_y = static_cast<int>(std::floor(continuum.grid.origin[1]));
    const int origin_z = static_cast<int>(std::floor(continuum.grid.origin[2]));

    if (continuum.base.thin_layer) {
        const auto clear_density = [&](const StructuredActiveBounds3D& bounds) {
            if (!bounds.valid) return;
            for (int y = bounds.y0; y < bounds.y1; ++y) {
                const std::size_t begin = index(bounds.x0, y, 0);
                const std::size_t end = index(bounds.x1 - 1, y, 0) + 1;
                for (auto& field : activation_density_) {
                    std::fill(field.begin() + begin, field.begin() + end, 0.0);
                }
            }
        };
        clear_density(activation_density_bounds_);
        activation_density_bounds_ = {};
        if (!population_bounds_.valid) return;

        // The ABM performs one query per 32x32 anchor block. Only blocks that
        // contain population can activate r mass, so use an exact local
        // integral image covering those blocks and their 70x70 query windows.
        const int first_qx = floor_div(
            origin_x + population_bounds_.x0, query);
        const int last_qx = floor_div(
            origin_x + population_bounds_.x1 - 1, query);
        const int first_qy = floor_div(
            origin_y + population_bounds_.y0, query);
        const int last_qy = floor_div(
            origin_y + population_bounds_.y1 - 1, query);
        const int source_x0 = std::clamp(
            first_qx * query + query / 2 - origin_x - lower, 0, nx);
        const int source_y0 = std::clamp(
            first_qy * query + query / 2 - origin_y - lower, 0, ny);
        const int source_x1 = std::clamp(
            last_qx * query + query / 2 - origin_x + upper + 1, 0, nx);
        const int source_y1 = std::clamp(
            last_qy * query + query / 2 - origin_y + upper + 1, 0, ny);
        const int width = source_x1 - source_x0;
        const int height = source_y1 - source_y0;
        const int pitch = width + 1;
        std::vector<double> prefix(
            static_cast<std::size_t>(pitch) * (height + 1), 0.0);
        for (int y = 1; y <= height; ++y) {
            double row = 0.0;
            for (int x = 1; x <= width; ++x) {
                const std::size_t source =
                    index(source_x0 + x - 1, source_y0 + y - 1, 0);
                row += (r_normal_[0][source] + r_normal_[1][source] +
                        r_active(StructuredStage3D::small, source) +
                        r_active(StructuredStage3D::large, source) +
                        K_[0][source] + K_[1][source]) * voxel_measure_;
                prefix[static_cast<std::size_t>(y) * pitch + x] = row;
            }
        }
        for (int x = 1; x <= width; ++x) {
            for (int y = 1; y <= height; ++y) {
                prefix[static_cast<std::size_t>(y) * pitch + x] +=
                    prefix[static_cast<std::size_t>(y - 1) * pitch + x];
            }
        }
        const auto sum = [&](int x0, int y0, int x1, int y1) {
            x0 -= source_x0;
            y0 -= source_y0;
            x1 -= source_x0;
            y1 -= source_y0;
            return prefix[static_cast<std::size_t>(y1) * pitch + x1]
                - prefix[static_cast<std::size_t>(y0) * pitch + x1]
                - prefix[static_cast<std::size_t>(y1) * pitch + x0]
                + prefix[static_cast<std::size_t>(y0) * pitch + x0];
        };
        for (int qy = first_qy; qy <= last_qy; ++qy) {
            for (int qx = first_qx; qx <= last_qx; ++qx) {
                const int center_x = qx * query + query / 2 - origin_x;
                const int center_y = qy * query + query / 2 - origin_y;
                const int x0 = std::clamp(center_x - lower, 0, nx);
                const int y0 = std::clamp(center_y - lower, 0, ny);
                const int x1 = std::clamp(center_x + upper + 1, 0, nx);
                const int y1 = std::clamp(center_y + upper + 1, 0, ny);
                const double count = x1 > x0 && y1 > y0
                    ? sum(x0, y0, x1, y1) : 0.0;
                const double small_density = count / static_cast<double>(edge * edge);
                const double large_density = count /
                    (static_cast<double>(edge * edge) / large_cell_volume_);
                const int bx0 = std::clamp(qx * query - origin_x, 0, nx);
                const int by0 = std::clamp(qy * query - origin_y, 0, ny);
                const int bx1 = std::clamp((qx + 1) * query - origin_x, 0, nx);
                const int by1 = std::clamp((qy + 1) * query - origin_y, 0, ny);
                if (bx1 <= bx0 || by1 <= by0) continue;
                if (!activation_density_bounds_.valid) {
                    activation_density_bounds_ =
                        {bx0, by0, 0, bx1, by1, 1, true};
                } else {
                    activation_density_bounds_.x0 = std::min(
                        activation_density_bounds_.x0, bx0);
                    activation_density_bounds_.y0 = std::min(
                        activation_density_bounds_.y0, by0);
                    activation_density_bounds_.x1 = std::max(
                        activation_density_bounds_.x1, bx1);
                    activation_density_bounds_.y1 = std::max(
                        activation_density_bounds_.y1, by1);
                }
                for (int y = by0; y < by1; ++y) {
                    for (int x = bx0; x < bx1; ++x) {
                        const std::size_t here = index(x, y, 0);
                        activation_density_[0][here] = small_density;
                        activation_density_[1][here] = large_density;
                    }
                }
            }
        }
        return;
    }

    activation_density_bounds_ = {0, 0, 0, nx, ny, nz, true};

    const int px = nx + 1;
    const int py = ny + 1;
    const int pz = nz + 1;
    const auto pindex = [px, py](int x, int y, int z) {
        return (static_cast<std::size_t>(z) * py + y) * px + x;
    };
    std::vector<double> prefix(static_cast<std::size_t>(px) * py * pz, 0.0);
    for (int z = 1; z <= nz; ++z) {
        for (int y = 1; y <= ny; ++y) {
            for (int x = 1; x <= nx; ++x) {
                const std::size_t source = index(x - 1, y - 1, z - 1);
                const double value = (r_normal_[0][source] + r_normal_[1][source] +
                    r_active(StructuredStage3D::small, source) +
                    r_active(StructuredStage3D::large, source) +
                    K_[0][source] + K_[1][source]) * voxel_measure_;
                prefix[pindex(x, y, z)] = value
                    + prefix[pindex(x - 1, y, z)]
                    + prefix[pindex(x, y - 1, z)]
                    + prefix[pindex(x, y, z - 1)]
                    - prefix[pindex(x - 1, y - 1, z)]
                    - prefix[pindex(x - 1, y, z - 1)]
                    - prefix[pindex(x, y - 1, z - 1)]
                    + prefix[pindex(x - 1, y - 1, z - 1)];
            }
        }
    }
    const auto sum = [&](int x0, int y0, int z0, int x1, int y1, int z1) {
        return prefix[pindex(x1, y1, z1)]
            - prefix[pindex(x0, y1, z1)] - prefix[pindex(x1, y0, z1)]
            - prefix[pindex(x1, y1, z0)] + prefix[pindex(x0, y0, z1)]
            + prefix[pindex(x0, y1, z0)] + prefix[pindex(x1, y0, z0)]
            - prefix[pindex(x0, y0, z0)];
    };
    const int first_qx = floor_div(origin_x, query);
    const int last_qx = floor_div(origin_x + nx - 1, query);
    const int first_qy = floor_div(origin_y, query);
    const int last_qy = floor_div(origin_y + ny - 1, query);
    const int first_qz = floor_div(origin_z, query);
    const int last_qz = floor_div(origin_z + nz - 1, query);
    for (int qz = first_qz; qz <= last_qz; ++qz) {
        for (int qy = first_qy; qy <= last_qy; ++qy) {
            for (int qx = first_qx; qx <= last_qx; ++qx) {
                const int cx = qx * query + query / 2 - origin_x;
                const int cy = qy * query + query / 2 - origin_y;
                const int cz = qz * query + query / 2 - origin_z;
                const int x0 = std::clamp(cx - lower, 0, nx);
                const int y0 = std::clamp(cy - lower, 0, ny);
                const int z0 = std::clamp(cz - lower, 0, nz);
                const int x1 = std::clamp(cx + upper + 1, 0, nx);
                const int y1 = std::clamp(cy + upper + 1, 0, ny);
                const int z1 = std::clamp(cz + upper + 1, 0, nz);
                const double count = x1 > x0 && y1 > y0 && z1 > z0
                    ? sum(x0, y0, z0, x1, y1, z1) : 0.0;
                const double capacity = static_cast<double>(edge) * edge * edge;
                const double small_density = count / capacity;
                const double large_density = count / (capacity / large_cell_volume_);
                const int bx0 = std::clamp(qx * query - origin_x, 0, nx);
                const int by0 = std::clamp(qy * query - origin_y, 0, ny);
                const int bz0 = std::clamp(qz * query - origin_z, 0, nz);
                const int bx1 = std::clamp((qx + 1) * query - origin_x, 0, nx);
                const int by1 = std::clamp((qy + 1) * query - origin_y, 0, ny);
                const int bz1 = std::clamp((qz + 1) * query - origin_z, 0, nz);
                for (int z = bz0; z < bz1; ++z) {
                    for (int y = by0; y < by1; ++y) {
                        for (int x = bx0; x < bx1; ++x) {
                            const std::size_t here = index(x, y, z);
                            activation_density_[0][here] = small_density;
                            activation_density_[1][here] = large_density;
                        }
                    }
                }
            }
        }
    }
}

void StructuredPdeModel3D::refresh_activation(double dt) {
    build_activation_density();
    if (!population_bounds_.valid) return;
    const double threshold = config_.continuum.base.migration_activation_threshold;
    const double r_inherent = mean_growth_rate(CellType::r);
    const double full_cycle =
        config_.continuum.base.division_timing.base_cycle_hours /
        std::max(1.0e-12, r_inherent);
    const double mean_fraction =
        config_.continuum.base.migration_activation_duration_mean_fraction;
    const double mean_duration = mean_fraction * full_cycle;
    const auto bounds = population_bounds_;
    for (std::size_t stage = 0; stage < kStructuredStageCount3D; ++stage) {
        for (int z = bounds.z0; z < bounds.z1; ++z) {
            for (int y = bounds.y0; y < bounds.y1; ++y) {
                for (int x = bounds.x0; x < bounds.x1; ++x) {
                    const std::size_t location = index(x, y, z);
                    double mass = r_normal_[stage][location];
                    if (config_.schema_version >= 7) {
                        auto& refractory = r_refractory_[stage][location];
                        auto& clock = refractory_clock_[stage][location];
                        clock = std::max(0.0, clock - refractory * dt);
                        if (clock <= 1.0e-12 * std::max(1.0, refractory) &&
                            activation_density_[stage][location] <=
                                config_.migration.reactivation_density_threshold) {
                            refractory = 0.0;
                            clock = 0.0;
                        }
                        mass = std::max(0.0, mass - refractory);
                    } else if (config_.schema_version >= 5) {
                        auto& cooldown = activation_cooldown_[stage][location];
                        auto& armed = activation_armed_[stage][location];
                        if (active_total_[stage][location] >=
                            config_.migration.minimum_density) {
                            armed = 0U;
                            cooldown = static_cast<float>(
                                config_.migration.reactivation_cooldown_hours);
                        } else {
                            cooldown = static_cast<float>(std::max(
                                0.0, static_cast<double>(cooldown) - dt));
                            if (cooldown <= 0.0F &&
                                activation_density_[stage][location] <=
                                    config_.migration
                                        .reactivation_density_threshold) {
                                armed = 1U;
                            }
                        }
                        if (armed == 0U) continue;
                    }
                    if (mass < config_.migration.minimum_density ||
                        activation_density_[stage][location] < threshold) continue;
                    r_normal_[stage][location] -= mass;
                    active_direction_[stage][0][location] +=
                        static_cast<float>(mass);
                    active_clock_[stage][0][location] +=
                        static_cast<float>(mass * mean_duration);
                    active_total_[stage][location] += static_cast<float>(mass);
                    if (config_.schema_version >= 5) {
                        activation_armed_[stage][location] = 0U;
                        activation_cooldown_[stage][location] =
                            static_cast<float>(
                                config_.migration.reactivation_cooldown_hours);
                    }
                    include_active_location(stage, x, y, z);
                }
            }
        }
    }
}

double StructuredPdeModel3D::mean_growth_rate(CellType type) const noexcept {
    const auto& base = config_.continuum.base;
    if (base.initial_growth_rate_model == "fixed") {
        return type == CellType::r ? base.initial_r_growth_rate
                                   : base.initial_K_growth_rate;
    }
    const auto& distribution = type == CellType::r
        ? base.initial_r_growth_truncated_normal
        : base.initial_K_growth_truncated_normal;
    return std::clamp(distribution.mean, distribution.minimum,
                      distribution.maximum);
}

double StructuredPdeModel3D::normal_diffusion(
    StructuredStage3D stage, CellType type) const noexcept {
    const auto& base = config_.continuum.base;
    const double rate = type == CellType::r
        ? beta_mean(base.normal_r_migration_beta)
        : (base.initial_K_migration_rate_model == "fixed"
               ? base.initial_K_migration_rate
               : beta_mean(base.initial_K_migration_beta));
    double result = (base.thin_layer &&
        config_.continuum.migration.mapping == "shared_fixed_lattice_means_v2"
        ? 3.0 / 8.0 : kFixed26DiffusionFactor) * rate *
        config_.continuum.migration.diffusion_scale;
    if (stage == StructuredStage3D::large) {
        result *= config_.continuum.migration.large_mobility_multiplier;
    }
    return result;
}

double StructuredPdeModel3D::active_rate(StructuredStage3D stage) const noexcept {
    double result = beta_mean(config_.continuum.base.normal_r_migration_beta) *
        config_.continuum.base.activated_r_normal_multiplier;
    if (stage == StructuredStage3D::large) {
        result *= config_.continuum.migration.large_mobility_multiplier;
    }
    return result;
}

void StructuredPdeModel3D::migrate_normal_and_K(double dt) {
    const auto& continuum = config_.continuum;
    if (!population_bounds_.valid) return;
    const int nx = continuum.grid.shape[0];
    const int ny = continuum.grid.shape[1];
    const int nz = continuum.grid.shape[2];
    const double inverse_h2 = 1.0 /
        (continuum.grid.spacing_voxels * continuum.grid.spacing_voxels);
    const int workers = std::max(
        1, std::min(continuum.base.threads, available_worker_threads()));
    const StructuredActiveBounds3D old_bounds = population_bounds_;
    const StructuredActiveBounds3D target_bounds = expanded_bounds(old_bounds);
    const auto clear = [&](const StructuredActiveBounds3D& bounds) {
        if (!bounds.valid) return;
        for (int z = bounds.z0; z < bounds.z1; ++z) {
            for (int y = bounds.y0; y < bounds.y1; ++y) {
                const std::size_t begin = index(bounds.x0, y, z);
                const std::size_t end = index(bounds.x1 - 1, y, z) + 1;
                for (std::size_t stage = 0; stage < 2; ++stage) {
                    std::fill(r_normal_work_[stage].begin() + begin,
                              r_normal_work_[stage].begin() + end, 0.0);
                    std::fill(K_work_[stage].begin() + begin,
                              K_work_[stage].begin() + end, 0.0);
                    if (config_.schema_version >= 7) {
                        std::fill(refractory_work_[stage].begin() + begin,
                                  refractory_work_[stage].begin() + end, 0.0);
                        std::fill(refractory_clock_work_[stage].begin() + begin,
                                  refractory_clock_work_[stage].begin() + end, 0.0);
                    }
                }
            }
        }
    };
    clear(normal_work_dirty_bounds_);
    clear(target_bounds);
    const std::size_t bx = static_cast<std::size_t>(
        target_bounds.x1 - target_bounds.x0);
    const std::size_t by = static_cast<std::size_t>(
        target_bounds.y1 - target_bounds.y0);
    const std::size_t bz = static_cast<std::size_t>(
        target_bounds.z1 - target_bounds.z0);
    const std::size_t box_size = bx * by * bz;

    const auto migrate_field = [&](const std::vector<double>& source,
                                   std::vector<double>& target,
                                   StructuredStage3D stage,
                                   CellType type) {
        const bool large = stage == StructuredStage3D::large;
        const double vacancy_exponent = large ? large_cell_volume_ : 1.0;
        const double diffusion = normal_diffusion(stage, type);
        std::atomic<bool> invalid_negative{false};
        deterministic_parallel_for(box_size, workers, [&](std::size_t offset) {
            const int x = target_bounds.x0 + static_cast<int>(offset % bx);
            const std::size_t yz = offset / bx;
            const int y = target_bounds.y0 + static_cast<int>(yz % by);
            const int z = target_bounds.z0 + static_cast<int>(yz / by);
            const std::size_t here = index(x, y, z);
            const double local_occupied = std::clamp(
                occupied_fraction(here) /
                    continuum.reaction.maximum_occupied_fraction,
                0.0, 1.0);
            const double vacancy = 1.0 - local_occupied;
            const double availability = vessel_blocks_cells(here)
                ? 0.0 : std::pow(vacancy, vacancy_exponent);
            double delta = 0.0;
            const auto exchange = [&](int ox, int oy, int oz) {
                if (ox < 0 || ox >= nx || oy < 0 || oy >= ny ||
                    oz < 0 || oz >= nz) return;
                const std::size_t other = index(ox, oy, oz);
                const double other_occupied = std::clamp(
                    occupied_fraction(other) /
                        continuum.reaction.maximum_occupied_fraction,
                    0.0, 1.0);
                const double other_vacancy = 1.0 - other_occupied;
                const double other_availability = vessel_blocks_cells(other)
                    ? 0.0 : std::pow(other_vacancy, vacancy_exponent);
                const double mobility = std::pow(
                    0.5 * (vacancy + other_vacancy),
                    continuum.migration.crowding_exponent);
                delta += diffusion * inverse_h2 * mobility *
                    (source[other] * availability -
                     source[here] * other_availability);
            };
            exchange(x - 1, y, z);
            exchange(x + 1, y, z);
            exchange(x, y - 1, z);
            exchange(x, y + 1, z);
            if (!continuum.base.thin_layer) {
                exchange(x, y, z - 1);
                exchange(x, y, z + 1);
            }
            const double value = source[here] + dt * delta;
            if (value < -1.0e-10) invalid_negative.store(true);
            target[here] = value >= config_.migration.minimum_density
                ? std::max(0.0, value) : 0.0;
        });
        if (invalid_negative.load()) {
            throw std::runtime_error(
                "structured normal migration produced negative density");
        }
    };

    for (std::size_t stage = 0; stage < kStructuredStageCount3D; ++stage) {
        const auto value = static_cast<StructuredStage3D>(stage);
        migrate_field(r_normal_[stage], r_normal_work_[stage], value, CellType::r);
        migrate_field(K_[stage], K_work_[stage], value, CellType::K);
        if (config_.schema_version >= 7) {
            migrate_field(r_refractory_[stage], refractory_work_[stage], value, CellType::r);
            migrate_field(refractory_clock_[stage], refractory_clock_work_[stage], value, CellType::r);
        }
    }
    r_normal_.swap(r_normal_work_);
    K_.swap(K_work_);
    if (config_.schema_version >= 7) {
        r_refractory_.swap(refractory_work_);
        refractory_clock_.swap(refractory_clock_work_);
        for (std::size_t stage = 0; stage < 2; ++stage) {
            for (int z = target_bounds.z0; z < target_bounds.z1; ++z) {
                for (int y = target_bounds.y0; y < target_bounds.y1; ++y) {
                    for (int x = target_bounds.x0; x < target_bounds.x1; ++x) {
                        const auto here = index(x, y, z);
                        const double mass = r_refractory_[stage][here];
                        if (mass > r_normal_[stage][here]) {
                            refractory_clock_[stage][here] *= r_normal_[stage][here] / mass;
                            r_refractory_[stage][here] = r_normal_[stage][here];
                        }
                        if (r_refractory_[stage][here] == 0.0) refractory_clock_[stage][here] = 0.0;
                    }
                }
            }
        }
    }
    normal_work_dirty_bounds_ = old_bounds;
    population_bounds_ = target_bounds;
}

std::vector<std::size_t> StructuredPdeModel3D::eligible_initial_directions(
    std::size_t location) const {
    const auto& continuum = config_.continuum;
    const int nx = continuum.grid.shape[0];
    const int ny = continuum.grid.shape[1];
    const int nz = continuum.grid.shape[2];
    const int x = static_cast<int>(location % static_cast<std::size_t>(nx));
    const std::size_t yz = location / static_cast<std::size_t>(nx);
    const int y = static_cast<int>(yz % static_cast<std::size_t>(ny));
    const int z = static_cast<int>(yz / static_cast<std::size_t>(ny));
    const int radius = continuum.base.direction_density_radius;
    const double minimum_cosine = std::cos(
        continuum.base.direction_density_half_angle_degrees *
        std::acos(-1.0) / 180.0);
    std::vector<std::size_t> result;
    for (std::size_t direction_index = 0;
         direction_index < direction_ids_.size(); ++direction_index) {
        const Vec3i forward = direction_vector(direction_ids_[direction_index]);
        const int target_x = x + forward.x;
        const int target_y = y + forward.y;
        const int target_z = z + forward.z;
        if (target_x < 0 || target_x >= nx || target_y < 0 || target_y >= ny ||
            target_z < 0 || target_z >= nz ||
            vessel_blocks_cells(index(target_x, target_y, target_z))) {
            continue;
        }
        const double forward_length =
            std::sqrt(static_cast<double>(squared_length(forward)));
        double count = 0.0;
        std::size_t sites = 0;
        for (int dz = continuum.base.thin_layer ? 0 : -radius;
             dz <= (continuum.base.thin_layer ? 0 : radius); ++dz) {
            for (int dy = -radius; dy <= radius; ++dy) {
                for (int dx = -radius; dx <= radius; ++dx) {
                    const int distance = std::max(
                        {std::abs(dx), std::abs(dy), std::abs(dz)});
                    if (distance == 0 || distance > radius) continue;
                    const Vec3i offset{dx, dy, dz};
                    const double offset_length =
                        std::sqrt(static_cast<double>(squared_length(offset)));
                    const double cosine = static_cast<double>(dot(offset, forward)) /
                        (offset_length * forward_length);
                    if (cosine + 1.0e-12 < minimum_cosine) continue;
                    ++sites;
                    const int ox = x + dx;
                    const int oy = y + dy;
                    const int oz = z + dz;
                    if (ox < 0 || ox >= nx || oy < 0 || oy >= ny ||
                        oz < 0 || oz >= nz) continue;
                    const std::size_t other = index(ox, oy, oz);
                    count += (r_normal_[0][other] + r_normal_[1][other] +
                              r_active(StructuredStage3D::small, other) +
                              r_active(StructuredStage3D::large, other) +
                              K_[0][other] + K_[1][other]) * voxel_measure_;
                }
            }
        }
        const double density = sites > 0 ? count / static_cast<double>(sites) : 1.0;
        if (density <= continuum.base.direction_density_threshold) {
            result.push_back(direction_index + 1);
        }
    }
    return result;
}

void StructuredPdeModel3D::build_guidance_prefix(double dt) {
    guidance_prefix_bounds_ = {};
    guidance_prefix_pitch_ = 0;
    if (config_.schema_version < 3 ||
        !config_.continuum.base.thin_layer) return;

    StructuredActiveBounds3D bounds;
    for (const auto& stage_bounds : active_bounds_) {
        if (!stage_bounds.valid) continue;
        if (!bounds.valid) {
            bounds = stage_bounds;
        } else {
            bounds.x0 = std::min(bounds.x0, stage_bounds.x0);
            bounds.y0 = std::min(bounds.y0, stage_bounds.y0);
            bounds.x1 = std::max(bounds.x1, stage_bounds.x1);
            bounds.y1 = std::max(bounds.y1, stage_bounds.y1);
        }
    }
    if (!bounds.valid) return;

    int maximum_substeps = 1;
    for (std::size_t stage = 0; stage < kStructuredStageCount3D; ++stage) {
        const double rate = active_rate(static_cast<StructuredStage3D>(stage));
        if (!(rate > 0.0)) continue;
        maximum_substeps = std::max(maximum_substeps,
            static_cast<int>(std::ceil(
                rate * dt /
                -std::log1p(-config_.migration.maximum_move_probability_per_substep))));
    }
    const int edge = config_.schema_version >= 5
        ? config_.migration.direction_nutrient_window_edge
        : config_.migration.direction_density_window_edge;
    const int lower = (edge - 1) / 2;
    const int upper = edge - lower - 1;
    const int margin = std::max(lower, upper) + maximum_substeps + 1;
    const int nx = config_.continuum.grid.shape[0];
    const int ny = config_.continuum.grid.shape[1];
    bounds.x0 = std::max(0, bounds.x0 - margin);
    bounds.y0 = std::max(0, bounds.y0 - margin);
    bounds.x1 = std::min(nx, bounds.x1 + margin);
    bounds.y1 = std::min(ny, bounds.y1 + margin);
    bounds.z0 = 0;
    bounds.z1 = 1;
    guidance_prefix_bounds_ = bounds;
    guidance_prefix_pitch_ = bounds.x1 - bounds.x0 + 1;
    const int height = bounds.y1 - bounds.y0;
    const std::size_t size = static_cast<std::size_t>(height) *
        static_cast<std::size_t>(guidance_prefix_pitch_);
    guidance_density_row_prefix_.assign(size, 0.0);
    guidance_resource_row_prefix_.assign(size, 0.0);
    const double vessel_value = config_.continuum.nutrient.vessel_value;
    for (int y = bounds.y0; y < bounds.y1; ++y) {
        const std::size_t row = static_cast<std::size_t>(y - bounds.y0) *
            static_cast<std::size_t>(guidance_prefix_pitch_);
        double density_sum = 0.0;
        double resource_sum = 0.0;
        for (int x = bounds.x0; x < bounds.x1; ++x) {
            const std::size_t location = index(x, y, 0);
            if (config_.schema_version < 5) {
                density_sum += (r_normal_[0][location] +
                    r_normal_[1][location] +
                    r_active(StructuredStage3D::small, location) +
                    r_active(StructuredStage3D::large, location) +
                    K_[0][location] + K_[1][location]) * voxel_measure_;
            }
            resource_sum += std::clamp(
                nutrient_[location] / vessel_value, 0.0, 1.0);
            const std::size_t column =
                static_cast<std::size_t>(x - bounds.x0 + 1);
            guidance_density_row_prefix_[row + column] = density_sum;
            guidance_resource_row_prefix_[row + column] = resource_sum;
        }
    }
}

std::vector<double> StructuredPdeModel3D::guided_direction_weights(
    std::size_t location,
    bool allow_crowded_sectors) const {
    const auto& continuum = config_.continuum;
    const auto& base = continuum.base;
    const int nx = continuum.grid.shape[0];
    const int ny = continuum.grid.shape[1];
    const int nz = continuum.grid.shape[2];
    const int x = static_cast<int>(location % static_cast<std::size_t>(nx));
    const std::size_t yz = location / static_cast<std::size_t>(nx);
    const int y = static_cast<int>(yz % static_cast<std::size_t>(ny));
    const int z = static_cast<int>(yz / static_cast<std::size_t>(ny));
    const int edge = config_.schema_version >= 5
        ? config_.migration.direction_nutrient_window_edge
        : config_.schema_version >= 3
        ? config_.migration.direction_density_window_edge
        : 2 * base.direction_density_radius + 1;
    const int lower = (edge - 1) / 2;
    const int upper = edge - lower - 1;
    const double minimum_cosine = std::cos(
        base.direction_density_half_angle_degrees * std::acos(-1.0) / 180.0);
    const double floor = base.direction_minimum_guidance_weight;
    std::vector<double> result(direction_ids_.size() + 1, 0.0);
    for (std::size_t direction_index = 0;
         direction_index < direction_ids_.size(); ++direction_index) {
        const Vec3i forward = direction_vector(direction_ids_[direction_index]);
        const int target_x = x + forward.x;
        const int target_y = y + forward.y;
        const int target_z = z + forward.z;
        if (target_x < 0 || target_x >= nx || target_y < 0 || target_y >= ny ||
            target_z < 0 || target_z >= nz ||
            vessel_blocks_cells(index(target_x, target_y, target_z))) {
            continue;
        }
        const double forward_length =
            std::sqrt(static_cast<double>(squared_length(forward)));
        double count = 0.0;
        double resource = 0.0;
        std::size_t sites = 0;
        std::size_t resource_sites = 0;
        if (config_.schema_version >= 3 && base.thin_layer &&
            guidance_prefix_bounds_.valid) {
            if (config_.schema_version < 7) {
                sites = direction_sector_site_counts_[direction_index];
            }
            for (const auto& span : direction_row_spans_[direction_index]) {
                const int oy = y + span.dy;
                if (oy < guidance_prefix_bounds_.y0 ||
                    oy >= guidance_prefix_bounds_.y1) continue;
                const int ox0 = std::max(
                    guidance_prefix_bounds_.x0, x + span.dx0);
                const int ox1 = std::min(
                    guidance_prefix_bounds_.x1, x + span.dx1 + 1);
                if (ox1 <= ox0) continue;
                const std::size_t row =
                    static_cast<std::size_t>(oy - guidance_prefix_bounds_.y0) *
                    static_cast<std::size_t>(guidance_prefix_pitch_);
                const std::size_t left =
                    static_cast<std::size_t>(ox0 - guidance_prefix_bounds_.x0);
                const std::size_t right =
                    static_cast<std::size_t>(ox1 - guidance_prefix_bounds_.x0);
                resource_sites += static_cast<std::size_t>(ox1 - ox0);
                if (config_.schema_version >= 7) sites += static_cast<std::size_t>(ox1 - ox0);
                count += guidance_density_row_prefix_[row + right] -
                    guidance_density_row_prefix_[row + left];
                resource += guidance_resource_row_prefix_[row + right] -
                    guidance_resource_row_prefix_[row + left];
            }
        } else {
            for (int dz = base.thin_layer ? 0 : -lower;
                 dz <= (base.thin_layer ? 0 : upper); ++dz) {
                for (int dy = -lower; dy <= upper; ++dy) {
                    for (int dx = -lower; dx <= upper; ++dx) {
                        if (dx == 0 && dy == 0 && dz == 0) continue;
                        const Vec3i offset{dx, dy, dz};
                        const double offset_length =
                            std::sqrt(static_cast<double>(squared_length(offset)));
                        const double cosine =
                            static_cast<double>(dot(offset, forward)) /
                            (offset_length * forward_length);
                        if (cosine + 1.0e-12 < minimum_cosine) continue;
                        if (config_.schema_version < 7) ++sites;
                        const int ox = x + dx;
                        const int oy = y + dy;
                        const int oz = z + dz;
                        if (ox < 0 || ox >= nx || oy < 0 || oy >= ny ||
                            oz < 0 || oz >= nz) continue;
                        if (config_.schema_version >= 7) ++sites;
                        const std::size_t other = index(ox, oy, oz);
                        ++resource_sites;
                        count += (r_normal_[0][other] + r_normal_[1][other] +
                                  r_active(StructuredStage3D::small, other) +
                                  r_active(StructuredStage3D::large, other) +
                                  K_[0][other] + K_[1][other]) * voxel_measure_;
                        resource += std::clamp(
                            nutrient_[other] / continuum.nutrient.vessel_value,
                            0.0, 1.0);
                    }
                }
            }
        }
        if (sites == 0) continue;
        if (config_.schema_version >= 5) {
            if (resource_sites == 0) continue;
            const double local_resource = std::clamp(
                nutrient_[location] / continuum.nutrient.vessel_value,
                0.0, 1.0);
            const double directional_resource =
                resource / static_cast<double>(resource_sites);
            double gradient = directional_resource - local_resource;
            if (base.direction_guidance_model == "nutrient_gradient_shared_resource_v3") {
                gradient *= continuum.nutrient.vessel_value;
            }
            if (std::abs(gradient) <
                config_.migration.zero_gradient_tolerance) {
                gradient = 0.0;
            }
            const double exponent = std::clamp(
                config_.migration.chemotaxis_strength * gradient,
                -40.0, 40.0);
            result[direction_index + 1] =
                std::pow(forward_length, -base.distance_weight_exponent) *
                std::exp(exponent);
            continue;
        }
        const double density = count / static_cast<double>(sites);
        if (!allow_crowded_sectors &&
            density > base.direction_density_threshold) continue;
        const double density_fraction = std::clamp(
            density / base.direction_density_threshold, 0.0, 1.0);
        const double density_score =
            floor + (1.0 - floor) * (1.0 - density_fraction);
        const double resource_score = floor + (1.0 - floor) *
            resource / static_cast<double>(sites);
        result[direction_index + 1] =
            std::pow(forward_length, -base.distance_weight_exponent) *
            std::pow(density_score,
                     base.direction_density_guidance_exponent) *
            std::pow(resource_score,
                     base.direction_resource_guidance_exponent);
    }
    return result;
}

void StructuredPdeModel3D::expire_active(std::size_t stage, double dt) {
    const double epsilon = config_.migration.minimum_density;
    const auto bounds = active_bounds_[stage];
    if (!bounds.valid) return;
    const auto expire_location = [&](std::size_t location) {
        double total = 0.0;
        for (std::size_t bucket = 0;
             bucket < active_direction_[stage].size(); ++bucket) {
            auto& mass = active_direction_[stage][bucket][location];
            auto& clock_value = active_clock_[stage][bucket][location];
            const double mass_value = mass;
            const double clock = std::max(
                0.0, static_cast<double>(clock_value));
            if (mass_value <= epsilon ||
                clock <= mass_value * dt +
                    1.0e-6 * std::max(1.0, clock)) {
                if (mass_value > 0.0) {
                    r_normal_[stage][location] += mass_value;
                    if (config_.schema_version >= 7) {
                        r_refractory_[stage][location] += mass_value;
                        refractory_clock_[stage][location] += mass_value *
                            config_.migration.reactivation_cooldown_hours;
                    } else if (config_.schema_version >= 5) {
                        activation_armed_[stage][location] = 0U;
                        activation_cooldown_[stage][location] =
                            static_cast<float>(std::max(
                                static_cast<double>(
                                    activation_cooldown_[stage][location]),
                                config_.migration
                                    .reactivation_cooldown_hours));
                    }
                }
                mass = 0.0F;
                clock_value = 0.0F;
            } else {
                clock_value = static_cast<float>(
                    clock - mass_value * dt);
                total += mass_value;
            }
        }
        active_total_[stage][location] = static_cast<float>(total);
    };

    if (config_.schema_version >= 4) {
        const std::size_t bx = static_cast<std::size_t>(
            bounds.x1 - bounds.x0);
        const std::size_t by = static_cast<std::size_t>(
            bounds.y1 - bounds.y0);
        const std::size_t bz = static_cast<std::size_t>(
            bounds.z1 - bounds.z0);
        const std::size_t box_size = bx * by * bz;
        const int workers = std::max(1, std::min(
            config_.continuum.base.threads, available_worker_threads()));
        deterministic_parallel_for(
            box_size, workers, [&](std::size_t offset) {
                const int x = bounds.x0 + static_cast<int>(offset % bx);
                const std::size_t yz = offset / bx;
                const int y = bounds.y0 + static_cast<int>(yz % by);
                const int z = bounds.z0 + static_cast<int>(yz / by);
                expire_location(index(x, y, z));
            });
        return;
    }

    for (int z = bounds.z0; z < bounds.z1; ++z) {
        for (int y = bounds.y0; y < bounds.y1; ++y) {
            for (int x = bounds.x0; x < bounds.x1; ++x) {
                expire_location(index(x, y, z));
            }
        }
    }
}

void StructuredPdeModel3D::migrate_active(double dt) {
    const auto& continuum = config_.continuum;
    const int nx = continuum.grid.shape[0];
    const int ny = continuum.grid.shape[1];
    const int nz = continuum.grid.shape[2];
    build_guidance_prefix(dt);
    if (config_.schema_version >= 4) {
        ++guidance_weight_cache_generation_;
        if (guidance_weight_cache_generation_ == 0U) {
            std::fill(guidance_weight_cache_stamp_.begin(),
                      guidance_weight_cache_stamp_.end(), 0U);
            guidance_weight_cache_generation_ = 1U;
        }
    }
    const auto ensure_guidance_cache = [&](std::size_t location) {
        if (config_.schema_version < 4 ||
            guidance_weight_cache_stamp_[location] ==
                guidance_weight_cache_generation_) {
            return;
        }
        const std::vector<double> weights = guided_direction_weights(location);
        for (std::size_t bucket = 0; bucket < weights.size(); ++bucket) {
            guidance_weight_cache_[bucket][location] =
                static_cast<float>(weights[bucket]);
        }
        guidance_weight_cache_stamp_[location] =
            guidance_weight_cache_generation_;
    };

    for (std::size_t stage = 0; stage < kStructuredStageCount3D; ++stage) {
        if (!active_bounds_[stage].valid) continue;
        const auto stage_value = static_cast<StructuredStage3D>(stage);
        const double rate = active_rate(stage_value);
        if (!(rate > 0.0)) continue;
        const bool resource_guided = config_.schema_version >= 5 ||
            (config_.schema_version >= 2 &&
             continuum.base.direction_guidance_model == "low_density_high_resource_v1");
        // Resolve the no-persistent-direction state before transport. V2/v3
        // rule uses the same low-density/high-resource directional cone score
        // as the nutrient-coupled ABM.
        const auto initial_bounds = active_bounds_[stage];
        for (int z = initial_bounds.z0; z < initial_bounds.z1; ++z) {
            for (int y = initial_bounds.y0; y < initial_bounds.y1; ++y) {
                for (int x = initial_bounds.x0; x < initial_bounds.x1; ++x) {
                    const std::size_t location = index(x, y, z);
                    const double mass = active_direction_[stage][0][location];
                    if (mass < config_.migration.minimum_density) continue;
                    std::vector<std::size_t> targets;
                    std::vector<double> weights(
                        active_direction_[stage].size(), 0.0);
                    double weight_sum = 0.0;
                    if (resource_guided) {
                        if (config_.schema_version >= 4) {
                            ensure_guidance_cache(location);
                            for (std::size_t target = 0;
                                 target < weights.size(); ++target) {
                                weights[target] =
                                    guidance_weight_cache_[target][location];
                            }
                        } else {
                            weights = guided_direction_weights(location);
                        }
                        for (std::size_t target = 1;
                             target < weights.size(); ++target) {
                            if (weights[target] > 0.0) {
                                targets.push_back(target);
                                weight_sum += weights[target];
                            }
                        }
                    } else {
                        targets = eligible_initial_directions(location);
                        for (const std::size_t target : targets) {
                            const double length = std::sqrt(static_cast<double>(
                                squared_length(direction_vector(
                                    direction_ids_[target - 1]))));
                            weights[target] = std::pow(
                                length,
                                -continuum.base.distance_weight_exponent);
                            weight_sum += weights[target];
                        }
                    }
                    if (targets.empty() || !(weight_sum > 0.0)) continue;
                    const double clock = active_clock_[stage][0][location];
                    for (const std::size_t target : targets) {
                        const double fraction = weights[target] / weight_sum;
                        active_direction_[stage][target][location] +=
                            static_cast<float>(mass * fraction);
                        active_clock_[stage][target][location] +=
                            static_cast<float>(clock * fraction);
                    }
                    active_direction_[stage][0][location] = 0.0F;
                    active_clock_[stage][0][location] = 0.0F;
                }
            }
        }
        const double maximum_probability =
            config_.migration.maximum_move_probability_per_substep;
        const int substeps = std::max(1, static_cast<int>(std::ceil(
            rate * dt / -std::log1p(-maximum_probability))));
        const double sub_dt = dt / substeps;
        const double move_probability = 1.0 - std::exp(-rate * sub_dt);
        const double vacancy_exponent = stage == 1 ? large_cell_volume_ : 1.0;
        // Active-r PDE tails occupy most sites inside their bounds, so a
        // contiguous dense traversal is faster than maintaining and sorting a
        // sparse index. Keep the sparse implementation available for future
        // genuinely sparse closures, but use the deterministic dense path for
        // this v4 model.
        constexpr bool sparse_transport = false;

        std::vector<std::size_t> sparse_locations;
        if (sparse_transport) {
            for (int z = initial_bounds.z0; z < initial_bounds.z1; ++z) {
                for (int y = initial_bounds.y0; y < initial_bounds.y1; ++y) {
                    for (int x = initial_bounds.x0; x < initial_bounds.x1; ++x) {
                        const std::size_t location = index(x, y, z);
                        if (active_total_[stage][location] >=
                            config_.migration.minimum_density) {
                            sparse_locations.push_back(location);
                        }
                    }
                }
            }
        }

        for (int substep = 0; substep < substeps; ++substep) {
            const StructuredActiveBounds3D old_bounds = active_bounds_[stage];
            const StructuredActiveBounds3D target_bounds =
                expanded_bounds(old_bounds);
            StructuredActiveBounds3D clear_bounds = target_bounds;
            if (work_dirty_bounds_.valid) {
                clear_bounds.x0 = std::min(clear_bounds.x0, work_dirty_bounds_.x0);
                clear_bounds.y0 = std::min(clear_bounds.y0, work_dirty_bounds_.y0);
                clear_bounds.z0 = std::min(clear_bounds.z0, work_dirty_bounds_.z0);
                clear_bounds.x1 = std::max(clear_bounds.x1, work_dirty_bounds_.x1);
                clear_bounds.y1 = std::max(clear_bounds.y1, work_dirty_bounds_.y1);
                clear_bounds.z1 = std::max(clear_bounds.z1, work_dirty_bounds_.z1);
            }
            // V4 clears the complete prior dirty region once. After each
            // sparse swap below, only the actual source locations in the old
            // buffer need clearing. Legacy schemas retain their dense path.
            if (!sparse_transport || substep == 0) {
                clear_active_work(clear_bounds);
            }

            std::vector<std::size_t> sparse_targets;
            if (sparse_transport) {
                ++active_location_generation_;
                if (active_location_generation_ == 0U) {
                    std::fill(active_location_stamp_.begin(),
                              active_location_stamp_.end(), 0U);
                    active_location_generation_ = 1U;
                }
                sparse_targets.reserve(sparse_locations.size() * 2);
            }
            const auto mark_sparse_target = [&](std::size_t location) {
                if (!sparse_transport ||
                    active_location_stamp_[location] ==
                        active_location_generation_) {
                    return;
                }
                active_location_stamp_[location] =
                    active_location_generation_;
                sparse_targets.push_back(location);
            };
            const auto transport_location = [&](int x, int y, int z,
                                                std::size_t location) {
                if (active_total_[stage][location] <
                    config_.migration.minimum_density) {
                    return;
                }
                if (resource_guided && config_.schema_version >= 4) {
                    ensure_guidance_cache(location);
                }
                const std::vector<double> guidance =
                    resource_guided && config_.schema_version < 4
                    ? guided_direction_weights(location)
                    : std::vector<double>{};
                const auto guidance_value = [&](std::size_t bucket) {
                    return config_.schema_version >= 4
                        ? static_cast<double>(
                              guidance_weight_cache_[bucket][location])
                        : guidance[bucket];
                };
                for (std::size_t bucket = 0;
                     bucket < active_direction_[stage].size(); ++bucket) {
                    const double mass =
                        active_direction_[stage][bucket][location];
                    if (!(mass > 0.0)) continue;
                    const double clock = std::max(0.0, static_cast<double>(
                        active_clock_[stage][bucket][location]));
                    if (mass < config_.migration.minimum_density) {
                        active_work_[bucket][location] +=
                            static_cast<float>(mass);
                        clock_work_[bucket][location] +=
                            static_cast<float>(clock);
                        mark_sparse_target(location);
                        continue;
                    }
                    const double clock_per_mass = clock / mass;
                    double accepted_total = 0.0;
                    const auto available = [&](std::size_t target_bucket) {
                        const Vec3i step = direction_vector(
                            direction_ids_[target_bucket - 1]);
                        if (!(x + step.x >= 0 && x + step.x < nx &&
                            y + step.y >= 0 && y + step.y < ny &&
                            z + step.z >= 0 && z + step.z < nz)) {
                            return false;
                        }
                        return !vessel_blocks_cells(index(
                            x + step.x, y + step.y, z + step.z));
                    };
                    const auto move = [&](std::size_t target_bucket,
                                          double direction_probability) {
                        const Vec3i step = direction_vector(
                            direction_ids_[target_bucket - 1]);
                        const int ox = x + step.x;
                        const int oy = y + step.y;
                        const int oz = z + step.z;
                        const std::size_t other = index(ox, oy, oz);
                        if (vessel_blocks_cells(other)) return;
                        const double occupied = std::clamp(
                            occupied_fraction(other) /
                                continuum.reaction.maximum_occupied_fraction,
                            0.0, 1.0);
                        const double availability = std::pow(
                            std::max(0.0, 1.0 - occupied),
                            vacancy_exponent);
                        const double moved = mass * move_probability *
                            direction_probability * availability;
                        if (!(moved > 0.0)) return;
                        accepted_total += moved;
                        active_work_[target_bucket][other] +=
                            static_cast<float>(moved);
                        clock_work_[target_bucket][other] +=
                            static_cast<float>(moved * clock_per_mass);
                        mark_sparse_target(other);
                    };
                    if (bucket != 0) {
                        if (resource_guided) {
                            const bool forward_available = available(bucket) &&
                                guidance_value(bucket) > 0.0;
                            std::array<std::size_t, 26> valid_turns{};
                            std::size_t valid_turn_count = 0;
                            for (const std::size_t turn :
                                 turn_buckets_[bucket]) {
                                if (available(turn) &&
                                    guidance_value(turn) > 0.0) {
                                    valid_turns[valid_turn_count++] = turn;
                                }
                            }
                            double total_weight = 0.0;
                            double forward_weight = 0.0;
                            if (forward_available) {
                                forward_weight = valid_turn_count == 0
                                    ? guidance_value(bucket)
                                    : continuum.base.continue_probability *
                                        guidance_value(bucket);
                                total_weight += forward_weight;
                            }
                            const double turn_prior = valid_turn_count == 0
                                ? 0.0
                                : (forward_available
                                       ? 1.0 - continuum.base.continue_probability
                                       : 1.0) /
                                    static_cast<double>(valid_turn_count);
                            for (std::size_t turn_index = 0;
                                 turn_index < valid_turn_count; ++turn_index) {
                                const std::size_t turn = valid_turns[turn_index];
                                total_weight +=
                                    turn_prior * guidance_value(turn);
                            }
                            if (total_weight > 0.0) {
                                if (forward_weight > 0.0) {
                                    move(bucket, forward_weight / total_weight);
                                }
                                for (std::size_t turn_index = 0;
                                     turn_index < valid_turn_count; ++turn_index) {
                                    const std::size_t turn =
                                        valid_turns[turn_index];
                                    move(turn,
                                         turn_prior * guidance_value(turn) /
                                             total_weight);
                                }
                            }
                        } else {
                            const bool forward_available = available(bucket);
                            std::size_t valid_turns = 0;
                            for (const std::size_t turn :
                                 turn_buckets_[bucket]) {
                                if (available(turn)) ++valid_turns;
                            }
                            if (forward_available) {
                                move(bucket, valid_turns == 0 ? 1.0
                                    : continuum.base.continue_probability);
                            }
                            if (valid_turns > 0) {
                                const double total_turn_probability =
                                    forward_available
                                    ? 1.0 - continuum.base.continue_probability
                                    : 1.0;
                                const double turn_probability =
                                    total_turn_probability / valid_turns;
                                for (const std::size_t turn :
                                     turn_buckets_[bucket]) {
                                    if (available(turn)) {
                                        move(turn, turn_probability);
                                    }
                                }
                            }
                        }
                    }
                    const double retained =
                        std::max(0.0, mass - accepted_total);
                    active_work_[bucket][location] +=
                        static_cast<float>(retained);
                    clock_work_[bucket][location] +=
                        static_cast<float>(retained * clock_per_mass);
                    if (retained > 0.0) mark_sparse_target(location);
                }
            };

            if (sparse_transport) {
                for (const std::size_t location : sparse_locations) {
                    const int x = static_cast<int>(
                        location % static_cast<std::size_t>(nx));
                    const std::size_t yz =
                        location / static_cast<std::size_t>(nx);
                    const int y = static_cast<int>(
                        yz % static_cast<std::size_t>(ny));
                    const int z = static_cast<int>(
                        yz / static_cast<std::size_t>(ny));
                    transport_location(x, y, z, location);
                }
            } else if (config_.schema_version >= 4 &&
                       continuum.base.thin_layer) {
                // Sources on rows separated by three voxels have disjoint
                // Moore-neighbourhood write regions. Process each of the
                // three row colours in parallel without atomics; x remains
                // ordered within a row, preserving the local accumulation
                // order. This is algebraically the same transport operator.
                const int workers = std::max(1, std::min(
                    continuum.base.threads, available_worker_threads()));
                for (int row_colour = 0; row_colour < 3; ++row_colour) {
                    const int first_y = old_bounds.y0 + row_colour;
                    if (first_y >= old_bounds.y1) continue;
                    const std::size_t row_count = static_cast<std::size_t>(
                        (old_bounds.y1 - first_y + 2) / 3);
                    const auto transport_row = [&](std::size_t row_offset) {
                        const int y = first_y +
                            3 * static_cast<int>(row_offset);
                        const int z = old_bounds.z0;
                        for (int x = old_bounds.x0;
                             x < old_bounds.x1; ++x) {
                            transport_location(x, y, z, index(x, y, z));
                        }
                    };
                    // Interleave rows across workers. Central rows contain
                    // more active mass and therefore more directional work;
                    // contiguous static chunks leave several workers idle at
                    // the barrier even though the row write sets are disjoint.
#ifdef _OPENMP
#pragma omp parallel for schedule(static, 1) num_threads(workers)
                    for (std::int64_t row_offset = 0;
                         row_offset < static_cast<std::int64_t>(row_count);
                         ++row_offset) {
                        transport_row(static_cast<std::size_t>(row_offset));
                    }
#else
                    for (std::size_t row_offset = 0;
                         row_offset < row_count; ++row_offset) {
                        transport_row(row_offset);
                    }
#endif
                }
            } else {
                for (int z = old_bounds.z0; z < old_bounds.z1; ++z) {
                    for (int y = old_bounds.y0; y < old_bounds.y1; ++y) {
                        for (int x = old_bounds.x0; x < old_bounds.x1; ++x) {
                            transport_location(x, y, z, index(x, y, z));
                        }
                    }
                }
            }
            active_direction_[stage].swap(active_work_);
            active_clock_[stage].swap(clock_work_);
            if (!population_bounds_.valid) {
                population_bounds_ = target_bounds;
            } else {
                population_bounds_.x0 = std::min(
                    population_bounds_.x0, target_bounds.x0);
                population_bounds_.y0 = std::min(
                    population_bounds_.y0, target_bounds.y0);
                population_bounds_.z0 = std::min(
                    population_bounds_.z0, target_bounds.z0);
                population_bounds_.x1 = std::max(
                    population_bounds_.x1, target_bounds.x1);
                population_bounds_.y1 = std::max(
                    population_bounds_.y1, target_bounds.y1);
                population_bounds_.z1 = std::max(
                    population_bounds_.z1, target_bounds.z1);
            }
            if (sparse_transport) {
                // After the swap, active_work_/clock_work_ contain the old
                // input. Clear exactly those source sites so the work fields
                // are globally clean for the next microstep.
                for (const std::size_t location : sparse_locations) {
                    active_total_[stage][location] = 0.0F;
                    for (std::size_t bucket = 0;
                         bucket < active_work_.size(); ++bucket) {
                        active_work_[bucket][location] = 0.0F;
                        clock_work_[bucket][location] = 0.0F;
                    }
                }
                work_dirty_bounds_ = {};
                std::sort(sparse_targets.begin(), sparse_targets.end());

                std::vector<std::size_t> remaining_locations;
                remaining_locations.reserve(sparse_targets.size());
                StructuredActiveBounds3D next_bounds;
                const double epsilon = config_.migration.minimum_density;
                for (const std::size_t location : sparse_targets) {
                    double total = 0.0;
                    for (std::size_t bucket = 0;
                         bucket < active_direction_[stage].size(); ++bucket) {
                        auto& mass =
                            active_direction_[stage][bucket][location];
                        auto& clock_value =
                            active_clock_[stage][bucket][location];
                        const double mass_value = mass;
                        const double clock = std::max(
                            0.0, static_cast<double>(clock_value));
                        if (mass_value <= epsilon ||
                            clock <= mass_value * sub_dt +
                                1.0e-6 * std::max(1.0, clock)) {
                            if (mass_value > 0.0) {
                                r_normal_[stage][location] += mass_value;
                            }
                            mass = 0.0F;
                            clock_value = 0.0F;
                        } else {
                            clock_value = static_cast<float>(
                                clock - mass_value * sub_dt);
                            total += mass_value;
                        }
                    }
                    active_total_[stage][location] =
                        static_cast<float>(total);
                    if (total < epsilon) continue;
                    remaining_locations.push_back(location);
                    const int x = static_cast<int>(
                        location % static_cast<std::size_t>(nx));
                    const std::size_t yz =
                        location / static_cast<std::size_t>(nx);
                    const int y = static_cast<int>(
                        yz % static_cast<std::size_t>(ny));
                    const int z = static_cast<int>(
                        yz / static_cast<std::size_t>(ny));
                    if (!next_bounds.valid) {
                        next_bounds = {x, y, z, x + 1, y + 1, z + 1, true};
                    } else {
                        next_bounds.x0 = std::min(next_bounds.x0, x);
                        next_bounds.y0 = std::min(next_bounds.y0, y);
                        next_bounds.z0 = std::min(next_bounds.z0, z);
                        next_bounds.x1 = std::max(next_bounds.x1, x + 1);
                        next_bounds.y1 = std::max(next_bounds.y1, y + 1);
                        next_bounds.z1 = std::max(next_bounds.z1, z + 1);
                    }
                }
                sparse_locations.swap(remaining_locations);
                active_bounds_[stage] = next_bounds;
            } else {
                work_dirty_bounds_ = old_bounds;
                active_bounds_[stage] = target_bounds;
                expire_active(stage, sub_dt);
                shrink_active_bounds(stage);
            }
            if (!active_bounds_[stage].valid) break;
        }
    }
}

void StructuredPdeModel3D::exchange_active_r_with_K(double dt) {
    if (config_.schema_version < 3 ||
        config_.migration.crowding_exchange !=
            "active_r_K_stage1_conservative_v1" ||
        !active_bounds_[0].valid) {
        return;
    }

    // This is the structured-PDE closure of the ABM singleton swap. Only
    // small active-r and small K participate. A swap moves equal mass in
    // opposite directions, so both global species mass and the occupied
    // fraction at each endpoint are conserved by this operator.
    struct SwapProposal {
        std::size_t source{};
        std::size_t target{};
        std::size_t direction_bucket{};
        double requested{};
        double clock_per_mass{};
    };

    const auto& continuum = config_.continuum;
    const int nx = continuum.grid.shape[0];
    const int ny = continuum.grid.shape[1];
    const int nz = continuum.grid.shape[2];
    const auto bounds = active_bounds_[0];
    const double rate = active_rate(StructuredStage3D::small);
    const double move_probability = 1.0 - std::exp(-rate * dt);
    const double maximum = continuum.reaction.maximum_occupied_fraction;
    const double epsilon = config_.migration.minimum_density;
    std::vector<SwapProposal> proposals;
    std::unordered_map<std::size_t, double> demand_by_target;

    for (int z = bounds.z0; z < bounds.z1; ++z) {
        for (int y = bounds.y0; y < bounds.y1; ++y) {
            for (int x = bounds.x0; x < bounds.x1; ++x) {
                const std::size_t source = index(x, y, z);
                const double active = active_total_[0][source];
                if (active < epsilon || vessel_blocks_cells(source)) continue;

                // The ABM only enters its swap path when no empty direction is
                // available. In a density closure, the minimum neighbour
                // occupancy is a continuous approximation of that condition:
                // one empty neighbour makes blockedness zero; a fully packed
                // neighbourhood makes it one.
                double maximum_vacancy = 0.0;
                bool has_cell_site_neighbor = false;
                for (const DirectionId direction : direction_ids_) {
                    const Vec3i step = direction_vector(direction);
                    const int tx = x + step.x;
                    const int ty = y + step.y;
                    const int tz = z + step.z;
                    if (tx < 0 || tx >= nx || ty < 0 || ty >= ny ||
                        tz < 0 || tz >= nz) continue;
                    const std::size_t target = index(tx, ty, tz);
                    if (vessel_blocks_cells(target)) continue;
                    has_cell_site_neighbor = true;
                    const double occupied = std::clamp(
                        occupied_fraction(target) / maximum, 0.0, 1.0);
                    maximum_vacancy = std::max(maximum_vacancy, 1.0 - occupied);
                }
                if (!has_cell_site_neighbor) continue;
                const double blockedness = 1.0 - maximum_vacancy;
                if (blockedness <= epsilon) continue;

                const std::vector<double> weights =
                    guided_direction_weights(source, true);
                double weight_sum = 0.0;
                for (std::size_t bucket = 1; bucket < weights.size(); ++bucket) {
                    if (!(weights[bucket] > 0.0)) continue;
                    const Vec3i step = direction_vector(direction_ids_[bucket - 1]);
                    const std::size_t target = index(
                        x + step.x, y + step.y, z + step.z);
                    if (K_[0][target] >= epsilon) weight_sum += weights[bucket];
                }
                if (!(weight_sum > 0.0)) continue;

                double clock = 0.0;
                for (const auto& field : active_clock_[0]) clock += field[source];
                const double clock_per_mass = clock / active;
                for (std::size_t bucket = 1; bucket < weights.size(); ++bucket) {
                    if (!(weights[bucket] > 0.0)) continue;
                    const Vec3i step = direction_vector(direction_ids_[bucket - 1]);
                    const std::size_t target = index(
                        x + step.x, y + step.y, z + step.z);
                    const double partner = std::clamp(K_[0][target], 0.0, 1.0);
                    if (partner < epsilon) continue;
                    const double requested = active * move_probability *
                        blockedness * weights[bucket] / weight_sum * partner;
                    if (requested < epsilon) continue;
                    proposals.push_back(
                        {source, target, bucket, requested, clock_per_mass});
                    demand_by_target[target] += requested;
                }
            }
        }
    }
    if (proposals.empty()) return;

    std::unordered_map<std::size_t, double> outgoing_by_source;
    std::unordered_map<std::size_t, double> K_delta;
    StructuredActiveBounds3D clear_bounds = expanded_bounds(bounds);
    if (work_dirty_bounds_.valid) {
        clear_bounds.x0 = std::min(clear_bounds.x0, work_dirty_bounds_.x0);
        clear_bounds.y0 = std::min(clear_bounds.y0, work_dirty_bounds_.y0);
        clear_bounds.z0 = std::min(clear_bounds.z0, work_dirty_bounds_.z0);
        clear_bounds.x1 = std::max(clear_bounds.x1, work_dirty_bounds_.x1);
        clear_bounds.y1 = std::max(clear_bounds.y1, work_dirty_bounds_.y1);
        clear_bounds.z1 = std::max(clear_bounds.z1, work_dirty_bounds_.z1);
    }
    clear_active_work(clear_bounds);
    for (const SwapProposal& proposal : proposals) {
        const double demand = demand_by_target.at(proposal.target);
        const double target_scale = demand > K_[0][proposal.target]
            ? K_[0][proposal.target] / demand : 1.0;
        const double accepted = proposal.requested * target_scale;
        if (!(accepted > 0.0)) continue;
        outgoing_by_source[proposal.source] += accepted;
        K_delta[proposal.source] += accepted;
        K_delta[proposal.target] -= accepted;
        active_work_[proposal.direction_bucket][proposal.target] +=
            static_cast<float>(accepted);
        clock_work_[proposal.direction_bucket][proposal.target] +=
            static_cast<float>(accepted * proposal.clock_per_mass);
    }

    const auto target_bounds = expanded_bounds(bounds);
    for (int z = target_bounds.z0; z < target_bounds.z1; ++z) {
        for (int y = target_bounds.y0; y < target_bounds.y1; ++y) {
            for (int x = target_bounds.x0; x < target_bounds.x1; ++x) {
                const std::size_t location = index(x, y, z);
                const auto found = outgoing_by_source.find(location);
                const double active = active_total_[0][location];
                const double fraction = found == outgoing_by_source.end() ||
                    !(active > 0.0) ? 0.0 :
                    std::clamp(found->second / active, 0.0, 1.0);
                double total = 0.0;
                for (std::size_t bucket = 0;
                     bucket < active_direction_[0].size(); ++bucket) {
                    const double mass =
                        active_direction_[0][bucket][location] * (1.0 - fraction) +
                        active_work_[bucket][location];
                    const double clock =
                        active_clock_[0][bucket][location] * (1.0 - fraction) +
                        clock_work_[bucket][location];
                    active_direction_[0][bucket][location] =
                        static_cast<float>(std::max(0.0, mass));
                    active_clock_[0][bucket][location] =
                        static_cast<float>(std::max(0.0, clock));
                    total += std::max(0.0, mass);
                }
                active_total_[0][location] = static_cast<float>(total);
            }
        }
    }
    for (const auto& [location, delta] : K_delta) {
        K_[0][location] = std::max(0.0, K_[0][location] + delta);
    }
    work_dirty_bounds_ = clear_bounds;
    active_bounds_[0] = target_bounds;
    shrink_active_bounds(0);
}

void StructuredPdeModel3D::build_local_counts(
    std::vector<double>& r_counts,
    std::vector<double>& K_counts) const {
    const auto& continuum = config_.continuum;
    if (!population_bounds_.valid) {
        r_counts.clear();
        K_counts.clear();
        return;
    }
    const auto bounds = population_bounds_;
    const int width = bounds.x1 - bounds.x0;
    const int height = bounds.y1 - bounds.y0;
    const int depth = bounds.z1 - bounds.z0;
    const std::size_t box_size = static_cast<std::size_t>(width) * height * depth;
    const int radius = std::max(0, static_cast<int>(std::floor(
        0.5 * continuum.base.growth_density_window_edge /
        continuum.grid.spacing_voxels)));
    const int lower = config_.schema_version >= 7
        ? static_cast<int>(std::floor(((continuum.base.growth_density_window_edge - 1) / 2) / continuum.grid.spacing_voxels))
        : radius;
    const int upper = config_.schema_version >= 7
        ? static_cast<int>(std::floor((continuum.base.growth_density_window_edge -
            (continuum.base.growth_density_window_edge - 1) / 2 - 1) / continuum.grid.spacing_voxels))
        : radius;
    if (continuum.base.thin_layer) {
        const int pitch = width + 1;
        const std::size_t prefix_size =
            static_cast<std::size_t>(pitch) * (height + 1);
        std::vector<double> r_prefix(prefix_size, 0.0);
        std::vector<double> K_prefix(prefix_size, 0.0);
        for (int y = 1; y <= height; ++y) {
            double r_row = 0.0;
            double K_row = 0.0;
            for (int x = 1; x <= width; ++x) {
                const std::size_t source =
                    index(bounds.x0 + x - 1, bounds.y0 + y - 1, bounds.z0);
                r_row += (r_normal_[0][source] + r_normal_[1][source] +
                    r_active(StructuredStage3D::small, source) +
                    r_active(StructuredStage3D::large, source)) * voxel_measure_;
                K_row += (K_[0][source] + K_[1][source]) * voxel_measure_;
                const std::size_t here = static_cast<std::size_t>(y) * pitch + x;
                r_prefix[here] = r_row;
                K_prefix[here] = K_row;
            }
        }
        for (int x = 1; x <= width; ++x) {
            for (int y = 1; y <= height; ++y) {
                const std::size_t here = static_cast<std::size_t>(y) * pitch + x;
                r_prefix[here] += r_prefix[here - pitch];
                K_prefix[here] += K_prefix[here - pitch];
            }
        }
        const auto sum = [pitch](const std::vector<double>& prefix,
                                 int x0, int y0, int x1, int y1) {
            return prefix[static_cast<std::size_t>(y1) * pitch + x1]
                - prefix[static_cast<std::size_t>(y0) * pitch + x1]
                - prefix[static_cast<std::size_t>(y1) * pitch + x0]
                + prefix[static_cast<std::size_t>(y0) * pitch + x0];
        };
        r_counts.assign(box_size, 0.0);
        K_counts.assign(box_size, 0.0);
        const int workers = std::max(
            1, std::min(continuum.base.threads, available_worker_threads()));
        deterministic_parallel_for(box_size, workers, [&](std::size_t offset) {
            const int local_x = static_cast<int>(offset %
                static_cast<std::size_t>(width));
            const int local_y = static_cast<int>(offset /
                static_cast<std::size_t>(width));
            const std::size_t here = index(
                bounds.x0 + local_x, bounds.y0 + local_y, bounds.z0);
            if (config_.schema_version < 7 && occupied_fraction(here) <= 0.0) return;
            const int x0 = std::max(0, local_x - lower);
            const int y0 = std::max(0, local_y - lower);
            const int x1 = std::min(width, local_x + upper + 1);
            const int y1 = std::min(height, local_y + upper + 1);
            r_counts[offset] = sum(r_prefix, x0, y0, x1, y1);
            K_counts[offset] = sum(K_prefix, x0, y0, x1, y1);
        });
        return;
    }

    const int px = width + 1;
    const int py = height + 1;
    const int pz = depth + 1;
    const auto pindex = [px, py](int x, int y, int z) {
        return (static_cast<std::size_t>(z) * py + y) * px + x;
    };
    std::vector<double> r_prefix(static_cast<std::size_t>(px) * py * pz, 0.0);
    std::vector<double> K_prefix(static_cast<std::size_t>(px) * py * pz, 0.0);
    for (int z = 1; z <= depth; ++z) {
        for (int y = 1; y <= height; ++y) {
            for (int x = 1; x <= width; ++x) {
                const std::size_t source = index(
                    bounds.x0 + x - 1, bounds.y0 + y - 1,
                    bounds.z0 + z - 1);
                const double r_value = (r_normal_[0][source] + r_normal_[1][source] +
                    r_active(StructuredStage3D::small, source) +
                    r_active(StructuredStage3D::large, source)) * voxel_measure_;
                const double K_value = (K_[0][source] + K_[1][source]) * voxel_measure_;
                const auto update = [&](std::vector<double>& prefix, double value) {
                    prefix[pindex(x, y, z)] = value
                        + prefix[pindex(x - 1, y, z)]
                        + prefix[pindex(x, y - 1, z)]
                        + prefix[pindex(x, y, z - 1)]
                        - prefix[pindex(x - 1, y - 1, z)]
                        - prefix[pindex(x - 1, y, z - 1)]
                        - prefix[pindex(x, y - 1, z - 1)]
                        + prefix[pindex(x - 1, y - 1, z - 1)];
                };
                update(r_prefix, r_value);
                update(K_prefix, K_value);
            }
        }
    }
    const auto sum = [&](const std::vector<double>& prefix,
                         int x0, int y0, int z0,
                         int x1, int y1, int z1) {
        return prefix[pindex(x1, y1, z1)]
            - prefix[pindex(x0, y1, z1)] - prefix[pindex(x1, y0, z1)]
            - prefix[pindex(x1, y1, z0)] + prefix[pindex(x0, y0, z1)]
            + prefix[pindex(x0, y1, z0)] + prefix[pindex(x1, y0, z0)]
            - prefix[pindex(x0, y0, z0)];
    };
    r_counts.assign(box_size, 0.0);
    K_counts.assign(box_size, 0.0);
    for (int z = 0; z < depth; ++z) {
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                const int x0 = std::max(0, x - lower);
                const int y0 = std::max(0, y - lower);
                const int z0 = std::max(0, z - lower);
                const int x1 = std::min(width, x + upper + 1);
                const int y1 = std::min(height, y + upper + 1);
                const int z1 = std::min(depth, z + upper + 1);
                const std::size_t offset =
                    (static_cast<std::size_t>(z) * height + y) * width + x;
                r_counts[offset] = sum(r_prefix, x0, y0, z0, x1, y1, z1);
                K_counts[offset] = sum(K_prefix, x0, y0, z0, x1, y1, z1);
            }
        }
    }
}

void StructuredPdeModel3D::react(double dt) {
    if (!population_bounds_.valid) return;
    std::vector<double> r_counts;
    std::vector<double> K_counts;
    build_local_counts(r_counts, K_counts);
    const auto& continuum = config_.continuum;
    const double r_inherent = mean_growth_rate(CellType::r);
    const double K_inherent = mean_growth_rate(CellType::K);
    const double legacy_r_limit = continuum.base.thin_layer
        ? continuum.base.legacy_mapping.source_r_limit : continuum.base.r_limit;
    const double legacy_K_limit = continuum.base.thin_layer
        ? continuum.base.legacy_mapping.source_K_limit : continuum.base.K_limit;
    const double legacy_r_capacity = continuum.base.thin_layer
        ? continuum.base.legacy_mapping.source_carrying_capacity_r
        : continuum.base.carrying_capacity_r;
    const double legacy_K_capacity = continuum.base.thin_layer
        ? continuum.base.legacy_mapping.source_carrying_capacity_K
        : continuum.base.carrying_capacity_K;
    const double r_limit = config_.schema_version >= 5
        ? continuum.nutrient.common_density_limit : legacy_r_limit;
    const double K_limit = config_.schema_version >= 5
        ? continuum.nutrient.common_density_limit : legacy_K_limit;
    const double r_capacity = config_.schema_version >= 5
        ? continuum.nutrient.common_carrying_capacity : legacy_r_capacity;
    const double K_capacity = config_.schema_version >= 5
        ? continuum.nutrient.common_carrying_capacity : legacy_K_capacity;
    const int workers = std::max(
        1, std::min(continuum.base.threads, available_worker_threads()));
    const auto bounds = population_bounds_;
    const std::size_t bx = static_cast<std::size_t>(bounds.x1 - bounds.x0);
    const std::size_t by = static_cast<std::size_t>(bounds.y1 - bounds.y0);
    const std::size_t bz = static_cast<std::size_t>(bounds.z1 - bounds.z0);
    deterministic_parallel_for(bx * by * bz, workers, [&](std::size_t offset) {
        const int x = bounds.x0 + static_cast<int>(offset % bx);
        const std::size_t yz = offset / bx;
        const int y = bounds.y0 + static_cast<int>(yz % by);
        const int z = bounds.z0 + static_cast<int>(yz / by);
        const std::size_t location = index(x, y, z);
        if (vessel_blocks_cells(location)) return;
        const double occupied = occupied_fraction(location);
        if (occupied <= 0.0) return;
        const double multiplier = capacity_multiplier(nutrient_[location]);
        const double r_count = r_counts[offset] / multiplier;
        const double K_count = K_counts[offset] / multiplier;
        const double total_count = r_count + K_count;
        const double raw_r_growth = calculate_density_growth_rate_continuous(
            static_cast<int>(CellType::r), r_inherent,
            r_count, K_count, total_count, r_limit, K_limit,
            continuum.base.alpha, continuum.base.beta, r_capacity, K_capacity);
        const double raw_K_growth = calculate_density_growth_rate_continuous(
            static_cast<int>(CellType::K), K_inherent,
            r_count, K_count, total_count, r_limit, K_limit,
            continuum.base.alpha, continuum.base.beta, r_capacity, K_capacity);
        const double nutrient_factor = config_.schema_version >= 5
            ? nutrient_[location] /
                (continuum.nutrient.growth_half_saturation +
                 nutrient_[location])
            : 1.0;
        const auto resource_limited = [&](double growth) {
            return growth > 0.0 ? growth * nutrient_factor : growth;
        };
        const double r_growth = resource_limited(raw_r_growth);
        const double K_growth = resource_limited(raw_K_growth);
        const bool abm_work_clock = continuum.reaction.model ==
            "abm_work_clock_neighbor_availability_v2";
        const double r_division_rate = positive_part(r_growth) /
            (continuum.base.division_timing.base_cycle_hours *
             (abm_work_clock ? 1.0 : std::max(1.0e-12, r_inherent)));
        const double K_division_rate = positive_part(K_growth) /
            (continuum.base.division_timing.base_cycle_hours *
             (abm_work_clock ? 1.0 : std::max(1.0e-12, K_inherent)));
        const double r_death_rate =
            r_growth <= continuum.base.death_growth_rate_threshold
            ? 1.0 / continuum.base.r_death_delay_hours : 0.0;
        const double K_death_rate =
            K_growth <= continuum.base.death_growth_rate_threshold
            ? 1.0 / continuum.base.K_death_delay_hours : 0.0;
        const double vacancy = std::clamp(
            1.0 - occupied / continuum.reaction.maximum_occupied_fraction,
            0.0, 1.0);
        const double large_success = std::pow(
            vacancy, continuum.reaction.large_daughter_vacancy_exponent);
        const double small_success = abm_work_clock
            ? 1.0 - std::pow(
                  1.0 - vacancy,
                  continuum.reaction.small_daughter_vacancy_exponent)
            : std::pow(
                  vacancy,
                  continuum.reaction.small_daughter_vacancy_exponent);
        const auto conversion_probability = [&](std::size_t stage) {
            if (!continuum.base.r_to_K_conversion.enabled) return 0.0;
            const double density = config_.schema_version >= 3
                ? activation_density_[stage][location] : occupied;
            return density >= continuum.base.r_to_K_conversion.density_threshold
                ? continuum.base.r_to_K_conversion.probability_per_division
                : 0.0;
        };

        const std::array<double, 2> old_r{
            r_normal_[0][location] + r_active(StructuredStage3D::small, location),
            r_normal_[1][location] + r_active(StructuredStage3D::large, location)};
        const std::array<double, 2> old_K{K_[0][location], K_[1][location]};
        std::array<double, 2> delta_r{
            -r_death_rate * old_r[0], -r_death_rate * old_r[1]};
        std::array<double, 2> delta_K{
            -K_death_rate * old_K[0], -K_death_rate * old_K[1]};

        const double r_large_events = r_division_rate * old_r[1];
        const double r_large_daughters = large_success * r_large_events;
        const double r_shape_reductions =
            (1.0 - large_success) * r_large_events;
        const double large_conversion = conversion_probability(1);
        delta_r[1] += (1.0 - large_conversion) * r_large_daughters
            - r_shape_reductions;
        delta_K[1] += large_conversion * r_large_daughters;
        delta_r[0] += (2.0 - large_conversion) * r_shape_reductions;
        delta_K[0] += large_conversion * r_shape_reductions;

        const double K_large_events = K_division_rate * old_K[1];
        const double K_large_daughters = large_success * K_large_events;
        const double K_shape_reductions =
            (1.0 - large_success) * K_large_events;
        delta_K[1] += K_large_daughters - K_shape_reductions;
        delta_K[0] += 2.0 * K_shape_reductions;

        const double r_small_events = r_division_rate * old_r[0];
        const double r_small_births = small_success * r_small_events;
        const double small_conversion = conversion_probability(0);
        delta_r[0] += (1.0 - small_conversion) * r_small_births
            - (1.0 - small_success) * r_small_events *
                continuum.reaction.failed_r_division_death_fraction;
        delta_K[0] += small_conversion * r_small_births;
        delta_K[0] += small_success * K_division_rate * old_K[0];

        std::array<double, 2> positive_r{};
        std::array<double, 2> positive_K{};
        double base_occupied = 0.0;
        double positive_occupied = 0.0;
        for (std::size_t stage = 0; stage < 2; ++stage) {
            const double r_negative = dt * std::min(0.0, delta_r[stage]);
            const double r_factor = old_r[stage] > 0.0
                ? std::clamp((old_r[stage] + r_negative) / old_r[stage], 0.0, 1.0)
                : 0.0;
            r_normal_[stage][location] *= r_factor;
            if (config_.schema_version >= 7) {
                r_refractory_[stage][location] *= r_factor;
                refractory_clock_[stage][location] *= r_factor;
            }
            for (std::size_t bucket = 0;
                 bucket < active_direction_[stage].size(); ++bucket) {
                active_direction_[stage][bucket][location] = static_cast<float>(
                    active_direction_[stage][bucket][location] * r_factor);
                active_clock_[stage][bucket][location] = static_cast<float>(
                    active_clock_[stage][bucket][location] * r_factor);
            }
            if (config_.schema_version >= 7) {
                double active_sum = 0.0;
                for (const auto& bucket : active_direction_[stage]) {
                    active_sum += bucket[location];
                }
                active_total_[stage][location] = static_cast<float>(active_sum);
            } else {
                active_total_[stage][location] = static_cast<float>(
                    active_total_[stage][location] * r_factor);
            }
            K_[stage][location] = std::max(
                0.0, old_K[stage] + dt * std::min(0.0, delta_K[stage]));
            positive_r[stage] = dt * std::max(0.0, delta_r[stage]);
            positive_K[stage] = dt * std::max(0.0, delta_K[stage]);
            const double volume = stage == 0 ? 1.0 : large_cell_volume_;
            base_occupied += volume *
                (r_normal_[stage][location] +
                 r_active(static_cast<StructuredStage3D>(stage), location) +
                 K_[stage][location]);
            positive_occupied += volume *
                (positive_r[stage] + positive_K[stage]);
        }
        const double available = std::max(
            0.0, continuum.reaction.maximum_occupied_fraction - base_occupied);
        const double scale = positive_occupied > available && positive_occupied > 0.0
            ? available / positive_occupied : 1.0;
        for (std::size_t stage = 0; stage < 2; ++stage) {
            // All newly created r mass begins in the ordinary state, matching
            // the ABM division-cycle reset before a later density refresh.
            r_normal_[stage][location] += scale * positive_r[stage];
            K_[stage][location] += scale * positive_K[stage];
        }
    });
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
    return config_.schema_version >= 7 ? r_refractory_[static_cast<std::size_t>(stage)][location] : 0.0;
}

double StructuredPdeModel3D::refractory_mean_hours(
    StructuredStage3D stage, std::size_t location) const noexcept {
    const double mass = refractory_mass(stage, location);
    return mass > 0.0 ? refractory_clock_[static_cast<std::size_t>(stage)][location] / mass : 0.0;
}

bool StructuredPdeModel3D::step() {
    if (!initialized_) throw std::logic_error("structured PDE is not initialized");
    const auto& continuum = config_.continuum;
    if (time_hours_ >= continuum.end_time_hours ||
        same_time(time_hours_, continuum.end_time_hours)) return false;
    double target = std::min(
        time_hours_ + continuum.time_step_hours, continuum.end_time_hours);
    if (config_.schema_version < 5 &&
        next_nutrient_refresh_hours_ < target &&
        !same_time(next_nutrient_refresh_hours_, target)) {
        target = next_nutrient_refresh_hours_;
    }
    const double dt = target - time_hours_;
    // Density is only the trigger. Active cohorts keep their own remaining
    // clock and therefore do not deactivate when this field later falls.
    refresh_activation(dt);
    migrate_normal_and_K(dt);
    migrate_active(dt);
    exchange_active_r_with_K(dt);
    // V4 conversion is evaluated at the post-transport division location.
    // At 200x migration, reusing the pre-transport 70x70 field can otherwise
    // classify mass hundreds of voxels away using its source neighbourhood.
    if (config_.schema_version >= 4 &&
        continuum.base.r_to_K_conversion.enabled) {
        build_activation_density();
    }
    react(dt);
    // Diffusion expands the conservative work box by one voxel. Remove the
    // numerical halo after all operators have finished so subsequent steps
    // scale with the occupied support rather than the full production grid.
    shrink_population_bounds();
    time_hours_ = target;
    ++step_count_;
    if (config_.schema_version >= 5) {
        rebuild_moving_tumour_front();
        advance_transient_nutrient(dt);
        next_nutrient_refresh_hours_ =
            time_hours_ + continuum.time_step_hours;
    } else if (time_hours_ > next_nutrient_refresh_hours_ ||
               same_time(time_hours_, next_nutrient_refresh_hours_)) {
        solve_nutrient();
        do {
            next_nutrient_refresh_hours_ +=
                continuum.nutrient.refresh_every_hours;
        } while (time_hours_ > next_nutrient_refresh_hours_ ||
                 same_time(time_hours_, next_nutrient_refresh_hours_));
    }
    validate_state();
    return true;
}

StructuredPdeDiagnostics3D StructuredPdeModel3D::diagnostics() const {
    StructuredPdeDiagnostics3D result;
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
        if (config_.schema_version >= 6 && tumour_mask_[location] != 0U) {
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
    const bool local_validation = initialized_ && config_.schema_version >= 6 &&
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
                        if (config_.schema_version >= 7 &&
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
    if (config_.schema_version >= 5) {
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
    const auto hash_double_region = [&](const std::vector<double>& values,
                                        const StructuredActiveBounds3D& bounds) {
        if (!bounds.valid) return;
        for (int z = bounds.z0; z < bounds.z1; ++z) {
            for (int y = bounds.y0; y < bounds.y1; ++y) {
                for (int x = bounds.x0; x < bounds.x1; ++x) {
                    state = hash_mix(state, std::bit_cast<std::uint64_t>(
                        values[index(x, y, z)]));
                }
            }
        }
    };
    hash_bounds(population_bounds_);
    for (std::size_t stage = 0; stage < 2; ++stage) {
        hash_double_region(r_normal_[stage], population_bounds_);
        hash_double_region(K_[stage], population_bounds_);
        if (config_.schema_version >= 7) {
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
    if (config_.schema_version >= 5) {
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
    return state;
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
    write_pod(stream, config_.schema_version >= 7 ? kCohortCheckpointVersion : config_.schema_version >= 5
        ? kRefractoryCheckpointVersion : kLegacyCheckpointVersion);
    write_pod(stream, config_.dynamics_fingerprint());
    write_pod(stream, static_cast<std::uint64_t>(voxel_count_));
    write_pod(stream, static_cast<std::uint64_t>(direction_ids_.size()));
    write_pod(stream, time_hours_);
    write_pod(stream, next_nutrient_refresh_hours_);
    write_pod(stream, step_count_);
    write_pod(stream, nutrient_solve_count_);
    const int nx = config_.continuum.grid.shape[0];
    const int ny = config_.continuum.grid.shape[1];
    write_bounds(stream, population_bounds_);
    for (std::size_t stage = 0; stage < 2; ++stage) {
        write_region(stream, r_normal_[stage], population_bounds_, nx, ny);
        write_region(stream, K_[stage], population_bounds_, nx, ny);
        write_bounds(stream, active_bounds_[stage]);
        for (const auto& bucket : active_direction_[stage]) {
            write_region(stream, bucket, active_bounds_[stage], nx, ny);
        }
        for (const auto& bucket : active_clock_[stage]) {
            write_region(stream, bucket, active_bounds_[stage], nx, ny);
        }
    }
    if (config_.schema_version >= 5) {
        for (std::size_t stage = 0; stage < 2; ++stage) {
            write_vector(stream, activation_cooldown_[stage]);
            write_vector(stream, activation_armed_[stage]);
        }
    }
    if (config_.schema_version >= 7) {
        for (std::size_t stage = 0; stage < 2; ++stage) {
            write_region(stream, r_refractory_[stage], population_bounds_, nx, ny);
            write_region(stream, refractory_clock_[stage], population_bounds_, nx, ny);
        }
    }
    write_vector(stream, nutrient_);
    write_vector(stream, vessel_);
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
    const std::uint32_t checkpoint_version = read_pod<std::uint32_t>(stream);
    const auto expected_version = config_.schema_version >= 7 ? kCohortCheckpointVersion :
        config_.schema_version >= 5 ? kRefractoryCheckpointVersion : kLegacyCheckpointVersion;
    if (magic != kCheckpointMagic || checkpoint_version != expected_version) {
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
    population_bounds_ = read_bounds(stream, nx, ny, nz);
    for (std::size_t stage = 0; stage < 2; ++stage) {
        read_region(stream, r_normal_[stage], population_bounds_, nx, ny);
        read_region(stream, K_[stage], population_bounds_, nx, ny);
        active_bounds_[stage] = read_bounds(stream, nx, ny, nz);
        for (auto& bucket : active_direction_[stage]) {
            read_region(stream, bucket, active_bounds_[stage], nx, ny);
        }
        for (auto& bucket : active_clock_[stage]) {
            read_region(stream, bucket, active_bounds_[stage], nx, ny);
        }
        std::fill(active_total_[stage].begin(), active_total_[stage].end(), 0.0F);
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
    if (checkpoint_version >= kRefractoryCheckpointVersion) {
        for (std::size_t stage = 0; stage < 2; ++stage) {
            read_vector(stream, activation_cooldown_[stage], voxel_count_);
            read_vector(stream, activation_armed_[stage], voxel_count_);
        }
    }
    if (config_.schema_version >= 7) {
        for (std::size_t stage = 0; stage < 2; ++stage) {
            read_region(stream, r_refractory_[stage], population_bounds_, nx, ny);
            read_region(stream, refractory_clock_[stage], population_bounds_, nx, ny);
        }
    }
    read_vector(stream, nutrient_, voxel_count_);
    read_vector(stream, vessel_, voxel_count_);
    const std::uint64_t expected_checksum = read_pod<std::uint64_t>(stream);
    if (stream.peek() != std::char_traits<char>::eof()) {
        throw std::runtime_error("structured checkpoint has trailing data");
    }
    rebuild_moving_tumour_front();
    nutrient_next_ = nutrient_;
    validate_resources();
    initialized_ = true;
    validate_state();
    if (state_checksum() != expected_checksum) {
        throw std::runtime_error("structured checkpoint checksum mismatch");
    }
}

}  // namespace atcg3d::structured_pde
