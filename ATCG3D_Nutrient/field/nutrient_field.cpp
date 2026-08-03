#include "field/nutrient_field.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <fstream>
#include <limits>
#include <set>
#include <stdexcept>
#include <type_traits>

#include "core/cell_store.hpp"
#include "geometry/footprint.hpp"
#include "vasculature/vessel_grid.hpp"

namespace atcg3d::nutrient {
namespace {

constexpr std::array<char, 8> kCheckpointMagic{{'A', 'T', 'C', 'G', 'N', 'U', 'T', '1'}};
constexpr std::uint32_t kCheckpointVersion = 1;

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
    if (!stream) throw std::runtime_error("unable to write nutrient checkpoint");
}

template <class T>
T read_pod(std::istream& stream) {
    static_assert(std::is_trivially_copyable_v<T>);
    T value{};
    stream.read(reinterpret_cast<char*>(&value), sizeof(value));
    if (!stream) throw std::runtime_error("truncated nutrient checkpoint");
    return value;
}

}  // namespace

NutrientEnvironment3D::NutrientEnvironment3D(
    NutrientFieldConfig3D config, bool thin_layer)
    : config_(std::move(config)),
      thin_layer_(thin_layer),
      block_voxels_(static_cast<std::size_t>(config_.block_edge) *
                    static_cast<std::size_t>(config_.block_edge) *
                    static_cast<std::size_t>(config_.block_edge)) {
    config_.validate();
}

int NutrientEnvironment3D::floor_div(int value, int divisor) noexcept {
    int quotient = value / divisor;
    const int remainder = value % divisor;
    if (remainder != 0 && ((remainder < 0) != (divisor < 0))) --quotient;
    return quotient;
}

NutrientEnvironment3D::Address NutrientEnvironment3D::address(
    Vec3i site) const noexcept {
    const Vec3i block{floor_div(site.x, config_.block_edge),
                      floor_div(site.y, config_.block_edge),
                      floor_div(site.z, config_.block_edge)};
    const int local_x = site.x - block.x * config_.block_edge;
    const int local_y = site.y - block.y * config_.block_edge;
    const int local_z = site.z - block.z * config_.block_edge;
    return {block, static_cast<std::uint32_t>(
                       (local_z * config_.block_edge + local_y) *
                           config_.block_edge + local_x)};
}

Vec3i NutrientEnvironment3D::site_from(Vec3i block,
                                       std::uint32_t index) const noexcept {
    const auto edge = static_cast<std::uint32_t>(config_.block_edge);
    const int local_x = static_cast<int>(index % edge);
    const std::uint32_t yz = index / edge;
    const int local_y = static_cast<int>(yz % edge);
    const int local_z = static_cast<int>(yz / edge);
    return {block.x * config_.block_edge + local_x,
            block.y * config_.block_edge + local_y,
            block.z * config_.block_edge + local_z};
}

NutrientEnvironment3D::Block* NutrientEnvironment3D::find_block(
    Vec3i coordinate) noexcept {
    const auto found = blocks_.find(coordinate);
    return found == blocks_.end() ? nullptr : &found->second;
}

const NutrientEnvironment3D::Block* NutrientEnvironment3D::find_block(
    Vec3i coordinate) const noexcept {
    const auto found = blocks_.find(coordinate);
    return found == blocks_.end() ? nullptr : &found->second;
}

bool NutrientEnvironment3D::active_site(Vec3i site) const noexcept {
    return !thin_layer_ || site.z == 0;
}

float NutrientEnvironment3D::value(Vec3i site) const noexcept {
    if (!active_site(site)) return 0.0F;
    const Address location = address(site);
    const Block* block = find_block(location.block);
    return block == nullptr ? 0.0F : block->value[location.index];
}

double NutrientEnvironment3D::capacity_multiplier(Vec3i site) const noexcept {
    const double local = std::clamp(
        static_cast<double>(value(site)), 0.0, config_.vessel_value);
    const double raw = local / (config_.capacity_half_saturation + local);
    const double at_vessel = config_.vessel_value /
        (config_.capacity_half_saturation + config_.vessel_value);
    const double saturation = at_vessel > 0.0
        ? std::clamp(raw / at_vessel, 0.0, 1.0)
        : 0.0;
    return 1.0 +
        (config_.maximum_capacity_multiplier - 1.0) * saturation;
}

double NutrientEnvironment3D::retained_density(Vec3i site) const noexcept {
    return 1.0 / capacity_multiplier(site);
}

void NutrientEnvironment3D::rebuild_sources_and_sinks(
    const CellStore3D& cells,
    const SparseVesselGrid3D& vessels) {
    std::set<Vec3i> active_blocks;
    std::set<Vec3i> source_blocks;
    const int block_radius =
        (config_.halo_voxels + config_.block_edge - 1) /
            config_.block_edge +
        1;
    const auto add_halo = [&](Vec3i center) {
        for (int dx = -block_radius; dx <= block_radius; ++dx) {
            for (int dy = -block_radius; dy <= block_radius; ++dy) {
                const int minimum_z = thin_layer_ ? 0 : -block_radius;
                const int maximum_z = thin_layer_ ? 0 : block_radius;
                for (int dz = minimum_z; dz <= maximum_z; ++dz) {
                    const std::int64_t x = static_cast<std::int64_t>(center.x) + dx;
                    const std::int64_t y = static_cast<std::int64_t>(center.y) + dy;
                    const std::int64_t z = thin_layer_ ? 0 :
                        static_cast<std::int64_t>(center.z) + dz;
                    if (x < std::numeric_limits<std::int32_t>::min() ||
                        x > std::numeric_limits<std::int32_t>::max() ||
                        y < std::numeric_limits<std::int32_t>::min() ||
                        y > std::numeric_limits<std::int32_t>::max() ||
                        z < std::numeric_limits<std::int32_t>::min() ||
                        z > std::numeric_limits<std::int32_t>::max()) {
                        throw std::overflow_error("nutrient active block coordinate overflow");
                    }
                    active_blocks.insert({static_cast<std::int32_t>(x),
                                          static_cast<std::int32_t>(y),
                                          static_cast<std::int32_t>(z)});
                }
            }
        }
    };

    for (const Slot slot : cells.alive_slots()) {
        source_blocks.insert(address(cells.anchor(slot)).block);
    }
    const std::vector<Vec3i> vessel_sites = vessels.occupied_sites();
    for (const Vec3i site : vessel_sites) {
        if (vessels.perfused(site)) source_blocks.insert(address(site).block);
    }
    for (const Vec3i source_block : source_blocks) add_halo(source_block);

    std::map<Vec3i, Block> rebuilt;
    for (const Vec3i coordinate : active_blocks) {
        Block block(block_voxels_);
        const auto previous = blocks_.find(coordinate);
        if (previous != blocks_.end()) block.value = previous->second.value;
        rebuilt.emplace(coordinate, std::move(block));
    }
    blocks_ = std::move(rebuilt);

    const auto add_consumption = [&](Vec3i site, CellType type) {
        if (!active_site(site)) return;
        const Address location = address(site);
        Block* block = find_block(location.block);
        if (block == nullptr) return;
        if (type == CellType::r) {
            block->r_consumption[location.index] += static_cast<float>(
                config_.r_consumption_per_voxel_hour);
        } else {
            block->K_consumption[location.index] += static_cast<float>(
                config_.K_consumption_per_voxel_hour);
        }
    };
    for (const Slot slot : cells.alive_slots()) {
        if (cells.stage(slot) == CellStage::large) {
            for (const Vec3i site : large_footprint(cells.anchor(slot))) {
                add_consumption(site, cells.type(slot));
            }
        } else {
            add_consumption(cells.anchor(slot), cells.type(slot));
        }
    }
    for (const Vec3i site : vessel_sites) {
        if (!vessels.perfused(site) || !active_site(site)) continue;
        const Address location = address(site);
        if (Block* block = find_block(location.block)) {
            block->perfused[location.index] = 1U;
        }
    }
}

void NutrientEnvironment3D::solve() {
    const double diffusion = config_.diffusion_voxels2_per_hour;
    const int dimensions = thin_layer_ ? 2 : 3;
    const double laplacian_diagonal = 2.0 * dimensions * diffusion;
    const std::array<Vec3i, 6> offsets{{
        {-1, 0, 0}, {1, 0, 0}, {0, -1, 0},
        {0, 1, 0}, {0, 0, -1}, {0, 0, 1}}};
    double maximum_update = 0.0;
    for (int iteration = 0; iteration < config_.solver_iterations; ++iteration) {
        maximum_update = 0.0;
        for (auto& [coordinate, block] : blocks_) {
            for (std::uint32_t index = 0; index < block.value.size(); ++index) {
                const Vec3i site = site_from(coordinate, index);
                if (!active_site(site)) {
                    block.next[index] = 0.0F;
                    continue;
                }
                const double old = block.value[index];
                double neighbor_sum = 0.0;
                const int neighbor_count = thin_layer_ ? 4 : 6;
                for (int direction = 0; direction < neighbor_count; ++direction) {
                    neighbor_sum += value(site + offsets[direction]);
                }
                const double r_sink = block.r_consumption[index] /
                    (config_.r_consumption_half_saturation + old);
                const double K_sink = block.K_consumption[index] /
                    (config_.K_consumption_half_saturation + old);
                const double vessel_exchange = block.perfused[index] != 0
                    ? config_.vessel_exchange_per_hour : 0.0;
                const double denominator = laplacian_diagonal +
                    config_.decay_per_hour + vessel_exchange + r_sink + K_sink;
                const double right_hand_side = diffusion * neighbor_sum +
                    vessel_exchange * config_.vessel_value;
                const double candidate = denominator > 0.0
                    ? right_hand_side / denominator : 0.0;
                const double relaxed = std::clamp(
                    old + config_.relaxation * (candidate - old),
                    0.0, config_.vessel_value);
                block.next[index] = static_cast<float>(relaxed);
                maximum_update = std::max(maximum_update, std::abs(relaxed - old));
            }
        }
        for (auto& [coordinate, block] : blocks_) {
            (void)coordinate;
            block.value.swap(block.next);
        }
    }
    diagnostics_.last_max_update = maximum_update;
}

void NutrientEnvironment3D::update_diagnostics() {
    diagnostics_.block_count = blocks_.size();
    diagnostics_.active_voxel_count = 0;
    diagnostics_.perfused_source_voxels = 0;
    diagnostics_.consuming_voxels = 0;
    diagnostics_.minimum = std::numeric_limits<double>::infinity();
    diagnostics_.maximum = 0.0;
    long double sum = 0.0L;
    for (const auto& [coordinate, block] : blocks_) {
        for (std::uint32_t index = 0; index < block.value.size(); ++index) {
            const Vec3i site = site_from(coordinate, index);
            if (!active_site(site)) continue;
            ++diagnostics_.active_voxel_count;
            const double current = block.value[index];
            diagnostics_.minimum = std::min(diagnostics_.minimum, current);
            diagnostics_.maximum = std::max(diagnostics_.maximum, current);
            sum += current;
            if (block.perfused[index] != 0) ++diagnostics_.perfused_source_voxels;
            if (block.r_consumption[index] > 0.0F ||
                block.K_consumption[index] > 0.0F) {
                ++diagnostics_.consuming_voxels;
            }
        }
    }
    if (diagnostics_.active_voxel_count == 0) {
        diagnostics_.minimum = 0.0;
        diagnostics_.mean = 0.0;
    } else {
        diagnostics_.mean = static_cast<double>(
            sum / diagnostics_.active_voxel_count);
    }
}

EnvironmentInitializationResult3D NutrientEnvironment3D::initialize(
    double now_hours,
    const CellStore3D& cells,
    const SparseVesselGrid3D& vessels) {
    if (!std::isfinite(now_hours) || now_hours < 0.0) {
        throw std::invalid_argument("nutrient initialization time is invalid");
    }
    if (loaded_checkpoint_) {
        if (!same_time(now_hours, loaded_checkpoint_time_)) {
            throw std::runtime_error(
                "nutrient checkpoint time does not match simulation restore time");
        }
        loaded_checkpoint_ = false;
        update_diagnostics();
        return {false};
    }
    rebuild_sources_and_sinks(cells, vessels);
    solve();
    last_refresh_time_hours_ = now_hours;
    next_refresh_time_hours_ = now_hours + config_.refresh_every_hours;
    schedule_generation_ = 1;
    refresh_count_ = 1;
    update_diagnostics();
    return {true};
}

void NutrientEnvironment3D::refresh(
    double now_hours,
    const CellStore3D& cells,
    const SparseVesselGrid3D& vessels) {
    if (!same_time(now_hours, next_refresh_time_hours_)) {
        throw std::logic_error("nutrient refresh occurred outside its scheduled time");
    }
    if (schedule_generation_ == std::numeric_limits<std::uint32_t>::max() ||
        refresh_count_ == std::numeric_limits<std::uint64_t>::max()) {
        throw std::overflow_error("nutrient refresh counter overflow");
    }
    rebuild_sources_and_sinks(cells, vessels);
    solve();
    last_refresh_time_hours_ = now_hours;
    next_refresh_time_hours_ = now_hours + config_.refresh_every_hours;
    ++schedule_generation_;
    ++refresh_count_;
    update_diagnostics();
}

std::size_t NutrientEnvironment3D::allocated_bytes() const noexcept {
    std::size_t result = blocks_.size() * sizeof(Block);
    for (const auto& [coordinate, block] : blocks_) {
        (void)coordinate;
        result += block.value.capacity() * sizeof(float) +
            block.next.capacity() * sizeof(float) +
            block.r_consumption.capacity() * sizeof(float) +
            block.K_consumption.capacity() * sizeof(float) +
            block.perfused.capacity() * sizeof(std::uint8_t);
    }
    return result;
}

std::uint64_t NutrientEnvironment3D::field_checksum() const noexcept {
    std::uint64_t result = config_.fingerprint();
    result = hash_mix(result, std::bit_cast<std::uint64_t>(last_refresh_time_hours_));
    result = hash_mix(result, std::bit_cast<std::uint64_t>(next_refresh_time_hours_));
    result = hash_mix(result, schedule_generation_);
    result = hash_mix(result, refresh_count_);
    for (const auto& [coordinate, block] : blocks_) {
        result = hash_mix(result, static_cast<std::uint32_t>(coordinate.x));
        result = hash_mix(result, static_cast<std::uint32_t>(coordinate.y));
        result = hash_mix(result, static_cast<std::uint32_t>(coordinate.z));
        for (const float current : block.value) {
            result = hash_mix(result, std::bit_cast<std::uint32_t>(current));
        }
    }
    return result;
}

std::vector<NutrientVoxelSample3D> NutrientEnvironment3D::nonzero_voxels() const {
    std::vector<NutrientVoxelSample3D> result;
    for (const auto& [coordinate, block] : blocks_) {
        for (std::uint32_t index = 0; index < block.value.size(); ++index) {
            if (block.value[index] <= 0.0F) continue;
            const Vec3i site = site_from(coordinate, index);
            if (active_site(site)) result.push_back({site, block.value[index]});
        }
    }
    return result;
}

void NutrientEnvironment3D::save_checkpoint(
    const std::filesystem::path& path,
    std::uint64_t base_state_checksum,
    double time_hours,
    std::uint64_t completed_events) const {
    if (!same_time(time_hours, last_refresh_time_hours_) &&
        time_hours < last_refresh_time_hours_) {
        throw std::logic_error("nutrient checkpoint precedes field state");
    }
    if (std::filesystem::exists(path)) {
        throw std::runtime_error("refusing to overwrite nutrient checkpoint: " +
                                 path.string());
    }
    std::filesystem::create_directories(path.parent_path());
    const std::filesystem::path temporary = path.string() + ".tmp";
    std::filesystem::remove(temporary);
    std::ofstream stream(temporary, std::ios::binary | std::ios::trunc);
    if (!stream) {
        throw std::runtime_error("unable to create nutrient checkpoint: " +
                                 temporary.string());
    }
    stream.write(kCheckpointMagic.data(), kCheckpointMagic.size());
    write_pod(stream, kCheckpointVersion);
    write_pod(stream, config_.fingerprint());
    write_pod(stream, base_state_checksum);
    write_pod(stream, std::bit_cast<std::uint64_t>(time_hours));
    write_pod(stream, completed_events);
    write_pod(stream, std::bit_cast<std::uint64_t>(last_refresh_time_hours_));
    write_pod(stream, std::bit_cast<std::uint64_t>(next_refresh_time_hours_));
    write_pod(stream, schedule_generation_);
    write_pod(stream, refresh_count_);
    write_pod(stream, static_cast<std::uint32_t>(config_.block_edge));
    write_pod(stream, static_cast<std::uint64_t>(blocks_.size()));
    for (const auto& [coordinate, block] : blocks_) {
        write_pod(stream, coordinate.x);
        write_pod(stream, coordinate.y);
        write_pod(stream, coordinate.z);
        stream.write(reinterpret_cast<const char*>(block.value.data()),
                     static_cast<std::streamsize>(block.value.size() * sizeof(float)));
        if (!stream) throw std::runtime_error("unable to write nutrient field values");
    }
    stream.flush();
    if (!stream) throw std::runtime_error("unable to finish nutrient checkpoint");
    stream.close();
    std::filesystem::rename(temporary, path);
}

void NutrientEnvironment3D::load_checkpoint(
    const std::filesystem::path& path,
    std::uint64_t expected_base_state_checksum,
    double expected_time_hours,
    std::uint64_t expected_completed_events) {
    std::ifstream stream(path, std::ios::binary);
    if (!stream) {
        throw std::runtime_error("unable to open nutrient checkpoint: " +
                                 path.string());
    }
    std::array<char, 8> magic{};
    stream.read(magic.data(), magic.size());
    if (!stream || magic != kCheckpointMagic ||
        read_pod<std::uint32_t>(stream) != kCheckpointVersion) {
        throw std::runtime_error("unsupported nutrient checkpoint format");
    }
    if (read_pod<std::uint64_t>(stream) != config_.fingerprint()) {
        throw std::runtime_error("nutrient checkpoint configuration mismatch");
    }
    if (read_pod<std::uint64_t>(stream) != expected_base_state_checksum) {
        throw std::runtime_error("nutrient/base checkpoint checksum mismatch");
    }
    const double time_hours = std::bit_cast<double>(read_pod<std::uint64_t>(stream));
    if (!same_time(time_hours, expected_time_hours) ||
        read_pod<std::uint64_t>(stream) != expected_completed_events) {
        throw std::runtime_error("nutrient/base checkpoint clock mismatch");
    }
    last_refresh_time_hours_ =
        std::bit_cast<double>(read_pod<std::uint64_t>(stream));
    next_refresh_time_hours_ =
        std::bit_cast<double>(read_pod<std::uint64_t>(stream));
    schedule_generation_ = read_pod<std::uint32_t>(stream);
    refresh_count_ = read_pod<std::uint64_t>(stream);
    if (!std::isfinite(last_refresh_time_hours_) ||
        !std::isfinite(next_refresh_time_hours_) ||
        last_refresh_time_hours_ < 0.0 ||
        next_refresh_time_hours_ <= expected_time_hours ||
        last_refresh_time_hours_ > expected_time_hours ||
        schedule_generation_ == 0 || refresh_count_ == 0) {
        throw std::runtime_error("nutrient checkpoint schedule is invalid");
    }
    const std::uint32_t block_edge = read_pod<std::uint32_t>(stream);
    const std::uint64_t block_count = read_pod<std::uint64_t>(stream);
    if (block_edge != static_cast<std::uint32_t>(config_.block_edge) ||
        block_count > 100000000ULL) {
        throw std::runtime_error("nutrient checkpoint field dimensions are invalid");
    }
    blocks_.clear();
    for (std::uint64_t block_index = 0; block_index < block_count; ++block_index) {
        const Vec3i coordinate{read_pod<std::int32_t>(stream),
                               read_pod<std::int32_t>(stream),
                               read_pod<std::int32_t>(stream)};
        Block block(block_voxels_);
        stream.read(reinterpret_cast<char*>(block.value.data()),
                    static_cast<std::streamsize>(block.value.size() * sizeof(float)));
        if (!stream || !blocks_.emplace(coordinate, std::move(block)).second) {
            throw std::runtime_error("nutrient checkpoint blocks are invalid");
        }
    }
    if (stream.peek() != std::char_traits<char>::eof()) {
        throw std::runtime_error("nutrient checkpoint has trailing data");
    }
    for (const auto& [coordinate, block] : blocks_) {
        (void)coordinate;
        for (const float current : block.value) {
            if (!std::isfinite(current) || current < 0.0F ||
                current > static_cast<float>(config_.vessel_value)) {
                throw std::runtime_error("nutrient checkpoint contains invalid values");
            }
        }
    }
    loaded_checkpoint_ = true;
    loaded_checkpoint_time_ = expected_time_hours;
    update_diagnostics();
}

}  // namespace atcg3d::nutrient
