#include "rules/initialization.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "core/stateless_rng.hpp"
#include "rules/density.hpp"
#include "rules/migration.hpp"

namespace atcg3d {
namespace {

std::int64_t radius_squared(Vec3i site) {
    return static_cast<std::int64_t>(site.x) * site.x +
           static_cast<std::int64_t>(site.y) * site.y +
           static_cast<std::int64_t>(site.z) * site.z;
}

std::vector<Vec3i> sphere_candidates(int radius, bool thin_layer) {
    std::vector<Vec3i> result;
    const std::int64_t limit = static_cast<std::int64_t>(radius) * radius;
    const int minimum_z = thin_layer ? 0 : -radius;
    const int maximum_z = thin_layer ? 0 : radius;
    for (int z = minimum_z; z <= maximum_z; ++z) {
        const std::int64_t z_squared = static_cast<std::int64_t>(z) * z;
        const int maximum_y = static_cast<int>(std::floor(std::sqrt(
            static_cast<double>(limit - z_squared))));
        for (int y = -maximum_y; y <= maximum_y; ++y) {
            const std::int64_t remaining = limit - z_squared -
                static_cast<std::int64_t>(y) * y;
            const int maximum_x = static_cast<int>(std::floor(std::sqrt(
                static_cast<double>(remaining))));
            for (int x = -maximum_x; x <= maximum_x; ++x) {
                result.push_back({x, y, z});
            }
        }
    }
    return result;
}

std::vector<Vec3i> shell_anchor_candidates(int outer_radius,
                                           int inner_radius,
                                           bool thin_layer) {
    std::vector<Vec3i> result;
    const std::int64_t outer_squared = static_cast<std::int64_t>(outer_radius) * outer_radius;
    const std::int64_t inner_squared = static_cast<std::int64_t>(inner_radius) * inner_radius;
    const int minimum_z = thin_layer ? 0 : -outer_radius;
    const int maximum_z = thin_layer ? 0 : outer_radius;
    for (int z = minimum_z; z <= maximum_z; ++z) {
        if (z % 2 != 0) continue;
        const std::int64_t z_squared = static_cast<std::int64_t>(z) * z;
        const int maximum_y = static_cast<int>(std::floor(std::sqrt(
            static_cast<double>(outer_squared - z_squared))));
        for (int y = -maximum_y; y <= maximum_y; ++y) {
            if (y % 2 != 0) continue;
            const std::int64_t remaining = outer_squared - z_squared -
                static_cast<std::int64_t>(y) * y;
            const int maximum_x = static_cast<int>(std::floor(std::sqrt(
                static_cast<double>(remaining))));
            for (int x = -maximum_x; x <= maximum_x; ++x) {
                if (x % 2 != 0) continue;
                const Vec3i candidate{x, y, z};
                if (radius_squared(candidate) >= inner_squared) {
                    result.push_back(candidate);
                }
            }
        }
    }
    return result;
}

CellInit initial_cell(CellUid uid, CellType type, CellStage stage, Vec3i anchor,
                      const Model3DConfig& config) {
    CellInit cell;
    cell.anchor = anchor;
    cell.uid = uid;
    cell.parent_uid = 0;
    cell.clone_id = static_cast<std::uint32_t>(uid);
    cell.type = type;
    cell.stage = stage;
    cell.inherent_growth_rate = static_cast<float>(type == CellType::r
        ? config.initial_r_growth_rate : config.initial_K_growth_rate);
    cell.density_growth_rate = cell.inherent_growth_rate;
    cell.migration_rate = static_cast<float>(type == CellType::r
        ? config.initial_r_migration_rate : config.initial_K_migration_rate);
    cell.last_update_time = 0.0;
    return cell;
}

}  // namespace

InitializationResult initialize_sphere_and_shell(CellStore3D& cells,
                                                  SparseChunkGrid3D& grid,
                                                  BlockDensityIndex3D& density,
                                                  const Model3DConfig& config) {
    const std::uint64_t total = config.initial_r_cells + config.initial_K_cells;
    cells.reserve(static_cast<std::size_t>(total));

    const int radius_limit = config.initial_radius;
    const int inner_shell = std::max(0, radius_limit - config.initial_shell_thickness);
    std::vector<Vec3i> shell_candidates = shell_anchor_candidates(
        radius_limit, inner_shell, config.thin_layer);
    std::vector<Vec3i> volume_candidates = sphere_candidates(inner_shell, config.thin_layer);
    deterministic_shuffle(shell_candidates.begin(), shell_candidates.end(), config.seed, 1,
                          static_cast<std::uint64_t>(RngEventKind::initialization), 0);
    deterministic_shuffle(volume_candidates.begin(), volume_candidates.end(), config.seed, 2,
                          static_cast<std::uint64_t>(RngEventKind::initialization), 0);

    const std::uint64_t desired_large = total / 2;
    std::uint64_t r_remaining = config.initial_r_cells;
    std::uint64_t K_remaining = config.initial_K_cells;
    CellUid next_uid = 1;
    std::size_t shell_index = 0;
    std::size_t volume_index = 0;

    auto next_type = [&]() {
        if (r_remaining == 0) {
            --K_remaining;
            return CellType::K;
        }
        if (K_remaining == 0) {
            --r_remaining;
            return CellType::r;
        }
        const bool choose_r = rng_bounded(config.seed, next_uid,
                                          static_cast<std::uint64_t>(RngEventKind::initialization),
                                          1, r_remaining + K_remaining) < r_remaining;
        if (choose_r) {
            --r_remaining;
            return CellType::r;
        }
        --K_remaining;
        return CellType::K;
    };

    for (std::uint64_t index = 0; index < total; ++index) {
        const bool large = index < desired_large;
        Vec3i anchor{};
        bool found = false;
        auto& candidates = large ? shell_candidates : volume_candidates;
        std::size_t& candidate_index = large ? shell_index : volume_index;
        while (candidate_index < candidates.size()) {
            anchor = candidates[candidate_index++];
            if ((large && grid.can_place_large(anchor)) || (!large && grid.available(anchor))) {
                found = true;
                break;
            }
        }
        if (!found) {
            throw std::runtime_error("initial sphere/shell is too small for requested cell count");
        }
        const CellType type = next_type();
        CellInit cell = initial_cell(next_uid++, type,
                                     large ? CellStage::large : CellStage::small,
                                     anchor, config);
        const Slot slot = cells.create(cell);
        const bool placed = large ? grid.place_large(anchor, slot) : grid.place_single(anchor, slot);
        if (!placed) {
            throw std::logic_error("initialization selected an unavailable footprint");
        }
        density.add(anchor, type);
    }

    for (const Slot slot : cells.alive_slots()) {
        refresh_growth_state(slot, 0.0, cells, density, config);
        const double migration_interval = cells.migration_rate(slot) > 0.0F
            ? 1.0 / static_cast<double>(cells.migration_rate(slot)) : 0.0;
        cells.set_next_migration_time(slot, migration_interval);
    }

    return {next_uid, {}};
}

}  // namespace atcg3d
