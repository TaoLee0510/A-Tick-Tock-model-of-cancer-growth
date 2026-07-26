#include "rules/initialization.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "core/stateless_rng.hpp"
#include "rules/density.hpp"
#include "rules/initial_rates.hpp"
#include "rules/migration.hpp"

namespace atcg3d {
namespace {

std::int64_t squared(int value) noexcept {
    const std::int64_t wide = value;
    return wide * wide;
}

int integer_square_root(std::int64_t value) {
    if (value < 0) {
        throw std::invalid_argument("cannot take the square root of a negative integer");
    }
    std::int64_t root = static_cast<std::int64_t>(
        std::floor(std::sqrt(static_cast<long double>(value))));
    while ((root + 1) <= value / (root + 1)) {
        ++root;
    }
    while (root > 0 && root > value / root) {
        --root;
    }
    if (root > std::numeric_limits<int>::max()) {
        throw std::overflow_error("initial geometry radius exceeds integer coordinates");
    }
    return static_cast<int>(root);
}

template <class Visitor>
void for_each_sphere_site(int radius, bool thin_layer, Visitor&& visitor) {
    if (radius < 0) {
        throw std::invalid_argument("initial sphere radius must be non-negative");
    }
    const std::int64_t radius_squared = squared(radius);
    const int minimum_z = thin_layer ? 0 : -radius;
    const int maximum_z = thin_layer ? 0 : radius;
    for (int z = minimum_z; z <= maximum_z; ++z) {
        const std::int64_t after_z = radius_squared - squared(z);
        const int maximum_y = integer_square_root(after_z);
        for (int y = -maximum_y; y <= maximum_y; ++y) {
            const int maximum_x = integer_square_root(after_z - squared(y));
            for (int x = -maximum_x; x <= maximum_x; ++x) {
                visitor(Vec3i{x, y, z});
            }
        }
    }
}

// A stage-0 cell represents one site of the legacy coarse lattice. Mapping a
// coarse coordinate c to anchor 2*c gives a disjoint tiling by 2x2x2
// footprints. The radial test is performed in the full-resolution coordinates,
// so the YAML radii retain their legacy physical meaning.
template <class Visitor>
void for_each_shell_anchor(int outer_radius,
                           int inner_radius,
                           bool thin_layer,
                           Visitor&& visitor) {
    if (outer_radius < 0 || inner_radius < 0 || inner_radius >= outer_radius) {
        throw std::invalid_argument("initial shell radii are inconsistent");
    }
    const std::int64_t outer_squared = squared(outer_radius);
    const std::int64_t inner_squared = squared(inner_radius);
    const int coarse_limit = outer_radius / 2;
    const int minimum_z = thin_layer ? 0 : -coarse_limit;
    const int maximum_z = thin_layer ? 0 : coarse_limit;
    for (int coarse_z = minimum_z; coarse_z <= maximum_z; ++coarse_z) {
        const std::int64_t after_z = outer_squared - 4 * squared(coarse_z);
        if (after_z < 0) {
            continue;
        }
        const int maximum_y = integer_square_root(after_z / 4);
        for (int coarse_y = -maximum_y; coarse_y <= maximum_y; ++coarse_y) {
            const std::int64_t after_y = after_z - 4 * squared(coarse_y);
            const int maximum_x = integer_square_root(after_y / 4);
            for (int coarse_x = -maximum_x; coarse_x <= maximum_x; ++coarse_x) {
                const std::int64_t coarse_squared =
                    squared(coarse_x) + squared(coarse_y) + squared(coarse_z);
                const std::int64_t anchor_squared = 4 * coarse_squared;
                if (anchor_squared < inner_squared) {
                    continue;
                }
                visitor(Vec3i{2 * coarse_x, 2 * coarse_y, 2 * coarse_z});
            }
        }
    }
}

std::vector<Vec3i> sphere_candidates(int radius, bool thin_layer) {
    std::vector<Vec3i> result;
    for_each_sphere_site(radius, thin_layer,
                         [&](Vec3i site) { result.push_back(site); });
    return result;
}

std::vector<Vec3i> shell_anchor_candidates(int outer_radius,
                                           int inner_radius,
                                           bool thin_layer) {
    std::vector<Vec3i> result;
    for_each_shell_anchor(outer_radius, inner_radius, thin_layer,
                          [&](Vec3i anchor) { result.push_back(anchor); });
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
    const InitialCellRates3D rates =
        sample_initial_cell_rates(config, type, config.seed, uid);
    cell.inherent_growth_rate = static_cast<float>(rates.inherent_growth_rate);
    cell.density_growth_rate = cell.inherent_growth_rate;
    cell.migration_rate = static_cast<float>(rates.migration_rate);
    cell.normal_migration_rate = sample_normal_migration_rate(
        type, cell.migration_rate, config, uid, 0);
    cell.last_update_time = 0.0;
    return cell;
}

void apply_initial_type(CellStore3D& cells,
                        Slot slot,
                        CellType type,
                        const Model3DConfig& config) {
    cells.set_type(slot, type);
    const InitialCellRates3D rates = sample_initial_cell_rates(
        config, type, config.seed, cells.uid(slot));
    const float growth = static_cast<float>(rates.inherent_growth_rate);
    cells.set_inherent_growth_rate(slot, growth);
    cells.set_density_growth_rate(slot, growth);
    cells.set_migration_rate(slot, static_cast<float>(rates.migration_rate));
    cells.set_normal_migration_rate(
        slot, sample_normal_migration_rate(
                  type, static_cast<float>(rates.migration_rate), config,
                  cells.uid(slot), 0));
    cells.set_migration_activation_end_time(slot, 0.0);
    cells.set_flags(slot, cells.flags(slot) &
        static_cast<std::uint8_t>(~kMigrationActive));
}

class DeterministicTypeAllocator {
public:
    DeterministicTypeAllocator(std::uint64_t r_count,
                               std::uint64_t K_count,
                               std::uint64_t seed,
                               std::uint64_t stream)
        : r_remaining_(r_count), K_remaining_(K_count), seed_(seed), stream_(stream) {}

    CellType next(CellUid uid) {
        if (r_remaining_ == 0) {
            if (K_remaining_ == 0) {
                throw std::logic_error("initial type allocation exhausted");
            }
            --K_remaining_;
            return CellType::K;
        }
        if (K_remaining_ == 0) {
            --r_remaining_;
            return CellType::r;
        }
        const std::uint64_t remaining = r_remaining_ + K_remaining_;
        const bool choose_r = rng_bounded(
            seed_, uid, static_cast<std::uint64_t>(RngEventKind::initialization),
            stream_, remaining) < r_remaining_;
        if (choose_r) {
            --r_remaining_;
            return CellType::r;
        }
        --K_remaining_;
        return CellType::K;
    }

private:
    std::uint64_t r_remaining_{};
    std::uint64_t K_remaining_{};
    std::uint64_t seed_{};
    std::uint64_t stream_{};
};

std::uint64_t fraction_count(std::uint64_t total, double fraction) {
    if (fraction <= 0.0) {
        return 0;
    }
    if (fraction >= 1.0) {
        return total;
    }
    // The 2D model converted the product to an integer, which truncated the
    // fractional remainder. Preserve that rule while evaluating in long double.
    return static_cast<std::uint64_t>(
        std::floor(static_cast<long double>(total) * fraction));
}

void require_cell_capacity(std::uint64_t count) {
    if (count >= static_cast<std::uint64_t>(kEmptySlot)) {
        throw std::overflow_error("initial geometry exceeds stable-slot capacity");
    }
}

CellUid initialize_explicit_counts(CellStore3D& cells,
                                   SparseChunkGrid3D& grid,
                                   BlockDensityIndex3D& density,
                                   const Model3DConfig& config) {
    if (config.initial_r_cells >
        std::numeric_limits<std::uint64_t>::max() - config.initial_K_cells) {
        throw std::overflow_error("explicit initial cell count overflow");
    }
    const std::uint64_t total = config.initial_r_cells + config.initial_K_cells;
    require_cell_capacity(total);
    cells.reserve(static_cast<std::size_t>(total));

    // Keep the original smoke-test geometry: shell thickness determines both
    // the stage-0 annulus and the stage-1 sampling sphere.
    const int radius_limit = config.initial_radius;
    const int inner_shell = std::max(0, radius_limit - config.initial_shell_thickness);
    std::vector<Vec3i> shell = shell_anchor_candidates(
        radius_limit, inner_shell, config.thin_layer);
    std::vector<Vec3i> volume = sphere_candidates(inner_shell, config.thin_layer);
    deterministic_shuffle(shell.begin(), shell.end(), config.seed, 1,
                          static_cast<std::uint64_t>(RngEventKind::initialization), 0);
    deterministic_shuffle(volume.begin(), volume.end(), config.seed, 2,
                          static_cast<std::uint64_t>(RngEventKind::initialization), 0);

    const std::uint64_t desired_large = fraction_count(total, config.initial_large_fraction);
    DeterministicTypeAllocator types(config.initial_r_cells, config.initial_K_cells,
                                     config.seed, 1);
    CellUid next_uid = 1;
    std::size_t shell_index = 0;
    std::size_t volume_index = 0;
    for (std::uint64_t index = 0; index < total; ++index) {
        const bool large = index < desired_large;
        Vec3i anchor{};
        bool found = false;
        auto& candidates = large ? shell : volume;
        std::size_t& candidate_index = large ? shell_index : volume_index;
        while (candidate_index < candidates.size()) {
            anchor = candidates[candidate_index++];
            if ((large && grid.can_place_large(anchor)) ||
                (!large && grid.available(anchor))) {
                found = true;
                break;
            }
        }
        if (!found) {
            throw std::runtime_error("initial sphere/shell is too small for requested cell count");
        }
        const CellType type = types.next(next_uid);
        const CellStage stage = large ? CellStage::large : CellStage::small;
        const Slot slot = cells.create(initial_cell(next_uid++, type, stage, anchor, config));
        const bool placed = large ? grid.place_large(anchor, slot) : grid.place_single(anchor, slot);
        if (!placed) {
            throw std::logic_error("initialization selected an unavailable footprint");
        }
        density.add(anchor, type, slot);
    }
    return next_uid;
}

CellUid initialize_legacy_geometry_fill(CellStore3D& cells,
                                        SparseChunkGrid3D& grid,
                                        BlockDensityIndex3D& density,
                                        const Model3DConfig& config) {
    std::vector<Vec3i> shell = shell_anchor_candidates(
        config.initial_radius, config.initial_shell_inner_radius, config.thin_layer);
    if (shell.size() >= static_cast<std::size_t>(kEmptySlot)) {
        throw std::overflow_error("initial geometry exceeds stable-slot capacity");
    }
    cells.reserve(shell.size());

    // Place the non-overlapping shell first. Type-dependent fields are assigned
    // after the inner-cell count is known so the requested fraction is exact.
    CellUid next_uid = 1;
    std::vector<Slot> large_slots;
    large_slots.reserve(shell.size());
    for (const Vec3i anchor : shell) {
        if (!grid.can_place_large(anchor)) {
            continue;
        }
        const Slot slot = cells.create(initial_cell(
            next_uid++, CellType::r, CellStage::large, anchor, config));
        if (!grid.place_large(anchor, slot)) {
            throw std::logic_error("geometry fill failed to place a feasible large footprint");
        }
        large_slots.push_back(slot);
    }

    std::uint64_t small_count = 0;
    for_each_sphere_site(config.initial_inner_small_radius, config.thin_layer,
                         [&](Vec3i site) {
                             if (grid.available(site)) {
                                 ++small_count;
                             }
                         });
    const std::uint64_t large_count = large_slots.size();
    if (large_count > std::numeric_limits<std::uint64_t>::max() - small_count) {
        throw std::overflow_error("initial cell count overflow");
    }
    const std::uint64_t total = large_count + small_count;
    require_cell_capacity(total);
    if (total == 0) {
        throw std::runtime_error("legacy geometry fill produced no placeable cells");
    }
    // Reserve only after overlap with the stage-0 shell has been accounted for.
    // Reserving the whole unfiltered inner sphere over-allocated roughly 31% in
    // the production 60/50/55 geometry.
    cells.reserve(static_cast<std::size_t>(total));

    const std::uint64_t r_count = fraction_count(total, config.initial_r_fraction);
    DeterministicTypeAllocator types(r_count, total - r_count, config.seed, 2);
    for (const Slot slot : large_slots) {
        const CellType type = types.next(cells.uid(slot));
        apply_initial_type(cells, slot, type, config);
        density.add(cells.anchor(slot), type, slot);
    }

    for_each_sphere_site(config.initial_inner_small_radius, config.thin_layer,
                         [&](Vec3i site) {
                             if (!grid.available(site)) {
                                 return;
                             }
                             const CellType type = types.next(next_uid);
                             const Slot slot = cells.create(initial_cell(
                                 next_uid++, type, CellStage::small, site, config));
                             if (!grid.place_single(site, slot)) {
                                 throw std::logic_error(
                                     "geometry fill failed to place a feasible small cell");
                             }
                             density.add(site, type, slot);
                         });
    if (cells.alive_count() != total) {
        throw std::logic_error("legacy geometry fill count changed between enumeration passes");
    }
    return next_uid;
}

void initialize_schedules(CellStore3D& cells,
                          BlockDensityIndex3D& density,
                          const Model3DConfig& config) {
    for (const Slot slot : cells.alive_slots()) {
        refresh_growth_state(slot, 0.0, cells, density, config);
        initialize_division_cycle(slot, 0.0, cells, config);
        (void)refresh_migration_activation_state(
            slot, 0.0, cells, density, config);
        const double migration_interval =
            migration_allowed_for_cell(slot, cells, config) &&
                    effective_migration_rate(slot, cells, config) > 0.0
            ? 1.0 / effective_migration_rate(slot, cells, config) : 0.0;
        cells.set_next_migration_time(slot, migration_interval);
    }
}

}  // namespace

InitializationResult initialize_sphere_and_shell(CellStore3D& cells,
                                                  SparseChunkGrid3D& grid,
                                                  BlockDensityIndex3D& density,
                                                  const Model3DConfig& config) {
    if (cells.alive_count() != 0) {
        throw std::invalid_argument("initialization requires an empty CellStore3D");
    }

    CellUid next_uid = 1;
    if (config.initialization_mode == "legacy_geometry_fill_v1") {
        next_uid = initialize_legacy_geometry_fill(cells, grid, density, config);
    } else if (config.initialization_mode == "explicit_counts") {
        next_uid = initialize_explicit_counts(cells, grid, density, config);
    } else {
        throw std::invalid_argument("unsupported initial initialization mode");
    }
    density.finish_local_window_bulk_load(cells.slot_count(), config.threads);
    initialize_schedules(cells, density, config);
    return {next_uid, {}};
}

}  // namespace atcg3d
