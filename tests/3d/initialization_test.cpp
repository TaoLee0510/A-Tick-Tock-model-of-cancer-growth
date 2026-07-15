#include <cassert>
#include <cmath>
#include <cstdint>
#include <tuple>
#include <unordered_set>
#include <vector>

#include "config/model_config.hpp"
#include "core/cell_store.hpp"
#include "geometry/footprint.hpp"
#include "rules/initialization.hpp"
#include "space/chunk_grid.hpp"
#include "space/density_index.hpp"

namespace {

using namespace atcg3d;

struct Signature {
    CellUid uid{};
    Vec3i anchor{};
    CellType type{};
    CellStage stage{};
    float inherent_growth_rate{};
    float migration_rate{};

    bool operator==(const Signature&) const = default;
};

struct InitializedModel {
    explicit InitializedModel(const Model3DConfig& config)
        : grid(config.chunk_edge, DomainPolicy(config)), density(config.density_block_edge) {
        result = initialize_sphere_and_shell(cells, grid, density, config);
    }

    CellStore3D cells;
    SparseChunkGrid3D grid;
    BlockDensityIndex3D density;
    InitializationResult result;
};

std::int64_t radius_squared(Vec3i value) {
    const std::int64_t x = value.x;
    const std::int64_t y = value.y;
    const std::int64_t z = value.z;
    return x * x + y * y + z * z;
}

std::vector<Signature> signatures(const CellStore3D& cells) {
    std::vector<Signature> result;
    result.reserve(cells.alive_count());
    for (const Slot slot : cells.alive_slots()) {
        result.push_back({cells.uid(slot), cells.anchor(slot), cells.type(slot),
                          cells.stage(slot), cells.inherent_growth_rate(slot),
                          cells.migration_rate(slot)});
    }
    return result;
}

void verify_geometry(const InitializedModel& model,
                     const Model3DConfig& config,
                     bool require_nonzero_z) {
    assert(model.cells.alive_count() > 0);
    assert(model.cells.stage_count(CellStage::large) > 0);
    assert(model.cells.stage_count(CellStage::small) > 0);
    assert(model.result.next_uid == model.cells.alive_count() + 1);

    std::unordered_set<Vec3i, Vec3iHash> occupied;
    occupied.reserve(model.cells.stage_count(CellStage::large) * 8 +
                     model.cells.stage_count(CellStage::small));
    std::uint64_t r_count = 0;
    std::uint64_t K_count = 0;
    bool saw_nonzero_z = false;
    for (const Slot slot : model.cells.alive_slots()) {
        const Vec3i anchor = model.cells.anchor(slot);
        saw_nonzero_z = saw_nonzero_z || anchor.z != 0;
        if (model.cells.type(slot) == CellType::r) {
            ++r_count;
        } else {
            ++K_count;
        }
        if (config.initial_growth_rate_model == "fixed") {
            const double expected_growth = model.cells.type(slot) == CellType::r
                ? config.initial_r_growth_rate : config.initial_K_growth_rate;
            assert(model.cells.inherent_growth_rate(slot) ==
                   static_cast<float>(expected_growth));
        }
        if (config.initial_migration_rate_model == "fixed") {
            const double expected_migration = model.cells.type(slot) == CellType::r
                ? config.initial_r_migration_rate : config.initial_K_migration_rate;
            assert(model.cells.migration_rate(slot) ==
                   static_cast<float>(expected_migration));
        }

        if (model.cells.stage(slot) == CellStage::large) {
            assert(anchor.x % 2 == 0 && anchor.y % 2 == 0 && anchor.z % 2 == 0);
            assert(radius_squared(anchor) >=
                   static_cast<std::int64_t>(config.initial_shell_inner_radius) *
                       config.initial_shell_inner_radius);
            assert(radius_squared(anchor) <=
                   static_cast<std::int64_t>(config.initial_radius) *
                       config.initial_radius);
            for (const Vec3i voxel : large_footprint(anchor)) {
                assert(model.grid.owner(voxel) == slot);
                assert(occupied.insert(voxel).second);
            }
        } else {
            assert(model.cells.stage(slot) == CellStage::small);
            assert(radius_squared(anchor) <=
                   static_cast<std::int64_t>(config.initial_inner_small_radius) *
                       config.initial_inner_small_radius);
            assert(model.grid.owner(anchor) == slot);
            assert(occupied.insert(anchor).second);
        }
    }
    assert(!require_nonzero_z || saw_nonzero_z);
    assert(r_count + K_count == model.cells.alive_count());
    assert(r_count == static_cast<std::uint64_t>(
        std::floor(static_cast<long double>(model.cells.alive_count()) *
                   config.initial_r_fraction)));

    const int bound = config.initial_radius + 1;
    const DensityCounts3D counts = model.density.estimate_box(
        {-bound, -bound, config.thin_layer ? 0 : -bound},
        {bound, bound, config.thin_layer ? 0 : bound});
    assert(counts.r == r_count && counts.K == K_count);
}

Model3DConfig geometry_config() {
    Model3DConfig config;
    config.output_enabled = false;
    config.initialization_mode = "legacy_geometry_fill_v1";
    config.initial_r_cells = 0;
    config.initial_K_cells = 0;
    config.initial_radius = 12;
    config.initial_shell_inner_radius = 8;
    config.initial_shell_thickness = 4;
    config.initial_inner_small_radius = 10;
    config.initial_r_fraction = 0.375;
    config.initial_r_growth_rate = 1.25;
    config.initial_K_growth_rate = 0.75;
    config.initial_r_migration_rate = 0.4;
    config.initial_K_migration_rate = 0.2;
    config.migration_activation_enabled = false;
    return config;
}

}  // namespace

int main() {
    using namespace atcg3d;

    // A true 3D geometry fill derives counts from a shell and inner sphere,
    // never overlaps footprints, and is reproducible for the same seed.
    {
        const Model3DConfig config = geometry_config();
        const InitializedModel first(config);
        const InitializedModel second(config);
        verify_geometry(first, config, true);
        verify_geometry(second, config, true);
        assert(signatures(first.cells) == signatures(second.cells));
    }

    // Changing the seed changes the spatial r/K assignment while retaining the
    // exact requested type counts and identical geometry.
    {
        Model3DConfig first_config = geometry_config();
        Model3DConfig second_config = first_config;
        second_config.seed = 9001;
        const InitializedModel first(first_config);
        const InitializedModel second(second_config);
        verify_geometry(first, first_config, true);
        verify_geometry(second, second_config, true);
        const auto lhs = signatures(first.cells);
        const auto rhs = signatures(second.cells);
        assert(lhs.size() == rhs.size());
        bool type_assignment_changed = false;
        for (std::size_t index = 0; index < lhs.size(); ++index) {
            assert(lhs[index].anchor == rhs[index].anchor);
            assert(lhs[index].stage == rhs[index].stage);
            type_assignment_changed = type_assignment_changed || lhs[index].type != rhs[index].type;
        }
        assert(type_assignment_changed);
    }

    // Thin-layer mode uses the same radial mapping in x/y, fixes all anchors at
    // z=0, and reserves z=1 only for the upper half of a stage-0 footprint.
    {
        Model3DConfig config = geometry_config();
        config.thin_layer = true;
        const InitializedModel model(config);
        verify_geometry(model, config, false);
        for (const Slot slot : model.cells.alive_slots()) {
            assert(model.cells.anchor(slot).z == 0);
            if (model.cells.stage(slot) == CellStage::large) {
                bool saw_upper_voxel = false;
                for (const Vec3i voxel : large_footprint(model.cells.anchor(slot))) {
                    saw_upper_voxel = saw_upper_voxel || voxel.z == 1;
                }
                assert(saw_upper_voxel);
            }
        }
    }

    // Explicit-count mode retains the small smoke-profile behavior and exact
    // requested r/K counts instead of filling the entire geometry.
    {
        Model3DConfig config;
        config.output_enabled = false;
        config.initialization_mode = "explicit_counts";
        config.initial_r_cells = 5;
        config.initial_K_cells = 7;
        config.initial_radius = 8;
        config.initial_shell_thickness = 2;
        config.initial_large_fraction = 0.5;
        config.migration_activation_enabled = false;
        const InitializedModel model(config);
        assert(model.cells.alive_count() == 12);
        assert(model.cells.stage_count(CellStage::large) == 6);
        assert(model.cells.stage_count(CellStage::small) == 6);
        std::size_t r_count = 0;
        for (const Slot slot : model.cells.alive_slots()) {
            r_count += model.cells.type(slot) == CellType::r ? 1 : 0;
        }
        assert(r_count == 5);
    }
}
