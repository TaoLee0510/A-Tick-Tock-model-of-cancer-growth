#include <algorithm>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <memory>

#include "engine/simulation.hpp"
#include "geometry/footprint.hpp"
#include "model/shared_resource_environment.hpp"
#include "model/structured_pde_model.hpp"
#include "rules/density.hpp"

int main() {
    using namespace atcg3d;
    using namespace atcg3d::shared_rules;
    using namespace atcg3d::structured_pde;
    auto config = StructuredPdeConfig3D::load(std::filesystem::path(ATCG_SOURCE_DIR) /
        "ATCG3D_SharedRules/config/validation_v7.yaml");
    config.continuum.grid.shape = {12, 12, 1};
    config.continuum.grid.origin = {-6.0, -6.0, -0.5};
    config.continuum.vascular.source_mode = "static_voxels";
    config.continuum.vascular.static_sources = {{3, 0, 0}};
    config.continuum.nutrient.boundary_mode = "vessels_dirichlet_v1";
    config.continuum.migration.diffusion_scale = 1.0e-12;
    config.continuum.base.division_timing.base_cycle_hours = 1.0e30;
    config.continuum.base.growth_density_window_edge = 1;
    SharedResourceEnvironment3D resource(config);
    CellStore3D cells;
    BlockDensityIndex3D density(1);
    SparseVesselGrid3D vessels(4, DomainPolicy(config.continuum.base));
    StructuredInitialFields3D fields;
    for (std::size_t stage = 0; stage < 2; ++stage) {
        fields.r_normal[stage].assign(144, 0.0);
        fields.r_active[stage].assign(144, 0.0);
        fields.K[stage].assign(144, 0.0);
        fields.active_remaining_hours[stage].assign(144, 0.0);
    }
    fields.vessel_fraction.assign(144, 0.0);
    for (int y = -6; y < 6; ++y) {
        for (int x = -6; x < 6; ++x) {
            if (x == 3 && y == 0) continue;
            CellInit cell;
            cell.uid = static_cast<CellUid>((y + 6) * 12 + x + 7);
            cell.type = CellType::K;
            cell.stage = CellStage::small;
            cell.anchor = {x, y, 0};
            const auto slot = cells.create(cell);
            density.add(cell.anchor, cell.type, slot);
            fields.K[0][static_cast<std::size_t>((y + 6) * 12 + x + 6)] = 1.0;
        }
    }
    resource.initialize(0.0, cells, vessels);
    StructuredPdeModel3D pde(config);
    pde.initialize_from_arrays(fields);
    for (int step = 1; step <= 3; ++step) {
        resource.refresh(step * config.continuum.time_step_hours, cells, vessels);
        assert(pde.step());
        assert(resource.nutrient() == pde.nutrient());
    }
    assert(resource.vessel_fraction() == pde.vessel_fraction());
    assert(resource.normalized_resource({3, 0, 0}) == 1.0);
    const auto slot = cells.alive_slots().front();
    const auto counts = growth_counts(density, cells.anchor(slot), 1, true);
    const double expected = calculate_density_growth_rate_continuous(static_cast<int>(CellType::K),
        cells.inherent_growth_rate(slot), counts.rc, counts.kc, counts.cells_number,
        resource.shared_density_limit(), resource.shared_density_limit(), config.continuum.base.alpha,
        config.continuum.base.beta, resource.shared_carrying_capacity(), resource.shared_carrying_capacity()) *
        resource.growth_resource_scale(cells.anchor(slot));
    assert(density_growth_rate_for_cell(cells, slot, density, config.continuum.base, &resource) == expected);

    // The cached row spans implement the literal bounded 70x70 sector sum.
    for (DirectionId direction = 1; direction <= 26; ++direction) {
        const auto forward = direction_vector(direction);
        if (forward.z != 0) continue;
        const Vec3i anchor{0, 0, 0};
        double sum = 0.0;
        std::size_t count = 0;
        const double cosine_limit = std::cos(config.continuum.base.direction_density_half_angle_degrees * std::acos(-1.0) / 180.0);
        for (int dy = -34; dy <= 35; ++dy) {
            for (int dx = -34; dx <= 35; ++dx) {
                const Vec3i offset{dx, dy, 0};
                if (squared_length(offset) == 0 || !resource.contains_resource_site(anchor + offset)) continue;
                const double cosine = static_cast<double>(dot(offset, forward)) /
                    std::sqrt(static_cast<double>(squared_length(offset) * squared_length(forward)));
                if (cosine + 1.0e-12 < cosine_limit) continue;
                sum += resource.normalized_resource(anchor + offset) * config.continuum.nutrient.vessel_value;
                ++count;
            }
        }
        const double gradient = sum / count - resource.normalized_resource(anchor) * config.continuum.nutrient.vessel_value;
        const double expected_weight = std::exp(config.migration.chemotaxis_strength * gradient);
        assert(std::abs(resource.nutrient_direction_weight(anchor, direction) - expected_weight) < 1.0e-12);
    }
    resource.activation_expired(41, 1.0);
    assert(!resource.activation_ready(41, 2.0, 0.1));
    assert(resource.activation_ready(42, 2.0, 0.95));
    assert(!resource.activation_ready(41, 25.0, 0.95));
    assert(resource.activation_ready(41, 25.0, 0.5));
    assert(resource.activation_ready(41, 25.0, 0.95));

    // Match the second moment of eight equiprobable unit lattice jumps:
    // E[dx^2 + dy^2] = 3/2, so D = 3 lambda / 8.
    auto diffusion_config = config;
    diffusion_config.continuum.migration.diffusion_scale = 1.0;
    diffusion_config.continuum.migration.crowding_exponent = 0.0;
    diffusion_config.continuum.vascular.source_mode = "abm_perfusion";
    diffusion_config.continuum.vascular.static_sources.clear();
    diffusion_config.continuum.nutrient.initial_value = 1.0;
    for (auto& field : fields.K) std::fill(field.begin(), field.end(), 0.0);
    fields.K[0][6 * 12 + 6] = 1.0;
    StructuredPdeModel3D diffusion(diffusion_config);
    diffusion.initialize_from_arrays(fields);
    assert(diffusion.step());
    double second_moment = 0.0;
    for (int y = 0; y < 12; ++y) {
        for (int x = 0; x < 12; ++x) {
            second_moment += ((x - 6) * (x - 6) + (y - 6) * (y - 6)) *
                diffusion.K(StructuredStage3D::small, y * 12 + x);
        }
    }
    const auto& beta = diffusion_config.continuum.base.initial_K_migration_beta;
    const double rate = beta.alpha / (beta.alpha + beta.beta) * beta.scale;
    assert(std::abs(second_moment - 1.5 * rate * diffusion_config.continuum.time_step_hours) < 1.0e-12);

    // Resource state and UID-specific cooldowns survive an event-time restart,
    // including a snapshot between nutrient ticks and a changed thread count.
    config = StructuredPdeConfig3D::load(std::filesystem::path(ATCG_SOURCE_DIR) /
        "ATCG3D_SharedRules/config/validation_v7.yaml");
    config.continuum.end_time_hours = 16.0;
    auto base = abm_config(config);
    auto environment = std::make_unique<SharedResourceEnvironment3D>(config);
    auto* original_resource = environment.get();
    Simulation3D original(base, std::move(environment));
    original.initialize();
    for (int event = 0; event < 50; ++event) assert(original.step());
    const auto checkpoint = std::filesystem::temp_directory_path() / "atcg_shared_rules.resource.bin";
    std::filesystem::remove(checkpoint);
    original_resource->save_checkpoint(checkpoint, original.state_checksum());
    auto parallel_config = config;
    parallel_config.continuum.base.threads = 4;
    auto restored_environment = std::make_unique<SharedResourceEnvironment3D>(parallel_config);
    auto* restored_resource = restored_environment.get();
    restored_resource->load_checkpoint(checkpoint, original.state_checksum());
    auto parallel_base = abm_config(parallel_config);
    Simulation3D resumed(parallel_base, std::move(restored_environment));
    resumed.restore(original.snapshot_cells(), original.next_uid(), original.clock(), original.stats(),
        original.lineage(), original.snapshot_vasculature(), original.cells().slot_count(),
        original.snapshot_cell_slots(), original.cells().free_slots());
    assert(original.state_checksum() == resumed.state_checksum());
    assert(original_resource->field_checksum() == restored_resource->field_checksum());
    for (int event = 0; event < 100; ++event) {
        assert(original.step() && resumed.step());
        assert(original.state_checksum() == resumed.state_checksum());
        assert(original_resource->field_checksum() == restored_resource->field_checksum());
    }
    std::filesystem::remove(checkpoint);

    // Three-dimensional sectors are built once, then reused by every event.
    parallel_config.continuum.base.thin_layer = false;
    parallel_config.continuum.reaction.small_daughter_vacancy_exponent = 26.0;
    parallel_config.continuum.grid.shape = {8, 8, 8};
    parallel_config.continuum.grid.origin = {-4.0, -4.0, -4.0};
    parallel_config.continuum.nutrient.boundary_mode = "vessels_dirichlet_v1";
    parallel_config.continuum.vascular.source_mode = "abm_perfusion";
    SharedResourceEnvironment3D spatial(parallel_config);
    assert(spatial.allocated_bytes() < 64U * 1024U * 1024U);
    for (DirectionId direction = 1; direction <= 26; ++direction) {
        assert(spatial.nutrient_direction_weight({0, 0, 0}, direction) == 1.0);
    }
}
