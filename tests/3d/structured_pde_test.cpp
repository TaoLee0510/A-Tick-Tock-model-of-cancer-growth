#include <cassert>
#include <cmath>
#include <filesystem>
#include <iostream>

#include "config/structured_config.hpp"
#include "engine/simulation.hpp"
#include "model/moving_tumor_front.hpp"
#include "model/structured_pde_model.hpp"

namespace {

bool close(double lhs, double rhs, double tolerance = 1.0e-6) {
    return std::abs(lhs - rhs) <= tolerance *
        std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

atcg3d::structured_pde::StructuredInitialFields3D empty_fields(
    std::size_t voxel_count) {
    atcg3d::structured_pde::StructuredInitialFields3D result;
    for (std::size_t stage = 0; stage < 2; ++stage) {
        result.r_normal[stage].assign(voxel_count, 0.0);
        result.r_active[stage].assign(voxel_count, 0.0);
        result.K[stage].assign(voxel_count, 0.0);
        result.active_remaining_hours[stage].assign(voxel_count, 0.0);
    }
    result.vessel_fraction.assign(voxel_count, 0.0);
    return result;
}

}  // namespace

int main() {
    using namespace atcg3d;
    using namespace atcg3d::structured_pde;

    const auto wrapper = std::filesystem::path(ATCG_SOURCE_DIR) /
        "ATCG3D_StructuredPDE/config/structured_legacy_2d_2000_r20_v1.yaml";
    StructuredPdeConfig3D loaded = StructuredPdeConfig3D::load(wrapper);
    assert(loaded.schema_version == 1);
    assert(loaded.continuum.base.thin_layer);
    assert(loaded.continuum.base.activated_r_migration_rate_model ==
           "normal_multiplier");
    assert(loaded.continuum.base.activated_r_normal_multiplier == 20.0);
    assert(loaded.continuum.base.migration_activation_window_edge == 70);
    assert(loaded.continuum.base.migration_activation_threshold == 0.90);
    assert(loaded.to_json().find("beta_mean_remaining_cycle_v1") !=
           std::string::npos);

    const auto wrapper_v2 = std::filesystem::path(ATCG_SOURCE_DIR) /
        "ATCG3D_StructuredPDE/config/"
        "structured_legacy_2d_2000_r20_resource_guided_v2.yaml";
    StructuredPdeConfig3D loaded_v2 = StructuredPdeConfig3D::load(wrapper_v2);
    assert(loaded_v2.schema_version == 2);
    assert(loaded_v2.continuum.nutrient.consumption_model ==
           "per_cell_ratio_v2");
    assert(close(
        loaded_v2.continuum.nutrient.r_consumption_rate_per_hour /
            loaded_v2.continuum.nutrient.K_consumption_rate_per_hour,
        1.2));
    assert(loaded_v2.continuum.base.direction_guidance_model ==
           "low_density_high_resource_v1");

    const auto wrapper_v3 = std::filesystem::path(ATCG_SOURCE_DIR) /
        "ATCG3D_StructuredPDE/config/"
        "structured_legacy_2d_2000_r20_exchange_v3.yaml";
    StructuredPdeConfig3D loaded_v3 = StructuredPdeConfig3D::load(wrapper_v3);
    assert(loaded_v3.schema_version == 3);
    assert(loaded_v3.migration.direction_density_window_edge == 70);
    assert(loaded_v3.migration.crowding_exchange ==
           "active_r_K_stage1_conservative_v1");
    assert(loaded_v3.migration.vessel_exclusion);

    const auto wrapper_v4 = std::filesystem::path(ATCG_SOURCE_DIR) /
        "ATCG3D_StructuredPDE/config/"
        "structured_legacy_2d_2000_r200_reaction_fixed_v4.yaml";
    StructuredPdeConfig3D loaded_v4 = StructuredPdeConfig3D::load(wrapper_v4);
    assert(loaded_v4.schema_version == 4);
    assert(loaded_v4.continuum.base.activated_r_normal_multiplier == 200.0);
    assert(loaded_v4.continuum.migration.activated_r_mobility_multiplier ==
           200.0);
    assert(loaded_v4.continuum.reaction.model ==
           "abm_work_clock_neighbor_availability_v2");
    assert(loaded_v4.continuum.reaction.small_daughter_vacancy_exponent == 8.0);

    const auto wrapper_v5 = std::filesystem::path(ATCG_SOURCE_DIR) /
        "ATCG3D_StructuredPDE_NutrientChemotaxis/config/"
        "structured_smoke_2d_256_v5.yaml";
    StructuredPdeConfig3D loaded_v5 = StructuredPdeConfig3D::load(wrapper_v5);
    assert(loaded_v5.schema_version == 5);
    assert(loaded_v5.continuum.schema_version == 3);
    assert(loaded_v5.migration.direction_nutrient_window_edge == 70);
    assert(loaded_v5.migration.direction_transport ==
           "nutrient_gradient_fixed_direction_jump_exchange_v4");
    assert(loaded_v5.continuum.nutrient.boundary_mode ==
           "planar_edges_and_vessels_dirichlet_v1");
    assert(close(loaded_v5.continuum.nutrient.initial_value, 1.0));
    assert(close(loaded_v5.continuum.nutrient.r_consumption_rate_per_hour,
                 loaded_v5.continuum.nutrient.K_consumption_rate_per_hour));
    assert(close(loaded_v5.continuum.nutrient.maximum_capacity_multiplier, 1.0));

    const auto wrapper_v6 = std::filesystem::path(ATCG_SOURCE_DIR) /
        "ATCG3D_StructuredPDE_NutrientChemotaxis/config/"
        "structured_smoke_2d_256_24h_v6.yaml";
    StructuredPdeConfig3D loaded_v6 = StructuredPdeConfig3D::load(wrapper_v6);
    assert(loaded_v6.schema_version == 6);
    assert(loaded_v6.continuum.schema_version == 4);
    assert(loaded_v6.continuum.nutrient.boundary_mode ==
           "moving_tumor_front_and_vessels_dirichlet_v2");
    assert(close(
        loaded_v6.continuum.nutrient.tumor_front_density_threshold, 0.05));
    assert(loaded_v6.continuum.nutrient.
        tumor_front_smoothing_radius_voxels == 4);

    // Validation supports large 2D domains without allocating their fields.
    // The aggregate 500-million-voxel safety guard remains in force.
    auto large_domain = loaded_v6.continuum;
    large_domain.grid.shape = {10000, 10000, 1};
    large_domain.grid.origin = {-5000.0, -5000.0, -0.5};
    large_domain.validate();

    // The moving-front geometry retains the main connected lesion, fills
    // enclosed low-density holes, and rejects a disconnected satellite.
    {
        constexpr int width = 41;
        constexpr int height = 41;
        std::vector<double> occupied(width * height, 0.0);
        const auto local_at = [](int x, int y) {
            return static_cast<std::size_t>(y) * width + x;
        };
        for (int y = 5; y <= 30; ++y) {
            for (int x = 5; x <= 30; ++x) occupied[local_at(x, y)] = 1.0;
        }
        for (int y = 14; y <= 18; ++y) {
            for (int x = 14; x <= 18; ++x) occupied[local_at(x, y)] = 0.0;
        }
        for (int y = 35; y <= 37; ++y) {
            for (int x = 35; x <= 37; ++x) occupied[local_at(x, y)] = 1.0;
        }
        continuum::MovingTumorFrontWorkspace2D workspace;
        std::vector<std::uint8_t> mask;
        const auto front = continuum::build_moving_tumor_front_mask_2d(
            occupied, width, height, 0, 0.5, workspace, mask);
        assert(mask[local_at(16, 16)] == 1U);
        assert(mask[local_at(36, 36)] == 0U);
        assert(front.tumour_voxels == 26U * 26U);
        assert(front.front_voxels > 0U);
        assert(front.front_voxels < front.tumour_voxels);
    }

    // The shared ABM contract is literal per-cell multiplication, not a second
    // unrelated Beta draw.
    Model3DConfig base = loaded.continuum.base;
    base.output_enabled = false;
    base.control_enabled = false;
    Simulation3D abm(base);
    abm.initialize();
    bool saw_r = false;
    for (const Slot slot : abm.cells().alive_slots()) {
        if (abm.cells().type(slot) != CellType::r) continue;
        saw_r = true;
        if (!close(abm.cells().migration_rate(slot),
                   20.0 * abm.cells().normal_migration_rate(slot))) {
            std::cerr << "r active=" << abm.cells().migration_rate(slot)
                      << " normal=" << abm.cells().normal_migration_rate(slot)
                      << " model=" << base.activated_r_migration_rate_model
                      << " multiplier=" << base.activated_r_normal_multiplier
                      << '\n';
            assert(false);
        }
    }
    assert(saw_r);

    StructuredPdeConfig3D small = loaded;
    small.continuum.output.enabled = false;
    small.continuum.grid.shape = {128, 128, 1};
    small.continuum.grid.origin = {-64.0, -64.0, -0.5};
    small.continuum.time_step_hours = 0.25;
    small.continuum.end_time_hours = 3.0;
    small.validate();
    const std::size_t voxel_count = 128U * 128U;
    const auto at = [](int x, int y) {
        return static_cast<std::size_t>(y) * 128U +
            static_cast<std::size_t>(x);
    };

    StructuredPdeConfig3D small_v2 = loaded_v2;
    small_v2.continuum.output.enabled = false;
    small_v2.continuum.grid.shape = {128, 128, 1};
    small_v2.continuum.grid.origin = {-64.0, -64.0, -0.5};
    small_v2.continuum.time_step_hours = 0.25;
    small_v2.continuum.end_time_hours = 3.0;
    small_v2.continuum.vascular.source_mode = "abm_perfusion";
    small_v2.validate();

    StructuredPdeConfig3D small_v3 = loaded_v3;
    small_v3.continuum.output.enabled = false;
    small_v3.continuum.grid.shape = {128, 128, 1};
    small_v3.continuum.grid.origin = {-64.0, -64.0, -0.5};
    small_v3.continuum.time_step_hours = 0.25;
    small_v3.continuum.end_time_hours = 3.0;
    small_v3.continuum.vascular.source_mode = "abm_perfusion";
    small_v3.validate();

    StructuredPdeConfig3D small_v4 = loaded_v4;
    small_v4.continuum.output.enabled = false;
    small_v4.continuum.grid.shape = {128, 128, 1};
    small_v4.continuum.grid.origin = {-64.0, -64.0, -0.5};
    small_v4.continuum.time_step_hours = 0.25;
    small_v4.continuum.end_time_hours = 0.25;
    small_v4.continuum.vascular.source_mode = "abm_perfusion";
    small_v4.continuum.migration.diffusion_scale = 1.0e-12;
    small_v4.continuum.base.r_to_K_conversion.enabled = false;
    small_v4.validate();

    StructuredPdeConfig3D small_v5 = loaded_v5;
    small_v5.continuum.output.enabled = false;
    small_v5.continuum.grid.shape = {128, 128, 1};
    small_v5.continuum.grid.origin = {-64.0, -64.0, -0.5};
    small_v5.continuum.time_step_hours = 0.25;
    small_v5.continuum.end_time_hours = 1.0;
    small_v5.continuum.vascular.source_mode = "abm_perfusion";
    small_v5.validate();

    StructuredPdeConfig3D small_v6 = loaded_v6;
    small_v6.continuum.output.enabled = false;
    small_v6.continuum.grid.shape = {128, 128, 1};
    small_v6.continuum.grid.origin = {-64.0, -64.0, -0.5};
    small_v6.continuum.time_step_hours = 0.25;
    small_v6.continuum.end_time_hours = 1.0;
    small_v6.continuum.vascular.source_mode = "abm_perfusion";
    small_v6.continuum.nutrient.boundary_mode =
        "moving_tumor_front_dirichlet_v2";
    small_v6.validate();

    // The ABM stores approximately one base-cycle of work and depletes it at
    // the density growth rate. The dilute PDE birth intensity must therefore
    // retain the different inherent r and K rates instead of cancelling them.
    const auto dilute_relative_growth = [&](CellType type) {
        auto fields = empty_fields(voxel_count);
        constexpr double initial = 1.0e-3;
        (type == CellType::r ? fields.r_normal[0] : fields.K[0])
            [at(64, 64)] = initial;
        StructuredPdeModel3D model(small_v4);
        model.initialize_from_arrays(std::move(fields));
        const auto before = model.diagnostics();
        assert(model.step());
        const auto after = model.diagnostics();
        return type == CellType::r
            ? after.r_total / before.r_total - 1.0
            : after.K_total / before.K_total - 1.0;
    };
    const double dilute_r_growth = dilute_relative_growth(CellType::r);
    const double dilute_K_growth = dilute_relative_growth(CellType::K);
    const double mean_r = loaded_v4.continuum.base
        .initial_r_growth_truncated_normal.mean;
    const double mean_K = loaded_v4.continuum.base
        .initial_K_growth_truncated_normal.mean;
    assert(close(dilute_r_growth,
                 0.25 * mean_r /
                     loaded_v4.continuum.base.division_timing.base_cycle_hours,
                 1.0e-5));
    assert(close(dilute_K_growth,
                 0.25 * mean_K /
                     loaded_v4.continuum.base.division_timing.base_cycle_hours,
                 1.0e-5));
    assert(close(dilute_r_growth / dilute_K_growth, mean_r / mean_K, 1.0e-5));

    // A thin-layer small cell fails only when all eight possible neighbours
    // are occupied. At phi=0.5 the mean-field success probability is
    // 1-0.5^8, not the legacy single-destination vacancy 1-phi.
    auto half_occupied = empty_fields(voxel_count);
    half_occupied.r_normal[0][at(64, 64)] = 0.5;
    StructuredPdeModel3D neighbor_closure(small_v4);
    neighbor_closure.initialize_from_arrays(std::move(half_occupied));
    const auto half_before = neighbor_closure.diagnostics();
    assert(neighbor_closure.step());
    const auto half_after = neighbor_closure.diagnostics();
    const double success = 1.0 - std::pow(0.5, 8.0);
    const double expected_relative = 0.25 * mean_r /
        loaded_v4.continuum.base.division_timing.base_cycle_hours *
        (2.0 * success - 1.0);
    assert(close(
        half_after.r_total / half_before.r_total - 1.0,
        expected_relative, 1.0e-5));

    const auto consumption_for = [&](bool large, CellType type) {
        auto fields = empty_fields(voxel_count);
        if (large) {
            for (int y = 64; y <= 65; ++y) {
                for (int x = 64; x <= 65; ++x) {
                    (type == CellType::r ? fields.r_normal[1] : fields.K[1])
                        [at(x, y)] = 0.025;
                }
            }
        } else {
            (type == CellType::r ? fields.r_normal[0] : fields.K[0])
                [at(64, 64)] = 0.1;
        }
        StructuredPdeModel3D model(small_v2);
        model.initialize_from_arrays(std::move(fields));
        return model.diagnostics();
    };
    const auto pde_small_r = consumption_for(false, CellType::r);
    const auto pde_large_r = consumption_for(true, CellType::r);
    const auto pde_small_K = consumption_for(false, CellType::K);
    assert(close(pde_small_r.assembled_r_consumption_per_hour,
                 pde_large_r.assembled_r_consumption_per_hour));
    assert(close(pde_small_r.assembled_r_consumption_per_hour /
                     pde_small_K.assembled_K_consumption_per_hour,
                 1.2));

    // An already-active cohort remains active in a dilute field after the
    // density trigger is absent, then returns to the ordinary compartment only
    // when its stored clock expires.
    auto pulse = empty_fields(voxel_count);
    pulse.r_active[0][at(64, 64)] = 0.1;
    pulse.active_remaining_hours[0][at(64, 64)] = 0.5;
    StructuredPdeModel3D clock_model(small);
    clock_model.initialize_from_arrays(std::move(pulse));
    assert(clock_model.diagnostics().r_active_total > 0.099);
    assert(clock_model.step());
    assert(clock_model.diagnostics().r_active_total > 0.09);
    assert(clock_model.step());
    assert(clock_model.diagnostics().r_active_total < 1.0e-7);
    assert(clock_model.diagnostics().r_normal_mass[0] > 0.09);

    // A quantized 70x70 ABM anchor window above 0.9 transfers ordinary r into
    // the active/no-direction-yet bucket on the next density refresh.
    StructuredPdeConfig3D dense_config = small;
    dense_config.continuum.end_time_hours = 0.25;
    auto dense = empty_fields(voxel_count);
    for (int y = 0; y < 128; ++y) {
        for (int x = 0; x < 128; ++x) {
            dense.K[0][at(x, y)] = 0.91;
        }
    }
    dense.r_normal[0][at(64, 64)] = 0.09;
    StructuredPdeModel3D triggered(dense_config);
    triggered.initialize_from_arrays(std::move(dense));
    assert(triggered.activation_density(
        StructuredStage3D::small, at(64, 64)) >= 0.90);
    assert(triggered.step());
    assert(triggered.diagnostics().r_active_total > 0.08);

    // The ABM directional-density filter suppresses motion into a dense
    // inward half-neighbourhood, so a boundary cohort has a positive outward
    // centroid shift rather than merely faster isotropic smoothing.
    auto boundary = empty_fields(voxel_count);
    boundary.r_active[0][at(80, 64)] = 0.1;
    boundary.active_remaining_hours[0][at(80, 64)] = 2.0;
    for (int y = 59; y <= 69; ++y) {
        for (int x = 75; x < 80; ++x) boundary.K[0][at(x, y)] = 0.8;
    }
    StructuredPdeModel3D directed(small);
    directed.initialize_from_arrays(std::move(boundary));
    const auto active_centroid_x = [&](const StructuredPdeModel3D& model) {
        long double mass = 0.0L;
        long double moment = 0.0L;
        for (std::size_t location = 0; location < model.voxel_count(); ++location) {
            const double active =
                model.r_active(StructuredStage3D::small, location);
            mass += active;
            moment += active * model.coordinate(location)[0];
        }
        return static_cast<double>(moment / mass);
    };
    const double before_x = active_centroid_x(directed);
    assert(directed.step());
    assert(active_centroid_x(directed) > before_x + 0.01);

    // With equal density in all directions, a supplied nutrient gradient is
    // sufficient to bias the v2 active-r flux toward the resource-rich side.
    StructuredPdeConfig3D guided_config = small_v2;
    guided_config.continuum.end_time_hours = 0.25;
    guided_config.continuum.nutrient.solver_iterations = 128;
    auto resource_gradient = empty_fields(voxel_count);
    resource_gradient.r_active[0][at(64, 64)] = 0.1;
    resource_gradient.active_remaining_hours[0][at(64, 64)] = 2.0;
    for (int y = 0; y < 128; ++y) {
        resource_gradient.vessel_fraction[at(68, y)] = 1.0;
    }
    StructuredPdeModel3D nutrient_directed(guided_config);
    nutrient_directed.initialize_from_arrays(std::move(resource_gradient));
    const double nutrient_before_x = active_centroid_x(nutrient_directed);
    assert(nutrient_directed.step());
    assert(active_centroid_x(nutrient_directed) > nutrient_before_x + 0.01);

    // V3 evaluates directional density and resource over an exact 70x70
    // boundary. A vessel 30 voxels away is outside the old radius-5 sector but
    // still biases the one-step active-r flux toward positive x.
    StructuredPdeConfig3D wide_guidance_config = small_v3;
    wide_guidance_config.continuum.end_time_hours = 0.25;
    auto distant_resource = empty_fields(voxel_count);
    distant_resource.r_active[0][at(64, 64)] = 0.1;
    distant_resource.active_remaining_hours[0][at(64, 64)] = 2.0;
    for (int y = 30; y < 99; ++y) {
        distant_resource.vessel_fraction[at(94, y)] = 1.0;
    }
    StructuredPdeModel3D wide_guided(wide_guidance_config);
    wide_guided.initialize_from_arrays(std::move(distant_resource));
    const double wide_before_x = active_centroid_x(wide_guided);
    assert(wide_guided.step());
    assert(active_centroid_x(wide_guided) > wide_before_x + 0.01);

    // R-to-K conversion uses the shared quantized 70x70 density rather than
    // the source voxel occupancy. The centre is locally dilute (<0.5), while
    // its 70x70 block window is dense enough (>0.5) to convert daughter mass.
    StructuredPdeConfig3D conversion_config = small_v3;
    conversion_config.continuum.end_time_hours = 0.25;
    auto conversion_fields = empty_fields(voxel_count);
    for (int y = 46; y < 116; ++y) {
        for (int x = 46; x < 116; ++x) {
            if (x >= 61 && x <= 67 && y >= 61 && y <= 67) continue;
            conversion_fields.K[0][at(x, y)] = 0.6;
        }
    }
    conversion_fields.r_normal[0][at(64, 64)] = 0.1;
    auto conversion_control_fields = conversion_fields;
    StructuredPdeModel3D conversion_model(conversion_config);
    conversion_model.initialize_from_arrays(std::move(conversion_fields));
    const double conversion_density = conversion_model.activation_density(
        StructuredStage3D::small, at(64, 64));
    assert(conversion_model.occupied_fraction(at(64, 64)) < 0.5);
    assert(conversion_density >= 0.5 && conversion_density < 0.9);
    StructuredPdeConfig3D conversion_control_config = conversion_config;
    conversion_control_config.continuum.base.r_to_K_conversion.enabled = false;
    conversion_control_config.validate();
    StructuredPdeModel3D conversion_control(conversion_control_config);
    conversion_control.initialize_from_arrays(
        std::move(conversion_control_fields));
    assert(conversion_model.step());
    assert(conversion_control.step());
    assert(conversion_model.diagnostics().K_total >
           conversion_control.diagnostics().K_total + 1.0e-7);

    // A packed small active-r cell can exchange with surrounding small K. The
    // operator moves equal masses in opposite directions and therefore keeps
    // total r+K mass constant when reaction death is made negligible.
    StructuredPdeConfig3D exchange_config = small_v3;
    exchange_config.continuum.end_time_hours = 0.25;
    exchange_config.continuum.base.r_death_delay_hours = 1.0e30;
    exchange_config.continuum.base.K_death_delay_hours = 1.0e30;
    auto packed = empty_fields(voxel_count);
    for (int y = 0; y < 128; ++y) {
        for (int x = 0; x < 128; ++x) packed.K[0][at(x, y)] = 1.0;
    }
    packed.K[0][at(64, 64)] = 0.0;
    packed.r_active[0][at(64, 64)] = 1.0;
    packed.active_remaining_hours[0][at(64, 64)] = 2.0;
    StructuredPdeModel3D exchanged(exchange_config);
    exchanged.initialize_from_arrays(std::move(packed));
    const auto packed_before = exchanged.diagnostics();
    assert(exchanged.step());
    const auto packed_after = exchanged.diagnostics();
    assert(exchanged.K(StructuredStage3D::small, at(64, 64)) > 0.1);
    assert(exchanged.r_active(StructuredStage3D::small, at(64, 64)) < 0.9);
    assert(close(packed_before.r_total + packed_before.K_total,
                 packed_after.r_total + packed_after.K_total, 1.0e-5));

    // V3 clears any initial cell mass on a vessel and all transport operators
    // keep the vessel raster cell-free thereafter.
    auto vessel_block = empty_fields(voxel_count);
    vessel_block.vessel_fraction[at(64, 64)] = 1.0;
    vessel_block.r_normal[0][at(64, 64)] = 0.4;
    vessel_block.K[0][at(64, 64)] = 0.4;
    vessel_block.r_active[0][at(63, 64)] = 0.1;
    vessel_block.active_remaining_hours[0][at(63, 64)] = 2.0;
    StructuredPdeModel3D vessel_excluded(small_v3);
    vessel_excluded.initialize_from_arrays(std::move(vessel_block));
    assert(close(vessel_excluded.occupied_fraction(at(64, 64)), 0.0));
    assert(vessel_excluded.step());
    assert(close(vessel_excluded.occupied_fraction(at(64, 64)), 0.0));

    // Restart preserves directional population buckets and their activation
    // clocks bit-for-bit.
    auto restart_fields = empty_fields(voxel_count);
    restart_fields.r_active[0][at(64, 64)] = 0.05;
    restart_fields.active_remaining_hours[0][at(64, 64)] = 2.0;
    StructuredPdeModel3D uninterrupted(small);
    uninterrupted.initialize_from_arrays(std::move(restart_fields));
    assert(uninterrupted.step());
    const auto checkpoint = std::filesystem::temp_directory_path() /
        "atcg3d_structured_pde_roundtrip.bin";
    std::filesystem::remove(checkpoint);
    uninterrupted.save_checkpoint(checkpoint);
    assert(std::filesystem::file_size(checkpoint) < 1024U * 1024U);
    StructuredPdeConfig3D resume_config = small;
    resume_config.continuum.run_mode = "resume";
    resume_config.continuum.resume_checkpoint = checkpoint;
    resume_config.continuum.end_time_hours = 4.0;
    resume_config.validate();
    assert(resume_config.dynamics_fingerprint() ==
           small.dynamics_fingerprint());
    StructuredPdeModel3D resumed(resume_config);
    resumed.load_checkpoint(checkpoint);
    assert(resumed.state_checksum() == uninterrupted.state_checksum());
    assert(uninterrupted.step());
    assert(resumed.step());
    assert(resumed.state_checksum() == uninterrupted.state_checksum());
    std::filesystem::remove(checkpoint);

    // V5 uses exact fixed-value planar-edge and vessel sources. Nutrient is
    // advanced transiently rather than being replaced by a quasi-steady solve.
    auto v5_sources = empty_fields(voxel_count);
    for (int y = 50; y < 79; ++y) {
        v5_sources.vessel_fraction[at(94, y)] = 1.0;
    }
    StructuredPdeConfig3D source_config = small_v5;
    source_config.continuum.nutrient.initial_value = 0.0;
    source_config.validate();
    StructuredPdeModel3D source_model(source_config);
    source_model.initialize_from_arrays(std::move(v5_sources));
    assert(close(source_model.nutrient()[at(0, 64)], 1.0));
    assert(close(source_model.nutrient()[at(127, 64)], 1.0));
    assert(close(source_model.nutrient()[at(64, 0)], 1.0));
    assert(close(source_model.nutrient()[at(94, 64)], 1.0));
    assert(close(source_model.nutrient()[at(64, 64)], 0.0));
    assert(source_model.step());
    assert(close(source_model.nutrient()[at(0, 64)], 1.0));
    assert(close(source_model.nutrient()[at(94, 64)], 1.0));
    assert(source_model.nutrient()[at(1, 64)] > 0.0);
    assert(source_model.nutrient()[at(93, 64)] > 0.0);

    // Without planar Dirichlet sources, the outer boundary is zero-flux. The
    // vessel-only operator must stay in-domain and preserve bounded nutrient.
    StructuredPdeConfig3D vessel_only_config = source_config;
    vessel_only_config.continuum.nutrient.boundary_mode =
        "vessels_dirichlet_v1";
    vessel_only_config.validate();
    auto vessel_only_fields = empty_fields(voxel_count);
    vessel_only_fields.vessel_fraction[at(94, 64)] = 1.0;
    StructuredPdeModel3D vessel_only_model(vessel_only_config);
    vessel_only_model.initialize_from_arrays(std::move(vessel_only_fields));
    assert(vessel_only_model.step());
    for (const double value : vessel_only_model.nutrient()) {
        assert(std::isfinite(value));
        assert(value >= 0.0 && value <= 1.0);
    }
    assert(close(vessel_only_model.nutrient()[at(94, 64)], 1.0));
    const double damkohler =
        small_v5.continuum.nutrient.K_consumption_rate_per_hour *
        35.0 * 35.0 /
        (small_v5.continuum.nutrient.diffusion_voxels2_per_hour * 1.25);
    assert(damkohler > 1.0);

    // Every biological cell has the same integrated demand in v5, regardless
    // of phenotype or footprint.
    const auto v5_consumption_for = [&](bool large, CellType type) {
        auto fields = empty_fields(voxel_count);
        if (large) {
            for (int y = 64; y <= 65; ++y) {
                for (int x = 64; x <= 65; ++x) {
                    (type == CellType::r ? fields.r_normal[1] : fields.K[1])
                        [at(x, y)] = 0.025;
                }
            }
        } else {
            (type == CellType::r ? fields.r_normal[0] : fields.K[0])
                [at(64, 64)] = 0.1;
        }
        StructuredPdeModel3D model(small_v5);
        model.initialize_from_arrays(std::move(fields));
        return model.diagnostics();
    };
    const auto v5_small_r = v5_consumption_for(false, CellType::r);
    const auto v5_large_r = v5_consumption_for(true, CellType::r);
    const auto v5_small_K = v5_consumption_for(false, CellType::K);
    assert(close(v5_small_r.assembled_r_consumption_per_hour,
                 v5_large_r.assembled_r_consumption_per_hour));
    assert(close(v5_small_r.assembled_r_consumption_per_hour,
                 v5_small_K.assembled_K_consumption_per_hour));

    // V6 clamps the host exterior rather than the computational box edge.
    // A hole inside the connected lesion remains tumour, a sub-threshold
    // satellite does not redefine the front, and nutrient enters uniformly
    // from the host/tumour interface.
    StructuredPdeConfig3D moving_front_config = small_v6;
    moving_front_config.continuum.nutrient.initial_value = 0.0;
    moving_front_config.continuum.end_time_hours = 0.5;
    moving_front_config.continuum.base.r_to_K_conversion.enabled = false;
    moving_front_config.validate();
    auto moving_front_fields = empty_fields(voxel_count);
    for (int y = 42; y <= 86; ++y) {
        for (int x = 42; x <= 86; ++x) {
            const int dx = x - 64;
            const int dy = y - 64;
            if (dx * dx + dy * dy <= 20 * 20) {
                moving_front_fields.K[0][at(x, y)] = 0.5;
            }
        }
    }
    for (int y = 61; y <= 67; ++y) {
        for (int x = 61; x <= 67; ++x) {
            moving_front_fields.K[0][at(x, y)] = 0.0;
        }
    }
    moving_front_fields.K[0][at(110, 110)] = 0.01;
    StructuredPdeModel3D moving_front_model(moving_front_config);
    moving_front_model.initialize_from_arrays(std::move(moving_front_fields));
    assert(moving_front_model.tumour_mask()[at(64, 64)] == 1U);
    assert(moving_front_model.tumour_mask()[at(110, 110)] == 0U);
    assert(close(moving_front_model.nutrient()[at(64, 64)], 0.0));
    assert(close(moving_front_model.nutrient()[at(10, 10)], 1.0));
    std::size_t interface_site = voxel_count;
    for (int y = 1; y < 127 && interface_site == voxel_count; ++y) {
        for (int x = 1; x < 127; ++x) {
            const auto here = at(x, y);
            if (moving_front_model.tumour_mask()[here] == 0U) continue;
            if (moving_front_model.tumour_mask()[at(x - 1, y)] == 0U ||
                moving_front_model.tumour_mask()[at(x + 1, y)] == 0U ||
                moving_front_model.tumour_mask()[at(x, y - 1)] == 0U ||
                moving_front_model.tumour_mask()[at(x, y + 1)] == 0U) {
                interface_site = here;
                break;
            }
        }
    }
    assert(interface_site != voxel_count);
    assert(moving_front_model.step());
    assert(moving_front_model.nutrient()[interface_site] > 0.0);
    assert(close(moving_front_model.nutrient()[at(10, 10)], 1.0));
    assert(moving_front_model.diagnostics().tumour_front_mean_nutrient > 0.0);

    // The embedded biological boundary is invariant to unused host padding.
    // The same centred lesion in 128 and 192 grids produces the same local
    // nutrient update and tumour size.
    const auto domain_invariant_probe = [&](int edge) {
        StructuredPdeConfig3D config = small_v6;
        config.continuum.grid.shape = {edge, edge, 1};
        config.continuum.grid.origin = {
            -0.5 * edge, -0.5 * edge, -0.5};
        config.continuum.end_time_hours = 0.25;
        config.continuum.nutrient.initial_value = 0.0;
        config.continuum.base.r_to_K_conversion.enabled = false;
        config.validate();
        const std::size_t count = static_cast<std::size_t>(edge) * edge;
        auto fields = empty_fields(count);
        const int centre = edge / 2;
        const auto location = [edge](int x, int y) {
            return static_cast<std::size_t>(y) * edge + x;
        };
        for (int y = centre - 20; y <= centre + 20; ++y) {
            for (int x = centre - 20; x <= centre + 20; ++x) {
                const int dx = x - centre;
                const int dy = y - centre;
                if (dx * dx + dy * dy <= 20 * 20) {
                    fields.K[0][location(x, y)] = 0.5;
                }
            }
        }
        StructuredPdeModel3D model(config);
        model.initialize_from_arrays(std::move(fields));
        assert(model.step());
        const auto value = model.diagnostics();
        return std::array<double, 3>{
            model.nutrient()[location(centre, centre)],
            model.nutrient()[location(centre + 21, centre)],
            value.tumour_volume};
    };
    const auto domain_128 = domain_invariant_probe(128);
    const auto domain_192 = domain_invariant_probe(192);
    for (std::size_t field = 0; field < domain_128.size(); ++field) {
        assert(close(domain_128[field], domain_192[field], 1.0e-12));
    }

    // The front is derived deterministically from checkpointed population
    // fields, so no additional checkpoint payload is required.
    const auto v6_checkpoint = std::filesystem::temp_directory_path() /
        "atcg3d_structured_pde_v6_roundtrip.bin";
    std::filesystem::remove(v6_checkpoint);
    moving_front_model.save_checkpoint(v6_checkpoint);
    StructuredPdeConfig3D moving_resume_config = moving_front_config;
    moving_resume_config.continuum.run_mode = "resume";
    moving_resume_config.continuum.resume_checkpoint = v6_checkpoint;
    moving_resume_config.continuum.end_time_hours = 0.75;
    moving_resume_config.validate();
    StructuredPdeModel3D moving_resumed(moving_resume_config);
    moving_resumed.load_checkpoint(v6_checkpoint);
    assert(moving_resumed.state_checksum() ==
           moving_front_model.state_checksum());
    assert(moving_resumed.tumour_mask() == moving_front_model.tumour_mask());
    assert(moving_front_model.step());
    assert(moving_resumed.step());
    assert(moving_resumed.state_checksum() ==
           moving_front_model.state_checksum());
    std::filesystem::remove(v6_checkpoint);

    // A high-density environment no longer removes any direction. With a
    // vessel/resource slab on the right, active r moves right even though all
    // eight directional sectors are above the legacy 0.60 density gate.
    StructuredPdeConfig3D chemotaxis_config = small_v5;
    chemotaxis_config.continuum.end_time_hours = 0.25;
    chemotaxis_config.continuum.nutrient.initial_value = 0.0;
    chemotaxis_config.continuum.nutrient.boundary_mode =
        "vessels_dirichlet_v1";
    chemotaxis_config.validate();
    auto dense_chemotaxis = empty_fields(voxel_count);
    for (int y = 0; y < 128; ++y) {
        for (int x = 0; x < 128; ++x) {
            dense_chemotaxis.K[0][at(x, y)] = 0.8;
        }
        for (int x = 90; x <= 96; ++x) {
            dense_chemotaxis.vessel_fraction[at(x, y)] = 1.0;
        }
    }
    dense_chemotaxis.K[0][at(64, 64)] = 0.8;
    dense_chemotaxis.r_active[0][at(64, 64)] = 0.1;
    dense_chemotaxis.active_remaining_hours[0][at(64, 64)] = 2.0;
    StructuredPdeModel3D chemotaxis_model(chemotaxis_config);
    chemotaxis_model.initialize_from_arrays(std::move(dense_chemotaxis));
    const double chemotaxis_before_x = active_centroid_x(chemotaxis_model);
    assert(chemotaxis_model.activation_density(
        StructuredStage3D::small, at(64, 64)) > 0.60);
    assert(chemotaxis_model.step());
    assert(active_centroid_x(chemotaxis_model) > chemotaxis_before_x + 0.001);

    // Clock expiry is terminal for the current activation episode. Dense
    // surroundings cannot immediately re-trigger the just-expired cohort;
    // the local cooldown and off-threshold must both clear first.
    StructuredPdeConfig3D stop_config = small_v5;
    stop_config.continuum.end_time_hours = 0.5;
    stop_config.continuum.nutrient.boundary_mode =
        "planar_edges_dirichlet_v1";
    stop_config.continuum.base.r_to_K_conversion.enabled = false;
    stop_config.validate();
    auto stopping = empty_fields(voxel_count);
    for (int y = 0; y < 128; ++y) {
        for (int x = 0; x < 128; ++x) stopping.K[0][at(x, y)] = 0.89;
    }
    stopping.K[0][at(64, 64)] = 0.80;
    stopping.r_active[0][at(64, 64)] = 0.1;
    stopping.active_remaining_hours[0][at(64, 64)] = 0.25;
    StructuredPdeModel3D stop_model(stop_config);
    stop_model.initialize_from_arrays(std::move(stopping));
    assert(stop_model.step());
    assert(stop_model.diagnostics().r_active_total < 1.0e-6);
    assert(stop_model.step());
    assert(stop_model.diagnostics().r_active_total < 1.0e-6);

    // The v3 checkpoint payload includes refractory and re-arm fields and is
    // deterministic across a restart.
    auto v5_restart_fields = empty_fields(voxel_count);
    v5_restart_fields.r_active[0][at(64, 64)] = 0.05;
    v5_restart_fields.active_remaining_hours[0][at(64, 64)] = 0.25;
    StructuredPdeModel3D v5_uninterrupted(small_v5);
    v5_uninterrupted.initialize_from_arrays(std::move(v5_restart_fields));
    assert(v5_uninterrupted.step());
    const auto v5_checkpoint = std::filesystem::temp_directory_path() /
        "atcg3d_structured_pde_v5_roundtrip.bin";
    std::filesystem::remove(v5_checkpoint);
    v5_uninterrupted.save_checkpoint(v5_checkpoint);
    StructuredPdeConfig3D v5_resume_config = small_v5;
    v5_resume_config.continuum.run_mode = "resume";
    v5_resume_config.continuum.resume_checkpoint = v5_checkpoint;
    v5_resume_config.continuum.end_time_hours = 2.0;
    v5_resume_config.validate();
    StructuredPdeModel3D v5_resumed(v5_resume_config);
    v5_resumed.load_checkpoint(v5_checkpoint);
    assert(v5_resumed.state_checksum() == v5_uninterrupted.state_checksum());
    assert(v5_uninterrupted.step());
    assert(v5_resumed.step());
    assert(v5_resumed.state_checksum() == v5_uninterrupted.state_checksum());
    std::filesystem::remove(v5_checkpoint);

    assert(uninterrupted.diagnostics().maximum_occupied_fraction <=
           small.continuum.reaction.maximum_occupied_fraction + 5.0e-5);
    return 0;
}
