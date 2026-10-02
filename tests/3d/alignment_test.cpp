#include <algorithm>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>

#include "config/nutrient_config.hpp"
#include "engine/simulation.hpp"
#include "field/nutrient_field.hpp"
#include "model/continuum_model.hpp"
#include "model/structured_pde_model.hpp"
#include "rules/density.hpp"
#include "rules/migration.hpp"
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
#include "io/checkpoint_hdf5.hpp"
#endif

namespace {
using namespace atcg3d;
using namespace atcg3d::structured_pde;

bool close(double a, double b, double tolerance = 1.0e-9) {
    return std::abs(a - b) <= tolerance * std::max({1.0, std::abs(a), std::abs(b)});
}

StructuredInitialFields3D empty(std::size_t size) {
    StructuredInitialFields3D fields;
    for (std::size_t stage = 0; stage < 2; ++stage) {
        fields.r_normal[stage].assign(size, 0.0);
        fields.r_active[stage].assign(size, 0.0);
        fields.K[stage].assign(size, 0.0);
        fields.active_remaining_hours[stage].assign(size, 0.0);
        fields.r_refractory[stage].assign(size, 0.0);
        fields.refractory_remaining_hours[stage].assign(size, 0.0);
    }
    fields.vessel_fraction.assign(size, 0.0);
    return fields;
}

template <typename Action> bool rejects(Action action) {
    try { action(); } catch (const std::exception&) { return true; }
    return false;
}

struct BoundedResource final : LocalDensityModifier3D {
    double retained_density(Vec3i) const noexcept override { return 1.0; }
    double normalized_resource(Vec3i site) const noexcept override {
        return contains_resource_site(site) ? 1.0 : 0.0;
    }
    bool contains_resource_site(Vec3i site) const noexcept override {
        return site.x >= 0 && site.x <= 4 && site.y >= 0 && site.y <= 4 && site.z == 0;
    }
};
}  // namespace

int main() {
    using namespace atcg3d::continuum;
    const auto root = std::filesystem::path(ATCG_SOURCE_DIR);
    auto config = StructuredPdeConfig3D::load(root /
        "ATCG3D_StructuredPDE_NutrientChemotaxis/config/structured_smoke_2d_256_v7.yaml");
    assert(config.schema_version == 7 && config.continuum.schema_version == 5);
    config.continuum.output.enabled = false;
    config.continuum.grid.shape = {128, 128, 1};
    config.continuum.grid.origin = {-64.0, -64.0, -0.5};
    config.continuum.base.output_enabled = false;
    config.continuum.base.control_enabled = false;
    config.continuum.base.static_vasculature = config.continuum.shared_vascular_geometry();

    // Both finite-field imports use the same exclusion mask before creating
    // cells. No cells are lost when the ABM state is mapped into the PDE.
    Simulation3D abm(config.continuum.base);
    abm.initialize();
    StructuredPdeModel3D imported(config);
    imported.initialize_from_abm(abm);
    ContinuumModel3D imported_continuum(config.continuum);
    imported_continuum.initialize_from_abm(abm);
    assert(close(imported.diagnostics().r_total + imported.diagnostics().K_total,
        static_cast<double>(abm.cells().alive_slots().size())));
    const auto masses = imported_continuum.diagnostics().type_mass;
    assert(close(masses[0] + masses[1], static_cast<double>(abm.cells().alive_slots().size())));
    assert(imported.vessel_fraction() == imported_continuum.vessel_fraction());
    for (const Slot slot : abm.cells().alive_slots()) {
        assert(!config.continuum.base.static_vasculature.source(abm.cells().anchor(slot)));
    }
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
    const auto abm_checkpoint = std::filesystem::temp_directory_path() / "atcg_static_geometry.h5";
    std::filesystem::remove(abm_checkpoint);
    write_hdf5_checkpoint(abm_checkpoint, abm);
    const auto saved = read_hdf5_checkpoint(abm_checkpoint, config.continuum.base);
    Simulation3D restored_abm(config.continuum.base);
    restored_abm.restore(saved.cells, saved.next_uid, saved.clock, saved.stats,
        saved.lineage, saved.vasculature, saved.cell_slot_count, saved.cell_slots, saved.cell_free_slots);
    assert(abm.state_checksum() == restored_abm.state_checksum());
    for (int event = 0; event < 12; ++event) {
        assert(abm.step() && restored_abm.step());
        assert(abm.state_checksum() == restored_abm.state_checksum());
    }
    std::filesystem::remove(abm_checkpoint);
#endif

    // Even windows contain exactly edge sites; odd windows retain symmetry.
    // Test both endpoints and clipped domain corners against the ABM index.
    config.continuum.vascular.source_mode = "abm_perfusion";
    config.continuum.grid.shape = {80, 80, 1};
    config.continuum.grid.origin = {-40.0, -40.0, -0.5};
    const auto at = [](int x, int y) { return static_cast<std::size_t>((y + 40) * 80 + x + 40); };
    auto initial = empty(6400);
    BlockDensityIndex3D density(4);
    Slot next_slot = 0;
    for (const Vec3i site : {Vec3i{-35, 0, 0}, {-34, 0, 0}, {35, 0, 0},
                            {36, 0, 0}, {0, -34, 0}, {0, 35, 0}, {-40, -40, 0}}) {
        const CellType type = next_slot % 2 == 0 ? CellType::r : CellType::K;
        density.add(site, type, next_slot++);
        (type == CellType::r ? initial.r_normal[0] : initial.K[0])[at(site.x, site.y)] = 1.0;
    }
    for (const int edge : {70, 69, 1}) {
        config.continuum.base.growth_density_window_edge = edge;
        StructuredPdeModel3D structured(config);
        structured.initialize_from_arrays(initial);
        std::array<std::vector<double>, 4> fields{
            initial.r_normal[0], initial.r_normal[1], initial.K[0], initial.K[1]};
        ContinuumModel3D continuum(config.continuum);
        continuum.initialize_from_arrays(fields, initial.vessel_fraction);
        for (const Vec3i site : {Vec3i{0, 0, 0}, {-40, -40, 0}, {35, 0, 0}}) {
            const auto reference = growth_counts(density, site, edge, true);
            const auto sc = structured.growth_counts_at(site);
            const auto cc = continuum.growth_counts_at(site);
            assert(sc[0] == reference.rc && sc[1] == reference.kc);
            if (cc != sc) std::cerr << "edge=" << edge << " site=" << site.x << "," << site.y << " sc=" << sc[0] << "," << sc[1] << " cc=" << cc[0] << "," << cc[1] << std::endl;
            assert(cc == sc);
        }
    }

    // r200 runs independently with deterministic explicit migration substeps.
    // Negligible biological rates isolate conservation in the transport step.
    config.continuum.base.division_timing.base_cycle_hours = 1.0e30;
    config.continuum.base.migration_activation_threshold = 0.001;
    config.migration.reactivation_density_threshold = 0.0005;
    config.continuum.base.migration_activation_window_edge = 1;
    config.continuum.base.migration_activation_block_edge = 1;
    config.continuum.vascular.source_mode = "static_voxels";
    config.continuum.vascular.static_sources = {{1, 0, 0}};
    initial = empty(6400);
    initial.r_normal[0][at(0, 0)] = 0.5;
    std::array<std::vector<double>, 4> fields{
        initial.r_normal[0], initial.r_normal[1], initial.K[0], initial.K[1]};
    ContinuumModel3D fast(config.continuum);
    fast.initialize_from_arrays(fields, initial.vessel_fraction);
    assert(fast.migration_substeps() > 1);
    for (int step = 0; step < 8; ++step) assert(fast.step());
    assert(close(fast.diagnostics().type_mass[0], 0.5));
    assert(fast.occupied_fraction(at(1, 0)) == 0.0);
    const auto checkpoint = std::filesystem::temp_directory_path() / "atcg_alignment.continuum.bin";
    std::filesystem::remove(checkpoint);
    fast.save_checkpoint(checkpoint);
    auto parallel_continuum = config.continuum;
    parallel_continuum.base.threads = 4;
    ContinuumModel3D restarted(parallel_continuum);
    restarted.load_checkpoint(checkpoint);
    for (int step = 0; step < 4; ++step) {
        assert(fast.step() && restarted.step());
        assert(fast.state_checksum() == restarted.state_checksum());
    }
    std::filesystem::remove(checkpoint);
    auto invalid = config;
    invalid.continuum.time_step_hours = 4.0;
    assert(rejects([&] { StructuredPdeModel3D model(invalid); }));
    invalid = config;
    invalid.continuum.base.direction_guidance_model = "density_gate_uniform_v1";
    assert(rejects([&] { invalid.validate(); }));

    // A refractory subset does not suppress an unrelated eligible cohort.
    config.continuum.vascular.source_mode = "abm_perfusion";
    config.continuum.vascular.static_sources.clear();
    config.continuum.migration.diffusion_scale = 1.0e-12;
    initial.r_normal[0][at(0, 0)] = 0.6;
    initial.r_refractory[0][at(0, 0)] = 0.4;
    initial.refractory_remaining_hours[0][at(0, 0)] = 4.0;
    StructuredPdeModel3D mixed(config);
    mixed.initialize_from_arrays(initial);
    assert(mixed.step());
    assert(close(mixed.diagnostics().r_active_total, 0.2, 1.0e-6));
    assert(close(mixed.refractory_mass(StructuredStage3D::small, at(0, 0)), 0.4));

    // Cooldown clock mass moves with ordinary r and mixes by mass weighting.
    config.continuum.migration.diffusion_scale = 1.0;
    config.continuum.base.migration_activation_threshold = 0.9;
    config.migration.reactivation_density_threshold = 0.8;
    config.continuum.time_step_hours = 0.1;
    initial = empty(6400);
    for (const int x : {-1, 1}) {
        initial.r_normal[0][at(x, 0)] = x < 0 ? 0.2 : 0.4;
        initial.r_refractory[0][at(x, 0)] = initial.r_normal[0][at(x, 0)];
        initial.refractory_remaining_hours[0][at(x, 0)] = x < 0 ? 2.0 : 6.0;
    }
    StructuredPdeModel3D transported(config);
    transported.initialize_from_arrays(initial);
    assert(transported.step());
    const auto stage = StructuredStage3D::small;
    double mass = 0.0, clock_mass = 0.0;
    for (std::size_t i = 0; i < 6400; ++i) {
        const double refractory = transported.refractory_mass(stage, i);
        mass += refractory;
        clock_mass += refractory * transported.refractory_mean_hours(stage, i);
    }
    assert(close(mass, 0.6));
    assert(close(clock_mass, 0.2 * 1.9 + 0.4 * 5.9));
    assert(transported.refractory_mass(stage, at(0, 0)) > 0.0);
    const double mean = transported.refractory_mean_hours(stage, at(0, 0));
    assert(mean > 1.9 && mean < 5.9);
    const auto structured_checkpoint = std::filesystem::temp_directory_path() / "atcg_alignment.structured.bin";
    std::filesystem::remove(structured_checkpoint);

    // Mortality rounds individual direction buckets. Their cached total must
    // be rebuilt in the same order as restart, rather than scaled separately.
    auto mortality_config = config;
    mortality_config.continuum.base.division_timing.base_cycle_hours = 24.0;
    mortality_config.continuum.base.growth_density_window_edge = 70;
    initial = empty(6400);
    for (int y = -4; y < 4; ++y) {
        for (int x = -4; x < 4; ++x) {
            initial.r_active[0][at(x, y)] = 0.3;
            initial.active_remaining_hours[0][at(x, y)] = 4.0;
            initial.K[0][at(x, y)] = 0.5;
        }
    }
    StructuredPdeModel3D mortality(mortality_config);
    mortality.initialize_from_arrays(initial);
    for (int step = 0; step < 4; ++step) assert(mortality.step());
    mortality.save_checkpoint(structured_checkpoint);
    StructuredPdeModel3D mortality_resumed(mortality_config);
    mortality_resumed.load_checkpoint(structured_checkpoint);
    for (int step = 0; step < 4; ++step) {
        assert(mortality.step() && mortality_resumed.step());
        assert(mortality.state_checksum() == mortality_resumed.state_checksum());
    }
    std::filesystem::remove(structured_checkpoint);
    transported.save_checkpoint(structured_checkpoint);
    config.continuum.base.threads = 4;
    StructuredPdeModel3D resumed(config);
    resumed.load_checkpoint(structured_checkpoint);
    for (int step = 0; step < 5; ++step) {
        assert(transported.step() && resumed.step());
        assert(transported.state_checksum() == resumed.state_checksum());
    }
    std::filesystem::remove(structured_checkpoint);

    // A zero clock releases the subset after local density falls below off.
    initial = empty(6400);
    initial.r_normal[0][at(0, 0)] = 0.2;
    initial.r_refractory[0][at(0, 0)] = 0.2;
    initial.refractory_remaining_hours[0][at(0, 0)] = 0.01;
    StructuredPdeModel3D released(config);
    released.initialize_from_arrays(initial);
    assert(released.step());
    assert(released.refractory_mass(stage, at(0, 0)) == 0.0);

    // Expired active mass starts its own cooldown wherever transport put it.
    initial = empty(6400);
    initial.r_active[0][at(0, 0)] = 0.1;
    initial.active_remaining_hours[0][at(0, 0)] = 0.01;
    StructuredPdeModel3D expired(config);
    expired.initialize_from_arrays(initial);
    assert(expired.step());
    mass = 0.0;
    clock_mass = 0.0;
    for (std::size_t i = 0; i < 6400; ++i) {
        const double refractory = expired.refractory_mass(stage, i);
        mass += refractory;
        clock_mass += refractory * expired.refractory_mean_hours(stage, i);
    }
    assert(close(mass, 0.1, 1.0e-6));
    assert(close(clock_mass, mass * config.migration.reactivation_cooldown_hours));
    assert(expired.diagnostics().r_active_total == 0.0);
    initial = empty(6400);
    initial.r_normal[0][at(0, 0)] = 0.9;
    initial.r_refractory[0][at(0, 0)] = 0.9;
    initial.refractory_remaining_hours[0][at(0, 0)] = 0.01;
    StructuredPdeModel3D hysteresis(config);
    hysteresis.initialize_from_arrays(initial);
    assert(hysteresis.step());
    assert(hysteresis.refractory_mass(stage, at(0, 0)) > 0.0);
    assert(hysteresis.refractory_mean_hours(stage, at(0, 0)) == 0.0);
    assert(hysteresis.diagnostics().r_active_total == 0.0);

    // The static nutrient mask uses the same voxel centre convention and
    // clamps both synthetic and explicitly supplied vessels to source value.
    auto nutrient_config = nutrient::NutrientModelConfig3D::load(root /
        "ATCG3D_Nutrient/config/nutrient_smoke_static_guided_v3.yaml");
    nutrient::NutrientEnvironment3D resource(nutrient_config.nutrient, true);
    resource.initialize(0.0, abm.cells(), abm.vessel_grid());
    assert(resource.value({-40, 0, 0}) == 1.0F);
    assert(resource.value({40, 0, 0}) == 1.0F);
    assert(!resource.contains_resource_site({256, 0, 0}));
    assert(resource.value({256, 0, 0}) == 0.0F);
    auto vascular_config = config.continuum;
    vascular_config.grid.shape = nutrient_config.nutrient.static_vasculature.shape;
    vascular_config.grid.origin = nutrient_config.nutrient.static_vasculature.origin;
    vascular_config.vascular.source_mode = "synthetic_central_line";
    vascular_config.vascular.static_sources = {{40, 0, 0}};
    const auto mask_initial = empty(256U * 256U);
    const std::array<std::vector<double>, 4> mask_fields{
        mask_initial.r_normal[0], mask_initial.r_normal[1], mask_initial.K[0], mask_initial.K[1]};
    ContinuumModel3D mask_model(vascular_config);
    mask_model.initialize_from_arrays(mask_fields, mask_initial.vessel_fraction);
    for (int y = -128; y < 128; ++y) {
        for (int x = -128; x < 128; ++x) {
            const auto location = static_cast<std::size_t>((y + 128) * 256 + x + 128);
            assert((mask_model.vessel_fraction()[location] > 0.0) ==
                nutrient_config.nutrient.static_vasculature.source({x, y, 0}));
        }
    }
    const auto bad_yaml = std::filesystem::temp_directory_path() / "atcg_bad_halo.yaml";
    std::ifstream source(root / "ATCG3D_Nutrient/config/nutrient_smoke_static_guided_v3.yaml");
    std::string text{std::istreambuf_iterator<char>(source), {}};
    text.replace(text.find("halo_voxels: 35"), std::string("halo_voxels: 35").size(), "halo_voxels: 1");
    const auto base_pos = text.find("../../configs/atcg2d_static_bounded_r20_v4.yaml");
    text.replace(base_pos, std::string("../../configs/atcg2d_static_bounded_r20_v4.yaml").size(), (root / "configs/atcg2d_static_bounded_r20_v4.yaml").string());
    std::ofstream(bad_yaml) << text;
    assert(rejects([&] { nutrient::NutrientModelConfig3D::load(bad_yaml); }));
    std::filesystem::remove(bad_yaml);

    // A persistent bounded guide ignores density only when explicitly told
    // to do so. Resource averaging excludes samples outside the finite grid.
    auto guide = config.continuum.base;
    guide.direction_guidance_model = "low_density_high_resource_bounded_v2";
    guide.static_vasculature = {};
    guide.direction_density_radius = 1;
    guide.direction_density_threshold = 0.01;
    guide.continue_probability = 1.0;
    CellStore3D cells;
    CellInit cell;
    cell.uid = 1;
    cell.anchor = {0, 0, 0};
    cell.type = CellType::r;
    cell.stage = CellStage::small;
    cell.flags = static_cast<std::uint8_t>(kMigrationActive);
    const Slot slot = cells.create(cell);
    DirectionId forward = kStayDirection;
    for (DirectionId direction = 1; direction <= 26; ++direction) {
        if (direction_vector(direction) == Vec3i{1, 0, 0}) forward = direction;
    }
    cells.set_last_direction(slot, forward);
    SparseChunkGrid3D grid(guide.chunk_edge, DomainPolicy(guide));
    assert(grid.place_single(cell.anchor, slot));
    BlockDensityIndex3D crowded(1);
    crowded.add({1, 0, 0}, CellType::K, 1);
    BoundedResource bounded;
    guide.direction_density_radius = 2;
    assert(migration_direction_resource({0, 0, 0}, forward, guide, &bounded) == 1.0);
    guide.direction_guidance_model = "low_density_high_resource_v1";
    assert(migration_direction_resource({0, 0, 0}, forward, guide, &bounded) < 1.0);
    guide.direction_guidance_model = "low_density_high_resource_bounded_v2";
    guide.direction_density_radius = 1;
    guide.persistence_uses_density = false;
    assert(select_migration_direction(slot, cells, grid, crowded, guide, 0, &bounded) == forward);
    guide.persistence_uses_density = true;
    assert(select_migration_direction(slot, cells, grid, crowded, guide, 0, &bounded) != forward);

    // The same edge contract applies in all three spatial dimensions.
    auto spatial_config = config.continuum;
    spatial_config.base.thin_layer = false;
    spatial_config.nutrient.boundary_mode = "planar_edges_and_vessels_dirichlet_v1";
    spatial_config.grid.shape = {8, 8, 8};
    spatial_config.grid.origin = {-4.0, -4.0, -4.0};
    spatial_config.vascular.static_sources.clear();
    auto spatial = empty(512);
    BlockDensityIndex3D spatial_density(1);
    next_slot = 0;
    for (const Vec3i site : {Vec3i{-2, 0, 0}, {-1, 0, 0}, {2, 0, 0},
                            {0, -1, 0}, {0, 2, 0}, {0, 0, -1}, {0, 0, 2}}) {
        const auto location = static_cast<std::size_t>(((site.z + 4) * 8 + site.y + 4) * 8 + site.x + 4);
        spatial.r_normal[0][location] = 1.0;
        spatial_density.add(site, CellType::r, next_slot++);
    }
    spatial_config.base.growth_density_window_edge = 4;
    ContinuumModel3D spatial_model(spatial_config);
    spatial_model.initialize_from_arrays({spatial.r_normal[0], spatial.r_normal[1], spatial.K[0], spatial.K[1]},
        spatial.vessel_fraction);
    assert(spatial_model.growth_counts_at({0, 0, 0})[0] == growth_counts(spatial_density, {0, 0, 0}, 4).rc);
}
