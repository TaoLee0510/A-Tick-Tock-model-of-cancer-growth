#include "engine/simulation.hpp"
#include "model/structured_pde_model.hpp"
#include <cassert>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <sys/resource.h>

int main(int argc, char **argv) {
    using namespace atcg3d;
    using namespace atcg3d::structured_pde;
    auto c = StructuredPdeConfig3D::load(
        std::filesystem::path(ATCG_SOURCE_DIR) /
        "ATCG3D_SharedRules/config/angiogenesis_v8.yaml");
    c.schema_version = 9;
    c.storage_model = "sparse_zero_pages_v1";
    c.continuum.angiogenesis.model = "disabled";
    c.continuum.base.angiogenesis.enabled = false;
    c.continuum.vascular.source_mode = "static_voxels";
    c.continuum.vascular.static_sources.clear();
    c.continuum.nutrient.initial_value = 1;
    c.continuum.nutrient.boundary_mode = "moving_tumor_front_dirichlet_v2";
    c.continuum.end_time_hours = 1;
    c.continuum.base.initial_r_cells = 16;
    c.continuum.base.initial_K_cells = 16;
    if (argc > 1 && std::string(argv[1]) == "--large") {
        c = StructuredPdeConfig3D::load(
            std::filesystem::path(ATCG_SOURCE_DIR) /
            "ATCG3D_StructuredPDE/config/structured_sparse_v9.yaml");
        Simulation3D source(c.continuum.base);
        source.initialize();
        StructuredPdeModel3D large(c);
        large.initialize_from_abm(source);
        while (large.step()) {
        }
        const auto diagnostics = large.diagnostics();
        const auto checksum = large.state_checksum();
        rusage usage{};
        assert(getrusage(RUSAGE_SELF, &usage) == 0);
#ifdef __APPLE__
        const std::uint64_t bytes = usage.ru_maxrss;
#else
        const std::uint64_t bytes = std::uint64_t(usage.ru_maxrss) * 1024;
#endif
        std::cout << "{\"grid\":10000,\"peak_resident_bytes\":" << bytes
                  << ",\"mass\":"
                  << diagnostics.r_total + diagnostics.K_total
                  << ",\"state_checksum\":" << checksum << "}\n";
        assert(bytes < 8000000000ULL);
        return 0;
    }
    Simulation3D source(c.continuum.base);
    source.initialize();
    StructuredPdeModel3D sparse(c);
    auto dense_config = c;
    dense_config.storage_model = "dense_v1";
    StructuredPdeModel3D dense(dense_config);
    StructuredInitialFields3D initial;
    for (int stage = 0; stage < 2; ++stage) {
        initial.r_normal[stage].resize(2304);
        initial.r_active[stage].resize(2304);
        initial.K[stage].resize(2304);
        initial.active_remaining_hours[stage].resize(2304);
        initial.r_refractory[stage].resize(2304);
        initial.refractory_remaining_hours[stage].resize(2304);
    }
    initial.vessel_fraction.resize(2304);
    for (int y = 20; y < 26; ++y)
        for (int x = 20; x < 26; ++x) {
            const auto i = std::size_t(y * 48 + x);
            initial.r_normal[0][i] = 0.3;
            initial.r_active[0][i] = 0.2;
            initial.K[0][i] = 0.2;
            initial.active_remaining_hours[0][i] = 6;
            initial.r_refractory[0][i] = 0.1;
            initial.refractory_remaining_hours[0][i] = 8;
            initial.r_normal[1][i] = 0.02;
            initial.r_active[1][i] = 0.01;
            initial.K[1][i] = 0.02;
            initial.active_remaining_hours[1][i] = 12;
            initial.r_refractory[1][i] = 0.01;
            initial.refractory_remaining_hours[1][i] = 4;
        }
    dense.initialize_from_arrays(initial);
    sparse.initialize_from_arrays(std::move(initial));
    for (int i = 0; i < 2; ++i) {
        assert(sparse.step());
        assert(dense.step());
    }
    for (std::size_t i = 0; i < sparse.voxel_count(); ++i)
        for (auto stage :
             {StructuredStage3D::small, StructuredStage3D::large}) {
            assert(sparse.r_normal(stage, i) == dense.r_normal(stage, i));
            assert(sparse.r_active(stage, i) == dense.r_active(stage, i));
            assert(sparse.K(stage, i) == dense.K(stage, i));
            assert(sparse.refractory_mass(stage, i) ==
                   dense.refractory_mass(stage, i));
        }
    assert(sparse.diagnostics().r_active_total > 0);
    assert(sparse.nutrient() == dense.nutrient());
    assert(sparse.vessel_fraction() == dense.vessel_fraction());
    const auto sd = sparse.diagnostics();
    const auto dd = dense.diagnostics();
    assert(sd.r_total == dd.r_total && sd.K_total == dd.K_total);
    assert(sd.r_active_total == dd.r_active_total);
    assert(sd.mean_nutrient == dd.mean_nutrient);
    assert(sd.vessel_volume == dd.vessel_volume);
    assert(sd.tumour_volume == dd.tumour_volume);
    assert(sd.tumour_front_volume == dd.tumour_front_volume);
    assert(sd.r_radius_50 == dd.r_radius_50);
    assert(sd.r_radius_90 == dd.r_radius_90);
    assert(sd.r_radius_99 == dd.r_radius_99);
    const auto path = std::filesystem::current_path() / "sparse-pde-test.bin";
    std::filesystem::remove(path);
    sparse.save_checkpoint(path);
    auto parallel = c;
    parallel.continuum.base.threads = 4;
    StructuredPdeModel3D resumed(parallel);
    resumed.load_checkpoint(path);
    assert(resumed.state_checksum() == sparse.state_checksum());
    while (sparse.step()) {
    }
    while (resumed.step()) {
    }
    assert(resumed.state_checksum() == sparse.state_checksum());
    std::filesystem::remove(path);
    auto budget = c;
    budget.maximum_active_voxels = 1;
    StructuredPdeModel3D limited(budget);
    limited.initialize_from_abm(source);
    bool rejected = false;
    try {
        limited.step();
    } catch (const std::runtime_error &) {
        rejected = true;
    }
    assert(rejected);
    // Mapping resets must remove existing dirty values rather than retaining
    // discarded pages, which differs across operating-system advice APIs.
    PagedField<double> field;
    field.set_sparse(true);
    field.assign(1000000, 0);
    field[5000] = 7;
    field.fill(0);
    assert(field[5000] == 0);
    field[6000] = 3;
    PagedField<double> moved(std::move(field));
    assert(moved[6000] == 3);
}
