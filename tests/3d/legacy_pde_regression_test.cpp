#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>

#include "model/continuum_model.hpp"
#include "model/structured_pde_model.hpp"

int main(int argc, char**) {
    using namespace atcg3d;
    using namespace atcg3d::structured_pde;
    const auto root = std::filesystem::path(ATCG_SOURCE_DIR);
    const char* wrappers[] = {
        "ATCG3D_StructuredPDE/config/structured_legacy_2d_2000_r20_v1.yaml",
        "ATCG3D_StructuredPDE/config/structured_legacy_2d_2000_r20_resource_guided_v2.yaml",
        "ATCG3D_StructuredPDE/config/structured_legacy_2d_2000_r20_exchange_v3.yaml",
        "ATCG3D_StructuredPDE/config/structured_legacy_2d_2000_r200_reaction_fixed_v4.yaml",
        "ATCG3D_StructuredPDE_NutrientChemotaxis/config/structured_smoke_2d_256_v5.yaml",
        "ATCG3D_StructuredPDE_NutrientChemotaxis/config/structured_smoke_2d_256_24h_v6.yaml"
    };
    std::ifstream fixture(root / "tests/3d/fixtures/legacy_pde_checksums.txt");
    if (argc == 1) assert(fixture);
    for (const char* wrapper : wrappers) {
        auto config = StructuredPdeConfig3D::load(root / wrapper);
        auto& c = config.continuum;
        c.output.enabled = false;
        c.base.threads = 1;
        c.grid.shape = {32, 32, 1};
        c.grid.origin = {-16.0, -16.0, -0.5};
        c.grid.spacing_voxels = 1.0;
        c.time_step_hours = 0.01;
        c.end_time_hours = 1.0;
        c.vascular.source_mode = "abm_perfusion";
        c.nutrient.solver_iterations = 3;
        StructuredInitialFields3D initial;
        for (std::size_t stage = 0; stage < 2; ++stage) {
            initial.r_normal[stage].assign(1024, 0.0);
            initial.r_active[stage].assign(1024, 0.0);
            initial.K[stage].assign(1024, 0.0);
            initial.active_remaining_hours[stage].assign(1024, 0.0);
        }
        initial.vessel_fraction.assign(1024, 0.0);
        initial.vessel_fraction[0] = 1.0;
        for (int y = 12; y < 20; ++y) {
            for (int x = 12; x < 20; ++x) {
                const auto location = static_cast<std::size_t>(y * 32 + x);
                initial.r_normal[0][location] = 0.2;
                initial.K[0][location] = 0.15;
                initial.r_active[1][location] = 0.05;
                initial.active_remaining_hours[1][location] = 0.035;
            }
        }
        StructuredPdeModel3D structured(config);
        structured.initialize_from_arrays(initial);
        for (int step = 0; step < 8; ++step) assert(structured.step());
        const auto structured_checksum = structured.state_checksum();
        std::array<std::vector<double>, 4> fields;
        fields[0] = initial.r_normal[0];
        fields[1] = initial.r_active[1];
        fields[2] = initial.K[0];
        fields[3] = initial.K[1];
        continuum::ContinuumModel3D continuum(c);
        continuum.initialize_from_arrays(fields, initial.vessel_fraction);
        for (int step = 0; step < 8; ++step) assert(continuum.step());
        const auto continuum_checksum = continuum.state_checksum();
        if (argc > 1) {
            std::cout << structured_checksum << ' ' << continuum_checksum << '\n';
        } else {
            std::uint64_t expected_structured{}, expected_continuum{};
            fixture >> expected_structured >> expected_continuum;
            assert(fixture);
            assert(structured_checksum == expected_structured);
            assert(continuum_checksum == expected_continuum);
        }
    }
}
