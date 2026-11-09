#include <cassert>
#include <cmath>
#include <filesystem>

#include "model/hybrid_model.hpp"

namespace {
using namespace atcg3d;
using namespace atcg3d::hybrid;

HybridConfig3D config(bool thin) {
    auto c = HybridConfig3D::load(std::filesystem::path(ATCG_SOURCE_DIR) /
        "ATCG3D_Hybrid/config/hybrid_invasion_r20_v4.yaml");
    c.rules.continuum.grid.shape = {32, 32, thin ? 1 : 32};
    c.rules.continuum.grid.origin = {-16, -16, thin ? -0.5 : -16};
    c.rules.continuum.base.thin_layer = thin;
    c.rules.continuum.reaction.small_daughter_vacancy_exponent = thin ? 8 : 26;
    c.rules.continuum.nutrient.boundary_mode = "planar_edges_and_vessels_dirichlet_v1";
    c.rules.continuum.vascular.source_mode = "static_voxels";
    c.rules.continuum.vascular.static_sources.clear();
    c.rules.continuum.end_time_hours = 2;
    c.rules.continuum.nutrient.K_consumption_rate_per_hour = 0;
    c.rules.continuum.nutrient.r_consumption_rate_per_hour = 0;
    c.rules.continuum.nutrient.decay_per_hour = 0;
    c.rules.continuum.base.initial_r_cells = 32;
    c.rules.continuum.base.initial_K_cells = 32;
    c.rules.continuum.base.initial_shell_inner_radius = 6;
    c.rules.continuum.base.initial_radius = 9;
    c.rules.continuum.base.initial_shell_thickness = 3;
    c.rules.continuum.base.initial_inner_small_radius = 5;
    c.front_band_voxels = 1;
    c.front_hysteresis_voxels = 1;
    c.active_guard_voxels = 1;
    c.smoothing_radius = 1;
    c.rules.continuum.nutrient.tumor_front_smoothing_radius_voxels = 1;
    c.rules.continuum.nutrient.tumor_front_density_threshold = 0.02;
    return c;
}

structured_pde::StructuredInitialFields3D fields(std::size_t size) {
    structured_pde::StructuredInitialFields3D f;
    for (int stage = 0; stage < 2; ++stage) {
        f.r_normal[stage].resize(size);
        f.r_active[stage].resize(size);
        f.K[stage].resize(size);
        f.active_remaining_hours[stage].resize(size);
    }
    f.vessel_fraction.resize(size);
    return f;
}

std::size_t index(int x, int y, int z, bool thin) {
    return (std::size_t(thin ? 0 : z + 16) * 32 + y + 16) * 32 + x + 16;
}

void check_conversion(bool thin) {
    auto c = config(thin);
    c.rules.continuum.base.migration_activation_threshold = 1;
    c.rules.continuum.nutrient.common_density_limit = 1000000.0;
    c.rules.continuum.nutrient.common_carrying_capacity = 1000000.0;
    c.core_on = 0.8;
    c.core_off = 0.6;
    auto f = fields(thin ? 1024 : 32768);
    for (int z = thin ? 0 : -7; z <= (thin ? 0 : 7); ++z)
        for (int y = -7; y <= 7; ++y)
            for (int x = -7; x <= 7; ++x)
                f.K[0][index(x, y, z, thin)] = 1.0;
    const double expected = thin ? 225 : 3375;
    HybridModel3D model(c);
    model.initialize_from_arrays(f);
    model.exchange();
    assert(model.is_core({0, 0, 0}));
    assert(!model.is_core({7, 0, 0}));
    assert(model.diagnostics().abm_mass > 0);
    assert(model.diagnostics().pde_mass > 0);
    assert(std::abs(model.diagnostics().total_mass - expected) < 1e-9);
    const auto conversions = model.diagnostics().to_abm;
    for (int repeat = 0; repeat < 3; ++repeat) {
        model.exchange();
        assert(std::abs(model.diagnostics().total_mass - expected) < 1e-9);
    }
    assert(model.diagnostics().to_abm == conversions);
    std::vector<double> front_mass(model.pde().voxel_count());
    for (std::size_t i = 0; i < front_mass.size(); ++i) {
        front_mass[i] = model.pde().K(structured_pde::StructuredStage3D::small, i);
    }
    // Between exchanges, reflecting core transport must not leave fractional
    // density on front sites reserved for whole individuals.
    assert(model.step());
    for (int z = thin ? 0 : -16; z < (thin ? 1 : 16); ++z) {
        for (int y = -16; y < 16; ++y) {
            for (int x = -16; x < 16; ++x) {
                const auto i = index(x, y, z, thin);
                if (!model.is_core({x, y, z})) {
                    assert(model.pde().K(structured_pde::StructuredStage3D::small, i)
                        <= front_mass[i] + 1e-14);
                }
            }
        }
    }
    assert(std::abs(model.diagnostics().total_mass - expected) < 1e-9);
    auto invalid = fields(thin ? 1024 : 32768);
    invalid.r_normal[0][index(0, 0, 0, thin)] = 1.0;
    bool rejected = false;
    try {
        HybridModel3D forbidden(c);
        forbidden.initialize_from_arrays(std::move(invalid));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    assert(rejected);
}

void clean(const std::filesystem::path& path) {
    for (auto suffix : {"", ".pde.bin", ".resource.bin"})
        std::filesystem::remove(path.string() + suffix);
}

void check_restart(bool thin, bool vascular) {
    auto c = config(thin);
    c.rules.continuum.base.migration_activation_threshold = 0.001;
    c.rules.continuum.base.migration_activation_window_edge = 4;
    c.rules.continuum.base.migration_activation_block_edge = 2;
    c.rules.migration.reactivation_density_threshold = 0.0001;
    if (vascular) {
        c.rules.schema_version = 15;
        auto& law = c.rules.continuum.angiogenesis;
        law.model = "shared_vegf_lattice_v2";
        law.seed_tips_per_hour = 0.1;
        c.rules.continuum.nutrient.initial_value = 0.05;
    }
    HybridModel3D original(c);
    original.initialize();
    for (auto slot : original.abm().cells().alive_slots())
        assert(!original.abm().grid().blocked_by_vessel(original.abm().cells().anchor(slot)));
    assert(original.diagnostics().active_mass > 0.0);
    assert(original.diagnostics().abm_mass >= c.rules.continuum.base.initial_r_cells);
    assert(original.step());
    const auto path = std::filesystem::current_path() / "hybrid-invasion.bin";
    clean(path);
    original.save_checkpoint(path);
    c.rules.continuum.base.threads = 4;
    HybridModel3D resumed(c);
    resumed.load_checkpoint(path);
    assert(original.state_checksum() == resumed.state_checksum());
    while (original.step()) {
        assert(resumed.step());
        assert(original.state_checksum() == resumed.state_checksum());
        assert(original.diagnostics().active_mass == resumed.diagnostics().active_mass);
        assert(original.pde().diagnostics().r_total == 0.0);
        for (auto slot : original.abm().cells().alive_slots())
            if (original.abm().cells().flags(slot) & kMigrationActive)
                assert(original.abm().cells().type(slot) == CellType::r);
    }
    clean(path);
}

void check_prefix(bool thin) {
    auto c = config(thin);
    c.rules.continuum.base.growth_density_window_edge = 6;
    c.rules.continuum.base.migration_activation_window_edge = 4;
    c.rules.continuum.base.migration_activation_block_edge = 2;
    auto f = fields(thin ? 1024 : 32768);
    for (int z = thin ? 0 : -2; z <= (thin ? 0 : 2); ++z)
        for (int y = -2; y <= 2; ++y)
            for (int x = -2; x <= 2; ++x)
                f.K[0][index(x, y, z, thin)] = double(x + 3) / 16.0;
    HybridModel3D model(c);
    model.initialize_from_arrays(f);
    for (auto point : {Vec3i{0, 0, 0}, Vec3i{-15, -15, thin ? 0 : -15},
                       Vec3i{15, 15, thin ? 0 : 15}}) {
        double expected = 0.0;
        for (int z = thin ? 0 : -2; z <= (thin ? 0 : 3); ++z)
            for (int y = -2; y <= 3; ++y)
                for (int x = -2; x <= 3; ++x) {
                    const auto q = point + Vec3i{x, y, z};
                    if (q.x >= -16 && q.x < 16 && q.y >= -16 && q.y < 16 &&
                        (thin || (q.z >= -16 && q.z < 16)))
                        expected += f.K[0][index(q.x, q.y, q.z, thin)];
                }
        const auto actual = model.abm().environment()->external_growth_counts(point);
        assert(actual[0] == 0.0);
        assert(actual[1] == expected);
    }
    double expected = 0.0;
    for (int z = 0; z <= (thin ? 0 : 3); ++z)
        for (int y = 0; y <= 3; ++y)
            for (int x = 0; x <= 3; ++x)
                expected += f.K[0][index(x, y, z, thin)];
    const double denominator = thin ? 4.0 : 8.0;
    assert(model.abm().environment()->external_activation_density({0, 0, 0}, CellStage::large) ==
           expected / denominator);
}
}

int main() {
    for (bool thin : {true, false}) {
        check_conversion(thin);
        check_prefix(thin);
        check_restart(thin, false);
        check_restart(thin, true);
    }
}
