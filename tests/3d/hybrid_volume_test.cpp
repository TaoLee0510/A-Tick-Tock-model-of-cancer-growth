#include <cassert>
#include <cmath>
#include <filesystem>

#include "geometry/footprint.hpp"
#include "model/hybrid_model.hpp"
#include "model/shared_resource.hpp"
#include "rules/density.hpp"

namespace {
using namespace atcg3d;
using namespace atcg3d::hybrid;

structured_pde::StructuredInitialFields3D empty_fields(std::size_t size) {
    structured_pde::StructuredInitialFields3D fields;
    for (int stage = 0; stage < 2; ++stage) {
        fields.r_normal[stage].resize(size);
        fields.r_active[stage].resize(size);
        fields.K[stage].resize(size);
        fields.active_remaining_hours[stage].resize(size);
    }
    fields.vessel_fraction.resize(size);
    return fields;
}

HybridConfig3D config(bool thin) {
    auto c = HybridConfig3D::load(std::filesystem::path(ATCG_SOURCE_DIR) /
        "ATCG3D_Hybrid/config/hybrid_regular_cycle_v3.yaml");
    c.rules.continuum.base.thin_layer = thin;
    c.rules.continuum.grid.shape = {12, 12, thin ? 1 : 12};
    c.rules.continuum.grid.origin = {-6, -6, thin ? -0.5 : -6};
    c.rules.continuum.end_time_hours = 1;
    c.rules.continuum.base.growth_density_window_edge = 1;
    c.rules.continuum.base.migration_activation_window_edge = 4;
    c.rules.continuum.base.migration_activation_block_edge = 2;
    c.rules.continuum.base.migration_activation_threshold = 1;
    c.rules.continuum.migration.diffusion_scale = 1e-12;
    c.rules.continuum.reaction.small_daughter_vacancy_exponent = thin ? 8 : 26;
    c.rules.continuum.vascular.source_mode = "static_voxels";
    c.rules.continuum.vascular.static_sources.clear();
    c.rules.continuum.nutrient.boundary_mode = "vessels_dirichlet_v1";
    c.rules.continuum.nutrient.diffusion_voxels2_per_hour = 1e-12;
    c.rules.continuum.nutrient.decay_per_hour = 0;
    c.rules.continuum.nutrient.initial_value = 0.2;
    c.rules.continuum.nutrient.K_consumption_rate_per_hour = 0.01;
    c.rules.continuum.nutrient.r_consumption_rate_per_hour = 0.01;
    auto& vascular = c.rules.continuum.angiogenesis;
    vascular.model = "vegf_tip_density_v1";
    vascular.seed_tips_per_hour = 0;
    vascular.taf_diffusion_voxels2_per_hour = 0;
    vascular.taf_decay_per_hour = 0;
    c.exchange_every_hours = 2;
    c.smoothing_radius = 2;
    c.core_on = 0.99;
    c.core_off = 0.98;
    return c;
}

std::size_t location(Vec3i point, bool thin) {
    return (std::size_t(thin ? 0 : point.z + 6) * 12 + point.y + 6) * 12 +
        point.x + 6;
}

double activation(const HybridModel3D& model, Vec3i anchor) {
    const auto& c = model.abm().config();
    return migration_activation_density(model.abm().density(), anchor,
        CellStage::large, c.migration_activation_window_edge,
        c.migration_activation_block_edge, c.thin_layer) +
        model.abm().environment()->external_activation_density(anchor, CellStage::large);
}

structured_pde::StructuredInitialFields3D large_fields(bool thin, CellType type) {
    const double volume = thin ? 4 : 8;
    auto fields = empty_fields(thin ? 144 : 1728);
    for (auto point : large_footprint({0, 0, 0})) {
        if (thin && point.z != 0) continue;
        auto& population = type == CellType::r ? fields.r_normal[1] : fields.K[1];
        population[location(point, thin)] = 1 / volume;
    }
    return fields;
}

void check_large_cell(bool thin, CellType type) {
    auto c = config(thin);
    const double volume = thin ? 4 : 8;
    HybridModel3D model(c);
    model.initialize_from_arrays(large_fields(thin, type));
    const double density = activation(model, {0, 0, 0});
    assert(std::abs(density - volume / std::pow(4.0, thin ? 2 : 3)) < 1e-15);
    model.exchange();
    assert(model.diagnostics().abm_mass == 1);
    assert(model.diagnostics().pde_mass == 0);
    const auto slot = model.abm().cells().alive_slots().front();
    assert(model.abm().cells().stage(slot) == CellStage::large);
    const auto anchor = model.abm().cells().anchor(slot);
    assert(std::abs(activation(model, anchor) - density) < 1e-15);
    double consumers = 0;
    for (auto point : large_footprint(anchor)) {
        if (thin && point.z != 0) continue;
        const auto i = location(point, thin);
        assert(model.pde().external_consumers(i) == 1 / volume);
        assert(!model.abm().environment()->destination_available(point));
        consumers += model.pde().external_consumers(i);
    }
    assert(consumers == 1);
    assert(model.step());
    for (auto point : large_footprint(anchor)) {
        if (thin && point.z != 0) continue;
        const auto i = location(point, thin);
        const double expected = continuum::resource_after_uptake(0.2,
            1 / volume, 0.01, c.rules.continuum.nutrient.K_consumption_half_saturation,
            c.rules.continuum.time_step_hours, 1.0);
        assert(std::abs(model.pde().nutrient()[i] - expected) < 1e-15);
        const double hypoxia = 1 - 0.2 / c.rules.continuum.angiogenesis.hypoxia_threshold;
        const double taf = c.rules.continuum.time_step_hours * hypoxia / volume;
        assert(std::abs(model.pde().angiogenesis()->taf()[i] - taf) < 1e-15);
    }
    const auto path = std::filesystem::current_path() / "hybrid-volume.bin";
    for (auto suffix : {"", ".pde.bin", ".resource.bin"})
        std::filesystem::remove(path.string() + suffix);
    model.save_checkpoint(path);
    c.rules.continuum.base.threads = 4;
    HybridModel3D resumed(c);
    resumed.load_checkpoint(path);
    assert(model.state_checksum() == resumed.state_checksum());
    while (model.step()) {
        assert(resumed.step());
        assert(model.state_checksum() == resumed.state_checksum());
    }
    for (auto suffix : {"", ".pde.bin", ".resource.bin"})
        std::filesystem::remove(path.string() + suffix);
}

void check_agent_removal(bool thin, CellType type) {
    auto c = config(thin);
    c.rules.continuum.angiogenesis.seed_tips_per_hour = 100;
    c.rules.continuum.angiogenesis.tip_speed_voxels_per_hour = 4;
    c.rules.continuum.angiogenesis.exclusion_fraction = 0.001;
    HybridModel3D model(c);
    model.initialize_from_arrays(large_fields(thin, type));
    model.exchange();
    assert(model.diagnostics().abm_mass == 1);
    assert(model.step());
    assert(model.diagnostics().total_mass == 0);
    const auto removed = model.pde().diagnostics().vascular_removed_mass;
    assert(removed.total() == 1);
    assert((type == CellType::r ? removed.r_normal[1] : removed.K[1]) == 1);
    assert(model.step());
    assert(model.pde().diagnostics().vascular_removed_mass.total() == 1);
}

void check_capacity() {
    auto c = config(true);
    c.rules.migration.minimum_density = 1e-4;
    auto fields = empty_fields(144);
    fields.r_normal[0][location({0, 0, 0}, true)] = 1e-5;
    HybridModel3D volume(c);
    volume.initialize_from_arrays(fields);
    assert(!volume.abm().environment()->destination_available({0, 0, 0}));
    assert(volume.abm().environment()->destination_available({1, 0, 0}));
    c.schema_version = 2;
    c.model = "hybrid_distributions_v2";
    c.rules.schema_version = 12;
    HybridModel3D legacy(c);
    legacy.initialize_from_arrays(fields);
    assert(legacy.abm().environment()->destination_available({0, 0, 0}));
    c.schema_version = 3;
    c.model = "hybrid_volume_coupling_v3";
    c.rules.schema_version = 13;
    c.rules.continuum.reaction.maximum_occupied_fraction = 0.8;
    HybridModel3D limited(c);
    limited.initialize_from_arrays(empty_fields(144));
    assert(!limited.abm().environment()->destination_available({0, 0, 0}));
}
} // namespace

int main() {
    for (bool thin : {true, false}) {
        for (auto type : {atcg3d::CellType::r, atcg3d::CellType::K}) {
            check_large_cell(thin, type);
            check_agent_removal(thin, type);
        }
    }
    check_capacity();
}
