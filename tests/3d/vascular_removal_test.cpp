#include <algorithm>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "io/structured_pde_output.hpp"
#include "model/structured_pde_model.hpp"

namespace {
using namespace atcg3d;
using namespace atcg3d::structured_pde;

std::vector<std::string> columns(const std::string& line) {
    std::istringstream stream(line);
    std::vector<std::string> values;
    std::string value;
    while (std::getline(stream, value, ',')) values.push_back(value);
    return values;
}

void check_removal(bool thin) {
    auto c = StructuredPdeConfig3D::load(std::filesystem::path(ATCG_SOURCE_DIR) /
        "ATCG3D_SharedRules/config/angiogenesis_v8.yaml");
    c.schema_version = 13;
    c.continuum.base.thin_layer = thin;
    c.continuum.grid.shape = {12, 12, thin ? 1 : 12};
    c.continuum.grid.origin = {-6, -6, thin ? -0.5 : -6};
    c.continuum.nutrient.boundary_mode = "vessels_dirichlet_v1";
    c.continuum.nutrient.initial_value = 0.01;
    c.continuum.vascular.source_mode = "static_voxels";
    c.continuum.vascular.static_sources.clear();
    c.continuum.end_time_hours = 1;
    c.continuum.time_step_hours = 0.25;
    c.continuum.reaction.small_daughter_vacancy_exponent = thin ? 8 : 26;
    c.continuum.angiogenesis.seed_tips_per_hour = 100;
    c.continuum.angiogenesis.tip_speed_voxels_per_hour = 4;
    c.continuum.angiogenesis.exclusion_fraction = 0.001;
    StructuredInitialFields3D fields;
    const std::size_t size = thin ? 144 : 1728;
    const std::size_t here = thin ? 78 : 942;
    for (int stage = 0; stage < 2; ++stage) {
        fields.r_normal[stage].resize(size);
        fields.r_active[stage].resize(size);
        fields.K[stage].resize(size);
        fields.active_remaining_hours[stage].resize(size);
        fields.r_normal[stage][here] = stage == 0 ? 0.05 : 0.02;
        fields.r_active[stage][here] = stage == 0 ? 0.04 : 0.01;
        fields.K[stage][here] = stage == 0 ? 0.03 : 0.02;
        fields.active_remaining_hours[stage][here] = 4;
    }
    fields.vessel_fraction.resize(size);
    StructuredPdeModel3D model(c);
    model.initialize_from_arrays(fields);
    const auto initial = model.diagnostics();
    assert(initial.vascular_removed_mass.total() == 0);
    assert(model.step());
    const auto after = model.diagnostics();
    assert(after.r_total == 0 && after.K_total == 0);
    for (int stage = 0; stage < 2; ++stage) {
        assert(after.vascular_removed_mass.r_normal[stage] == initial.r_normal_mass[stage]);
        assert(after.vascular_removed_mass.r_active[stage] == initial.r_active_mass[stage]);
        assert(after.vascular_removed_mass.K[stage] == initial.K_mass[stage]);
    }
    assert(std::abs(after.vascular_removed_mass.total() -
        initial.r_total - initial.K_total) < 1e-15);
    const auto directory = std::filesystem::current_path() /
        (thin ? "vascular-removal-2d" : "vascular-removal-3d");
    std::filesystem::remove_all(directory);
    std::filesystem::create_directory(directory);
    const auto checkpoint = directory / "checkpoint.bin";
    model.save_checkpoint(checkpoint);
    auto parallel = c;
    parallel.continuum.base.threads = 4;
    StructuredPdeModel3D resumed(parallel);
    resumed.load_checkpoint(checkpoint);
    assert(model.state_checksum() == resumed.state_checksum());
    while (model.step()) {
        assert(resumed.step());
        assert(model.state_checksum() == resumed.state_checksum());
        assert(model.diagnostics().vascular_removed_mass.total() ==
            after.vascular_removed_mass.total());
    }
    c.continuum.output.directory = directory / "output";
    c.continuum.output.field_every_hours = 0;
    c.continuum.output.checkpoint_every_hours = 0;
    {
        StructuredPdeOutput3D output(c, model.time_hours());
        output.observe(model);
        output.finalize(model);
    }
    std::ifstream metrics(c.continuum.output.directory / "metrics.csv");
    std::string header, row;
    std::getline(metrics, header);
    std::getline(metrics, row);
    const auto names = columns(header);
    const auto values = columns(row);
    assert(names.size() == values.size());
    const std::vector<std::string> removed_names{
        "vascular_removed_r_normal_small", "vascular_removed_r_normal_large",
        "vascular_removed_r_active_small", "vascular_removed_r_active_large",
        "vascular_removed_K_small", "vascular_removed_K_large", "vascular_removed_mass"};
    std::vector<double> expected;
    for (const auto* masses : {&after.vascular_removed_mass.r_normal,
                               &after.vascular_removed_mass.r_active,
                               &after.vascular_removed_mass.K})
        for (double mass : *masses) expected.push_back(mass);
    expected.push_back(after.vascular_removed_mass.total());
    for (std::size_t index = 0; index < removed_names.size(); ++index) {
        const auto found = std::find(names.begin(), names.end(), removed_names[index]);
        assert(found != names.end());
        assert(std::stod(values[std::size_t(found - names.begin())]) == expected[index]);
    }
    c.continuum.run_mode = "resume";
    {
        StructuredPdeOutput3D output(c, model.time_hours());
    }
    {
        std::ofstream incompatible(c.continuum.output.directory / "metrics.csv");
        incompatible << "time_hours,state_checksum\n";
    }
    bool rejected = false;
    try {
        StructuredPdeOutput3D output(c, model.time_hours());
    } catch (const std::runtime_error&) {
        rejected = true;
    }
    assert(rejected);
    std::filesystem::remove_all(directory);

    // Published schemas still delete mass without new state or output columns.
    c.schema_version = 12;
    c.continuum.run_mode = "new";
    StructuredPdeModel3D legacy(c);
    legacy.initialize_from_arrays(std::move(fields));
    assert(legacy.step());
    assert(legacy.diagnostics().r_total == 0 && legacy.diagnostics().K_total == 0);
    assert(legacy.diagnostics().vascular_removed_mass.total() == 0);
}
} // namespace

int main() {
    check_removal(true);
    check_removal(false);
}
