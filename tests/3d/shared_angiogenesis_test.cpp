#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <sstream>
#include <vector>

#include "model/shared_angiogenesis.hpp"

namespace {
using atcg3d::continuum::AngiogenesisFieldConfig3D;
using atcg3d::continuum::SharedAngiogenesis3D;

AngiogenesisFieldConfig3D law() {
    AngiogenesisFieldConfig3D config;
    config.model = "shared_vegf_lattice_v2";
    config.seed_tips_per_hour = 20.0;
    config.taf_production_per_cell_hour = 0.03;
    config.taf_diffusion_voxels2_per_hour = 0.2;
    config.tip_chemotaxis = 0.5;
    config.tip_branching_per_hour = 0.1;
    config.tip_anastomosis_per_hour = 0.1;
    return config;
}

void check_parallel_restart(bool individual, bool thin) {
    const std::array<int, 3> shape{17, 17, thin ? 1 : 17};
    const std::size_t size = static_cast<std::size_t>(shape[0]) * shape[1] * shape[2];
    std::vector<double> cells(size, 0.0), nutrient(size, 0.0), vessels(size, 0.0);
    const std::size_t middle = size / 2;
    cells[middle] = cells[middle + 1] = cells[middle + shape[0]] = 1.0;
    vessels[0] = 1.0;
    const auto center_z = thin ? 0 : 8;
    const atcg3d::continuum::VascularConsumerBounds3D range{
        {8, 8, center_z}, {9, 9, center_z}, true};
    SharedAngiogenesis3D one(law(), shape, 1.0, thin, 1);
    SharedAngiogenesis3D eight(law(), shape, 1.0, thin, 8);
    if (individual) {
        one.use_individual_tips(42);
        eight.use_individual_tips(42);
    }
    one.initialize(vessels);
    eight.initialize(vessels);
    for (int step = 0; step < 48; ++step) {
        one.advance(0.25, cells, nutrient, 1.0, true);
        eight.advance(0.25, cells, nutrient, 1.0, true, &range);
        assert(one.checksum() == eight.checksum());
    }
    const auto& diagnostic = one.diagnostics();
    const double tip_mass = std::accumulate(one.tips().begin(), one.tips().end(), 0.0);
    const double counted = diagnostic.seeded_tips + diagnostic.branches - diagnostic.anastomoses -
        diagnostic.discarded_tip_mass;
    assert(std::abs(tip_mass - counted) < 1.0e-10);
    assert(diagnostic.centerline_growth > 0.0);
    if (individual) assert(diagnostic.centerline_growth == std::floor(diagnostic.centerline_growth));
    assert(diagnostic.initial_perfused_volume == 1.0);
    assert(one.vessels()[0] == 1.0);
    for (const double value : one.vessels()) assert(value >= 0.0 && value <= 1.0);
    assert(diagnostic.active_voxels <= size);
    std::stringstream checkpoint(std::ios::in | std::ios::out | std::ios::binary);
    one.save(checkpoint);
    SharedAngiogenesis3D restored(law(), shape, 1.0, thin, 8);
    if (individual) restored.use_individual_tips(42);
    restored.load(checkpoint);
    assert(restored.checksum() == one.checksum());
    for (int step = 0; step < 12; ++step) {
        one.advance(0.25, cells, nutrient, 1.0, true);
        restored.advance(0.25, cells, nutrient, 1.0, true, &range);
        assert(restored.checksum() == one.checksum());
    }
}

void check_gradient_flux() {
    auto config = law();
    config.taf_production_per_cell_hour = 1.0;
    config.taf_diffusion_voxels2_per_hour = config.taf_decay_per_hour = 0.0;
    config.tip_diffusion_voxels2_per_hour = 0.0;
    config.tip_chemotaxis = 1.0;
    config.tip_branching_per_hour = config.tip_anastomosis_per_hour = config.seed_tips_per_hour = 0.0;
    SharedAngiogenesis3D field(config, {5, 5, 1}, 1.0, true, 4);
    std::vector<double> cells(25, 0.0), nutrient(25, 0.0), empty(25, 0.0), tips(25, 0.0);
    cells[13] = 1.0;
    tips[12] = 1.0;
    field.initialize(empty);
    field.advance(0.25, cells, nutrient, 1.0, false);
    field.initialize_tip_density(tips);
    field.advance(0.5, cells, nutrient, 1.0, true);
    assert(field.tips()[12] == 0.875);
    assert(field.tips()[13] == 0.125);
    assert(field.tips()[11] == 0.0 && field.tips()[7] == 0.0);
    assert(std::abs(field.diagnostics().centerline_growth - (1.0 - std::exp(-0.125))) < 1.0e-15);
    const auto length = field.diagnostics().centerline_growth;
    field.advance(0.0, cells, nutrient, 1.0, true);
    assert(field.diagnostics().centerline_growth == length);
}

void check_cutoff_and_validation() {
    auto config = law();
    config.taf_diffusion_voxels2_per_hour = config.taf_decay_per_hour = config.seed_tips_per_hour = 0.0;
    config.taf_production_per_cell_hour = 1.0;
    SharedAngiogenesis3D field(config, {5, 5, 1}, 1.0, true, 1);
    std::vector<double> cells(25, 0.0), nutrient(25, 0.0), empty(25, 0.0);
    cells[12] = 1.0e-15;
    field.initialize(empty);
    field.advance(0.25, cells, nutrient, 1.0, true);
    assert(field.diagnostics().discarded_taf_mass == 0.25e-15);
    assert(std::accumulate(field.taf().begin(), field.taf().end(), 0.0) == 0.0);
    cells[12] = 0.0;
    cells[1] = 1.0;
    field.advance(0.25, cells, nutrient, 1.0, true);
    assert(field.taf()[12] == 0.0);
    assert(field.taf()[1] == 0.25);
    SharedAngiogenesis3D individual(config, {5, 5, 1}, 1.0, true, 1);
    individual.use_individual_tips(1);
    auto fractional = empty;
    fractional[12] = 0.5;
    bool rejected = false;
    try { individual.initialize_tip_density(fractional); }
    catch (const std::invalid_argument&) { rejected = true; }
    assert(rejected);
}

void check_branch_and_connection_moments() {
    constexpr double initial_count = 200000.0;
    for (const bool branching : {true, false}) {
        auto config = law();
        config.taf_production_per_cell_hour = 1.0;
        config.taf_diffusion_voxels2_per_hour = config.taf_decay_per_hour = 0.0;
        config.tip_diffusion_voxels2_per_hour = config.tip_chemotaxis = config.seed_tips_per_hour = 0.0;
        config.tip_branching_per_hour = branching ? 1.0 : 0.0;
        config.tip_anastomosis_per_hour = branching ? 0.0 : 0.2 / initial_count;
        SharedAngiogenesis3D field(config, {5, 5, 1}, 1.0, true, 1);
        field.use_individual_tips(2026);
        std::vector<double> cells(25, 0.0), nutrient(25, 0.0), empty(25, 0.0), tips(25, 0.0);
        cells[12] = 1.0;
        tips[12] = initial_count;
        field.initialize(empty);
        field.advance(0.25, cells, nutrient, 1.0, false);
        field.initialize_tip_density(tips);
        field.advance(0.5, cells, nutrient, 1.0, true);
        const double factor = std::exp(branching ? 0.125 : -0.1);
        const double expected = initial_count * factor;
        const double variance = initial_count * factor * (branching ? factor - 1.0 : 1.0 - factor);
        assert(std::abs(field.tips()[12] - expected) < 6.0 * std::sqrt(variance));
        assert(field.tips()[12] == initial_count + field.diagnostics().branches - field.diagnostics().anastomoses);
    }
}
}  // namespace

int main() {
    for (const bool thin : {true, false}) {
        for (const bool individual : {true, false}) check_parallel_restart(individual, thin);
    }
    check_gradient_flux();
    check_cutoff_and_validation();
    check_branch_and_connection_moments();
}
