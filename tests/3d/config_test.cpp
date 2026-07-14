#include <cassert>
#include <limits>
#include <stdexcept>

#include "config/model_config.hpp"

int main() {
    using namespace atcg3d;
    Model3DConfig config;
    config.validate();
    config.apply_override("direction.continue_probability=0.75");
    config.apply_override("space.thin_layer=true");
    config.apply_override("output.directory='test-output'");
    assert(config.continue_probability == 0.75);
    assert(config.thin_layer);
    assert(config.output_directory == "test-output");
    assert(config.to_json().find("legacy_like_v1") != std::string::npos);

    bool rejected = false;
    try {
        config.apply_override("unknown.parameter=1");
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    assert(rejected);
    Model3DConfig contradictory;
    contradictory.domain_policy = "bounded";
    contradictory.bounded_domain = false;
    rejected = false;
    try {
        contradictory.validate();
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    assert(rejected);
    rejected = false;
    try {
        config.apply_override("direction.continue_probability=2");
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    assert(rejected);
    assert(config.continue_probability == 0.75);
    rejected = false;
    try {
        config.apply_override("simulation.end_time_hours=nan");
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    assert(rejected);
    rejected = false;
    try {
        config.apply_override("simulation.end_time_hours=-1");
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    assert(rejected);
}
