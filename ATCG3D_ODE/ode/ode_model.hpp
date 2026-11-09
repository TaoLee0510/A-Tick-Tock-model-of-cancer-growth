#pragma once

#include <array>
#include <filesystem>
#include <vector>

#include "config/structured_config.hpp"

namespace atcg3d::ode {

// Cell-number densities, active time-mass, nutrient, vascular capacity,
// refractory ordinary-r subsets and their time-mass, respectively.
using State = std::array<double, 14>;
inline constexpr std::size_t nutrient_index = 8;
inline constexpr std::size_t vessel_index = 9;

struct SolverConfig {
    double relative_tolerance{1.0e-8};
    double absolute_tolerance{1.0e-11};
    double maximum_step_hours{0.1};
};
struct OdeConfig3D {
    structured_pde::StructuredPdeConfig3D rules;
    SolverConfig solver;
    State initial{};
    std::filesystem::path output_directory{"atcg3d_ode_run"};
    static OdeConfig3D load(const std::filesystem::path& path);
    void validate() const;
};

class Reaction3D {
public:
    explicit Reaction3D(structured_pde::StructuredPdeConfig3D rules);
    State derivative(const State& state, double r_count, double K_count) const;
    void advance(State& state, double duration, const SolverConfig& solver,
                 double r_count = -1.0, double K_count = -1.0) const;
    double total_mass(const State& state) const noexcept;
    double occupied(const State& state) const noexcept;
    double window_measure() const noexcept { return window_measure_; }
private:
    void transitions(State& state) const;
    structured_pde::StructuredPdeConfig3D rules_;
    double window_measure_{}, large_volume_{}, r_inherent_{}, K_inherent_{};
};

class OdeModel3D {
public:
    explicit OdeModel3D(OdeConfig3D config);
    bool step();
    void run();
    void save_checkpoint(const std::filesystem::path& path) const;
    void load_checkpoint(const std::filesystem::path& path);
    std::uint64_t state_checksum() const;
    const State& state() const noexcept { return state_; }
    double time_hours() const noexcept { return time_; }
private:
    OdeConfig3D config_;
    Reaction3D reaction_;
    State state_;
    double time_{};
};

// Periodic finite-volume reference for the same reaction contract. All
// compartments, including clocks, share their corresponding mass flux.
class PeriodicReactionPde3D {
public:
    explicit PeriodicReactionPde3D(OdeConfig3D config);
    void initialize(std::vector<State> fields);
    bool step();
    const std::vector<State>& fields() const noexcept { return fields_; }
    State mean() const;
    double time_hours() const noexcept { return time_; }
private:
    std::size_t index(int x, int y, int z) const noexcept;
    OdeConfig3D config_;
    Reaction3D reaction_;
    std::vector<State> fields_, work_;
    double time_{};
};
}  // namespace atcg3d::ode
