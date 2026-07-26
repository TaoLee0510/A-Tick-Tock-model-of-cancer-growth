#pragma once

#include <cstdint>

namespace atcg3d {

struct AngiogenesisProcessState3D {
    bool eligible{};
    double next_seed_time_hours{};
    double eligibility_started_hours{};
    double accumulated_eligible_hours{};
    // Unit-exponential hazard remaining until the next site arrival. Keeping
    // this value, rather than resampling an absolute time, makes hourly
    // density-dependent rate changes statistically exact and restart-safe.
    double remaining_hazard{};
    double hazard_last_update_hours{};
    double hazard_not_before_hours{};
    double current_rate_sites_per_30_days{};
    double current_density_stress{};
    std::uint64_t event_sequence{};
    std::uint32_t schedule_generation{};
    std::uint64_t attempted_events{};
    std::uint64_t committed_roots{};
    std::uint64_t rejected_events{};
};

// Event-driven homogeneous Poisson process for one spatial lesion. The
// configured intensity is the expected number of attempted surface-root sites
// per 30 eligible days for that lesion, not a per-surface-voxel rate. One
// arrival attempts at most one root.
class AngiogenesisProcess3D {
public:
    explicit AngiogenesisProcess3D(std::uint64_t seed = 1,
                                   std::uint64_t process_uid = 0) noexcept;

    const AngiogenesisProcessState3D& state() const noexcept { return state_; }
    void restore(AngiogenesisProcessState3D state);

    // Returns true when eligibility or the pending event changed.
    bool update_volume(double now_hours,
                       double biological_volume_voxels3,
                       double activation_volume_voxels3,
                       double deactivation_volume_voxels3,
                       double activation_delay_hours,
                       double rate_sites_per_30_days);

    // Changes the piecewise-constant Poisson intensity after first consuming
    // hazard accumulated under the previous rate. Returns true if the pending
    // event time/generation changed.
    bool update_rate(double now_hours,
                     double rate_sites_per_30_days,
                     double density_stress);

    bool event_current(double time_hours, std::uint32_t generation) const noexcept;

    // Consumes exactly one current Poisson arrival and schedules the following
    // arrival. A rejected surface sample still consumes one site arrival.
    void consume_event(double now_hours,
                       bool committed,
                       double rate_sites_per_30_days);

    void stop(double now_hours);

    static double rate_per_hour(double rate_sites_per_30_days);
    static double sample_waiting_hours(std::uint64_t seed,
                                       std::uint64_t event_sequence,
                                       double rate_sites_per_30_days);
    static double sample_unit_hazard(std::uint64_t seed,
                                     std::uint64_t process_uid,
                                     std::uint64_t event_sequence);
    static double sample_waiting_hours(std::uint64_t seed,
                                       std::uint64_t process_uid,
                                       std::uint64_t event_sequence,
                                       double rate_sites_per_30_days);

private:
    void schedule_next(double now_hours,
                       double delay_hours,
                       double rate_sites_per_30_days);
    void consume_hazard_until(double now_hours);
    void recompute_next_time(double now_hours);
    void accumulate_eligible_time(double now_hours);

    std::uint64_t seed_{};
    std::uint64_t process_uid_{};
    AngiogenesisProcessState3D state_{};
};

}  // namespace atcg3d
