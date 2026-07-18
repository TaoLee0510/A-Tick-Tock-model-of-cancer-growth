#include "vasculature/angiogenesis_process.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include "core/stateless_rng.hpp"

namespace atcg3d {
namespace {

constexpr double kHoursPerThirtyDays = 30.0 * 24.0;
constexpr std::uint64_t kAngiogenesisProcessUid = 0x414e47494f47454eULL;
constexpr std::uint64_t kSeedWaitingEventKind = 1001;

bool same_time(double lhs, double rhs) noexcept {
    return std::abs(lhs - rhs) <=
           1e-10 * std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

void finite_nonnegative(double value, const char* name) {
    if (!std::isfinite(value) || value < 0.0) {
        throw std::invalid_argument(std::string(name) + " must be finite and nonnegative");
    }
}

}  // namespace

AngiogenesisProcess3D::AngiogenesisProcess3D(
    std::uint64_t seed, std::uint64_t process_uid) noexcept
    : seed_(seed),
      process_uid_(process_uid == 0 ? kAngiogenesisProcessUid : process_uid) {}

void AngiogenesisProcess3D::restore(AngiogenesisProcessState3D state) {
    finite_nonnegative(state.next_seed_time_hours, "next seed time");
    finite_nonnegative(state.eligibility_started_hours, "eligibility start time");
    finite_nonnegative(state.accumulated_eligible_hours, "eligible time");
    if (!state.eligible && state.next_seed_time_hours != 0.0) {
        throw std::invalid_argument("ineligible angiogenesis state has a pending event");
    }
    if (state.committed_roots > state.attempted_events ||
        state.rejected_events != state.attempted_events - state.committed_roots) {
        throw std::invalid_argument(
            "angiogenesis site-arrival counters are inconsistent");
    }
    state_ = state;
}

double AngiogenesisProcess3D::rate_per_hour(double rate_sites_per_30_days) {
    finite_nonnegative(rate_sites_per_30_days, "angiogenesis rate");
    return rate_sites_per_30_days / kHoursPerThirtyDays;
}

double AngiogenesisProcess3D::sample_waiting_hours(
    std::uint64_t seed,
    std::uint64_t event_sequence,
    double rate_sites_per_30_days) {
    return sample_waiting_hours(seed, kAngiogenesisProcessUid, event_sequence,
                                rate_sites_per_30_days);
}

double AngiogenesisProcess3D::sample_waiting_hours(
    std::uint64_t seed,
    std::uint64_t process_uid,
    std::uint64_t event_sequence,
    double rate_sites_per_30_days) {
    const double rate = rate_per_hour(rate_sites_per_30_days);
    if (!(rate > 0.0)) {
        return std::numeric_limits<double>::infinity();
    }
    const double uniform = std::clamp(
        rng_unit(seed, process_uid, kSeedWaitingEventKind, event_sequence),
        1e-15, 1.0 - 1e-15);
    return -std::log1p(-uniform) / rate;
}

bool AngiogenesisProcess3D::update_volume(
    double now_hours,
    double biological_volume_voxels3,
    double activation_volume_voxels3,
    double deactivation_volume_voxels3,
    double activation_delay_hours,
    double rate_sites_per_30_days) {
    finite_nonnegative(now_hours, "current time");
    finite_nonnegative(biological_volume_voxels3, "biological tumour volume");
    finite_nonnegative(activation_volume_voxels3, "activation volume");
    finite_nonnegative(deactivation_volume_voxels3, "deactivation volume");
    finite_nonnegative(activation_delay_hours, "activation delay");
    if (deactivation_volume_voxels3 > activation_volume_voxels3) {
        throw std::invalid_argument("angiogenesis deactivation volume exceeds activation volume");
    }
    (void)rate_per_hour(rate_sites_per_30_days);

    if (!state_.eligible && biological_volume_voxels3 >= activation_volume_voxels3) {
        state_.eligible = true;
        state_.eligibility_started_hours = now_hours;
        schedule_next(now_hours, activation_delay_hours, rate_sites_per_30_days);
        return true;
    }
    if (state_.eligible && biological_volume_voxels3 < deactivation_volume_voxels3) {
        stop(now_hours);
        return true;
    }
    return false;
}

bool AngiogenesisProcess3D::event_current(double time_hours,
                                          std::uint32_t generation) const noexcept {
    return state_.eligible && state_.next_seed_time_hours > 0.0 &&
           generation == state_.schedule_generation &&
           same_time(time_hours, state_.next_seed_time_hours);
}

void AngiogenesisProcess3D::consume_event(double now_hours,
                                          bool committed,
                                          double rate_sites_per_30_days) {
    finite_nonnegative(now_hours, "current time");
    if (!state_.eligible || !same_time(now_hours, state_.next_seed_time_hours)) {
        throw std::logic_error("attempted to consume a stale angiogenesis event");
    }
    ++state_.attempted_events;
    if (committed) {
        ++state_.committed_roots;
    } else {
        ++state_.rejected_events;
    }
    schedule_next(now_hours, 0.0, rate_sites_per_30_days);
}

void AngiogenesisProcess3D::stop(double now_hours) {
    finite_nonnegative(now_hours, "current time");
    if (state_.eligible) {
        accumulate_eligible_time(now_hours);
    }
    state_.eligible = false;
    state_.next_seed_time_hours = 0.0;
    ++state_.schedule_generation;
}

void AngiogenesisProcess3D::schedule_next(double now_hours,
                                          double delay_hours,
                                          double rate_sites_per_30_days) {
    const double wait = sample_waiting_hours(seed_, process_uid_,
                                             state_.event_sequence++,
                                             rate_sites_per_30_days);
    if (!std::isfinite(wait)) {
        state_.next_seed_time_hours = 0.0;
    } else {
        state_.next_seed_time_hours = now_hours + delay_hours + wait;
    }
    ++state_.schedule_generation;
}

void AngiogenesisProcess3D::accumulate_eligible_time(double now_hours) {
    if (now_hours < state_.eligibility_started_hours) {
        throw std::logic_error("angiogenesis time moved backwards");
    }
    state_.accumulated_eligible_hours += now_hours - state_.eligibility_started_hours;
    state_.eligibility_started_hours = now_hours;
}

}  // namespace atcg3d
