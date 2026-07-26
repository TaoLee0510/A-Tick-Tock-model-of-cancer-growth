#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

#include "engine/simulation.hpp"
#include "vasculature/angiogenesis_process.hpp"

int main() {
    using namespace atcg3d;

    assert(std::abs(AngiogenesisProcess3D::rate_per_hour(1.0) - 1.0 / 720.0) < 1e-15);
    const double waiting = AngiogenesisProcess3D::sample_waiting_hours(7, 0, 2.0);
    assert(waiting > 0.0 && std::isfinite(waiting));
    assert(waiting == AngiogenesisProcess3D::sample_waiting_hours(7, 0, 2.0));
    const double lesion_one_wait =
        AngiogenesisProcess3D::sample_waiting_hours(7, 101, 0, 2.0);
    const double lesion_two_wait =
        AngiogenesisProcess3D::sample_waiting_hours(7, 102, 0, 2.0);
    assert(lesion_one_wait ==
           AngiogenesisProcess3D::sample_waiting_hours(7, 101, 0, 2.0));
    assert(lesion_one_wait != lesion_two_wait);

    double accumulated_wait = 0.0;
    constexpr std::uint64_t sample_count = 20000;
    for (std::uint64_t sequence = 0; sequence < sample_count; ++sequence) {
        accumulated_wait += AngiogenesisProcess3D::sample_waiting_hours(19, sequence, 1.0);
    }
    const double empirical_mean_hours = accumulated_wait / static_cast<double>(sample_count);
    assert(std::abs(empirical_mean_hours - 720.0) / 720.0 < 0.03);

    AngiogenesisProcess3D process(7);
    assert(!process.update_volume(0.0, 99.0, 100.0, 80.0, 0.0, 2.0));
    assert(process.update_volume(1.0, 100.0, 100.0, 80.0, 0.0, 2.0));
    const auto activated = process.state();
    assert(activated.eligible);
    assert(activated.next_seed_time_hours > 1.0);
    assert(process.event_current(activated.next_seed_time_hours,
                                 activated.schedule_generation));

    process.consume_event(activated.next_seed_time_hours, true, 2.0);
    assert(process.state().attempted_events == 1);
    assert(process.state().committed_roots == 1);
    assert(process.state().rejected_events == 0);
    assert(process.state().next_seed_time_hours > activated.next_seed_time_hours);

    const auto second_arrival = process.state();
    process.consume_event(second_arrival.next_seed_time_hours, false, 2.0);
    assert(process.state().attempted_events == 2);
    assert(process.state().committed_roots == 1);
    assert(process.state().rejected_events == 1);

    const auto before_pause = process.state();
    assert(process.update_volume(before_pause.next_seed_time_hours - 0.1,
                                 79.0, 100.0, 80.0, 0.0, 2.0));
    assert(!process.state().eligible);
    assert(process.state().next_seed_time_hours == 0.0);
    assert(!process.event_current(before_pause.next_seed_time_hours,
                                  before_pause.schedule_generation));

    const std::uint64_t sequence_before_reactivation = process.state().event_sequence;
    assert(process.update_volume(before_pause.next_seed_time_hours,
                                 101.0, 100.0, 80.0, 0.0, 2.0));
    assert(process.state().event_sequence == sequence_before_reactivation + 1);

    AngiogenesisProcess3D restored(7);
    restored.restore(process.state());
    assert(restored.state().next_seed_time_hours == process.state().next_seed_time_hours);
    assert(restored.event_current(process.state().next_seed_time_hours,
                                  process.state().schedule_generation));

    AngiogenesisProcessState3D invalid = process.state();
    invalid.rejected_events = invalid.attempted_events + 1;
    bool rejected_invalid = false;
    try {
        restored.restore(invalid);
    } catch (const std::invalid_argument&) {
        rejected_invalid = true;
    }
    assert(rejected_invalid);

    invalid = process.state();
    invalid.committed_roots = invalid.attempted_events + 1;
    rejected_invalid = false;
    try {
        restored.restore(invalid);
    } catch (const std::invalid_argument&) {
        rejected_invalid = true;
    }
    assert(rejected_invalid);

    // A density refresh changes the slope of the integrated hazard, not the
    // sampled unit-exponential target. This permits a genuinely dynamic
    // Poisson process without inventing a new random wait at every refresh.
    AngiogenesisProcess3D modulated(41, 9001);
    assert(modulated.update_volume(0.0, 1.0, 1.0, 0.0, 0.0, 2.0));
    const auto rate_two = modulated.state();
    const double first_wait = rate_two.next_seed_time_hours;
    const double first_hazard = rate_two.remaining_hazard;
    assert(first_wait == first_hazard / AngiogenesisProcess3D::rate_per_hour(2.0));
    const double refresh_time = first_wait * 0.25;
    assert(modulated.update_rate(refresh_time, 4.0, 0.75));
    const auto rate_four = modulated.state();
    assert(std::abs(rate_four.remaining_hazard - first_hazard * 0.75) < 1e-12);
    assert(std::abs(rate_four.next_seed_time_hours - first_wait * 0.625) < 1e-10);
    assert(rate_four.current_density_stress == 0.75);
    const double unchanged_time = rate_four.next_seed_time_hours;
    const std::uint32_t unchanged_generation = rate_four.schedule_generation;
    assert(!modulated.update_rate(refresh_time + 0.01, 4.0, 0.70));
    assert(modulated.state().next_seed_time_hours == unchanged_time);
    assert(modulated.state().schedule_generation == unchanged_generation);
    assert(modulated.state().current_density_stress == 0.70);
    assert(modulated.update_rate(refresh_time + 0.02, 0.0, 0.0));
    assert(modulated.state().next_seed_time_hours == 0.0);
    const double paused_hazard = modulated.state().remaining_hazard;
    assert(modulated.update_rate(refresh_time + 100.0, 2.0, 0.5));
    assert(modulated.state().remaining_hazard == paused_hazard);
    assert(modulated.state().next_seed_time_hours > refresh_time + 100.0);

    // The compatibility aggregate is a deterministic snapshot, not a second
    // scheduler. Input order must not change any bit, and every active
    // process contributes its elapsed interval through the snapshot time.
    LesionAngiogenesisState3D lesion_ten;
    lesion_ten.lesion_id = 10;
    lesion_ten.process.eligible = true;
    lesion_ten.process.next_seed_time_hours = 20.0;
    lesion_ten.process.eligibility_started_hours = 2.0;
    lesion_ten.process.accumulated_eligible_hours = 4.0;
    lesion_ten.process.remaining_hazard = 18.0;
    lesion_ten.process.hazard_last_update_hours = 2.0;
    lesion_ten.process.hazard_not_before_hours = 2.0;
    lesion_ten.process.current_rate_sites_per_30_days = 720.0;
    lesion_ten.process.current_density_stress = 0.8;
    lesion_ten.process.event_sequence = 10;
    lesion_ten.process.schedule_generation = 3;
    lesion_ten.process.attempted_events = 3;
    lesion_ten.process.committed_roots = 1;
    lesion_ten.process.rejected_events = 2;

    LesionAngiogenesisState3D lesion_twenty;
    lesion_twenty.lesion_id = 20;
    lesion_twenty.process.eligible = true;
    lesion_twenty.process.next_seed_time_hours = 18.0;
    lesion_twenty.process.eligibility_started_hours = 6.0;
    lesion_twenty.process.accumulated_eligible_hours = 8.0;
    lesion_twenty.process.remaining_hazard = 12.0;
    lesion_twenty.process.hazard_last_update_hours = 6.0;
    lesion_twenty.process.hazard_not_before_hours = 6.0;
    lesion_twenty.process.current_rate_sites_per_30_days = 720.0;
    lesion_twenty.process.current_density_stress = 0.6;
    lesion_twenty.process.event_sequence = 20;
    lesion_twenty.process.schedule_generation = 7;
    lesion_twenty.process.attempted_events = 5;
    lesion_twenty.process.committed_roots = 4;
    lesion_twenty.process.rejected_events = 1;

    LesionAngiogenesisState3D lesion_thirty;
    lesion_thirty.lesion_id = 30;
    lesion_thirty.process.accumulated_eligible_hours = 1.0e16;
    lesion_thirty.process.event_sequence = 100;
    lesion_thirty.process.schedule_generation = 9;
    lesion_thirty.process.attempted_events = 7;
    lesion_thirty.process.committed_roots = 2;
    lesion_thirty.process.rejected_events = 5;

    const std::vector<LesionAngiogenesisState3D> unordered_processes{
        lesion_thirty, lesion_ten, lesion_twenty};
    const std::vector<LesionAngiogenesisState3D> differently_ordered_processes{
        lesion_twenty, lesion_thirty, lesion_ten};
    const AngiogenesisProcessState3D aggregate =
        aggregate_angiogenesis_process_states(unordered_processes, 10.0);
    const AngiogenesisProcessState3D reordered_aggregate =
        aggregate_angiogenesis_process_states(differently_ordered_processes,
                                               10.0);
    assert(aggregate.eligible);
    assert(aggregate.next_seed_time_hours == 18.0);
    assert(aggregate.eligibility_started_hours == 10.0);
    assert(aggregate.remaining_hazard == 8.0);
    assert(aggregate.hazard_last_update_hours == 10.0);
    assert(aggregate.hazard_not_before_hours == 10.0);
    assert(aggregate.current_rate_sites_per_30_days == 720.0);
    assert(aggregate.current_density_stress == 0.6);
    assert(aggregate.accumulated_eligible_hours == 1.0e16 + 24.0);
    assert(aggregate.event_sequence == 100);
    assert(aggregate.schedule_generation == 9);
    assert(aggregate.attempted_events == 15);
    assert(aggregate.committed_roots == 7);
    assert(aggregate.rejected_events == 8);
    assert(aggregate.eligible == reordered_aggregate.eligible);
    assert(aggregate.next_seed_time_hours ==
           reordered_aggregate.next_seed_time_hours);
    assert(aggregate.eligibility_started_hours ==
           reordered_aggregate.eligibility_started_hours);
    assert(aggregate.accumulated_eligible_hours ==
           reordered_aggregate.accumulated_eligible_hours);
    assert(aggregate.event_sequence == reordered_aggregate.event_sequence);
    assert(aggregate.schedule_generation ==
           reordered_aggregate.schedule_generation);
    assert(aggregate.attempted_events == reordered_aggregate.attempted_events);
    assert(aggregate.committed_roots == reordered_aggregate.committed_roots);
    assert(aggregate.rejected_events == reordered_aggregate.rejected_events);

    std::vector<LesionAngiogenesisState3D> overflowing_processes(2);
    overflowing_processes[0].lesion_id = 1;
    overflowing_processes[0].process.attempted_events =
        std::numeric_limits<std::uint64_t>::max();
    overflowing_processes[0].process.committed_roots =
        std::numeric_limits<std::uint64_t>::max();
    overflowing_processes[1].lesion_id = 2;
    overflowing_processes[1].process.attempted_events = 1;
    overflowing_processes[1].process.committed_roots = 1;
    bool rejected_overflow = false;
    try {
        (void)aggregate_angiogenesis_process_states(overflowing_processes, 0.0);
    } catch (const std::overflow_error&) {
        rejected_overflow = true;
    }
    assert(rejected_overflow);
}
