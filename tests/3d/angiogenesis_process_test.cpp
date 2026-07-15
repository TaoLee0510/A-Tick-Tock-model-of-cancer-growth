#include <cassert>
#include <cmath>
#include <stdexcept>

#include "vasculature/angiogenesis_process.hpp"

int main() {
    using namespace atcg3d;

    assert(std::abs(AngiogenesisProcess3D::rate_per_hour(1.0) - 1.0 / 720.0) < 1e-15);
    const double waiting = AngiogenesisProcess3D::sample_waiting_hours(7, 0, 2.0);
    assert(waiting > 0.0 && std::isfinite(waiting));
    assert(waiting == AngiogenesisProcess3D::sample_waiting_hours(7, 0, 2.0));

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
}
