#include "core/cell_store.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>

namespace atcg3d {
namespace {

template <class T>
std::size_t vector_bytes(const std::vector<T>& values) {
    return values.capacity() * sizeof(T);
}

std::size_t stage_index(CellStage stage) {
    const auto index = static_cast<std::size_t>(stage);
    if (index >= 3) {
        throw std::invalid_argument("invalid cell stage");
    }
    return index;
}

enum class TimeRounding {
    nearest,
    upward,
};

float quantize_time(double value, TimeRounding rounding,
                    const char* field_name) {
    if (!std::isfinite(value) || value < 0.0 ||
        value > static_cast<double>(std::numeric_limits<float>::max())) {
        throw std::invalid_argument(
            std::string(field_name) +
            " must be finite, nonnegative, and representable as float");
    }
    float stored = static_cast<float>(value);
    if (rounding == TimeRounding::upward &&
        static_cast<double>(stored) < value) {
        stored = std::nextafter(
            stored, std::numeric_limits<float>::infinity());
    }
    if (!std::isfinite(stored)) {
        throw std::invalid_argument(
            std::string(field_name) + " cannot be represented as float");
    }
    return stored;
}

void validate_cell_times(const CellInit& cell) {
    (void)quantize_time(cell.migration_activation_end_time,
                        TimeRounding::upward,
                        "migration_activation_end_time");
    (void)quantize_time(cell.next_migration_time, TimeRounding::upward,
                        "next_migration_time");
    (void)quantize_time(cell.next_division_time, TimeRounding::upward,
                        "next_division_time");
    (void)quantize_time(cell.death_deadline, TimeRounding::upward,
                        "death_deadline");
    (void)quantize_time(cell.last_update_time, TimeRounding::nearest,
                        "last_update_time");
    (void)quantize_time(cell.swap_ready_time, TimeRounding::upward,
                        "swap_ready_time");
    if (cell.swap_wait_state > 1 ||
        cell.pending_swap_direction > 26) {
        throw std::invalid_argument("invalid crowding-swap state");
    }
    if (cell.swap_wait_state == 0 &&
        (cell.swap_ready_time != 0.0 ||
         cell.pending_swap_direction != kStayDirection)) {
        throw std::invalid_argument(
            "inactive crowding-swap state must have zero time and direction");
    }
    if (cell.swap_wait_state == 1 &&
        (!(cell.swap_ready_time > 0.0) ||
         cell.pending_swap_direction == kStayDirection)) {
        throw std::invalid_argument(
            "active crowding-swap state requires a ready time and direction");
    }
}

bool has_event_specific_generations(const CellInit& cell) noexcept {
    return cell.migration_schedule_generation != 0 ||
           cell.division_schedule_generation != 0 ||
           cell.death_schedule_generation != 0;
}

std::uint32_t migration_generation(const CellInit& cell) noexcept {
    return has_event_specific_generations(cell)
        ? cell.migration_schedule_generation : cell.schedule_generation;
}

std::uint32_t division_generation(const CellInit& cell) noexcept {
    return has_event_specific_generations(cell)
        ? cell.division_schedule_generation : cell.schedule_generation;
}

std::uint32_t death_generation(const CellInit& cell) noexcept {
    return has_event_specific_generations(cell)
        ? cell.death_schedule_generation : cell.schedule_generation;
}

}  // namespace

void CellStore3D::reserve(std::size_t count) {
    x_.reserve(count); y_.reserve(count); z_.reserve(count);
    uid_.reserve(count); parent_uid_.reserve(count); clone_id_.reserve(count);
    type_.reserve(count); stage_.reserve(count); viability_.reserve(count);
    flags_.reserve(count); last_direction_.reserve(count); alive_.reserve(count);
    inherent_growth_rate_.reserve(count); density_growth_rate_.reserve(count); migration_rate_.reserve(count);
    normal_migration_rate_.reserve(count); migration_activation_end_time_.reserve(count);
    division_work_remaining_.reserve(count);
    next_migration_time_.reserve(count); next_division_time_.reserve(count);
    death_deadline_.reserve(count); last_update_time_.reserve(count);
    pending_swap_direction_.reserve(count);
    event_sequence_.reserve(count);
    migration_schedule_generation_.reserve(count);
    division_schedule_generation_.reserve(count);
    death_schedule_generation_.reserve(count);
    checkpoint_dirty_.reserve(count);
    checkpoint_dirty_slots_.reserve(std::min<std::size_t>(count, 1U << 20U));
}

Slot CellStore3D::create(const CellInit& cell) {
    if (cell.uid == 0) {
        throw std::invalid_argument("cell uid must be nonzero");
    }
    const std::size_t new_stage = stage_index(cell.stage);
    // Validate before consuming a free slot so a rejected record cannot alter
    // the stable free-list order.
    validate_cell_times(cell);
    Slot slot = kEmptySlot;
    if (!free_slots_.empty()) {
        slot = free_slots_.back();
        free_slots_.pop_back();
        checkpoint_free_list_mutations_.push_back(
            {FreeListMutationKind3D::pop, slot});
        assign(slot, cell);
    } else {
        if (slot_count() >= static_cast<std::size_t>(kEmptySlot)) {
            throw std::overflow_error("CellStore3D slot space exhausted");
        }
        slot = static_cast<Slot>(slot_count());
        append(cell);
    }
    mark_checkpoint_dirty(slot);
    ++stage_counts_[new_stage];
    ++alive_count_;
    return slot;
}

void CellStore3D::erase(Slot slot) {
    if (!valid(slot)) {
        return;
    }
    const std::size_t old_stage = stage_index(static_cast<CellStage>(stage_[slot]));
    if (stage_counts_[old_stage] == 0) {
        throw std::logic_error("cell stage count underflow");
    }
    --stage_counts_[old_stage];
    alive_[slot] = 0;
    viability_[slot] = 0;
    free_slots_.push_back(slot);
    mark_checkpoint_dirty(slot);
    checkpoint_free_list_mutations_.push_back(
        {FreeListMutationKind3D::push, slot});
    --alive_count_;
}

bool CellStore3D::valid(Slot slot) const noexcept {
    return slot < alive_.size() && alive_[slot] != 0;
}

std::size_t CellStore3D::stage_count(CellStage stage) const {
    return stage_counts_[stage_index(stage)];
}

std::vector<Slot> CellStore3D::alive_slots() const {
    std::vector<Slot> result;
    result.reserve(alive_count_);
    for (std::size_t index = 0; index < alive_.size(); ++index) {
        if (alive_[index] != 0) {
            result.push_back(static_cast<Slot>(index));
        }
    }
    return result;
}

CellInit CellStore3D::snapshot(Slot slot) const {
    if (!valid(slot)) {
        throw std::out_of_range("cannot snapshot an inactive cell slot");
    }
    CellInit result;
    result.anchor = anchor(slot);
    result.uid = uid_[slot];
    result.parent_uid = parent_uid_[slot];
    result.clone_id = clone_id_[slot];
    result.type = static_cast<CellType>(type_[slot]);
    result.stage = static_cast<CellStage>(stage_[slot]);
    result.viability = viability_[slot];
    result.flags = flags_[slot];
    result.last_direction = last_direction_[slot];
    result.inherent_growth_rate = inherent_growth_rate_[slot];
    result.density_growth_rate = density_growth_rate_[slot];
    result.migration_rate = migration_rate_[slot];
    result.normal_migration_rate = normal_migration_rate_[slot];
    result.migration_activation_end_time = migration_activation_end_time_[slot];
    result.division_work_remaining = division_work_remaining_[slot];
    result.next_migration_time = next_migration_time_[slot];
    result.next_division_time = next_division_time_[slot];
    result.death_deadline = death_deadline_[slot];
    result.last_update_time = last_update_time_[slot];
    result.swap_ready_time = swap_ready_time(slot);
    result.swap_wait_state = swap_wait_state(slot);
    result.pending_swap_direction = pending_swap_direction_[slot];
    result.event_sequence = event_sequence_[slot];
    result.migration_schedule_generation = migration_schedule_generation_[slot];
    result.division_schedule_generation = division_schedule_generation_[slot];
    result.death_schedule_generation = death_schedule_generation_[slot];
    result.schedule_generation = std::max({result.migration_schedule_generation,
                                           result.division_schedule_generation,
                                           result.death_schedule_generation});
    return result;
}

void CellStore3D::mark_checkpoint_dirty(Slot slot) const {
    if (static_cast<std::size_t>(slot) >= checkpoint_dirty_.size()) {
        throw std::out_of_range("checkpoint dirty slot is out of range");
    }
    if (checkpoint_dirty_[slot] == 0) {
        checkpoint_dirty_[slot] = 1;
        checkpoint_dirty_slots_.push_back(slot);
    }
}

CheckpointCellJournal3D CellStore3D::take_checkpoint_journal() const {
    CheckpointCellJournal3D journal;
    journal.slot_count = slot_count();
    journal.mutations.reserve(checkpoint_dirty_slots_.size());
    for (const Slot slot : checkpoint_dirty_slots_) {
        CheckpointCellMutation3D mutation;
        mutation.slot = slot;
        mutation.alive = valid(slot);
        if (mutation.alive) mutation.cell = snapshot(slot);
        journal.mutations.push_back(std::move(mutation));
        checkpoint_dirty_[slot] = 0;
    }
    journal.free_list_mutations =
        std::move(checkpoint_free_list_mutations_);
    checkpoint_dirty_slots_.clear();
    checkpoint_free_list_mutations_.clear();
    return journal;
}

void CellStore3D::reset_checkpoint_journal() const {
    for (const Slot slot : checkpoint_dirty_slots_) {
        checkpoint_dirty_[slot] = 0;
    }
    checkpoint_dirty_slots_.clear();
    checkpoint_free_list_mutations_.clear();
}

double CellStore3D::set_migration_activation_end_time(Slot slot,
                                                       double value) {
    const float stored = quantize_time(
        value, TimeRounding::upward, "migration_activation_end_time");
    mark_checkpoint_dirty(slot);
    migration_activation_end_time_.at(slot) = stored;
    return static_cast<double>(stored);
}

double CellStore3D::set_next_migration_time(Slot slot, double value) {
    const float stored = quantize_time(
        value, TimeRounding::upward, "next_migration_time");
    mark_checkpoint_dirty(slot);
    next_migration_time_.at(slot) = stored;
    return static_cast<double>(stored);
}

double CellStore3D::set_next_division_time(Slot slot, double value) {
    const float stored = quantize_time(
        value, TimeRounding::upward, "next_division_time");
    mark_checkpoint_dirty(slot);
    next_division_time_.at(slot) = stored;
    return static_cast<double>(stored);
}

double CellStore3D::set_death_deadline(Slot slot, double value) {
    const float stored = quantize_time(
        value, TimeRounding::upward, "death_deadline");
    mark_checkpoint_dirty(slot);
    death_deadline_.at(slot) = stored;
    return static_cast<double>(stored);
}

double CellStore3D::set_last_update_time(Slot slot, double value) {
    const float stored = quantize_time(
        value, TimeRounding::nearest, "last_update_time");
    mark_checkpoint_dirty(slot);
    last_update_time_.at(slot) = stored;
    return static_cast<double>(stored);
}

double CellStore3D::set_swap_ready_time(Slot slot, double value) {
    const float stored = quantize_time(
        value, TimeRounding::upward, "swap_ready_time");
    mark_checkpoint_dirty(slot);
    next_migration_time_.at(slot) = stored;
    return static_cast<double>(stored);
}

void CellStore3D::restore_layout(std::size_t allocated_slot_count,
                                 const std::vector<Slot>& alive_slots,
                                 const std::vector<CellInit>& alive_cells,
                                 const std::vector<Slot>& free_slots) {
    if (slot_count() != 0 || alive_count_ != 0 || !free_slots_.empty()) {
        throw std::logic_error("CellStore3D layout restore requires a fresh store");
    }
    if (allocated_slot_count > static_cast<std::size_t>(kEmptySlot)) {
        throw std::overflow_error("CellStore3D restored slot space is too large");
    }
    if (alive_slots.size() != alive_cells.size()) {
        throw std::invalid_argument("restored alive slot and cell counts differ");
    }
    if (alive_slots.size() > allocated_slot_count ||
        free_slots.size() != allocated_slot_count - alive_slots.size()) {
        throw std::invalid_argument("restored alive/free slots do not cover the allocated layout");
    }

    // Validate the complete partition before mutating the store. A byte per
    // allocated slot keeps this check linear and bounded for checkpoint-scale
    // restores while detecting both duplicate and missing slots.
    std::vector<std::uint8_t> slot_state(allocated_slot_count, 0);
    for (std::size_t index = 0; index < alive_slots.size(); ++index) {
        const Slot slot = alive_slots[index];
        if (static_cast<std::size_t>(slot) >= allocated_slot_count) {
            throw std::invalid_argument("restored alive slot is out of range");
        }
        if (slot_state[slot] != 0) {
            throw std::invalid_argument("restored layout contains a duplicate slot");
        }
        if (alive_cells[index].uid == 0) {
            throw std::invalid_argument("restored cell uid must be nonzero");
        }
        (void)stage_index(alive_cells[index].stage);
        slot_state[slot] = 1;
    }
    for (const Slot slot : free_slots) {
        if (static_cast<std::size_t>(slot) >= allocated_slot_count) {
            throw std::invalid_argument("restored free slot is out of range");
        }
        if (slot_state[slot] != 0) {
            throw std::invalid_argument("restored layout contains a duplicate slot");
        }
        slot_state[slot] = 2;
    }
    if (std::find(slot_state.begin(), slot_state.end(), 0) != slot_state.end()) {
        throw std::invalid_argument("restored layout is missing an allocated slot");
    }

    // Build in a temporary store so validation/allocation failures leave this
    // fresh destination unchanged. Inactive payload bytes are intentionally
    // unspecified; only the slot topology and live-cell state are persistent.
    CellStore3D restored;
    restored.reserve(allocated_slot_count);
    CellInit placeholder;
    placeholder.uid = 1;
    for (std::size_t index = 0; index < allocated_slot_count; ++index) {
        restored.append(placeholder);
    }
    std::fill(restored.alive_.begin(), restored.alive_.end(), 0);

    for (std::size_t index = 0; index < alive_slots.size(); ++index) {
        const Slot slot = alive_slots[index];
        restored.assign(slot, alive_cells[index]);
        ++restored.stage_counts_[stage_index(alive_cells[index].stage)];
    }
    restored.alive_count_ = alive_cells.size();
    restored.free_slots_ = free_slots;
    restored.reset_checkpoint_journal();
    *this = std::move(restored);
}

Vec3i CellStore3D::anchor(Slot slot) const {
    return {x_.at(slot), y_.at(slot), z_.at(slot)};
}

void CellStore3D::set_anchor(Slot slot, Vec3i value) {
    mark_checkpoint_dirty(slot);
    x_.at(slot) = value.x;
    y_.at(slot) = value.y;
    z_.at(slot) = value.z;
}

void CellStore3D::set_stage(Slot slot, CellStage value) {
    if (!valid(slot)) {
        throw std::out_of_range("cannot change stage of an inactive cell slot");
    }
    const std::size_t old_stage = stage_index(static_cast<CellStage>(stage_[slot]));
    const std::size_t new_stage = stage_index(value);
    if (old_stage == new_stage) {
        return;
    }
    mark_checkpoint_dirty(slot);
    if (stage_counts_[old_stage] == 0) {
        throw std::logic_error("cell stage count underflow");
    }
    --stage_counts_[old_stage];
    ++stage_counts_[new_stage];
    stage_[slot] = static_cast<std::uint8_t>(value);
}

void CellStore3D::append(const CellInit& cell) {
    const float migration_activation_end_time = quantize_time(
        cell.migration_activation_end_time, TimeRounding::upward,
        "migration_activation_end_time");
    const float next_migration_time = quantize_time(
        cell.next_migration_time, TimeRounding::upward,
        "next_migration_time");
    const float next_division_time = quantize_time(
        cell.next_division_time, TimeRounding::upward,
        "next_division_time");
    const float death_deadline = quantize_time(
        cell.death_deadline, TimeRounding::upward, "death_deadline");
    const float last_update_time = quantize_time(
        cell.last_update_time, TimeRounding::nearest, "last_update_time");
    const float swap_ready_time = quantize_time(
        cell.swap_ready_time, TimeRounding::upward, "swap_ready_time");
    if (cell.swap_wait_state != 0 &&
        swap_ready_time != next_migration_time) {
        throw std::invalid_argument(
            "crowding-swap ready time must equal next migration time");
    }
    x_.push_back(cell.anchor.x); y_.push_back(cell.anchor.y); z_.push_back(cell.anchor.z);
    uid_.push_back(cell.uid); parent_uid_.push_back(cell.parent_uid); clone_id_.push_back(cell.clone_id);
    type_.push_back(static_cast<std::uint8_t>(cell.type)); stage_.push_back(static_cast<std::uint8_t>(cell.stage));
    viability_.push_back(cell.viability); flags_.push_back(cell.flags); last_direction_.push_back(cell.last_direction);
    alive_.push_back(1);
    checkpoint_dirty_.push_back(0);
    inherent_growth_rate_.push_back(cell.inherent_growth_rate);
    density_growth_rate_.push_back(cell.density_growth_rate);
    migration_rate_.push_back(cell.migration_rate);
    normal_migration_rate_.push_back(cell.normal_migration_rate);
    migration_activation_end_time_.push_back(migration_activation_end_time);
    division_work_remaining_.push_back(cell.division_work_remaining);
    next_migration_time_.push_back(next_migration_time);
    next_division_time_.push_back(next_division_time);
    death_deadline_.push_back(death_deadline);
    last_update_time_.push_back(last_update_time);
    pending_swap_direction_.push_back(cell.pending_swap_direction);
    event_sequence_.push_back(cell.event_sequence);
    migration_schedule_generation_.push_back(migration_generation(cell));
    division_schedule_generation_.push_back(division_generation(cell));
    death_schedule_generation_.push_back(death_generation(cell));
}

void CellStore3D::assign(Slot slot, const CellInit& cell) {
    const float migration_activation_end_time = quantize_time(
        cell.migration_activation_end_time, TimeRounding::upward,
        "migration_activation_end_time");
    const float next_migration_time = quantize_time(
        cell.next_migration_time, TimeRounding::upward,
        "next_migration_time");
    const float next_division_time = quantize_time(
        cell.next_division_time, TimeRounding::upward,
        "next_division_time");
    const float death_deadline = quantize_time(
        cell.death_deadline, TimeRounding::upward, "death_deadline");
    const float last_update_time = quantize_time(
        cell.last_update_time, TimeRounding::nearest, "last_update_time");
    const float swap_ready_time = quantize_time(
        cell.swap_ready_time, TimeRounding::upward, "swap_ready_time");
    if (cell.swap_wait_state != 0 &&
        swap_ready_time != next_migration_time) {
        throw std::invalid_argument(
            "crowding-swap ready time must equal next migration time");
    }
    x_[slot] = cell.anchor.x; y_[slot] = cell.anchor.y; z_[slot] = cell.anchor.z;
    uid_[slot] = cell.uid; parent_uid_[slot] = cell.parent_uid; clone_id_[slot] = cell.clone_id;
    type_[slot] = static_cast<std::uint8_t>(cell.type); stage_[slot] = static_cast<std::uint8_t>(cell.stage);
    viability_[slot] = cell.viability; flags_[slot] = cell.flags; last_direction_[slot] = cell.last_direction;
    alive_[slot] = 1;
    inherent_growth_rate_[slot] = cell.inherent_growth_rate;
    density_growth_rate_[slot] = cell.density_growth_rate;
    migration_rate_[slot] = cell.migration_rate;
    normal_migration_rate_[slot] = cell.normal_migration_rate;
    migration_activation_end_time_[slot] = migration_activation_end_time;
    division_work_remaining_[slot] = cell.division_work_remaining;
    next_migration_time_[slot] = next_migration_time;
    next_division_time_[slot] = next_division_time;
    death_deadline_[slot] = death_deadline;
    last_update_time_[slot] = last_update_time;
    pending_swap_direction_[slot] = cell.pending_swap_direction;
    event_sequence_[slot] = cell.event_sequence;
    migration_schedule_generation_[slot] = migration_generation(cell);
    division_schedule_generation_[slot] = division_generation(cell);
    death_schedule_generation_[slot] = death_generation(cell);
}

std::size_t CellStore3D::allocated_bytes() const noexcept {
    return vector_bytes(x_) + vector_bytes(y_) + vector_bytes(z_) +
           vector_bytes(uid_) + vector_bytes(parent_uid_) + vector_bytes(clone_id_) +
           vector_bytes(type_) + vector_bytes(stage_) + vector_bytes(viability_) + vector_bytes(flags_) +
           vector_bytes(last_direction_) + vector_bytes(alive_) + vector_bytes(inherent_growth_rate_) +
           vector_bytes(density_growth_rate_) + vector_bytes(migration_rate_) +
           vector_bytes(normal_migration_rate_) +
           vector_bytes(migration_activation_end_time_) +
           vector_bytes(division_work_remaining_) +
           vector_bytes(next_migration_time_) + vector_bytes(next_division_time_) +
           vector_bytes(death_deadline_) + vector_bytes(last_update_time_) +
           vector_bytes(pending_swap_direction_) +
           vector_bytes(event_sequence_) + vector_bytes(migration_schedule_generation_) +
           vector_bytes(division_schedule_generation_) +
           vector_bytes(death_schedule_generation_) + vector_bytes(free_slots_) +
           vector_bytes(checkpoint_dirty_) +
           vector_bytes(checkpoint_dirty_slots_) +
           vector_bytes(checkpoint_free_list_mutations_);
}

}  // namespace atcg3d
