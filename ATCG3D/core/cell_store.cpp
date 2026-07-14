#include "core/cell_store.hpp"

#include <algorithm>

namespace atcg3d {
namespace {

template <class T>
std::size_t vector_bytes(const std::vector<T>& values) {
    return values.capacity() * sizeof(T);
}

}  // namespace

void CellStore3D::reserve(std::size_t count) {
    x_.reserve(count); y_.reserve(count); z_.reserve(count);
    uid_.reserve(count); parent_uid_.reserve(count); clone_id_.reserve(count);
    type_.reserve(count); stage_.reserve(count); viability_.reserve(count);
    flags_.reserve(count); last_direction_.reserve(count); alive_.reserve(count);
    inherent_growth_rate_.reserve(count); density_growth_rate_.reserve(count); migration_rate_.reserve(count);
    next_migration_time_.reserve(count); next_division_time_.reserve(count);
    death_deadline_.reserve(count); last_update_time_.reserve(count);
    event_sequence_.reserve(count); schedule_generation_.reserve(count);
}

Slot CellStore3D::create(const CellInit& cell) {
    if (cell.uid == 0) {
        throw std::invalid_argument("cell uid must be nonzero");
    }
    Slot slot = kEmptySlot;
    if (!free_slots_.empty()) {
        slot = free_slots_.back();
        free_slots_.pop_back();
        assign(slot, cell);
    } else {
        if (slot_count() >= static_cast<std::size_t>(kEmptySlot)) {
            throw std::overflow_error("CellStore3D slot space exhausted");
        }
        slot = static_cast<Slot>(slot_count());
        append(cell);
    }
    ++alive_count_;
    return slot;
}

void CellStore3D::erase(Slot slot) {
    if (!valid(slot)) {
        return;
    }
    alive_[slot] = 0;
    viability_[slot] = 0;
    free_slots_.push_back(slot);
    --alive_count_;
}

bool CellStore3D::valid(Slot slot) const noexcept {
    return slot < alive_.size() && alive_[slot] != 0;
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
    result.next_migration_time = next_migration_time_[slot];
    result.next_division_time = next_division_time_[slot];
    result.death_deadline = death_deadline_[slot];
    result.last_update_time = last_update_time_[slot];
    result.event_sequence = event_sequence_[slot];
    result.schedule_generation = schedule_generation_[slot];
    return result;
}

Vec3i CellStore3D::anchor(Slot slot) const {
    return {x_.at(slot), y_.at(slot), z_.at(slot)};
}

void CellStore3D::set_anchor(Slot slot, Vec3i value) {
    x_.at(slot) = value.x;
    y_.at(slot) = value.y;
    z_.at(slot) = value.z;
}

void CellStore3D::append(const CellInit& cell) {
    x_.push_back(cell.anchor.x); y_.push_back(cell.anchor.y); z_.push_back(cell.anchor.z);
    uid_.push_back(cell.uid); parent_uid_.push_back(cell.parent_uid); clone_id_.push_back(cell.clone_id);
    type_.push_back(static_cast<std::uint8_t>(cell.type)); stage_.push_back(static_cast<std::uint8_t>(cell.stage));
    viability_.push_back(cell.viability); flags_.push_back(cell.flags); last_direction_.push_back(cell.last_direction);
    alive_.push_back(1);
    inherent_growth_rate_.push_back(cell.inherent_growth_rate);
    density_growth_rate_.push_back(cell.density_growth_rate);
    migration_rate_.push_back(cell.migration_rate);
    next_migration_time_.push_back(cell.next_migration_time);
    next_division_time_.push_back(cell.next_division_time);
    death_deadline_.push_back(static_cast<float>(cell.death_deadline));
    last_update_time_.push_back(cell.last_update_time);
    event_sequence_.push_back(cell.event_sequence);
    schedule_generation_.push_back(cell.schedule_generation);
}

void CellStore3D::assign(Slot slot, const CellInit& cell) {
    x_[slot] = cell.anchor.x; y_[slot] = cell.anchor.y; z_[slot] = cell.anchor.z;
    uid_[slot] = cell.uid; parent_uid_[slot] = cell.parent_uid; clone_id_[slot] = cell.clone_id;
    type_[slot] = static_cast<std::uint8_t>(cell.type); stage_[slot] = static_cast<std::uint8_t>(cell.stage);
    viability_[slot] = cell.viability; flags_[slot] = cell.flags; last_direction_[slot] = cell.last_direction;
    alive_[slot] = 1;
    inherent_growth_rate_[slot] = cell.inherent_growth_rate;
    density_growth_rate_[slot] = cell.density_growth_rate;
    migration_rate_[slot] = cell.migration_rate;
    next_migration_time_[slot] = cell.next_migration_time;
    next_division_time_[slot] = cell.next_division_time;
    death_deadline_[slot] = static_cast<float>(cell.death_deadline);
    last_update_time_[slot] = cell.last_update_time;
    event_sequence_[slot] = cell.event_sequence;
    schedule_generation_[slot] = cell.schedule_generation;
}

std::size_t CellStore3D::allocated_bytes() const noexcept {
    return vector_bytes(x_) + vector_bytes(y_) + vector_bytes(z_) +
           vector_bytes(uid_) + vector_bytes(parent_uid_) + vector_bytes(clone_id_) +
           vector_bytes(type_) + vector_bytes(stage_) + vector_bytes(viability_) + vector_bytes(flags_) +
           vector_bytes(last_direction_) + vector_bytes(alive_) + vector_bytes(inherent_growth_rate_) +
           vector_bytes(density_growth_rate_) + vector_bytes(migration_rate_) +
           vector_bytes(next_migration_time_) + vector_bytes(next_division_time_) +
           vector_bytes(death_deadline_) + vector_bytes(last_update_time_) +
           vector_bytes(event_sequence_) + vector_bytes(schedule_generation_) + vector_bytes(free_slots_);
}

}  // namespace atcg3d
