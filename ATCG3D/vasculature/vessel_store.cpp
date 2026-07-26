#include "vasculature/vessel_store.hpp"

#include <cmath>
#include <stdexcept>

namespace atcg3d {
namespace {

template <class T>
std::size_t vector_bytes(const std::vector<T>& values) noexcept {
    return values.capacity() * sizeof(T);
}

bool valid_role(VesselBranchRole role) noexcept {
    return role == VesselBranchRole::root || role == VesselBranchRole::inward ||
           role == VesselBranchRole::outward;
}

bool valid_status(VesselTipStatus status) noexcept {
    return status == VesselTipStatus::dormant || status == VesselTipStatus::active ||
           status == VesselTipStatus::blocked || status == VesselTipStatus::merged ||
           status == VesselTipStatus::reached_target ||
           status == VesselTipStatus::max_length ||
           status == VesselTipStatus::boundary_stop ||
           status == VesselTipStatus::complete ||
           status == VesselTipStatus::transiting;
}

void validate_node(const VesselNodeInit3D& node) {
    if (node.uid == 0 || node.vessel_id == 0) {
        throw std::invalid_argument("vessel node uid and vessel id must be nonzero");
    }
    if (!(node.diameter_voxels > 0.0F) || !std::isfinite(node.diameter_voxels)) {
        throw std::invalid_argument("vessel node diameter must be finite and positive");
    }
    if (node.created_time_hours < 0.0 || !std::isfinite(node.created_time_hours)) {
        throw std::invalid_argument("vessel node creation time must be finite and nonnegative");
    }
    if (!valid_role(node.role)) {
        throw std::invalid_argument("invalid vessel node branch role");
    }
    if (node.role == VesselBranchRole::root) {
        if (node.parent_uid != 0 || node.parent_node_slot != kEmptyVesselNodeSlot) {
            throw std::invalid_argument("a root vessel node cannot have a parent");
        }
    } else if (node.parent_uid == 0 || node.parent_node_slot == kEmptyVesselNodeSlot) {
        throw std::invalid_argument("a non-root vessel node requires parent uid and slot");
    }
}

void validate_tip(const VesselTipInit3D& tip) {
    if (tip.uid == 0 || tip.vessel_id == 0 || tip.current_node_uid == 0 ||
        tip.current_node_slot == kEmptyVesselNodeSlot) {
        throw std::invalid_argument("vessel tip ids and current node slot must be valid");
    }
    if (!(tip.diameter_voxels > 0.0F) || !std::isfinite(tip.diameter_voxels)) {
        throw std::invalid_argument("vessel tip diameter must be finite and positive");
    }
    if (tip.speed_voxels_per_hour < 0.0F || !std::isfinite(tip.speed_voxels_per_hour) ||
        tip.max_length_voxels < 0.0F || !std::isfinite(tip.max_length_voxels) ||
        tip.grown_length_voxels < 0.0F || !std::isfinite(tip.grown_length_voxels) ||
        tip.grown_length_voxels > tip.max_length_voxels) {
        throw std::invalid_argument("vessel tip lengths and speed are invalid");
    }
    if (tip.next_growth_time < 0.0 || !std::isfinite(tip.next_growth_time)) {
        throw std::invalid_argument("vessel tip next growth time must be finite and nonnegative");
    }
    if (!valid_role(tip.role) || !valid_status(tip.status) || tip.last_direction > 26 ||
        tip.pending_direction > 26) {
        throw std::invalid_argument("invalid vessel tip enum or direction");
    }
    if (vessel_tip_growing(tip.status) &&
        (!(tip.speed_voxels_per_hour > 0.0F) || !(tip.max_length_voxels > 0.0F) ||
         tip.grown_length_voxels >= tip.max_length_voxels)) {
        throw std::invalid_argument("an active vessel tip must have positive speed and remaining length");
    }
}

}  // namespace

void VesselNodeStore3D::reserve(std::size_t count) {
    x_.reserve(count);
    y_.reserve(count);
    z_.reserve(count);
    uid_.reserve(count);
    parent_uid_.reserve(count);
    parent_node_slot_.reserve(count);
    vessel_id_.reserve(count);
    source_lesion_id_.reserve(count);
    role_.reserve(count);
    perfused_.reserve(count);
    alive_.reserve(count);
    diameter_voxels_.reserve(count);
    created_time_hours_.reserve(count);
}

VesselNodeSlot VesselNodeStore3D::create(const VesselNodeInit3D& node) {
    validate_node(node);
    VesselNodeSlot slot{};
    if (free_slots_.empty()) {
        slot = static_cast<VesselNodeSlot>(alive_.size());
        append(node);
    } else {
        slot = free_slots_.back();
        free_slots_.pop_back();
        assign(slot, node);
    }
    ++alive_count_;
    return slot;
}

void VesselNodeStore3D::erase(VesselNodeSlot slot) {
    if (!valid(slot)) {
        throw std::out_of_range("invalid vessel node slot");
    }
    alive_[slot] = 0;
    free_slots_.push_back(slot);
    --alive_count_;
}

bool VesselNodeStore3D::valid(VesselNodeSlot slot) const noexcept {
    return slot < alive_.size() && alive_[slot] != 0;
}

std::vector<VesselNodeSlot> VesselNodeStore3D::alive_slots() const {
    std::vector<VesselNodeSlot> result;
    result.reserve(alive_count_);
    for (VesselNodeSlot slot = 0; slot < alive_.size(); ++slot) {
        if (alive_[slot] != 0) result.push_back(slot);
    }
    return result;
}

VesselNodeInit3D VesselNodeStore3D::snapshot(VesselNodeSlot slot) const {
    if (!valid(slot)) throw std::out_of_range("invalid vessel node slot");
    return {position(slot), uid(slot), parent_uid(slot), parent_node_slot(slot), vessel_id(slot),
            source_lesion_id(slot), role(slot), perfused(slot), diameter_voxels(slot),
            created_time_hours(slot)};
}

Vec3i VesselNodeStore3D::position(VesselNodeSlot slot) const {
    if (!valid(slot)) throw std::out_of_range("invalid vessel node slot");
    return {x_.at(slot), y_.at(slot), z_.at(slot)};
}

void VesselNodeStore3D::append(const VesselNodeInit3D& node) {
    x_.push_back(node.position.x);
    y_.push_back(node.position.y);
    z_.push_back(node.position.z);
    uid_.push_back(node.uid);
    parent_uid_.push_back(node.parent_uid);
    parent_node_slot_.push_back(node.parent_node_slot);
    vessel_id_.push_back(node.vessel_id);
    source_lesion_id_.push_back(node.source_lesion_id);
    role_.push_back(static_cast<std::uint8_t>(node.role));
    perfused_.push_back(node.perfused ? 1U : 0U);
    alive_.push_back(1U);
    diameter_voxels_.push_back(node.diameter_voxels);
    created_time_hours_.push_back(node.created_time_hours);
}

void VesselNodeStore3D::assign(VesselNodeSlot slot, const VesselNodeInit3D& node) {
    x_.at(slot) = node.position.x;
    y_.at(slot) = node.position.y;
    z_.at(slot) = node.position.z;
    uid_.at(slot) = node.uid;
    parent_uid_.at(slot) = node.parent_uid;
    parent_node_slot_.at(slot) = node.parent_node_slot;
    vessel_id_.at(slot) = node.vessel_id;
    source_lesion_id_.at(slot) = node.source_lesion_id;
    role_.at(slot) = static_cast<std::uint8_t>(node.role);
    perfused_.at(slot) = node.perfused ? 1U : 0U;
    alive_.at(slot) = 1U;
    diameter_voxels_.at(slot) = node.diameter_voxels;
    created_time_hours_.at(slot) = node.created_time_hours;
}

std::size_t VesselNodeStore3D::allocated_bytes() const noexcept {
    return vector_bytes(x_) + vector_bytes(y_) + vector_bytes(z_) + vector_bytes(uid_) +
           vector_bytes(parent_uid_) + vector_bytes(parent_node_slot_) + vector_bytes(vessel_id_) +
           vector_bytes(source_lesion_id_) +
           vector_bytes(role_) + vector_bytes(perfused_) + vector_bytes(alive_) +
           vector_bytes(diameter_voxels_) + vector_bytes(created_time_hours_) +
           vector_bytes(free_slots_);
}

void VesselTipStore3D::reserve(std::size_t count) {
    x_.reserve(count);
    y_.reserve(count);
    z_.reserve(count);
    bias_x_.reserve(count);
    bias_y_.reserve(count);
    bias_z_.reserve(count);
    target_x_.reserve(count);
    target_y_.reserve(count);
    target_z_.reserve(count);
    uid_.reserve(count);
    vessel_id_.reserve(count);
    source_lesion_id_.reserve(count);
    current_node_uid_.reserve(count);
    current_node_slot_.reserve(count);
    role_.reserve(count);
    status_.reserve(count);
    perfused_.reserve(count);
    last_direction_.reserve(count);
    pending_direction_.reserve(count);
    alive_.reserve(count);
    diameter_voxels_.reserve(count);
    speed_voxels_per_hour_.reserve(count);
    max_length_voxels_.reserve(count);
    grown_length_voxels_.reserve(count);
    next_growth_time_.reserve(count);
    event_sequence_.reserve(count);
    schedule_generation_.reserve(count);
}

VesselTipSlot VesselTipStore3D::create(const VesselTipInit3D& tip) {
    validate_tip(tip);
    VesselTipSlot slot{};
    if (free_slots_.empty()) {
        slot = static_cast<VesselTipSlot>(alive_.size());
        append(tip);
    } else {
        slot = free_slots_.back();
        free_slots_.pop_back();
        assign(slot, tip);
    }
    ++alive_count_;
    return slot;
}

void VesselTipStore3D::erase(VesselTipSlot slot) {
    if (!valid(slot)) throw std::out_of_range("invalid vessel tip slot");
    alive_[slot] = 0;
    free_slots_.push_back(slot);
    --alive_count_;
}

bool VesselTipStore3D::valid(VesselTipSlot slot) const noexcept {
    return slot < alive_.size() && alive_[slot] != 0;
}

std::vector<VesselTipSlot> VesselTipStore3D::alive_slots() const {
    std::vector<VesselTipSlot> result;
    result.reserve(alive_count_);
    for (VesselTipSlot slot = 0; slot < alive_.size(); ++slot) {
        if (alive_[slot] != 0) result.push_back(slot);
    }
    return result;
}

VesselTipInit3D VesselTipStore3D::snapshot(VesselTipSlot slot) const {
    if (!valid(slot)) throw std::out_of_range("invalid vessel tip slot");
    VesselTipInit3D tip;
    tip.position = position(slot);
    tip.bias_axis = bias_axis(slot);
    tip.target = target(slot);
    tip.uid = uid(slot);
    tip.vessel_id = vessel_id(slot);
    tip.source_lesion_id = source_lesion_id(slot);
    tip.current_node_uid = current_node_uid(slot);
    tip.current_node_slot = current_node_slot(slot);
    tip.role = role(slot);
    tip.status = status(slot);
    tip.perfused = perfused(slot);
    tip.last_direction = last_direction(slot);
    tip.pending_direction = pending_direction(slot);
    tip.diameter_voxels = diameter_voxels(slot);
    tip.speed_voxels_per_hour = speed_voxels_per_hour(slot);
    tip.max_length_voxels = max_length_voxels(slot);
    tip.grown_length_voxels = grown_length_voxels(slot);
    tip.next_growth_time = next_growth_time(slot);
    tip.event_sequence = event_sequence(slot);
    tip.schedule_generation = schedule_generation(slot);
    return tip;
}

Vec3i VesselTipStore3D::position(VesselTipSlot slot) const {
    if (!valid(slot)) throw std::out_of_range("invalid vessel tip slot");
    return {x_.at(slot), y_.at(slot), z_.at(slot)};
}

Vec3i VesselTipStore3D::bias_axis(VesselTipSlot slot) const {
    if (!valid(slot)) throw std::out_of_range("invalid vessel tip slot");
    return {bias_x_.at(slot), bias_y_.at(slot), bias_z_.at(slot)};
}

Vec3i VesselTipStore3D::target(VesselTipSlot slot) const {
    if (!valid(slot)) throw std::out_of_range("invalid vessel tip slot");
    return {target_x_.at(slot), target_y_.at(slot), target_z_.at(slot)};
}

void VesselTipStore3D::set_position(VesselTipSlot slot, Vec3i value) {
    if (!valid(slot)) throw std::out_of_range("invalid vessel tip slot");
    x_.at(slot) = value.x;
    y_.at(slot) = value.y;
    z_.at(slot) = value.z;
}

void VesselTipStore3D::set_bias_axis(VesselTipSlot slot, Vec3i value) {
    if (!valid(slot)) throw std::out_of_range("invalid vessel tip slot");
    bias_x_.at(slot) = value.x;
    bias_y_.at(slot) = value.y;
    bias_z_.at(slot) = value.z;
}

void VesselTipStore3D::set_target(VesselTipSlot slot, Vec3i value) {
    if (!valid(slot)) throw std::out_of_range("invalid vessel tip slot");
    target_x_.at(slot) = value.x;
    target_y_.at(slot) = value.y;
    target_z_.at(slot) = value.z;
}

void VesselTipStore3D::set_grown_length_voxels(VesselTipSlot slot, float value) {
    if (!std::isfinite(value) || value < 0.0F || value > max_length_voxels(slot)) {
        throw std::invalid_argument("invalid vessel tip grown length");
    }
    grown_length_voxels_.at(slot) = value;
}

void VesselTipStore3D::set_next_growth_time(VesselTipSlot slot, double value) {
    if (!std::isfinite(value) || value < 0.0) {
        throw std::invalid_argument("invalid vessel tip growth time");
    }
    next_growth_time_.at(slot) = value;
}

void VesselTipStore3D::append(const VesselTipInit3D& tip) {
    x_.push_back(tip.position.x);
    y_.push_back(tip.position.y);
    z_.push_back(tip.position.z);
    bias_x_.push_back(tip.bias_axis.x);
    bias_y_.push_back(tip.bias_axis.y);
    bias_z_.push_back(tip.bias_axis.z);
    target_x_.push_back(tip.target.x);
    target_y_.push_back(tip.target.y);
    target_z_.push_back(tip.target.z);
    uid_.push_back(tip.uid);
    vessel_id_.push_back(tip.vessel_id);
    source_lesion_id_.push_back(tip.source_lesion_id);
    current_node_uid_.push_back(tip.current_node_uid);
    current_node_slot_.push_back(tip.current_node_slot);
    role_.push_back(static_cast<std::uint8_t>(tip.role));
    status_.push_back(static_cast<std::uint8_t>(tip.status));
    perfused_.push_back(tip.perfused ? 1U : 0U);
    last_direction_.push_back(tip.last_direction);
    pending_direction_.push_back(tip.pending_direction);
    alive_.push_back(1U);
    diameter_voxels_.push_back(tip.diameter_voxels);
    speed_voxels_per_hour_.push_back(tip.speed_voxels_per_hour);
    max_length_voxels_.push_back(tip.max_length_voxels);
    grown_length_voxels_.push_back(tip.grown_length_voxels);
    next_growth_time_.push_back(tip.next_growth_time);
    event_sequence_.push_back(tip.event_sequence);
    schedule_generation_.push_back(tip.schedule_generation);
}

void VesselTipStore3D::assign(VesselTipSlot slot, const VesselTipInit3D& tip) {
    x_.at(slot) = tip.position.x;
    y_.at(slot) = tip.position.y;
    z_.at(slot) = tip.position.z;
    bias_x_.at(slot) = tip.bias_axis.x;
    bias_y_.at(slot) = tip.bias_axis.y;
    bias_z_.at(slot) = tip.bias_axis.z;
    target_x_.at(slot) = tip.target.x;
    target_y_.at(slot) = tip.target.y;
    target_z_.at(slot) = tip.target.z;
    uid_.at(slot) = tip.uid;
    vessel_id_.at(slot) = tip.vessel_id;
    source_lesion_id_.at(slot) = tip.source_lesion_id;
    current_node_uid_.at(slot) = tip.current_node_uid;
    current_node_slot_.at(slot) = tip.current_node_slot;
    role_.at(slot) = static_cast<std::uint8_t>(tip.role);
    status_.at(slot) = static_cast<std::uint8_t>(tip.status);
    perfused_.at(slot) = tip.perfused ? 1U : 0U;
    last_direction_.at(slot) = tip.last_direction;
    pending_direction_.at(slot) = tip.pending_direction;
    alive_.at(slot) = 1U;
    diameter_voxels_.at(slot) = tip.diameter_voxels;
    speed_voxels_per_hour_.at(slot) = tip.speed_voxels_per_hour;
    max_length_voxels_.at(slot) = tip.max_length_voxels;
    grown_length_voxels_.at(slot) = tip.grown_length_voxels;
    next_growth_time_.at(slot) = tip.next_growth_time;
    event_sequence_.at(slot) = tip.event_sequence;
    schedule_generation_.at(slot) = tip.schedule_generation;
}

std::size_t VesselTipStore3D::allocated_bytes() const noexcept {
    return vector_bytes(x_) + vector_bytes(y_) + vector_bytes(z_) + vector_bytes(bias_x_) +
           vector_bytes(bias_y_) + vector_bytes(bias_z_) + vector_bytes(target_x_) +
           vector_bytes(target_y_) + vector_bytes(target_z_) + vector_bytes(uid_) +
           vector_bytes(vessel_id_) + vector_bytes(current_node_uid_) +
           vector_bytes(source_lesion_id_) +
           vector_bytes(current_node_slot_) + vector_bytes(role_) +
           vector_bytes(status_) + vector_bytes(perfused_) + vector_bytes(last_direction_) +
           vector_bytes(pending_direction_) + vector_bytes(alive_) +
           vector_bytes(diameter_voxels_) + vector_bytes(speed_voxels_per_hour_) +
           vector_bytes(max_length_voxels_) + vector_bytes(grown_length_voxels_) +
           vector_bytes(next_growth_time_) + vector_bytes(event_sequence_) +
           vector_bytes(schedule_generation_) + vector_bytes(free_slots_);
}

}  // namespace atcg3d
