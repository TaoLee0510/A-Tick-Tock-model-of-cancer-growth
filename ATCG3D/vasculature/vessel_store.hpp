#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "vasculature/types.hpp"

namespace atcg3d {

struct VesselNodeInit3D {
    Vec3i position{};
    VesselNodeUid uid{};
    VesselNodeUid parent_uid{};
    VesselNodeSlot parent_node_slot{kEmptyVesselNodeSlot};
    VesselId vessel_id{};
    LesionId source_lesion_id{kNoLesionId};
    VesselBranchRole role{VesselBranchRole::root};
    bool perfused{};
    float diameter_voxels{1.0F};
    double created_time_hours{};
};

class VesselNodeStore3D {
public:
    void reserve(std::size_t count);
    VesselNodeSlot create(const VesselNodeInit3D& node);
    void erase(VesselNodeSlot slot);

    bool valid(VesselNodeSlot slot) const noexcept;
    std::size_t slot_count() const noexcept { return alive_.size(); }
    std::size_t alive_count() const noexcept { return alive_count_; }
    std::size_t free_count() const noexcept { return free_slots_.size(); }
    std::vector<VesselNodeSlot> alive_slots() const;
    VesselNodeInit3D snapshot(VesselNodeSlot slot) const;

    Vec3i position(VesselNodeSlot slot) const;
    VesselNodeUid uid(VesselNodeSlot slot) const { return uid_.at(slot); }
    VesselNodeUid parent_uid(VesselNodeSlot slot) const { return parent_uid_.at(slot); }
    VesselNodeSlot parent_node_slot(VesselNodeSlot slot) const {
        return parent_node_slot_.at(slot);
    }
    VesselId vessel_id(VesselNodeSlot slot) const { return vessel_id_.at(slot); }
    LesionId source_lesion_id(VesselNodeSlot slot) const {
        return source_lesion_id_.at(slot);
    }
    VesselBranchRole role(VesselNodeSlot slot) const {
        return static_cast<VesselBranchRole>(role_.at(slot));
    }
    bool perfused(VesselNodeSlot slot) const { return perfused_.at(slot) != 0; }
    float diameter_voxels(VesselNodeSlot slot) const { return diameter_voxels_.at(slot); }
    double created_time_hours(VesselNodeSlot slot) const { return created_time_hours_.at(slot); }

    void set_perfused(VesselNodeSlot slot, bool value) { perfused_.at(slot) = value ? 1U : 0U; }

    std::size_t allocated_bytes() const noexcept;

private:
    void append(const VesselNodeInit3D& node);
    void assign(VesselNodeSlot slot, const VesselNodeInit3D& node);

    std::vector<std::int32_t> x_;
    std::vector<std::int32_t> y_;
    std::vector<std::int32_t> z_;
    std::vector<VesselNodeUid> uid_;
    std::vector<VesselNodeUid> parent_uid_;
    std::vector<VesselNodeSlot> parent_node_slot_;
    std::vector<VesselId> vessel_id_;
    std::vector<LesionId> source_lesion_id_;
    std::vector<std::uint8_t> role_;
    std::vector<std::uint8_t> perfused_;
    std::vector<std::uint8_t> alive_;
    std::vector<float> diameter_voxels_;
    std::vector<double> created_time_hours_;
    std::vector<VesselNodeSlot> free_slots_;
    std::size_t alive_count_{};
};

struct VesselTipInit3D {
    Vec3i position{};
    Vec3i bias_axis{};
    Vec3i target{};
    VesselTipUid uid{};
    VesselId vessel_id{};
    LesionId source_lesion_id{kNoLesionId};
    VesselNodeUid current_node_uid{};
    VesselNodeSlot current_node_slot{kEmptyVesselNodeSlot};
    VesselBranchRole role{VesselBranchRole::inward};
    VesselTipStatus status{VesselTipStatus::dormant};
    bool perfused{};
    DirectionId last_direction{kStayDirection};
    DirectionId pending_direction{kStayDirection};
    float diameter_voxels{1.0F};
    float speed_voxels_per_hour{};
    float max_length_voxels{};
    float grown_length_voxels{};
    double next_growth_time{};
    std::uint64_t event_sequence{};
    std::uint32_t schedule_generation{};
};

class VesselTipStore3D {
public:
    void reserve(std::size_t count);
    VesselTipSlot create(const VesselTipInit3D& tip);
    void erase(VesselTipSlot slot);

    bool valid(VesselTipSlot slot) const noexcept;
    std::size_t slot_count() const noexcept { return alive_.size(); }
    std::size_t alive_count() const noexcept { return alive_count_; }
    std::size_t free_count() const noexcept { return free_slots_.size(); }
    std::vector<VesselTipSlot> alive_slots() const;
    VesselTipInit3D snapshot(VesselTipSlot slot) const;

    Vec3i position(VesselTipSlot slot) const;
    Vec3i bias_axis(VesselTipSlot slot) const;
    Vec3i target(VesselTipSlot slot) const;
    VesselTipUid uid(VesselTipSlot slot) const { return uid_.at(slot); }
    VesselId vessel_id(VesselTipSlot slot) const { return vessel_id_.at(slot); }
    LesionId source_lesion_id(VesselTipSlot slot) const {
        return source_lesion_id_.at(slot);
    }
    VesselNodeUid current_node_uid(VesselTipSlot slot) const { return current_node_uid_.at(slot); }
    VesselNodeSlot current_node_slot(VesselTipSlot slot) const {
        return current_node_slot_.at(slot);
    }
    VesselBranchRole role(VesselTipSlot slot) const {
        return static_cast<VesselBranchRole>(role_.at(slot));
    }
    VesselTipStatus status(VesselTipSlot slot) const {
        return static_cast<VesselTipStatus>(status_.at(slot));
    }
    bool perfused(VesselTipSlot slot) const { return perfused_.at(slot) != 0; }
    DirectionId last_direction(VesselTipSlot slot) const { return last_direction_.at(slot); }
    DirectionId pending_direction(VesselTipSlot slot) const { return pending_direction_.at(slot); }
    float diameter_voxels(VesselTipSlot slot) const { return diameter_voxels_.at(slot); }
    float speed_voxels_per_hour(VesselTipSlot slot) const { return speed_voxels_per_hour_.at(slot); }
    float max_length_voxels(VesselTipSlot slot) const { return max_length_voxels_.at(slot); }
    float grown_length_voxels(VesselTipSlot slot) const { return grown_length_voxels_.at(slot); }
    double next_growth_time(VesselTipSlot slot) const { return next_growth_time_.at(slot); }
    std::uint64_t event_sequence(VesselTipSlot slot) const { return event_sequence_.at(slot); }
    std::uint32_t schedule_generation(VesselTipSlot slot) const {
        return schedule_generation_.at(slot);
    }

    void set_position(VesselTipSlot slot, Vec3i value);
    void set_bias_axis(VesselTipSlot slot, Vec3i value);
    void set_target(VesselTipSlot slot, Vec3i value);
    void set_current_node_uid(VesselTipSlot slot, VesselNodeUid value) {
        current_node_uid_.at(slot) = value;
    }
    void set_current_node_slot(VesselTipSlot slot, VesselNodeSlot value) {
        current_node_slot_.at(slot) = value;
    }
    void set_status(VesselTipSlot slot, VesselTipStatus value) {
        status_.at(slot) = static_cast<std::uint8_t>(value);
    }
    void set_perfused(VesselTipSlot slot, bool value) { perfused_.at(slot) = value ? 1U : 0U; }
    void set_last_direction(VesselTipSlot slot, DirectionId value) { last_direction_.at(slot) = value; }
    void set_pending_direction(VesselTipSlot slot, DirectionId value) {
        pending_direction_.at(slot) = value;
    }
    void set_grown_length_voxels(VesselTipSlot slot, float value);
    void set_next_growth_time(VesselTipSlot slot, double value);
    void set_event_sequence(VesselTipSlot slot, std::uint64_t value) {
        event_sequence_.at(slot) = value;
    }
    std::uint64_t consume_event_sequence(VesselTipSlot slot) {
        return event_sequence_.at(slot)++;
    }
    void set_schedule_generation(VesselTipSlot slot, std::uint32_t value) {
        schedule_generation_.at(slot) = value;
    }
    std::uint32_t bump_schedule_generation(VesselTipSlot slot) {
        return ++schedule_generation_.at(slot);
    }

    std::size_t allocated_bytes() const noexcept;

private:
    void append(const VesselTipInit3D& tip);
    void assign(VesselTipSlot slot, const VesselTipInit3D& tip);

    std::vector<std::int32_t> x_;
    std::vector<std::int32_t> y_;
    std::vector<std::int32_t> z_;
    std::vector<std::int32_t> bias_x_;
    std::vector<std::int32_t> bias_y_;
    std::vector<std::int32_t> bias_z_;
    std::vector<std::int32_t> target_x_;
    std::vector<std::int32_t> target_y_;
    std::vector<std::int32_t> target_z_;
    std::vector<VesselTipUid> uid_;
    std::vector<VesselId> vessel_id_;
    std::vector<LesionId> source_lesion_id_;
    std::vector<VesselNodeUid> current_node_uid_;
    std::vector<VesselNodeSlot> current_node_slot_;
    std::vector<std::uint8_t> role_;
    std::vector<std::uint8_t> status_;
    std::vector<std::uint8_t> perfused_;
    std::vector<DirectionId> last_direction_;
    std::vector<DirectionId> pending_direction_;
    std::vector<std::uint8_t> alive_;
    std::vector<float> diameter_voxels_;
    std::vector<float> speed_voxels_per_hour_;
    std::vector<float> max_length_voxels_;
    std::vector<float> grown_length_voxels_;
    std::vector<double> next_growth_time_;
    std::vector<std::uint64_t> event_sequence_;
    std::vector<std::uint32_t> schedule_generation_;
    std::vector<VesselTipSlot> free_slots_;
    std::size_t alive_count_{};
};

}  // namespace atcg3d
