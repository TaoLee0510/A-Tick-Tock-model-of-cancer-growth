#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "core/types.hpp"

namespace atcg3d {

struct CellInit {
    Vec3i anchor{};
    CellUid uid{};
    CellUid parent_uid{};
    std::uint32_t clone_id{};
    CellType type{CellType::r};
    CellStage stage{CellStage::small};
    std::uint8_t viability{1};
    std::uint8_t flags{kDirtyDensity};
    DirectionId last_direction{kStayDirection};
    float inherent_growth_rate{1.0F};
    float density_growth_rate{1.0F};
    // Inherent/base rate used only while density-triggered migration is active.
    float migration_rate{0.0F};
    // Ordinary low-density rate. For K cells this equals migration_rate; for r
    // cells it is a separate configured Beta draw.
    float normal_migration_rate{0.0F};
    // Zero in normal state, otherwise the exact event time at which the finite
    // activation interval ends.
    double migration_activation_end_time{};
    // Remaining growth work for the current cell cycle, in unit-rate hours.
    // Density refreshes integrate this value lazily instead of drawing a new
    // division delay.
    float division_work_remaining{};
    double next_migration_time{};
    double next_division_time{24.0};
    double death_deadline{};
    double last_update_time{};
    std::uint64_t event_sequence{};
    // Kept as a checkpoint-v1 compatibility value. New in-memory scheduling
    // uses the three event-specific generations below.
    std::uint32_t schedule_generation{};
    std::uint32_t migration_schedule_generation{};
    std::uint32_t division_schedule_generation{};
    std::uint32_t death_schedule_generation{};
};

class CellStore3D {
public:
    void reserve(std::size_t count);
    Slot create(const CellInit& cell);
    void erase(Slot slot);

    bool valid(Slot slot) const noexcept;
    std::size_t slot_count() const noexcept { return alive_.size(); }
    std::size_t alive_count() const noexcept { return alive_count_; }
    std::size_t stage_count(CellStage stage) const;
    std::size_t free_count() const noexcept { return free_slots_.size(); }
    const std::vector<Slot>& free_slots() const noexcept { return free_slots_; }
    std::vector<Slot> alive_slots() const;
    CellInit snapshot(Slot slot) const;
    void restore_layout(std::size_t allocated_slot_count,
                        const std::vector<Slot>& alive_slots,
                        const std::vector<CellInit>& alive_cells,
                        const std::vector<Slot>& free_slots);

    Vec3i anchor(Slot slot) const;
    void set_anchor(Slot slot, Vec3i value);

    CellUid uid(Slot slot) const { return uid_.at(slot); }
    CellUid parent_uid(Slot slot) const { return parent_uid_.at(slot); }
    std::uint32_t clone_id(Slot slot) const { return clone_id_.at(slot); }
    CellType type(Slot slot) const { return static_cast<CellType>(type_.at(slot)); }
    CellStage stage(Slot slot) const { return static_cast<CellStage>(stage_.at(slot)); }
    std::uint8_t viability(Slot slot) const { return viability_.at(slot); }
    std::uint8_t flags(Slot slot) const { return flags_.at(slot); }
    DirectionId last_direction(Slot slot) const { return last_direction_.at(slot); }
    float inherent_growth_rate(Slot slot) const { return inherent_growth_rate_.at(slot); }
    float density_growth_rate(Slot slot) const { return density_growth_rate_.at(slot); }
    float migration_rate(Slot slot) const { return migration_rate_.at(slot); }
    float normal_migration_rate(Slot slot) const {
        return normal_migration_rate_.at(slot);
    }
    double migration_activation_end_time(Slot slot) const {
        return static_cast<double>(migration_activation_end_time_.at(slot));
    }
    float division_work_remaining(Slot slot) const {
        return division_work_remaining_.at(slot);
    }
    double next_migration_time(Slot slot) const {
        return static_cast<double>(next_migration_time_.at(slot));
    }
    double next_division_time(Slot slot) const {
        return static_cast<double>(next_division_time_.at(slot));
    }
    double death_deadline(Slot slot) const { return static_cast<double>(death_deadline_.at(slot)); }
    double last_update_time(Slot slot) const {
        return static_cast<double>(last_update_time_.at(slot));
    }
    std::uint64_t event_sequence(Slot slot) const { return event_sequence_.at(slot); }
    std::uint32_t migration_schedule_generation(Slot slot) const {
        return migration_schedule_generation_.at(slot);
    }
    std::uint32_t division_schedule_generation(Slot slot) const {
        return division_schedule_generation_.at(slot);
    }
    std::uint32_t death_schedule_generation(Slot slot) const {
        return death_schedule_generation_.at(slot);
    }

    void set_parent_uid(Slot slot, CellUid value) { parent_uid_.at(slot) = value; }
    void set_clone_id(Slot slot, std::uint32_t value) { clone_id_.at(slot) = value; }
    void set_type(Slot slot, CellType value) { type_.at(slot) = static_cast<std::uint8_t>(value); }
    void set_stage(Slot slot, CellStage value);
    void set_viability(Slot slot, std::uint8_t value) { viability_.at(slot) = value; }
    void set_flags(Slot slot, std::uint8_t value) { flags_.at(slot) = value; }
    void set_last_direction(Slot slot, DirectionId value) { last_direction_.at(slot) = value; }
    void set_inherent_growth_rate(Slot slot, float value) { inherent_growth_rate_.at(slot) = value; }
    void set_density_growth_rate(Slot slot, float value) { density_growth_rate_.at(slot) = value; }
    void set_migration_rate(Slot slot, float value) { migration_rate_.at(slot) = value; }
    void set_normal_migration_rate(Slot slot, float value) {
        normal_migration_rate_.at(slot) = value;
    }
    double set_migration_activation_end_time(Slot slot, double value);
    void set_division_work_remaining(Slot slot, float value) {
        division_work_remaining_.at(slot) = value;
    }
    double set_next_migration_time(Slot slot, double value);
    double set_next_division_time(Slot slot, double value);
    double set_death_deadline(Slot slot, double value);
    double set_last_update_time(Slot slot, double value);
    void set_event_sequence(Slot slot, std::uint64_t value) { event_sequence_.at(slot) = value; }
    std::uint64_t consume_event_sequence(Slot slot) { return event_sequence_.at(slot)++; }
    void set_migration_schedule_generation(Slot slot, std::uint32_t value) {
        migration_schedule_generation_.at(slot) = value;
    }
    void set_division_schedule_generation(Slot slot, std::uint32_t value) {
        division_schedule_generation_.at(slot) = value;
    }
    void set_death_schedule_generation(Slot slot, std::uint32_t value) {
        death_schedule_generation_.at(slot) = value;
    }
    std::uint32_t bump_migration_schedule_generation(Slot slot) {
        return ++migration_schedule_generation_.at(slot);
    }
    std::uint32_t bump_division_schedule_generation(Slot slot) {
        return ++division_schedule_generation_.at(slot);
    }
    std::uint32_t bump_death_schedule_generation(Slot slot) {
        return ++death_schedule_generation_.at(slot);
    }

    std::size_t allocated_bytes() const noexcept;
    static constexpr std::size_t logical_bytes_per_slot() noexcept {
        return sizeof(std::int32_t) * 3 + sizeof(CellUid) * 2 + sizeof(std::uint32_t) * 4 +
               sizeof(std::uint8_t) * 6 + sizeof(float) * 10 +
               sizeof(std::uint64_t);
    }

private:
    void append(const CellInit& cell);
    void assign(Slot slot, const CellInit& cell);

    std::vector<std::int32_t> x_;
    std::vector<std::int32_t> y_;
    std::vector<std::int32_t> z_;
    std::vector<CellUid> uid_;
    std::vector<CellUid> parent_uid_;
    std::vector<std::uint32_t> clone_id_;
    std::vector<std::uint8_t> type_;
    std::vector<std::uint8_t> stage_;
    std::vector<std::uint8_t> viability_;
    std::vector<std::uint8_t> flags_;
    std::vector<DirectionId> last_direction_;
    std::vector<std::uint8_t> alive_;
    std::vector<float> inherent_growth_rate_;
    std::vector<float> density_growth_rate_;
    std::vector<float> migration_rate_;
    std::vector<float> normal_migration_rate_;
    std::vector<float> migration_activation_end_time_;
    std::vector<float> division_work_remaining_;
    std::vector<float> next_migration_time_;
    std::vector<float> next_division_time_;
    std::vector<float> death_deadline_;
    std::vector<float> last_update_time_;
    std::vector<std::uint64_t> event_sequence_;
    std::vector<std::uint32_t> migration_schedule_generation_;
    std::vector<std::uint32_t> division_schedule_generation_;
    std::vector<std::uint32_t> death_schedule_generation_;
    std::vector<Slot> free_slots_;
    std::array<std::size_t, 3> stage_counts_{};
    std::size_t alive_count_{};
};

}  // namespace atcg3d
