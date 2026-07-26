#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <stdexcept>

#include "config/model_config.hpp"
#include "core/cell_store.hpp"
#include "space/chunk_grid.hpp"
#include "space/density_index.hpp"

namespace atcg3d {

template <class T, std::size_t Capacity>
class FixedList3D {
public:
    using iterator = typename std::array<T, Capacity>::iterator;
    using const_iterator = typename std::array<T, Capacity>::const_iterator;

    bool empty() const noexcept { return size_ == 0; }
    std::size_t size() const noexcept { return size_; }
    T& front() noexcept { return values_.front(); }
    const T& front() const noexcept { return values_.front(); }
    T& operator[](std::size_t index) noexcept { return values_[index]; }
    const T& operator[](std::size_t index) const noexcept {
        return values_[index];
    }
    iterator begin() noexcept { return values_.begin(); }
    const_iterator begin() const noexcept { return values_.begin(); }
    iterator end() noexcept {
        return values_.begin() + static_cast<std::ptrdiff_t>(size_);
    }
    const_iterator end() const noexcept {
        return values_.begin() + static_cast<std::ptrdiff_t>(size_);
    }
    void clear() noexcept { size_ = 0; }
    void push_back(const T& value) {
        if (size_ == Capacity) {
            throw std::length_error("fixed migration list capacity exceeded");
        }
        values_[size_++] = value;
    }
    template <class Iterator>
    void assign(Iterator first, Iterator last) {
        clear();
        for (; first != last; ++first) push_back(*first);
    }
    FixedList3D& operator=(std::initializer_list<T> values) {
        assign(values.begin(), values.end());
        return *this;
    }

private:
    std::array<T, Capacity> values_{};
    std::size_t size_{};
};

using DirectionCandidates3D = FixedList3D<DirectionId, 26>;
using MigrationReservedSites3D = FixedList3D<Vec3i, 7>;

enum class RngEventKind : std::uint64_t {
    migration_direction = 1,
    division_location = 2,
    division_timing = 3,
    conflict_priority = 4,
    initialization = 5,
    death_timing = 6,
    division_type_conversion = 7,
    division_conflict_priority = 8,
    stage_recovery_conflict_priority = 9,
};

DirectionCandidates3D feasible_directions(
    Slot slot,
    const CellStore3D& cells,
    const SparseChunkGrid3D& grid,
    bool thin_layer);

DirectionId select_migration_direction(Slot slot,
                                       const CellStore3D& cells,
                                       const SparseChunkGrid3D& grid,
                                       const BlockDensityIndex3D& density,
                                       const Model3DConfig& config,
                                       std::uint64_t event_sequence);

// Crowding exchange is deliberately limited to two singleton stage-1 cells.
// Large footprints and co-location groups require a different many-site
// transaction and are therefore excluded from the first version.
DirectionCandidates3D feasible_crowding_swap_directions(
    Slot slot,
    const CellStore3D& cells,
    const SparseChunkGrid3D& grid,
    bool thin_layer);

DirectionId select_crowding_swap_direction(
    Slot slot,
    const CellStore3D& cells,
    const SparseChunkGrid3D& grid,
    const Model3DConfig& config,
    std::uint64_t event_sequence);

struct MoveProposal {
    Slot slot{kEmptySlot};
    CellUid uid{};
    DirectionId direction{};
    Vec3i from{};
    Vec3i to{};
    MigrationReservedSites3D reserved_sites;
    std::uint64_t priority{};
    bool swaps_anchors{};
    Slot swap_partner{kEmptySlot};
    CellUid swap_partner_uid{};
};

MoveProposal make_move_proposal(Slot slot,
                                const CellStore3D& cells,
                                const SparseChunkGrid3D& grid,
                                const BlockDensityIndex3D& density,
                                const Model3DConfig& config,
                                std::uint64_t event_sequence,
                                std::uint64_t time_bucket);

MoveProposal make_crowding_swap_proposal(
    Slot slot,
    DirectionId direction,
    const CellStore3D& cells,
    const SparseChunkGrid3D& grid,
    const Model3DConfig& config,
    std::uint64_t time_bucket);

bool commit_move(const MoveProposal& proposal,
                 CellStore3D& cells,
                 SparseChunkGrid3D& grid,
                 BlockDensityIndex3D& density);

bool commit_crowding_swap(const MoveProposal& proposal,
                          CellStore3D& cells,
                          SparseChunkGrid3D& grid,
                          BlockDensityIndex3D& density);

}  // namespace atcg3d
