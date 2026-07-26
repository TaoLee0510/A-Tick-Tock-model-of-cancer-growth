#include "io/preview_sampler.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>

#include "core/stateless_rng.hpp"

namespace atcg3d {

std::vector<Slot> stable_preview_sample(const CellStore3D& cells,
                                        std::size_t maximum_cells,
                                        std::uint64_t seed) {
    if (maximum_cells == 0) {
        throw std::invalid_argument("preview maximum must be positive");
    }
    if (cells.alive_count() <= maximum_cells) {
        std::vector<Slot> result;
        result.reserve(cells.alive_count());
        for (std::size_t raw_slot = 0; raw_slot < cells.slot_count();
             ++raw_slot) {
            const Slot slot = static_cast<Slot>(raw_slot);
            if (cells.valid(slot)) result.push_back(slot);
        }
        std::sort(result.begin(), result.end(), [&cells](Slot lhs, Slot rhs) {
            return cells.uid(lhs) < cells.uid(rhs);
        });
        return result;
    }

    struct Candidate {
        std::uint64_t hash{};
        Slot slot{};
    };
    std::vector<Candidate> candidates;
    candidates.reserve(cells.alive_count());
    for (std::size_t raw_slot = 0; raw_slot < cells.slot_count(); ++raw_slot) {
        const Slot slot = static_cast<Slot>(raw_slot);
        if (!cells.valid(slot)) continue;
        const CellUid uid = cells.uid(slot);
        candidates.push_back({splitmix64(uid ^ splitmix64(seed)), slot});
    }
    const auto order = [&cells](const Candidate& lhs, const Candidate& rhs) {
        if (lhs.hash != rhs.hash) return lhs.hash < rhs.hash;
        return cells.uid(lhs.slot) < cells.uid(rhs.slot);
    };
    std::nth_element(candidates.begin(), candidates.begin() + static_cast<std::ptrdiff_t>(maximum_cells),
                     candidates.end(), order);
    candidates.resize(maximum_cells);
    std::sort(candidates.begin(), candidates.end(), [&cells](const Candidate& lhs, const Candidate& rhs) {
        return cells.uid(lhs.slot) < cells.uid(rhs.slot);
    });
    std::vector<Slot> result;
    result.reserve(candidates.size());
    for (const Candidate& candidate : candidates) {
        result.push_back(candidate.slot);
    }
    return result;
}

std::vector<std::size_t> stable_preview_sample(
    std::span<const CellInit> cells,
    std::size_t maximum_cells,
    std::uint64_t seed) {
    if (maximum_cells == 0) {
        throw std::invalid_argument("preview maximum must be positive");
    }
    struct Candidate {
        std::uint64_t hash{};
        std::size_t index{};
    };
    std::vector<Candidate> candidates;
    candidates.reserve(cells.size());
    for (std::size_t index = 0; index < cells.size(); ++index) {
        const CellUid uid = cells[index].uid;
        candidates.push_back({splitmix64(uid ^ splitmix64(seed)), index});
    }
    const auto hash_order = [&cells](const Candidate& lhs, const Candidate& rhs) {
        if (lhs.hash != rhs.hash) return lhs.hash < rhs.hash;
        return cells[lhs.index].uid < cells[rhs.index].uid;
    };
    if (candidates.size() > maximum_cells) {
        std::nth_element(
            candidates.begin(),
            candidates.begin() + static_cast<std::ptrdiff_t>(maximum_cells),
            candidates.end(), hash_order);
        candidates.resize(maximum_cells);
    }
    std::sort(candidates.begin(), candidates.end(),
              [&cells](const Candidate& lhs, const Candidate& rhs) {
                  return cells[lhs.index].uid < cells[rhs.index].uid;
              });
    std::vector<std::size_t> result;
    result.reserve(candidates.size());
    for (const Candidate& candidate : candidates) result.push_back(candidate.index);
    return result;
}

}  // namespace atcg3d
