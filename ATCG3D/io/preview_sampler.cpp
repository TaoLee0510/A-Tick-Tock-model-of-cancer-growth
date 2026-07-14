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
    const std::vector<Slot> alive = cells.alive_slots();
    if (alive.size() <= maximum_cells) {
        std::vector<Slot> result = alive;
        std::sort(result.begin(), result.end(), [&cells](Slot lhs, Slot rhs) {
            return cells.uid(lhs) < cells.uid(rhs);
        });
        return result;
    }

    struct Candidate {
        std::uint64_t hash{};
        CellUid uid{};
        Slot slot{};
    };
    std::vector<Candidate> candidates;
    candidates.reserve(alive.size());
    for (const Slot slot : alive) {
        const CellUid uid = cells.uid(slot);
        candidates.push_back({splitmix64(uid ^ splitmix64(seed)), uid, slot});
    }
    const auto order = [](const Candidate& lhs, const Candidate& rhs) {
        if (lhs.hash != rhs.hash) return lhs.hash < rhs.hash;
        return lhs.uid < rhs.uid;
    };
    std::nth_element(candidates.begin(), candidates.begin() + static_cast<std::ptrdiff_t>(maximum_cells),
                     candidates.end(), order);
    candidates.resize(maximum_cells);
    std::sort(candidates.begin(), candidates.end(), [](const Candidate& lhs, const Candidate& rhs) {
        return lhs.uid < rhs.uid;
    });
    std::vector<Slot> result;
    result.reserve(candidates.size());
    for (const Candidate& candidate : candidates) {
        result.push_back(candidate.slot);
    }
    return result;
}

}  // namespace atcg3d
