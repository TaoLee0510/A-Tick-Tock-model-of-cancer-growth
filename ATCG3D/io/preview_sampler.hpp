#pragma once

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "core/cell_store.hpp"

namespace atcg3d {

std::vector<Slot> stable_preview_sample(const CellStore3D& cells,
                                        std::size_t maximum_cells,
                                        std::uint64_t seed);

std::vector<std::size_t> stable_preview_sample(
    std::span<const CellInit> cells,
    std::size_t maximum_cells,
    std::uint64_t seed);

}  // namespace atcg3d
