#pragma once

#include <span>
#include <vector>

namespace atcg3d::structured_pde {

// Expected uniform choice among independent feasible directions. The
// unallocated probability is the event that no direction is feasible.
std::vector<double> uniform_feasible_jump_probabilities(std::span<const double> availability);

}  // namespace atcg3d::structured_pde
