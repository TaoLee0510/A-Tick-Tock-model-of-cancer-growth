#include "stateless_rng.hpp"

#include <algorithm>
#include <array>
#include <cassert>

int main()
{
    const double first = stateless_uniform(101, 25, 7);
    const double repeated = stateless_uniform(101, 25, 7);
    assert(first == repeated);

    assert(stateless_seed(101, 25, 7) != stateless_seed(102, 25, 7));
    assert(stateless_seed(101, 25, 7) != stateless_seed(101, 25, 8));

    std::array<int, 8> first_order{1, 2, 3, 4, 5, 6, 7, 8};
    std::array<int, 8> repeated_order = first_order;
    stateless_shuffle(first_order.data(), first_order.data() + first_order.size(), 101, 25, 9);
    stateless_shuffle(repeated_order.data(), repeated_order.data() + repeated_order.size(), 101, 25, 9);
    assert(first_order == repeated_order);

    std::sort(first_order.begin(), first_order.end());
    assert((first_order == std::array<int, 8>{1, 2, 3, 4, 5, 6, 7, 8}));
}
