#include "stateless_rng.hpp"

#include <algorithm>
#include <climits>
#include <limits>
#include <random>

std::uint64_t stateless_splitmix64(std::uint64_t x)
{
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

unsigned long stateless_seed(long cell_id, long time_step, long event_id)
{
    std::uint64_t seed = 0x6a09e667f3bcc909ULL;
    seed ^= stateless_splitmix64(static_cast<std::uint64_t>(cell_id) + 0x243f6a8885a308d3ULL);
    seed ^= stateless_splitmix64(static_cast<std::uint64_t>(time_step) + 0x13198a2e03707344ULL);
    seed ^= stateless_splitmix64(static_cast<std::uint64_t>(event_id) + 0xa4093822299f31d0ULL);

    unsigned long out = static_cast<unsigned long>(seed & static_cast<std::uint64_t>(ULONG_MAX));
    return out == 0 ? 1 : out;
}

double stateless_uniform(long cell_id, long time_step, long event_id)
{
    std::mt19937_64 engine(stateless_seed(cell_id, time_step, event_id));
    return std::generate_canonical<double, 53>(engine);
}

double stateless_beta(long cell_id, long time_step, long event_id, double alpha, double beta)
{
    if (alpha <= 0.0 || beta <= 0.0)
    {
        return 0.0;
    }

    std::mt19937_64 engine(stateless_seed(cell_id, time_step, event_id));
    std::gamma_distribution<double> gamma_alpha(alpha, 1.0);
    std::gamma_distribution<double> gamma_beta(beta, 1.0);
    double x = gamma_alpha(engine);
    double y = gamma_beta(engine);
    double sum = x + y;
    return sum > 0.0 ? x / sum : alpha / (alpha + beta);
}

unsigned int stateless_geometric(long cell_id, long time_step, long event_id, double p)
{
    if (p >= 1.0)
    {
        return 1;
    }
    if (p <= 0.0)
    {
        return std::numeric_limits<unsigned int>::max();
    }

    std::mt19937_64 engine(stateless_seed(cell_id, time_step, event_id));
    std::geometric_distribution<unsigned int> geometric(p);
    return geometric(engine) + 1;
}

void stateless_shuffle(int *first, int *last, long cell_id, long time_step, long event_id)
{
    std::mt19937_64 engine(stateless_seed(cell_id, time_step, event_id));
    std::shuffle(first, last, engine);
}
