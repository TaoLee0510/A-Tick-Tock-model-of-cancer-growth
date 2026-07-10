#ifndef stateless_rng_hpp
#define stateless_rng_hpp

#include <cstdint>

std::uint64_t stateless_splitmix64(std::uint64_t x);
unsigned long stateless_seed(long cell_id, long time_step, long event_id);
double stateless_uniform(long cell_id, long time_step, long event_id);
double stateless_beta(long cell_id, long time_step, long event_id, double alpha, double beta);
unsigned int stateless_geometric(long cell_id, long time_step, long event_id, double p);
void stateless_shuffle(int *first, int *last, long cell_id, long time_step, long event_id);

#endif
