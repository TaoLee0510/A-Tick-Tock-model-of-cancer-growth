#ifndef migration_hpp
#define migration_hpp

#include "cell_store.hpp"
#include "visual_range.hpp"

void migration(int i, double deltah, CellStore &cell_array, VisualRange &Visual_range, double &migration_judgement, long rng_time_step, long rng_event_base);

#endif
