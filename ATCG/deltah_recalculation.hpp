#ifndef deltah_recalculation_hpp
#define deltah_recalculation_hpp

#include "cell_store.hpp"

void set_deltah_from_max_migration_rate(double &deltah, int &MMR, double max_mig_r);
void deltah_recalculation(double &deltah, const CellStore &cells, int &MMR, int DDM);

#endif
