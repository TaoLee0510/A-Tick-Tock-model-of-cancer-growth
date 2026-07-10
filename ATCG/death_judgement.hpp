#ifndef death_judgement_hpp
#define death_judgement_hpp

#include "cell_store.hpp"
#include "visual_range.hpp"

void death_judgement(int Visual_range_x, int Visual_range_y, int N00, int N01, double r_limit, double K_limit, double lambda_r, double lambda_K, double alpha, double beta, double carrying_capacity_r, double carrying_capacity_K, double Cr, double CK, double death_time_range_r, double death_time_range_K, double deltah, double &h, CellStore &cell_array, VisualRange &Visual_range, double deathjudge, int Col, int nthreads, long rng_time_step);

#endif
