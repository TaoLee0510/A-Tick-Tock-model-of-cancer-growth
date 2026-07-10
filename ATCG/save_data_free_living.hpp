#ifndef save_data_free_living_hpp
#define save_data_free_living_hpp

#include "cell_store.hpp"
#include "cell_trace.hpp"
#include "color_space.hpp"

void save_data_free_living(int Visual_range_x, int Visual_range_y, int N0, int N00, int N01, int MMR, int H, int T, double alpha, double beta, const CellStore &cell_array, int migration_judgement, double deltah, const ColorSpace &colorspace, int DDM, int allpng, int Col, const CellTraceStore &cell_trace);

#endif
