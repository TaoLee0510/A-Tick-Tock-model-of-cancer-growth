#ifndef SaveClonePNGHR_hpp
#define SaveClonePNGHR_hpp

#include "cell_store.hpp"
#include "color_space.hpp"

void SaveClonePNGHR(int Visual_range_x, int Visual_range_y, const CellStore &cell_array, int H, int &T, double alpha, double beta, double deltah, const ColorSpace &colorspace);

#endif
