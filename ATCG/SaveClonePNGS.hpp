#ifndef SaveClonePNGS_hpp
#define SaveClonePNGS_hpp

#include "cell_store.hpp"
#include "color_space.hpp"

void SaveClonePNGS(int Visual_range_x, int Visual_range_y, int &T, double alpha, double beta, const CellStore &cell_array, const ColorSpace &colorspace);

#endif
