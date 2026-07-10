#include "color_space.hpp"

#include <algorithm>

ColorSpace::ColorSpace(int row_count, int column_count)
{
    resize(row_count, column_count);
}

void ColorSpace::resize(int row_count, int column_count)
{
    row_count_ = row_count;
    column_count_ = column_count;
    values_.assign((row_count_ + 1) * (column_count_ + 1), 0.0);
}

void ColorSpace::operator=(double value)
{
    std::fill(values_.begin(), values_.end(), value);
}

double &ColorSpace::operator()(int row, int col)
{
    return values_[row * (column_count_ + 1) + col];
}

double ColorSpace::operator()(int row, int col) const
{
    return values_[row * (column_count_ + 1) + col];
}
