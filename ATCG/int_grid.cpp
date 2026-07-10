#include "int_grid.hpp"

#include <algorithm>

IntGrid::IntGrid(int row_count, int column_count)
{
    resize(row_count, column_count);
}

void IntGrid::resize(int row_count, int column_count)
{
    row_count_ = row_count;
    column_count_ = column_count;
    values_.assign((row_count_ + 1) * (column_count_ + 1), 0);
}

void IntGrid::operator=(int value)
{
    std::fill(values_.begin(), values_.end(), value);
}

int &IntGrid::operator()(int row, int col)
{
    return values_[row * (column_count_ + 1) + col];
}

int IntGrid::operator()(int row, int col) const
{
    return values_[row * (column_count_ + 1) + col];
}
