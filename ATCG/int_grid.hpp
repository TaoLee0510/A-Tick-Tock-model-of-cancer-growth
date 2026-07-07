//
//  int_grid.hpp
//  ATCG
//

#ifndef int_grid_hpp
#define int_grid_hpp

#include <algorithm>
#include <vector>

class IntGrid
{
public:
    IntGrid(int row_count = 0, int column_count = 0)
    {
        resize(row_count, column_count);
    }

    void resize(int row_count, int column_count)
    {
        row_count_ = row_count;
        column_count_ = column_count;
        values_.assign((row_count_ + 1) * (column_count_ + 1), 0);
    }

    void operator=(int value)
    {
        std::fill(values_.begin(), values_.end(), value);
    }

    int &operator()(int row, int col)
    {
        return values_[row * (column_count_ + 1) + col];
    }

    int operator()(int row, int col) const
    {
        return values_[row * (column_count_ + 1) + col];
    }

private:
    int row_count_ = 0;
    int column_count_ = 0;
    std::vector<int> values_;
};

#endif /* int_grid_hpp */
