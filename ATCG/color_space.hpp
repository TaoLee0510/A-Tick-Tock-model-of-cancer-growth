//
//  color_space.hpp
//  ATCG
//

#ifndef color_space_hpp
#define color_space_hpp

#include <algorithm>
#include <vector>

class ColorSpace
{
public:
    ColorSpace(int row_count = 0, int column_count = 4)
    {
        resize(row_count, column_count);
    }

    void resize(int row_count, int column_count = 4)
    {
        row_count_ = row_count;
        column_count_ = column_count;
        values_.assign((row_count_ + 1) * (column_count_ + 1), 0.0);
    }

    void operator=(double value)
    {
        std::fill(values_.begin(), values_.end(), value);
    }

    double &operator()(int row, int col)
    {
        return values_[row * (column_count_ + 1) + col];
    }

    double operator()(int row, int col) const
    {
        return values_[row * (column_count_ + 1) + col];
    }

private:
    int row_count_ = 0;
    int column_count_ = 0;
    std::vector<double> values_;
};

#endif /* color_space_hpp */
