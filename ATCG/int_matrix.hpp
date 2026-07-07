//
//  int_matrix.hpp
//  ATCG
//

#ifndef int_matrix_hpp
#define int_matrix_hpp

#include <algorithm>
#include <vector>

class IntMatrix
{
public:
    IntMatrix() = default;

    IntMatrix(int row_count, int column_count)
    {
        resize(row_count, column_count);
    }

    void resize(int row_count, int column_count)
    {
        rows_ = row_count;
        columns_ = column_count;
        values_.assign(rows_ * columns_, 0);
    }

    IntMatrix &operator=(int value)
    {
        std::fill(values_.begin(), values_.end(), value);
        return *this;
    }

    int columns() const
    {
        return columns_;
    }

    int &operator()(int row, int column)
    {
        return values_[(row - 1) * columns_ + (column - 1)];
    }

    int operator()(int row, int column) const
    {
        return values_[(row - 1) * columns_ + (column - 1)];
    }

    void set_col(int column, int first, int second)
    {
        (*this)(1, column) = first;
        (*this)(2, column) = second;
    }

    void copy_col_from(int target_column, const IntMatrix &source, int source_column)
    {
        for (int row = 1; row <= rows_ && row <= source.rows_; ++row)
        {
            (*this)(row, target_column) = source(row, source_column);
        }
    }

private:
    int rows_ = 0;
    int columns_ = 0;
    std::vector<int> values_;
};

#endif /* int_matrix_hpp */
