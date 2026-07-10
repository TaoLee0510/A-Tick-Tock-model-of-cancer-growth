#include "int_matrix.hpp"

#include <algorithm>

IntMatrix::IntMatrix() = default;

IntMatrix::IntMatrix(int row_count, int column_count)
{
    resize(row_count, column_count);
}

void IntMatrix::resize(int row_count, int column_count)
{
    rows_ = row_count;
    columns_ = column_count;
    values_.assign(rows_ * columns_, 0);
}

IntMatrix &IntMatrix::operator=(int value)
{
    std::fill(values_.begin(), values_.end(), value);
    return *this;
}

int IntMatrix::columns() const
{
    return columns_;
}

int &IntMatrix::operator()(int row, int column)
{
    return values_[(row - 1) * columns_ + (column - 1)];
}

int IntMatrix::operator()(int row, int column) const
{
    return values_[(row - 1) * columns_ + (column - 1)];
}

void IntMatrix::set_col(int column, int first, int second)
{
    (*this)(1, column) = first;
    (*this)(2, column) = second;
}

void IntMatrix::copy_col_from(int target_column, const IntMatrix &source, int source_column)
{
    for (int row = 1; row <= rows_ && row <= source.rows_; ++row)
    {
        (*this)(row, target_column) = source(row, source_column);
    }
}
