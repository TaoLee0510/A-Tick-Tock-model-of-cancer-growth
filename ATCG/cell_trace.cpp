#include "cell_trace.hpp"

#include <algorithm>

CellTraceStore::CellTraceStore(int row_count, int column_count)
{
    resize(row_count, column_count);
}

int CellTraceStore::rows() const
{
    return row_count_;
}

int CellTraceStore::column_count() const
{
    return column_count_;
}

void CellTraceStore::resize(int row_count, int column_count)
{
    row_count_ = row_count;
    column_count_ = column_count;
    values_.assign((row_count_ + 1) * (column_count_ + 1), 0L);
}

void CellTraceStore::resizeAndPreserve(int row_count, int column_count)
{
    std::vector<long> resized((row_count + 1) * (column_count + 1), 0L);
    int copied_rows = std::min(row_count_, row_count);
    int copied_cols = std::min(column_count_, column_count);
    for (int row = 1; row <= copied_rows; ++row)
    {
        for (int col = 1; col <= copied_cols; ++col)
        {
            resized[row * (column_count + 1) + col] = (*this)(row, col);
        }
    }
    row_count_ = row_count;
    column_count_ = column_count;
    values_.swap(resized);
}

void CellTraceStore::operator=(long value)
{
    std::fill(values_.begin(), values_.end(), value);
}

long &CellTraceStore::operator()(int row, int col)
{
    return values_[row * (column_count_ + 1) + col];
}

long CellTraceStore::operator()(int row, int col) const
{
    return values_[row * (column_count_ + 1) + col];
}

void CellTraceStore::copy_row_from(const CellTraceStore &source, int source_row, int target_row)
{
    int copied_cols = std::min(column_count_, source.column_count());
    for (int col = 1; col <= copied_cols; ++col)
    {
        (*this)(target_row, col) = source(source_row, col);
    }
}

void CellTraceStore::copy_columns_from(const CellTraceStore &source, int source_row, int target_row, int first_col, int last_col)
{
    int copied_last_col = std::min(last_col, std::min(column_count_, source.column_count()));
    for (int col = first_col; col <= copied_last_col; ++col)
    {
        (*this)(target_row, col) = source(source_row, col);
    }
}

void CellTraceStore::append_row_from(const CellTraceStore &source, int source_row)
{
    int target_row = row_count_ + 1;
    resizeAndPreserve(target_row, column_count_);
    copy_row_from(source, source_row, target_row);
}

void CellTraceStore::append_from(const CellTraceStore &source)
{
    int current_rows = row_count_;
    resizeAndPreserve(row_count_ + source.rows(), column_count_);
    for (int row = 1; row <= source.rows(); ++row)
    {
        copy_row_from(source, row, current_rows + row);
    }
}
