//
//  cell_trace.hpp
//  ATCG
//

#ifndef cell_trace_hpp
#define cell_trace_hpp

#include <algorithm>
#include <vector>

class CellTraceStore
{
public:
    explicit CellTraceStore(int row_count = 0, int column_count = 150)
    {
        resize(row_count, column_count);
    }

    int rows() const
    {
        return row_count_;
    }

    int column_count() const
    {
        return column_count_;
    }

    void resize(int row_count, int column_count = 150)
    {
        row_count_ = row_count;
        column_count_ = column_count;
        values_.assign((row_count_ + 1) * (column_count_ + 1), 0L);
    }

    void resizeAndPreserve(int row_count, int column_count = 150)
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

    void operator=(long value)
    {
        std::fill(values_.begin(), values_.end(), value);
    }

    long &operator()(int row, int col)
    {
        return values_[row * (column_count_ + 1) + col];
    }

    long operator()(int row, int col) const
    {
        return values_[row * (column_count_ + 1) + col];
    }

    void copy_row_from(const CellTraceStore &source, int source_row, int target_row)
    {
        int copied_cols = std::min(column_count_, source.column_count());
        for (int col = 1; col <= copied_cols; ++col)
        {
            (*this)(target_row, col) = source(source_row, col);
        }
    }

    void copy_columns_from(const CellTraceStore &source, int source_row, int target_row, int first_col, int last_col)
    {
        int copied_last_col = std::min(last_col, std::min(column_count_, source.column_count()));
        for (int col = first_col; col <= copied_last_col; ++col)
        {
            (*this)(target_row, col) = source(source_row, col);
        }
    }

    void append_row_from(const CellTraceStore &source, int source_row)
    {
        int target_row = row_count_ + 1;
        resizeAndPreserve(target_row, column_count_);
        copy_row_from(source, source_row, target_row);
    }

    void append_from(const CellTraceStore &source)
    {
        int current_rows = row_count_;
        resizeAndPreserve(row_count_ + source.rows(), column_count_);
        for (int row = 1; row <= source.rows(); ++row)
        {
            copy_row_from(source, row, current_rows + row);
        }
    }

private:
    int row_count_ = 0;
    int column_count_ = 0;
    std::vector<long> values_;
};

#endif /* cell_trace_hpp */
