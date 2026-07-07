//
//  cell_store.hpp
//  ATCG
//

#ifndef cell_store_hpp
#define cell_store_hpp

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <vector>
#include "cell_columns.hpp"

class CellStore
{
public:
    using Column = std::vector<double>;

    explicit CellStore(int column_count = cell_col::kMaxColumnCount)
    : column_count_(column_count), row_count_(0), columns_(column_count + 1)
    {
    }

    int column_count() const
    {
        return column_count_;
    }

    int rows() const
    {
        return row_count_;
    }

    void resize(int row_count)
    {
        row_count_ = row_count;
        for (int col = 1; col <= column_count_; ++col)
        {
            columns_[col].resize(row_count);
        }
    }

    void reserve(int row_capacity)
    {
        for (int col = 1; col <= column_count_; ++col)
        {
            columns_[col].reserve(row_capacity);
        }
    }

    void push_empty()
    {
        ++row_count_;
        for (int col = 1; col <= column_count_; ++col)
        {
            columns_[col].push_back(0.0);
        }
    }

    void clear_row(int row)
    {
        for (int col = 1; col <= column_count_; ++col)
        {
            columns_[col][row - 1] = 0.0;
        }
    }

    void append_row_from(const CellStore &source, int source_row)
    {
        push_empty();
        int copied_cols = std::min(column_count_, source.column_count());
        for (int col = 1; col <= copied_cols; ++col)
        {
            columns_[col][row_count_ - 1] = source.column(col)[source_row - 1];
        }
    }

    void append_from(const CellStore &source)
    {
        reserve(row_count_ + source.rows());
        for (int row = 1; row <= source.rows(); ++row)
        {
            append_row_from(source, row);
        }
    }

    bool has_consistent_row_count() const
    {
        for (int col = 1; col <= column_count_; ++col)
        {
            if ((int)columns_[col].size() != row_count_)
            {
                return false;
            }
        }
        return true;
    }

    void validate_row_count() const
    {
        if (!has_consistent_row_count())
        {
            throw std::runtime_error("CellStore column row counts are inconsistent");
        }
    }

    void apply_permutation(const std::vector<int> &zero_based_order)
    {
        if ((int)zero_based_order.size() != row_count_)
        {
            throw std::runtime_error("CellStore permutation size does not match row count");
        }

        for (int col = 1; col <= column_count_; ++col)
        {
            Column reordered(row_count_);
            const Column &source = columns_[col];
            for (int new_row = 0; new_row < row_count_; ++new_row)
            {
                reordered[new_row] = source[zero_based_order[new_row]];
            }
            columns_[col].swap(reordered);
        }
    }

    void sort_by_column(int col, bool descending = false)
    {
        std::vector<int> order(row_count_);
        std::iota(order.begin(), order.end(), 0);

        const Column &sort_values = column(col);
        std::sort(order.begin(), order.end(), [&](int lhs, int rhs) {
            if (descending)
            {
                return sort_values[lhs] > sort_values[rhs];
            }
            return sort_values[lhs] < sort_values[rhs];
        });

        apply_permutation(order);
    }

    Column &column(int col)
    {
        return columns_[col];
    }

    const Column &column(int col) const
    {
        return columns_[col];
    }

    Column &x1() { return column(cell_col::kX1); }
    Column &x2() { return column(cell_col::kX2); }
    Column &x3() { return column(cell_col::kX3); }
    Column &x4() { return column(cell_col::kX4); }
    Column &y1() { return column(cell_col::kY1); }
    Column &y2() { return column(cell_col::kY2); }
    Column &y3() { return column(cell_col::kY3); }
    Column &y4() { return column(cell_col::kY4); }
    Column &type() { return column(cell_col::kType); }
    Column &growth_rate() { return column(cell_col::kGrowthRate); }
    Column &density_growth_rate() { return column(cell_col::kDensityGrowthRate); }
    Column &migration_rate_base() { return column(cell_col::kMigrationRateBase); }
    Column &random_label() { return column(cell_col::kRandomLabel); }
    Column &stage() { return column(cell_col::kStage); }
    Column &id() { return column(cell_col::kId); }
    Column &division_elapsed() { return column(cell_col::kDivisionElapsed); }
    Column &division_time() { return column(cell_col::kDivisionTime); }
    Column &death_time() { return column(cell_col::kDeathTime); }
    Column &death_elapsed() { return column(cell_col::kDeathElapsed); }
    Column &migration_elapsed() { return column(cell_col::kMigrationElapsed); }
    Column &migration_interval() { return column(cell_col::kMigrationInterval); }
    Column &viability() { return column(cell_col::kViability); }
    Column &migration_direction() { return column(cell_col::kMigrationDirection); }
    Column &migration_follow_flag() { return column(cell_col::kMigrationFollowFlag); }
    Column &migration_active() { return column(cell_col::kMigrationActive); }
    Column &migration_duration() { return column(cell_col::kMigrationDuration); }
    Column &migration_passed() { return column(cell_col::kMigrationPassed); }
    Column &migration_rate() { return column(cell_col::kMigrationRate); }
    Column &cell_trace_label() { return column(cell_col::kCellTraceLabel); }
    Column &parent_trace_label() { return column(cell_col::kParentTraceLabel); }
    Column &division_count() { return column(cell_col::kDivisionCount); }
    Column &division_marker() { return column(cell_col::kDivisionMarker); }

    const Column &x1() const { return column(cell_col::kX1); }
    const Column &x2() const { return column(cell_col::kX2); }
    const Column &x3() const { return column(cell_col::kX3); }
    const Column &x4() const { return column(cell_col::kX4); }
    const Column &y1() const { return column(cell_col::kY1); }
    const Column &y2() const { return column(cell_col::kY2); }
    const Column &y3() const { return column(cell_col::kY3); }
    const Column &y4() const { return column(cell_col::kY4); }
    const Column &type() const { return column(cell_col::kType); }
    const Column &growth_rate() const { return column(cell_col::kGrowthRate); }
    const Column &density_growth_rate() const { return column(cell_col::kDensityGrowthRate); }
    const Column &migration_rate_base() const { return column(cell_col::kMigrationRateBase); }
    const Column &random_label() const { return column(cell_col::kRandomLabel); }
    const Column &stage() const { return column(cell_col::kStage); }
    const Column &id() const { return column(cell_col::kId); }
    const Column &division_elapsed() const { return column(cell_col::kDivisionElapsed); }
    const Column &division_time() const { return column(cell_col::kDivisionTime); }
    const Column &death_time() const { return column(cell_col::kDeathTime); }
    const Column &death_elapsed() const { return column(cell_col::kDeathElapsed); }
    const Column &migration_elapsed() const { return column(cell_col::kMigrationElapsed); }
    const Column &migration_interval() const { return column(cell_col::kMigrationInterval); }
    const Column &viability() const { return column(cell_col::kViability); }
    const Column &migration_direction() const { return column(cell_col::kMigrationDirection); }
    const Column &migration_follow_flag() const { return column(cell_col::kMigrationFollowFlag); }
    const Column &migration_active() const { return column(cell_col::kMigrationActive); }
    const Column &migration_duration() const { return column(cell_col::kMigrationDuration); }
    const Column &migration_passed() const { return column(cell_col::kMigrationPassed); }
    const Column &migration_rate() const { return column(cell_col::kMigrationRate); }
    const Column &cell_trace_label() const { return column(cell_col::kCellTraceLabel); }
    const Column &parent_trace_label() const { return column(cell_col::kParentTraceLabel); }
    const Column &division_count() const { return column(cell_col::kDivisionCount); }
    const Column &division_marker() const { return column(cell_col::kDivisionMarker); }

private:
    int column_count_;
    int row_count_;
    std::vector<Column> columns_;
};

class CellRowBuffer
{
public:
    explicit CellRowBuffer(int row_count = 0, int column_count = cell_col::kMaxColumnCount)
    {
        resize(row_count, column_count);
    }

    void resize(int row_count, int column_count)
    {
        row_count_ = row_count;
        column_count_ = column_count;
        values_.assign((row_count_ + 1) * (column_count_ + 1), 0.0);
    }

    void resizeAndPreserve(int row_count, int column_count)
    {
        std::vector<double> resized((row_count + 1) * (column_count + 1), 0.0);
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

    int rows() const
    {
        return row_count_;
    }

    int column_count() const
    {
        return column_count_;
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

inline void cell_store_copy_row_to_array(const CellStore &source, int source_row, CellRowBuffer &target, int target_row, int column_count)
{
    int copied_cols = std::min(column_count, source.column_count());
    copied_cols = std::min(copied_cols, target.column_count());
    for (int col = 1; col <= copied_cols; ++col)
    {
        target(target_row, col) = source.column(col)[source_row - 1];
    }
}

inline void cell_store_assign_row_from_array(CellStore &target, int target_row, const CellRowBuffer &source, int source_row, int column_count)
{
    int copied_cols = std::min(column_count, target.column_count());
    copied_cols = std::min(copied_cols, source.column_count());
    for (int col = 1; col <= copied_cols; ++col)
    {
        target.column(col)[target_row - 1] = source(source_row, col);
    }
}

inline void cell_store_assign_row_from_array(CellRowBuffer &target, int target_row, const CellRowBuffer &source, int source_row, int column_count)
{
    int copied_cols = std::min(column_count, target.column_count());
    copied_cols = std::min(copied_cols, source.column_count());
    for (int col = 1; col <= copied_cols; ++col)
    {
        target(target_row, col) = source(source_row, col);
    }
}

inline void cell_store_append_row_from_array(CellStore &target, const CellRowBuffer &source, int source_row, int column_count)
{
    target.push_empty();
    cell_store_assign_row_from_array(target, target.rows(), source, source_row, column_count);
}

inline void cell_store_append_row_from_array(CellRowBuffer &target, const CellRowBuffer &source, int source_row, int column_count)
{
    int current_size = target.rows();
    target.resizeAndPreserve(current_size + 1, column_count);
    cell_store_assign_row_from_array(target, current_size + 1, source, source_row, column_count);
}

#endif /* cell_store_hpp */
