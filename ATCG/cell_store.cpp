#include "cell_store.hpp"

#include <algorithm>
#include <numeric>
#include <stdexcept>

CellStore::CellStore(int column_count)
    : column_count_(column_count), row_count_(0), columns_(column_count + 1)
{
}

int CellStore::column_count() const
{
    return column_count_;
}

int CellStore::rows() const
{
    return row_count_;
}

void CellStore::resize(int row_count)
{
    row_count_ = row_count;
    for (int col = 1; col <= column_count_; ++col)
    {
        columns_[col].resize(row_count);
    }
}

void CellStore::reserve(int row_capacity)
{
    for (int col = 1; col <= column_count_; ++col)
    {
        columns_[col].reserve(row_capacity);
    }
}

void CellStore::push_empty()
{
    ++row_count_;
    for (int col = 1; col <= column_count_; ++col)
    {
        columns_[col].push_back(0.0);
    }
}

void CellStore::clear_row(int row)
{
    for (int col = 1; col <= column_count_; ++col)
    {
        columns_[col][row - 1] = 0.0;
    }
}

void CellStore::append_row_from(const CellStore &source, int source_row)
{
    push_empty();
    int copied_cols = std::min(column_count_, source.column_count());
    for (int col = 1; col <= copied_cols; ++col)
    {
        columns_[col][row_count_ - 1] = source.column(col)[source_row - 1];
    }
}

void CellStore::append_from(const CellStore &source)
{
    reserve(row_count_ + source.rows());
    for (int row = 1; row <= source.rows(); ++row)
    {
        append_row_from(source, row);
    }
}

bool CellStore::has_consistent_row_count() const
{
    for (int col = 1; col <= column_count_; ++col)
    {
        if (static_cast<int>(columns_[col].size()) != row_count_)
        {
            return false;
        }
    }
    return true;
}

void CellStore::validate_row_count() const
{
    if (!has_consistent_row_count())
    {
        throw std::runtime_error("CellStore column row counts are inconsistent");
    }
}

void CellStore::apply_permutation(const std::vector<int> &zero_based_order)
{
    if (static_cast<int>(zero_based_order.size()) != row_count_)
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

void CellStore::sort_by_column(int col, bool descending)
{
    std::vector<int> order(row_count_);
    std::iota(order.begin(), order.end(), 0);

    const Column &sort_values = column(col);
    std::sort(order.begin(), order.end(), [&](int lhs, int rhs) {
        return descending ? sort_values[lhs] > sort_values[rhs] : sort_values[lhs] < sort_values[rhs];
    });
    apply_permutation(order);
}

CellStore::Column &CellStore::column(int col)
{
    return columns_[col];
}

const CellStore::Column &CellStore::column(int col) const
{
    return columns_[col];
}

#define ATCG_DEFINE_COLUMN_ACCESSOR(name, column_id)                  \
    CellStore::Column &CellStore::name() { return column(column_id); } \
    const CellStore::Column &CellStore::name() const { return column(column_id); }

ATCG_DEFINE_COLUMN_ACCESSOR(x1, cell_col::kX1)
ATCG_DEFINE_COLUMN_ACCESSOR(x2, cell_col::kX2)
ATCG_DEFINE_COLUMN_ACCESSOR(x3, cell_col::kX3)
ATCG_DEFINE_COLUMN_ACCESSOR(x4, cell_col::kX4)
ATCG_DEFINE_COLUMN_ACCESSOR(y1, cell_col::kY1)
ATCG_DEFINE_COLUMN_ACCESSOR(y2, cell_col::kY2)
ATCG_DEFINE_COLUMN_ACCESSOR(y3, cell_col::kY3)
ATCG_DEFINE_COLUMN_ACCESSOR(y4, cell_col::kY4)
ATCG_DEFINE_COLUMN_ACCESSOR(type, cell_col::kType)
ATCG_DEFINE_COLUMN_ACCESSOR(growth_rate, cell_col::kGrowthRate)
ATCG_DEFINE_COLUMN_ACCESSOR(density_growth_rate, cell_col::kDensityGrowthRate)
ATCG_DEFINE_COLUMN_ACCESSOR(migration_rate_base, cell_col::kMigrationRateBase)
ATCG_DEFINE_COLUMN_ACCESSOR(random_label, cell_col::kRandomLabel)
ATCG_DEFINE_COLUMN_ACCESSOR(stage, cell_col::kStage)
ATCG_DEFINE_COLUMN_ACCESSOR(id, cell_col::kId)
ATCG_DEFINE_COLUMN_ACCESSOR(division_elapsed, cell_col::kDivisionElapsed)
ATCG_DEFINE_COLUMN_ACCESSOR(division_time, cell_col::kDivisionTime)
ATCG_DEFINE_COLUMN_ACCESSOR(death_time, cell_col::kDeathTime)
ATCG_DEFINE_COLUMN_ACCESSOR(death_elapsed, cell_col::kDeathElapsed)
ATCG_DEFINE_COLUMN_ACCESSOR(migration_elapsed, cell_col::kMigrationElapsed)
ATCG_DEFINE_COLUMN_ACCESSOR(migration_interval, cell_col::kMigrationInterval)
ATCG_DEFINE_COLUMN_ACCESSOR(viability, cell_col::kViability)
ATCG_DEFINE_COLUMN_ACCESSOR(migration_direction, cell_col::kMigrationDirection)
ATCG_DEFINE_COLUMN_ACCESSOR(migration_follow_flag, cell_col::kMigrationFollowFlag)
ATCG_DEFINE_COLUMN_ACCESSOR(migration_active, cell_col::kMigrationActive)
ATCG_DEFINE_COLUMN_ACCESSOR(migration_duration, cell_col::kMigrationDuration)
ATCG_DEFINE_COLUMN_ACCESSOR(migration_passed, cell_col::kMigrationPassed)
ATCG_DEFINE_COLUMN_ACCESSOR(migration_rate, cell_col::kMigrationRate)
ATCG_DEFINE_COLUMN_ACCESSOR(cell_trace_label, cell_col::kCellTraceLabel)
ATCG_DEFINE_COLUMN_ACCESSOR(parent_trace_label, cell_col::kParentTraceLabel)
ATCG_DEFINE_COLUMN_ACCESSOR(division_count, cell_col::kDivisionCount)
ATCG_DEFINE_COLUMN_ACCESSOR(division_marker, cell_col::kDivisionMarker)

#undef ATCG_DEFINE_COLUMN_ACCESSOR

CellRowBuffer::CellRowBuffer(int row_count, int column_count)
{
    resize(row_count, column_count);
}

void CellRowBuffer::resize(int row_count, int column_count)
{
    row_count_ = row_count;
    column_count_ = column_count;
    values_.assign((row_count_ + 1) * (column_count_ + 1), 0.0);
}

void CellRowBuffer::resizeAndPreserve(int row_count, int column_count)
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

int CellRowBuffer::rows() const
{
    return row_count_;
}

int CellRowBuffer::column_count() const
{
    return column_count_;
}

void CellRowBuffer::operator=(double value)
{
    std::fill(values_.begin(), values_.end(), value);
}

double &CellRowBuffer::operator()(int row, int col)
{
    return values_[row * (column_count_ + 1) + col];
}

double CellRowBuffer::operator()(int row, int col) const
{
    return values_[row * (column_count_ + 1) + col];
}

void cell_store_copy_row_to_array(const CellStore &source, int source_row, CellRowBuffer &target, int target_row, int column_count)
{
    int copied_cols = std::min(column_count, source.column_count());
    copied_cols = std::min(copied_cols, target.column_count());
    for (int col = 1; col <= copied_cols; ++col)
    {
        target(target_row, col) = source.column(col)[source_row - 1];
    }
}

void cell_store_assign_row_from_array(CellStore &target, int target_row, const CellRowBuffer &source, int source_row, int column_count)
{
    int copied_cols = std::min(column_count, target.column_count());
    copied_cols = std::min(copied_cols, source.column_count());
    for (int col = 1; col <= copied_cols; ++col)
    {
        target.column(col)[target_row - 1] = source(source_row, col);
    }
}

void cell_store_assign_row_from_array(CellRowBuffer &target, int target_row, const CellRowBuffer &source, int source_row, int column_count)
{
    int copied_cols = std::min(column_count, target.column_count());
    copied_cols = std::min(copied_cols, source.column_count());
    for (int col = 1; col <= copied_cols; ++col)
    {
        target(target_row, col) = source(source_row, col);
    }
}

void cell_store_append_row_from_array(CellStore &target, const CellRowBuffer &source, int source_row, int column_count)
{
    target.push_empty();
    cell_store_assign_row_from_array(target, target.rows(), source, source_row, column_count);
}

void cell_store_append_row_from_array(CellRowBuffer &target, const CellRowBuffer &source, int source_row, int column_count)
{
    int current_size = target.rows();
    target.resizeAndPreserve(current_size + 1, column_count);
    cell_store_assign_row_from_array(target, current_size + 1, source, source_row, column_count);
}
