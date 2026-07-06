//
//  cell_store.hpp
//  ATCG
//

#ifndef cell_store_hpp
#define cell_store_hpp

#include <algorithm>
#include <cstddef>
#include <vector>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include "cell_columns.hpp"

using namespace blitz;

class CellStore
{
public:
    using Column = std::vector<double>;

    explicit CellStore(int column_count = cell_col::kMaxColumnCount)
    : column_count_(column_count), columns_(column_count + 1)
    {
    }

    int column_count() const
    {
        return column_count_;
    }

    int rows() const
    {
        return columns_.size() > 1 ? static_cast<int>(columns_[1].size()) : 0;
    }

    void resize(int row_count)
    {
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
        for (int col = 1; col <= column_count_; ++col)
        {
            columns_[col].push_back(0.0);
        }
    }

    double &operator()(int row, int col)
    {
        return columns_[col][row - 1];
    }

    double operator()(int row, int col) const
    {
        return columns_[col][row - 1];
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
    Column &y1() { return column(cell_col::kY1); }
    Column &type() { return column(cell_col::kType); }
    Column &growth_rate() { return column(cell_col::kGrowthRate); }
    Column &density_growth_rate() { return column(cell_col::kDensityGrowthRate); }
    Column &migration_rate_base() { return column(cell_col::kMigrationRateBase); }
    Column &stage() { return column(cell_col::kStage); }
    Column &id() { return column(cell_col::kId); }
    Column &division_elapsed() { return column(cell_col::kDivisionElapsed); }
    Column &division_time() { return column(cell_col::kDivisionTime); }
    Column &migration_elapsed() { return column(cell_col::kMigrationElapsed); }
    Column &migration_interval() { return column(cell_col::kMigrationInterval); }
    Column &viability() { return column(cell_col::kViability); }
    Column &migration_direction() { return column(cell_col::kMigrationDirection); }
    Column &migration_rate() { return column(cell_col::kMigrationRate); }

    const Column &x1() const { return column(cell_col::kX1); }
    const Column &y1() const { return column(cell_col::kY1); }
    const Column &type() const { return column(cell_col::kType); }
    const Column &growth_rate() const { return column(cell_col::kGrowthRate); }
    const Column &density_growth_rate() const { return column(cell_col::kDensityGrowthRate); }
    const Column &migration_rate_base() const { return column(cell_col::kMigrationRateBase); }
    const Column &stage() const { return column(cell_col::kStage); }
    const Column &id() const { return column(cell_col::kId); }
    const Column &division_elapsed() const { return column(cell_col::kDivisionElapsed); }
    const Column &division_time() const { return column(cell_col::kDivisionTime); }
    const Column &migration_elapsed() const { return column(cell_col::kMigrationElapsed); }
    const Column &migration_interval() const { return column(cell_col::kMigrationInterval); }
    const Column &viability() const { return column(cell_col::kViability); }
    const Column &migration_direction() const { return column(cell_col::kMigrationDirection); }
    const Column &migration_rate() const { return column(cell_col::kMigrationRate); }

private:
    int column_count_;
    std::vector<Column> columns_;
};

inline CellStore cell_store_from_array(const Array<double, 2> &cell_array, int column_count = 0)
{
    int cols = column_count > 0 ? column_count : cell_array.cols();
    CellStore cells(cols);
    int row_count = cell_array.rows();
    cells.resize(row_count);

    for (int col = 1; col <= cols; ++col)
    {
        CellStore::Column &target = cells.column(col);
        for (int row = 1; row <= row_count; ++row)
        {
            target[row - 1] = cell_array(row, col);
        }
    }

    return cells;
}

inline Array<double, 2> cell_array_from_store(const CellStore &cells, int column_count = 0)
{
    int cols = column_count > 0 ? column_count : cells.column_count();
    cols = std::min(cols, cells.column_count());

    Array<double, 2> cell_array(cells.rows(), cols, FortranArray<2>());
    for (int col = 1; col <= cols; ++col)
    {
        const CellStore::Column &source = cells.column(col);
        for (int row = 1; row <= cells.rows(); ++row)
        {
            cell_array(row, col) = source[row - 1];
        }
    }

    return cell_array;
}

inline void cell_store_to_array(const CellStore &cells, Array<double, 2> &cell_array, int column_count = 0)
{
    cell_array = cell_array_from_store(cells, column_count);
}

#endif /* cell_store_hpp */
