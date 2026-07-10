#ifndef cell_store_hpp
#define cell_store_hpp

#include <vector>

#include "cell_columns.hpp"

class CellStore
{
public:
    using Column = std::vector<double>;

    explicit CellStore(int column_count = cell_col::kMaxColumnCount);

    int column_count() const;
    int rows() const;
    void resize(int row_count);
    void reserve(int row_capacity);
    void push_empty();
    void clear_row(int row);
    void append_row_from(const CellStore &source, int source_row);
    void append_from(const CellStore &source);
    bool has_consistent_row_count() const;
    void validate_row_count() const;
    void apply_permutation(const std::vector<int> &zero_based_order);
    void sort_by_column(int col, bool descending = false);

    Column &column(int col);
    const Column &column(int col) const;

    Column &x1();
    Column &x2();
    Column &x3();
    Column &x4();
    Column &y1();
    Column &y2();
    Column &y3();
    Column &y4();
    Column &type();
    Column &growth_rate();
    Column &density_growth_rate();
    Column &migration_rate_base();
    Column &random_label();
    Column &stage();
    Column &id();
    Column &division_elapsed();
    Column &division_time();
    Column &death_time();
    Column &death_elapsed();
    Column &migration_elapsed();
    Column &migration_interval();
    Column &viability();
    Column &migration_direction();
    Column &migration_follow_flag();
    Column &migration_active();
    Column &migration_duration();
    Column &migration_passed();
    Column &migration_rate();
    Column &cell_trace_label();
    Column &parent_trace_label();
    Column &division_count();
    Column &division_marker();

    const Column &x1() const;
    const Column &x2() const;
    const Column &x3() const;
    const Column &x4() const;
    const Column &y1() const;
    const Column &y2() const;
    const Column &y3() const;
    const Column &y4() const;
    const Column &type() const;
    const Column &growth_rate() const;
    const Column &density_growth_rate() const;
    const Column &migration_rate_base() const;
    const Column &random_label() const;
    const Column &stage() const;
    const Column &id() const;
    const Column &division_elapsed() const;
    const Column &division_time() const;
    const Column &death_time() const;
    const Column &death_elapsed() const;
    const Column &migration_elapsed() const;
    const Column &migration_interval() const;
    const Column &viability() const;
    const Column &migration_direction() const;
    const Column &migration_follow_flag() const;
    const Column &migration_active() const;
    const Column &migration_duration() const;
    const Column &migration_passed() const;
    const Column &migration_rate() const;
    const Column &cell_trace_label() const;
    const Column &parent_trace_label() const;
    const Column &division_count() const;
    const Column &division_marker() const;

private:
    int column_count_;
    int row_count_;
    std::vector<Column> columns_;
};

class CellRowBuffer
{
public:
    explicit CellRowBuffer(int row_count = 0, int column_count = cell_col::kMaxColumnCount);

    void resize(int row_count, int column_count);
    void resizeAndPreserve(int row_count, int column_count);
    int rows() const;
    int column_count() const;
    void operator=(double value);
    double &operator()(int row, int col);
    double operator()(int row, int col) const;

private:
    int row_count_ = 0;
    int column_count_ = 0;
    std::vector<double> values_;
};

void cell_store_copy_row_to_array(const CellStore &source, int source_row, CellRowBuffer &target, int target_row, int column_count);
void cell_store_assign_row_from_array(CellStore &target, int target_row, const CellRowBuffer &source, int source_row, int column_count);
void cell_store_assign_row_from_array(CellRowBuffer &target, int target_row, const CellRowBuffer &source, int source_row, int column_count);
void cell_store_append_row_from_array(CellStore &target, const CellRowBuffer &source, int source_row, int column_count);
void cell_store_append_row_from_array(CellRowBuffer &target, const CellRowBuffer &source, int source_row, int column_count);

#endif
