#ifndef cell_trace_hpp
#define cell_trace_hpp

#include <vector>

class CellTraceStore
{
public:
    explicit CellTraceStore(int row_count = 0, int column_count = 150);

    int rows() const;
    int column_count() const;
    void resize(int row_count, int column_count = 150);
    void resizeAndPreserve(int row_count, int column_count = 150);
    void operator=(long value);
    long &operator()(int row, int col);
    long operator()(int row, int col) const;
    void copy_row_from(const CellTraceStore &source, int source_row, int target_row);
    void copy_columns_from(const CellTraceStore &source, int source_row, int target_row, int first_col, int last_col);
    void append_row_from(const CellTraceStore &source, int source_row);
    void append_from(const CellTraceStore &source);

private:
    int row_count_ = 0;
    int column_count_ = 0;
    std::vector<long> values_;
};

#endif
