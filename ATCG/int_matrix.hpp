#ifndef int_matrix_hpp
#define int_matrix_hpp

#include <vector>

class IntMatrix
{
public:
    IntMatrix();
    IntMatrix(int row_count, int column_count);

    void resize(int row_count, int column_count);
    IntMatrix &operator=(int value);
    int columns() const;
    int &operator()(int row, int column);
    int operator()(int row, int column) const;
    void set_col(int column, int first, int second);
    void copy_col_from(int target_column, const IntMatrix &source, int source_column);

private:
    int rows_ = 0;
    int columns_ = 0;
    std::vector<int> values_;
};

#endif
