#ifndef int_grid_hpp
#define int_grid_hpp

#include <vector>

class IntGrid
{
public:
    IntGrid(int row_count = 0, int column_count = 0);

    void resize(int row_count, int column_count);
    void operator=(int value);
    int &operator()(int row, int col);
    int operator()(int row, int col) const;

private:
    int row_count_ = 0;
    int column_count_ = 0;
    std::vector<int> values_;
};

#endif
