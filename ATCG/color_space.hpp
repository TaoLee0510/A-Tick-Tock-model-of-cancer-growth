#ifndef color_space_hpp
#define color_space_hpp

#include <vector>

class ColorSpace
{
public:
    ColorSpace(int row_count = 0, int column_count = 4);

    void resize(int row_count, int column_count = 4);
    void operator=(double value);
    double &operator()(int row, int col);
    double operator()(int row, int col) const;

private:
    int row_count_ = 0;
    int column_count_ = 0;
    std::vector<double> values_;
};

#endif
