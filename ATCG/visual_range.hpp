#ifndef visual_range_hpp
#define visual_range_hpp

#include <cstddef>
#include <vector>

class VisualRange
{
public:
    VisualRange();
    VisualRange(int width, int height);

    void resize(int width, int height);
    void clear();
    int width() const;
    int height() const;
    long &occupied(int x, int y);
    long occupied(int x, int y) const;
    long &density_label(int x, int y);
    long density_label(int x, int y) const;
    long &stage(int x, int y);
    long stage(int x, int y) const;
    long &cell_label(int x, int y);
    long cell_label(int x, int y) const;
    void clear_site(int x, int y);
    void clear_square(int x1, int y1);
    void write_site(int x, int y, long cell_id, long cell_stage, long label);
    void write_square(int x1, int y1, long cell_id, long cell_stage, long label);
    void set_square_occupied(int x1, int y1, long value);
    void set_square_density_label(int x1, int y1, long value);
    void set_square_stage(int x1, int y1, long value);
    void set_square_cell_label(int x1, int y1, long value);

private:
    std::size_t index(int x, int y) const;
    void set_layer(int x, int y, int layer, long value);
    void set_layers(int x, int y, int first_layer, int last_layer, long value);
    void set_rect(int first_x, int last_x, int first_y, int last_y, int first_layer, int last_layer, long value);

    int width_ = 0;
    int height_ = 0;
    std::vector<long> occupied_;
    std::vector<long> density_label_;
    std::vector<long> stage_;
    std::vector<long> cell_label_;
};

#endif
