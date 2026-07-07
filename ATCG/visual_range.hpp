//
//  visual_range.hpp
//  ATCG
//

#ifndef visual_range_hpp
#define visual_range_hpp

#include <algorithm>
#include <cstddef>
#include <vector>

class VisualRange
{
public:
    VisualRange() = default;

    VisualRange(int width, int height)
    {
        resize(width, height);
    }

    void resize(int width, int height)
    {
        width_ = width;
        height_ = height;
        std::size_t count = static_cast<std::size_t>(width_) * static_cast<std::size_t>(height_);
        occupied_.assign(count, 0);
        density_label_.assign(count, 0);
        stage_.assign(count, 0);
        cell_label_.assign(count, 0);
    }

    void clear()
    {
        std::fill(occupied_.begin(), occupied_.end(), 0);
        std::fill(density_label_.begin(), density_label_.end(), 0);
        std::fill(stage_.begin(), stage_.end(), 0);
        std::fill(cell_label_.begin(), cell_label_.end(), 0);
    }

    int width() const { return width_; }
    int height() const { return height_; }

    long &occupied(int x, int y) { return occupied_[index(x, y)]; }
    long occupied(int x, int y) const { return occupied_[index(x, y)]; }

    long &density_label(int x, int y) { return density_label_[index(x, y)]; }
    long density_label(int x, int y) const { return density_label_[index(x, y)]; }

    long &stage(int x, int y) { return stage_[index(x, y)]; }
    long stage(int x, int y) const { return stage_[index(x, y)]; }

    long &cell_label(int x, int y) { return cell_label_[index(x, y)]; }
    long cell_label(int x, int y) const { return cell_label_[index(x, y)]; }

    void clear_site(int x, int y)
    {
        set_layers(x, y, 1, 4, 0);
    }

    void clear_square(int x1, int y1)
    {
        set_rect(x1, x1 + 1, y1, y1 + 1, 1, 4, 0);
    }

    void write_site(int x, int y, long cell_id, long cell_stage, long label)
    {
        occupied(x, y) = 1;
        density_label(x, y) = cell_id;
        stage(x, y) = cell_stage;
        cell_label(x, y) = label;
    }

    void write_square(int x1, int y1, long cell_id, long cell_stage, long label)
    {
        set_rect(x1, x1 + 1, y1, y1 + 1, 1, 1, 1);
        set_rect(x1, x1 + 1, y1, y1 + 1, 2, 2, cell_id);
        set_rect(x1, x1 + 1, y1, y1 + 1, 3, 3, cell_stage);
        set_rect(x1, x1 + 1, y1, y1 + 1, 4, 4, label);
    }

    void set_square_occupied(int x1, int y1, long value)
    {
        set_rect(x1, x1 + 1, y1, y1 + 1, 1, 1, value);
    }

    void set_square_density_label(int x1, int y1, long value)
    {
        set_rect(x1, x1 + 1, y1, y1 + 1, 2, 2, value);
    }

    void set_square_stage(int x1, int y1, long value)
    {
        set_rect(x1, x1 + 1, y1, y1 + 1, 3, 3, value);
    }

    void set_square_cell_label(int x1, int y1, long value)
    {
        set_rect(x1, x1 + 1, y1, y1 + 1, 4, 4, value);
    }

private:
    std::size_t index(int x, int y) const
    {
        return static_cast<std::size_t>(x - 1) * static_cast<std::size_t>(height_) + static_cast<std::size_t>(y - 1);
    }

    void set_layer(int x, int y, int layer, long value)
    {
        std::size_t site = index(x, y);
        switch (layer)
        {
            case 1:
                occupied_[site] = value;
                break;
            case 2:
                density_label_[site] = value;
                break;
            case 3:
                stage_[site] = value;
                break;
            default:
                cell_label_[site] = value;
                break;
        }
    }

    void set_layers(int x, int y, int first_layer, int last_layer, long value)
    {
        for (int layer = first_layer; layer <= last_layer; ++layer)
        {
            set_layer(x, y, layer, value);
        }
    }

    void set_rect(int first_x, int last_x, int first_y, int last_y, int first_layer, int last_layer, long value)
    {
        for (int x = first_x; x <= last_x; ++x)
        {
            for (int y = first_y; y <= last_y; ++y)
            {
                set_layers(x, y, first_layer, last_layer, value);
            }
        }
    }

    int width_ = 0;
    int height_ = 0;
    std::vector<long> occupied_;
    std::vector<long> density_label_;
    std::vector<long> stage_;
    std::vector<long> cell_label_;
};

#endif /* visual_range_hpp */
