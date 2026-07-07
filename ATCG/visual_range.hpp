//
//  visual_range.hpp
//  ATCG
//

#ifndef visual_range_hpp
#define visual_range_hpp

#include <algorithm>
#include <cstddef>
#include <vector>
#include <blitz/range.h>

class VisualRange
{
public:
    VisualRange() = default;

    VisualRange(int width, int height)
    {
        resize(width, height);
    }

    template <typename Storage>
    VisualRange(int width, int height, int, const Storage &)
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

    long &operator()(int x, int y, int layer)
    {
        return layer_ref(x, y, layer);
    }

    long operator()(int x, int y, int layer) const
    {
        return layer_value(x, y, layer);
    }

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
        (*this)(x, y, 1) = 1;
        (*this)(x, y, 2) = cell_id;
        (*this)(x, y, 3) = cell_stage;
        (*this)(x, y, 4) = label;
    }

    void write_square(int x1, int y1, long cell_id, long cell_stage, long label)
    {
        set_rect(x1, x1 + 1, y1, y1 + 1, 1, 1, 1);
        set_rect(x1, x1 + 1, y1, y1 + 1, 2, 2, cell_id);
        set_rect(x1, x1 + 1, y1, y1 + 1, 3, 3, cell_stage);
        set_rect(x1, x1 + 1, y1, y1 + 1, 4, 4, label);
    }

    class SiteLayersProxy
    {
    public:
        SiteLayersProxy(VisualRange &range, int x, int y, int first_layer, int last_layer)
            : range_(range), x_(x), y_(y), first_layer_(first_layer), last_layer_(last_layer)
        {
        }

        SiteLayersProxy &operator=(long value)
        {
            range_.set_layers(x_, y_, first_layer_, last_layer_, value);
            return *this;
        }

    private:
        VisualRange &range_;
        int x_;
        int y_;
        int first_layer_;
        int last_layer_;
    };

    class RectLayersProxy
    {
    public:
        RectLayersProxy(VisualRange &range, int first_x, int last_x, int first_y, int last_y, int first_layer, int last_layer)
            : range_(range), first_x_(first_x), last_x_(last_x), first_y_(first_y), last_y_(last_y), first_layer_(first_layer), last_layer_(last_layer)
        {
        }

        RectLayersProxy &operator=(long value)
        {
            range_.set_rect(first_x_, last_x_, first_y_, last_y_, first_layer_, last_layer_, value);
            return *this;
        }

    private:
        VisualRange &range_;
        int first_x_;
        int last_x_;
        int first_y_;
        int last_y_;
        int first_layer_;
        int last_layer_;
    };

    SiteLayersProxy operator()(int x, int y, const blitz::Range &layers)
    {
        return SiteLayersProxy(*this, x, y, layers.first(1), layers.last(4));
    }

    RectLayersProxy operator()(const blitz::Range &xs, const blitz::Range &ys, int layer)
    {
        return RectLayersProxy(*this, xs.first(1), xs.last(width_), ys.first(1), ys.last(height_), layer, layer);
    }

    RectLayersProxy operator()(const blitz::Range &xs, const blitz::Range &ys, const blitz::Range &layers)
    {
        return RectLayersProxy(*this, xs.first(1), xs.last(width_), ys.first(1), ys.last(height_), layers.first(1), layers.last(4));
    }

private:
    std::size_t index(int x, int y) const
    {
        return static_cast<std::size_t>(x - 1) * static_cast<std::size_t>(height_) + static_cast<std::size_t>(y - 1);
    }

    long &layer_ref(int x, int y, int layer)
    {
        std::size_t site = index(x, y);
        switch (layer)
        {
            case 1:
                return occupied_[site];
            case 2:
                return density_label_[site];
            case 3:
                return stage_[site];
            default:
                return cell_label_[site];
        }
    }

    long layer_value(int x, int y, int layer) const
    {
        std::size_t site = index(x, y);
        switch (layer)
        {
            case 1:
                return occupied_[site];
            case 2:
                return density_label_[site];
            case 3:
                return stage_[site];
            default:
                return cell_label_[site];
        }
    }

    void set_layers(int x, int y, int first_layer, int last_layer, long value)
    {
        for (int layer = first_layer; layer <= last_layer; ++layer)
        {
            (*this)(x, y, layer) = value;
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
