#include "model/sector_mean_field.hpp"

#include "engine/parallelism.hpp"
#include "geometry/directions.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>

namespace atcg3d::continuum {
SectorMeanField3D::SectorMeanField3D(std::array<int, 3> shape, bool thin,
                                   int window_edge, double half_angle_degrees, int threads)
    : shape_(shape), thin_(thin), edge_(window_edge), lower_((window_edge - 1) / 2) {
    if (window_edge < 1 || window_edge > 70 ||
        !std::isfinite(half_angle_degrees) || half_angle_degrees <= 0.0 || half_angle_degrees > 180.0 ||
        shape[0] < 1 || shape[1] < 1 || shape[2] < 1 || thin != (shape[2] == 1))
        throw std::invalid_argument("invalid cached sector geometry");
    int effective_window = 1;
    for (int axis = 0; axis < (thin ? 2 : 3); ++axis) {
        tile_edge_ = std::min(tile_edge_, shape[axis]);
        effective_window = std::max(effective_window,
            std::min(window_edge, 2 * shape[axis] - 1));
    }
    fft_edge_ = 1;
    while (fft_edge_ < tile_edge_ + effective_window)
        fft_edge_ *= 2;
    tile_volume_ = std::size_t(tile_edge_) * tile_edge_ * (thin ? 1 : tile_edge_);
    fft_volume_ = std::size_t(fft_edge_) * fft_edge_ * (thin ? 1 : fft_edge_);
    reversed_.resize(fft_edge_);
    for (std::size_t index = 0; index < reversed_.size(); ++index) {
        std::size_t value = index, reversed = 0;
        for (int size = fft_edge_; size > 1; size /= 2) {
            reversed = (reversed << 1) | (value & 1);
            value >>= 1;
        }
        reversed_[index] = reversed;
    }
    roots_.resize(fft_edge_ / 2);
    for (std::size_t index = 0; index < roots_.size(); ++index) {
        const double angle = -2.0 * std::acos(-1.0) * index / fft_edge_;
        roots_[index] = {std::cos(angle), std::sin(angle)};
    }
    set_threads(threads);
    build_rows(half_angle_degrees);
    build_kernels();
}

void SectorMeanField3D::set_threads(int threads) {
    if (threads < 1)
        throw std::invalid_argument("cached sector threads must be positive");
    threads_ = threads;
}

std::size_t SectorMeanField3D::fft_index(int x, int y, int z) const noexcept {
    return (std::size_t(z) * fft_edge_ + y) * fft_edge_ + x;
}

std::size_t SectorMeanField3D::tile_key(int x, int y, int z) const noexcept {
    const int nx = (shape_[0] + tile_edge_ - 1) / tile_edge_;
    const int ny = (shape_[1] + tile_edge_ - 1) / tile_edge_;
    return (std::size_t(z) * ny + y) * nx + x;
}

void SectorMeanField3D::build_rows(double half_angle_degrees) {
    const double minimum_cosine = std::cos(half_angle_degrees * std::acos(-1.0) / 180.0);
    const int upper = edge_ - lower_ - 1;
    for (DirectionId direction = 1; direction <= 26; ++direction) {
        const auto forward = direction_vector(direction);
        if (thin_ && forward.z != 0)
            continue;
        const double forward_length = std::sqrt(double(squared_length(forward)));
        const int first_z = thin_ ? 0 : std::max(-lower_, 1 - shape_[2]);
        const int last_z = thin_ ? 0 : std::min(upper, shape_[2] - 1);
        const int first_y = std::max(-lower_, 1 - shape_[1]);
        const int last_y = std::min(upper, shape_[1] - 1);
        const int first_x = std::max(-lower_, 1 - shape_[0]);
        const int last_x = std::min(upper, shape_[0] - 1);
        for (int z = first_z; z <= last_z; ++z)
            for (int y = first_y; y <= last_y; ++y) {
                int first = 0;
                bool in_row = false;
                for (int x = first_x; x <= last_x + 1; ++x) {
                    const Vec3i offset{x, y, z};
                    const int squared = squared_length(offset);
                    const bool included = x <= last_x && squared > 0 &&
                        double(dot(offset, forward)) / (forward_length * std::sqrt(double(squared))) +
                        1e-12 >= minimum_cosine;
                    if (included && !in_row) {
                        first = x;
                        in_row = true;
                    }
                    if (!included && in_row) {
                        rows_[direction].push_back({y, z, first, x - 1});
                        in_row = false;
                    }
                }
            }
    }
}

void SectorMeanField3D::transform(std::vector<Complex>& values, bool inverse) const {
    const auto line = [&](std::size_t start, std::size_t stride) {
        std::array<Complex, 128> work;
        for (int index = 0; index < fft_edge_; ++index)
            work[reversed_[index]] = values[start + std::size_t(index) * stride];
        for (int length = 2; length <= fft_edge_; length *= 2)
            for (int begin = 0; begin < fft_edge_; begin += length)
                for (int offset = 0; offset < length / 2; ++offset) {
                    const auto root = roots_[std::size_t(offset) * fft_edge_ / length];
                    const auto odd = work[begin + offset + length / 2] *
                        (inverse ? std::conj(root) : root);
                    const auto even = work[begin + offset];
                    work[begin + offset] = even + odd;
                    work[begin + offset + length / 2] = even - odd;
                }
        const double scale = inverse ? 1.0 / fft_edge_ : 1.0;
        for (int index = 0; index < fft_edge_; ++index)
            values[start + std::size_t(index) * stride] = work[index] * scale;
    };
    const int nz = thin_ ? 1 : fft_edge_;
    for (int z = 0; z < nz; ++z)
        for (int y = 0; y < fft_edge_; ++y)
            line(fft_index(0, y, z), 1);
    for (int z = 0; z < nz; ++z)
        for (int x = 0; x < fft_edge_; ++x)
            line(fft_index(x, 0, z), fft_edge_);
    if (!thin_)
        for (int y = 0; y < fft_edge_; ++y)
            for (int x = 0; x < fft_edge_; ++x)
                line(fft_index(x, y, 0), std::size_t(fft_edge_) * fft_edge_);
}

void SectorMeanField3D::prefix(std::vector<Complex>& values) const {
    const int nz = thin_ ? 1 : fft_edge_;
    for (int z = 0; z < nz; ++z)
        for (int y = 0; y < fft_edge_; ++y) {
            double sum = 0.0;
            for (int x = 0; x < fft_edge_; ++x) {
                const auto i = fft_index(x, y, z);
                const double old = values[i].real();
                values[i] = sum;
                sum += old;
            }
        }
    for (int z = 0; z < nz; ++z)
        for (int x = 0; x < fft_edge_; ++x) {
            double sum = 0.0;
            for (int y = 0; y < fft_edge_; ++y) {
                const auto i = fft_index(x, y, z);
                const double old = values[i].real();
                values[i] = sum;
                sum += old;
            }
        }
    if (!thin_)
        for (int y = 0; y < fft_edge_; ++y)
            for (int x = 0; x < fft_edge_; ++x) {
                double sum = 0.0;
                for (int z = 0; z < fft_edge_; ++z) {
                    const auto i = fft_index(x, y, z);
                    const double old = values[i].real();
                    values[i] = sum;
                    sum += old;
                }
            }
}

void SectorMeanField3D::build_kernels() {
    deterministic_parallel_for(26, threads_, [&](std::size_t index) {
        const auto direction = DirectionId(index + 1);
        if (rows_[direction].empty())
            return;
        auto& kernel = kernels_[direction];
        kernel.assign(fft_volume_, {});
        const auto wrap = [&](int value) {
            return (fft_edge_ - value % fft_edge_) % fft_edge_;
        };
        for (const auto row : rows_[direction])
            for (int high_z = 0; high_z < (thin_ ? 1 : 2); ++high_z)
                for (int high_y = 0; high_y < 2; ++high_y)
                    for (int high_x = 0; high_x < 2; ++high_x) {
                        const int x = high_x ? row.last + 1 : row.first;
                        const int y = row.dy + high_y;
                        const int z = thin_ ? 0 : row.dz + high_z;
                        const int low_count = (1 - high_x) + (1 - high_y) +
                            (thin_ ? 0 : 1 - high_z);
                        kernel[fft_index(wrap(x), wrap(y), wrap(z))] +=
                            low_count % 2 ? -1.0 : 1.0;
                    }
        transform(kernel, false);
    });
}

void SectorMeanField3D::prepare(std::span<const double> resource, double maximum,
                               std::span<const SectorQueryBox3D> boxes, std::uint64_t field_epoch) {
    if (resource.size() != std::size_t(shape_[0]) * shape_[1] * shape_[2] ||
        !std::isfinite(maximum) || !(maximum > 0.0))
        throw std::invalid_argument("invalid cached sector resource field");
    std::set<std::array<int, 3>> requested;
    for (const auto& box : boxes) {
        std::array<int, 3> lo{}, hi{};
        bool empty = false;
        for (int axis = 0; axis < 3; ++axis) {
            lo[axis] = std::clamp(box.lower[axis], 0, shape_[axis]);
            hi[axis] = std::clamp(box.upper[axis], 0, shape_[axis]);
            empty = empty || lo[axis] >= hi[axis];
        }
        if (empty)
            continue;
        for (int z = lo[2] / tile_edge_; z <= (hi[2] - 1) / tile_edge_; ++z)
            for (int y = lo[1] / tile_edge_; y <= (hi[1] - 1) / tile_edge_; ++y)
                for (int x = lo[0] / tile_edge_; x <= (hi[0] - 1) / tile_edge_; ++x)
                    requested.insert({x, y, z});
    }
    if (field_epoch != 0 && field_epoch == cached_epoch_ && requested.size() == tiles_.size()) {
        bool same = true;
        std::size_t index = 0;
        for (const auto coordinate : requested)
            same = same && tiles_[index++].coordinate == coordinate;
        if (same)
            return;
    }
    tiles_.resize(requested.size());
    tile_lookup_.clear();
    std::size_t index = 0;
    for (const auto coordinate : requested) {
        if (tiles_[index].coordinate != coordinate)
            tiles_[index].counts_ready = false;
        tiles_[index].coordinate = coordinate;
        tile_lookup_.emplace(tile_key(coordinate[0], coordinate[1], coordinate[2]), index++);
    }
    if (tiles_.empty()) {
        cached_epoch_ = field_epoch;
        return;
    }
    const int workers = std::max(1, std::min(threads_, int(tiles_.size())));
    workspaces_.resize(workers);
    for (auto& workspace : workspaces_) {
        workspace.field.resize(fft_volume_);
        workspace.mask.resize(fft_volume_);
        workspace.product.resize(fft_volume_);
    }
    deterministic_parallel_for(workers, workers, [&](std::size_t worker) {
        for (std::size_t tile = worker; tile < tiles_.size(); tile += workers)
            prepare_tile(tiles_[tile], workspaces_[worker], resource, maximum);
    });
    cached_epoch_ = field_epoch;
}

void SectorMeanField3D::prepare_tile(Tile& tile, Workspace& workspace,
                                    std::span<const double> resource, double maximum) {
    const int nz = thin_ ? 1 : fft_edge_;
    const int upper = edge_ - lower_ - 1;
    bool interior = true;
    for (int axis = 0; axis < (thin_ ? 2 : 3); ++axis) {
        const int begin = tile.coordinate[axis] * tile_edge_;
        interior = interior && begin - lower_ >= 0 &&
            begin + tile_edge_ - 1 + upper < shape_[axis];
    }
    const bool prepare_counts = !tile.counts_ready && !interior;
    for (int z = 0; z < nz; ++z)
        for (int y = 0; y < fft_edge_; ++y)
            for (int x = 0; x < fft_edge_; ++x) {
                const int shift_x = std::min(lower_, shape_[0] - 1);
                const int shift_y = std::min(lower_, shape_[1] - 1);
                const int shift_z = thin_ ? 0 : std::min(lower_, shape_[2] - 1);
                const int gx = tile.coordinate[0] * tile_edge_ + x - shift_x;
                const int gy = tile.coordinate[1] * tile_edge_ + y - shift_y;
                const int gz = tile.coordinate[2] * tile_edge_ + z - shift_z;
                const bool inside = gx >= 0 && gx < shape_[0] && gy >= 0 && gy < shape_[1] &&
                    gz >= 0 && gz < shape_[2];
                const auto i = fft_index(x, y, z);
                workspace.field[i] = inside ? std::clamp(resource[
                    (std::size_t(gz) * shape_[1] + gy) * shape_[0] + gx] / maximum, 0.0, 1.0) : 0.0;
                if (prepare_counts) workspace.mask[i] = inside ? 1.0 : 0.0;
            }
    prefix(workspace.field);
    transform(workspace.field, false);
    if (prepare_counts) {
        prefix(workspace.mask);
        transform(workspace.mask, false);
    }
    const int tile_z = thin_ ? 1 : tile_edge_;
    const int shift_x = std::min(lower_, shape_[0] - 1);
    const int shift_y = std::min(lower_, shape_[1] - 1);
    const int shift_z = thin_ ? 0 : std::min(lower_, shape_[2] - 1);
    for (DirectionId direction = 1; direction <= 26; ++direction) {
        if (kernels_[direction].empty())
            continue;
        tile.means[direction].resize(tile_volume_);
        tile.counts[direction].resize(tile_volume_);
        if (!tile.counts_ready && interior) {
            std::uint32_t count = 0;
            for (const auto row : rows_[direction])
                count += row.last - row.first + 1;
            std::fill(tile.counts[direction].begin(), tile.counts[direction].end(), count);
        }
        if (prepare_counts) {
            for (std::size_t i = 0; i < fft_volume_; ++i)
                workspace.product[i] = workspace.mask[i] * kernels_[direction][i];
            transform(workspace.product, true);
            for (int z = 0; z < tile_z; ++z)
            for (int y = 0; y < tile_edge_; ++y)
                for (int x = 0; x < tile_edge_; ++x) {
                    const auto local = (std::size_t(z) * tile_edge_ + y) * tile_edge_ + x;
                    const auto i = fft_index(x + shift_x, y + shift_y, z + shift_z);
                    const auto count = std::llround(workspace.product[i].real());
                    if (count < 0 || count > std::int64_t(edge_) * edge_ * edge_)
                        throw std::logic_error("invalid cached sector site count");
                    tile.counts[direction][local] = std::uint32_t(count);
                }
        }
        for (std::size_t i = 0; i < fft_volume_; ++i)
            workspace.product[i] = workspace.field[i] * kernels_[direction][i];
        transform(workspace.product, true);
        for (int z = 0; z < tile_z; ++z)
            for (int y = 0; y < tile_edge_; ++y)
                for (int x = 0; x < tile_edge_; ++x) {
                    const auto local = (std::size_t(z) * tile_edge_ + y) * tile_edge_ + x;
                    const auto count = tile.counts[direction][local];
                    const auto i = fft_index(x + shift_x, y + shift_y, z + shift_z);
                    tile.means[direction][local] = count
                        ? std::clamp(workspace.product[i].real() / count, 0.0, 1.0) : 0.0;
                }
    }
    tile.counts_ready = true;
}

std::pair<std::size_t, std::size_t> SectorMeanField3D::lookup(std::size_t location) const {
    const std::size_t size = std::size_t(shape_[0]) * shape_[1] * shape_[2];
    if (location >= size)
        throw std::out_of_range("cached sector query outside the resource grid");
    const int x = int(location % shape_[0]);
    const int y = int((location / shape_[0]) % shape_[1]);
    const int z = int(location / (std::size_t(shape_[0]) * shape_[1]));
    const auto tile = tile_lookup_.find(tile_key(x / tile_edge_, y / tile_edge_, z / tile_edge_));
    if (tile == tile_lookup_.end())
        throw std::logic_error("cached sector query outside its prepared transport envelope");
    const auto local = (std::size_t(z % tile_edge_) * tile_edge_ + y % tile_edge_) * tile_edge_ + x % tile_edge_;
    return {tile->second, local};
}

double SectorMeanField3D::mean(std::size_t location, DirectionId direction) const {
    if (direction == 0 || direction > 26 || kernels_[direction].empty())
        return 0.0;
    const auto [tile, local] = lookup(location);
    return tiles_[tile].means[direction][local];
}

std::uint32_t SectorMeanField3D::count(std::size_t location, DirectionId direction) const {
    if (direction == 0 || direction > 26 || kernels_[direction].empty())
        return 0;
    const auto [tile, local] = lookup(location);
    return tiles_[tile].counts[direction][local];
}

std::size_t SectorMeanField3D::allocated_bytes() const noexcept {
    std::size_t bytes = reversed_.capacity() * sizeof(std::size_t) + roots_.capacity() * sizeof(Complex);
    for (const auto& kernel : kernels_)
        bytes += kernel.capacity() * sizeof(Complex);
    for (const auto& rows : rows_)
        bytes += rows.capacity() * sizeof(ConeRow);
    for (const auto& tile : tiles_)
        for (DirectionId direction = 1; direction <= 26; ++direction)
            bytes += tile.means[direction].capacity() * sizeof(double) +
                tile.counts[direction].capacity() * sizeof(std::uint32_t);
    for (const auto& workspace : workspaces_)
        bytes += (workspace.field.capacity() + workspace.mask.capacity() + workspace.product.capacity()) * sizeof(Complex);
    return bytes;
}
} // namespace atcg3d::continuum
