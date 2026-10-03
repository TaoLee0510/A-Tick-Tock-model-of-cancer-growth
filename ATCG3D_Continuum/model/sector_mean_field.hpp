#pragma once

#include "core/types.hpp"
#include <array>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <span>
#include <unordered_map>
#include <vector>

namespace atcg3d::continuum {

struct SectorQueryBox3D {
    std::array<int, 3> lower{}, upper{};
};

// A field-epoch cache, using summed-volume prefixes and precomputed cone-row
// corner kernels. Overlap-save tiles materialize every requested mean before
// transport; migration queries perform a lookup independent of cone volume.
class SectorMeanField3D {
public:
    SectorMeanField3D(std::array<int, 3> shape, bool thin, int window_edge,
                      double half_angle_degrees, int threads);
    void prepare(std::span<const double> resource, double maximum,
                 std::span<const SectorQueryBox3D> boxes, std::uint64_t field_epoch = 0);
    double mean(std::size_t location, DirectionId direction) const;
    std::uint32_t count(std::size_t location, DirectionId direction) const;
    std::size_t allocated_bytes() const noexcept;
    std::size_t prepared_tiles() const noexcept { return tiles_.size(); }
    void set_threads(int threads);

private:
    using Complex = std::complex<double>;
    struct ConeRow {
        int dy{}, dz{}, first{}, last{};
    };
    struct Tile {
        std::array<int, 3> coordinate{};
        bool counts_ready{};
        std::array<std::vector<double>, 27> means;
        std::array<std::vector<std::uint32_t>, 27> counts;
    };
    struct Workspace {
        std::vector<Complex> field, mask, product;
    };
    void build_rows(double half_angle_degrees);
    void build_kernels();
    void transform(std::vector<Complex>& values, bool inverse) const;
    void prefix(std::vector<Complex>& values) const;
    void prepare_tile(Tile& tile, Workspace& workspace,
                      std::span<const double> resource, double maximum);
    std::pair<std::size_t, std::size_t> lookup(std::size_t location) const;
    std::size_t fft_index(int x, int y, int z) const noexcept;
    std::size_t tile_key(int x, int y, int z) const noexcept;
    std::array<int, 3> shape_{};
    bool thin_{};
    int edge_{}, lower_{}, fft_edge_{}, threads_{};
    int tile_edge_{32};
    std::size_t tile_volume_{}, fft_volume_{};
    std::array<std::vector<ConeRow>, 27> rows_;
    std::array<std::vector<Complex>, 27> kernels_;
    std::vector<std::size_t> reversed_;
    std::vector<Complex> roots_;
    std::vector<Tile> tiles_;
    std::unordered_map<std::size_t, std::size_t> tile_lookup_;
    std::vector<Workspace> workspaces_;
    std::uint64_t cached_epoch_{};
};
} // namespace atcg3d::continuum
