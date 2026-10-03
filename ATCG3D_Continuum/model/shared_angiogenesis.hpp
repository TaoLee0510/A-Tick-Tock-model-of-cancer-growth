#pragma once

#include <array>
#include <cstdint>
#include <iosfwd>
#include <vector>

#include "model/angiogenesis_field.hpp"

namespace atcg3d::continuum {

struct SharedVascularDiagnostics3D {
    double centerline_growth{};
    double initial_perfused_volume{};
    double seeded_tips{};
    double branches{};
    double anastomoses{};
    double discarded_taf_mass{};
    double discarded_tip_mass{};
    std::size_t active_voxels{};
};

// Versioned lattice-tip law shared by individual and density representations.
// The published field implementation does not call this class.
class SharedAngiogenesis3D {
public:
    SharedAngiogenesis3D(AngiogenesisFieldConfig3D config,
                        std::array<int, 3> shape,
                        double spacing, bool thin, int threads);
    void use_individual_tips(std::uint64_t seed);
    void initialize(std::vector<double> vessels);
    void initialize_tip_density(const std::vector<double>& tips);
    void advance(double dt, const std::vector<double>& cells,
                 const std::vector<double>& nutrient, double maximum,
                 bool grow, const VascularConsumerBounds3D* consumer_bounds = nullptr);
    const std::vector<double>& taf() const noexcept { return taf_; }
    const std::vector<double>& tips() const noexcept { return tips_; }
    const std::vector<double>& vessels() const noexcept { return vessels_; }
    const SharedVascularDiagnostics3D& diagnostics() const noexcept { return diagnostics_; }
    std::size_t allocated_bytes() const noexcept;
    std::uint64_t checksum() const;
    void save(std::ostream& out) const;
    void load(std::istream& in);

private:
    using Bounds = VascularConsumerBounds3D;
    struct Tip { std::uint64_t uid{}; std::size_t location{}; };
    std::size_t index(int x, int y, int z) const noexcept;
    std::array<int, 3> coordinate(std::size_t location) const noexcept;
    bool neighbor(std::size_t here, int face, std::size_t& other) const noexcept;
    void include(Bounds& bounds, std::size_t location) const noexcept;
    Bounds expanded(Bounds bounds) const noexcept;
    std::size_t count(const Bounds& bounds) const noexcept;
    std::size_t bounded_location(const Bounds& bounds, std::size_t offset) const noexcept;
    Bounds prepare_sources(const std::vector<double>& cells,
                           const std::vector<double>& nutrient, double maximum,
                           const Bounds* consumer_bounds);
    double jump_rate(std::size_t here, std::size_t other) const noexcept;
    double outgoing_rate(std::size_t here) const noexcept;
    double substep(double remaining, const Bounds& bounds);
    void advance_taf(double dt, const Bounds& bounds);
    void advance_density_tips(double dt, const Bounds& bounds);
    void advance_individual_tips(double dt, const Bounds& bounds);
    void deposit_centerlines(double dt, const Bounds& bounds);
    void reduce_step(const Bounds& bounds);
    void clear_previous_workspace(const Bounds& bounds);
    double uniform();
    std::uint64_t poisson(double mean);
    std::size_t seed_location();

    AngiogenesisFieldConfig3D config_;
    std::array<int, 3> shape_;
    double spacing_{}, measure_{}, cross_section_{};
    int dimensions_{}, threads_{};
    bool individual_{};
    std::uint64_t seed_{}, random_counter_{}, next_uid_{1};
    std::vector<Tip> individuals_, next_individuals_;
    std::vector<double> taf_, tips_, vessels_, work_taf_, work_tips_;
    std::vector<double> hypoxia_, surface_, rates_, moved_, branch_counts_, connection_counts_;
    std::vector<double> discarded_taf_, discarded_tips_;
    std::array<std::vector<double>, 3> edges_, work_edges_;
    std::vector<std::size_t> seed_sites_;
    std::vector<double> seed_prefix_;
    double surface_sum_{}, seed_rate_{};
    Bounds support_, workspace_bounds_, previous_sources_;
    SharedVascularDiagnostics3D diagnostics_;
};

}  // namespace atcg3d::continuum
