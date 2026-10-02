#pragma once

#include <array>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <vector>

namespace atcg3d::continuum {
struct AngiogenesisFieldConfig3D {
    std::string model{"disabled"};
    double hypoxia_threshold{0.30};
    double taf_production_per_cell_hour{1.0};
    double taf_diffusion_voxels2_per_hour{1.0};
    double taf_decay_per_hour{0.1};
    double tip_diffusion_voxels2_per_hour{0.05};
    double tip_chemotaxis{2.0};
    double tip_branching_per_hour{0.1};
    double tip_anastomosis_per_hour{0.2};
    double seed_tips_per_hour{0.02};
    double tip_speed_voxels_per_hour{0.5};
    double vessel_radius_voxels{1.0};
    double perfusion_exchange_per_hour{2.0};
    double exclusion_fraction{0.5};
    double maximum_tip_density{1.0};
    void validate() const;
    std::string to_json() const;
    std::uint64_t fingerprint() const;
};

// Reflecting finite-volume VEGF and upwind tip transport. Vessel density is
// union volume fraction; its deposited volume/cross-section estimates length.
class AngiogenesisField3D {
public:
    AngiogenesisField3D(AngiogenesisFieldConfig3D config,std::array<int,3> shape,double spacing,bool thin);
    void initialize(std::vector<double> vessels);
    void advance(double dt,const std::vector<double>& cells,const std::vector<double>& nutrient,double nutrient_maximum,bool grow_vessels=true);
    const std::vector<double>& taf() const noexcept { return taf_; }
    const std::vector<double>& tips() const noexcept { return tips_; }
    const std::vector<double>& vessels() const noexcept { return vessels_; }
    double vessel_length() const;
    double perfused_volume() const;
    double lesion_perfused_fraction(const std::vector<double>& cells) const;
    std::uint64_t checksum() const;
    void save(std::ostream& out) const;
    void load(std::istream& in);
private:
    std::size_t index(int x,int y,int z) const noexcept;
    AngiogenesisFieldConfig3D config_;
    std::array<int,3> shape_;
    double spacing_{},measure_{},cross_section_{};
    bool thin_{};
    std::vector<double> taf_,tips_,vessels_,work_taf_,work_tips_;
};
}  // namespace atcg3d::continuum
