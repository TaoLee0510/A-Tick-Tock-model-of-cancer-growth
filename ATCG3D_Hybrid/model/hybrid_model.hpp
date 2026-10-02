#pragma once

#include "engine/simulation.hpp"
#include "model/shared_resource_environment.hpp"
#include "model/structured_pde_model.hpp"
#include <memory>

namespace atcg3d::hybrid {

struct HybridConfig3D {
    int schema_version{1};
    std::string model{"hybrid_shared_grid_v1"}, mode{"adaptive"};
    structured_pde::StructuredPdeConfig3D rules;
    std::filesystem::path output_directory{"atcg3d_hybrid_run"};
    double exchange_every_hours{1.0};
    int smoothing_radius{3};
    double core_on{0.5}, core_off{0.3};
    static HybridConfig3D load(const std::filesystem::path &path);
    void validate() const;
    std::uint64_t fingerprint() const;
};

struct HybridDiagnostics3D {
    double abm_mass{}, pde_mass{}, total_mass{}, active_mass{}, mean_nutrient{};
    double r_mass{}, K_mass{};
    std::uint64_t to_pde{}, to_abm{}, exchanges{};
};

class HybridModel3D {
  public:
    explicit HybridModel3D(HybridConfig3D config);
    void initialize();
    void
    initialize_from_arrays(structured_pde::StructuredInitialFields3D fields);
    bool step();
    double time_hours() const noexcept;
    HybridDiagnostics3D diagnostics() const;
    std::vector<double> radial_mass() const;
    std::uint64_t state_checksum() const;
    const Simulation3D &abm() const noexcept { return *abm_; }
    const structured_pde::StructuredPdeModel3D &pde() const noexcept {
        return *pde_;
    }
    void exchange();
    void save_checkpoint(const std::filesystem::path &path) const;
    void load_checkpoint(const std::filesystem::path &path);

  private:
    std::size_t location(Vec3i site) const;
    Vec3i site(std::size_t index) const;
    std::vector<Vec3i> footprint(Vec3i anchor, CellStage stage) const;
    void attach_coupling();
    void assemble();
    void synchronize_environment(bool refresh_rates);
    void remove_agent(Slot slot);
    void classify_core();
    void convert_to_density();
    void convert_to_agents();
    void convert_region(const std::vector<std::size_t> &region);
    void canonicalize();
    bool available(Vec3i point) const;
    HybridConfig3D config_;
    std::unique_ptr<Simulation3D> abm_;
    std::unique_ptr<structured_pde::StructuredPdeModel3D> pde_;
    shared_rules::SharedResourceEnvironment3D *environment_{};
    std::vector<std::uint8_t> core_;
    std::vector<double> own_r_, own_K_, own_occupied_;
    std::uint64_t exchange_count_{}, to_pde_{}, to_abm_{};
    bool initialized_{};
};
} // namespace atcg3d::hybrid
