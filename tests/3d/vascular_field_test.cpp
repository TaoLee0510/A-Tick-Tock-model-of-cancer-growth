#include <algorithm>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <memory>
#include <sstream>

#include "engine/simulation.hpp"
#include "model/angiogenesis_field.hpp"
#include "model/continuum_model.hpp"
#include "model/shared_resource_environment.hpp"
#include "model/structured_pde_model.hpp"
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
#include "io/checkpoint_hdf5.hpp"
#endif

namespace {
class UniformResource final : public atcg3d::EnvironmentCoupling3D {
public:
    explicit UniformResource(double value):value_(value) {}
    double retained_density(atcg3d::Vec3i) const noexcept override { return 1; }
    double normalized_resource(atcg3d::Vec3i) const noexcept override { return value_; }
    atcg3d::EnvironmentInitializationResult3D initialize(double,const atcg3d::CellStore3D&,const atcg3d::SparseVesselGrid3D&) override { return {}; }
    void refresh(double,const atcg3d::CellStore3D&,const atcg3d::SparseVesselGrid3D&) override {}
    double next_refresh_time_hours() const noexcept override { return 1.0e30; }
    std::uint32_t schedule_generation() const noexcept override { return 1; }
    std::uint64_t refresh_count() const noexcept override { return 0; }
    std::size_t allocated_bytes() const noexcept override { return 0; }
    std::uint64_t field_checksum() const noexcept override { return 0; }
private: double value_;
};
}
int main() {
    using namespace atcg3d;
    using namespace atcg3d::continuum;
    auto config=Model3DConfig::load(std::filesystem::path(ATCG_SOURCE_DIR)/"configs/shared_angiogenesis_r200_v1.yaml");
    config.angiogenesis.outward_speed_policy="strict_supremum_v1";
    bool rejected=false;
    try { config.validate(); } catch(const std::invalid_argument&) { rejected=true; }
    assert(rejected);
    config.angiogenesis.outward_speed_policy="mean_path_speed_v1";
    config.angiogenesis.outward_speed_voxels_per_hour=70;
    config.validate(); // Less than the supremum but above the nominal mean.
    config.angiogenesis.outward_speed_policy="disabled";
    config.angiogenesis.outward_speed_voxels_per_hour=2;
    config.validate();
    config.activated_r_migration_rate_model="beta";
    config.activated_r_migration_beta={1,1,1,true,0.5,0};
    config.angiogenesis.outward_speed_policy="mean_path_speed_v1";
    config.angiogenesis.inward_speed_voxels_per_hour=0.1;
    config.angiogenesis.outward_speed_voxels_per_hour=0.453;
    config.validate(); // Clamped uniform mean 3/8, step length (1+sqrt(2))/2.
    config.angiogenesis.outward_speed_voxels_per_hour=0.452;
    rejected=false; try { config.validate(); } catch(const std::invalid_argument&) { rejected=true; }
    assert(rejected);

    config=Model3DConfig::load(std::filesystem::path(ATCG_SOURCE_DIR)/"configs/shared_angiogenesis_r20_v1.yaml");
    config.angiogenesis.seed_rate_sites_per_30_days=720000;
    config.angiogenesis.seed_rate_sites_per_hour=1000;
    config.angiogenesis.max_total_roots=1;
    config.angiogenesis.surface_min_separation_voxels=0;
    config.output_enabled=false;config.end_time_hours=1;
    Simulation3D hypoxic(config,std::make_unique<UniformResource>(0));hypoxic.initialize();hypoxic.run();
    Simulation3D oxic(config,std::make_unique<UniformResource>(1));oxic.initialize();oxic.run();
    assert(hypoxic.stats().angiogenesis_roots==1);
    assert(oxic.stats().angiogenesis_roots==0);
    for(const auto slot:hypoxic.vessel_tips().alive_slots()) assert(hypoxic.vessel_tips().bias_axis(slot).z==0);

    AngiogenesisFieldConfig3D law;law.model="vegf_tip_density_v1";
    law.tip_branching_per_hour=law.tip_anastomosis_per_hour=0;law.seed_tips_per_hour=1;
    AngiogenesisField3D field(law,{9,9,1},1,true),warm(law,{9,9,1},1,true);
    std::vector<double> cells(81,0),N(81,0),empty(81,0),full(81,1);cells[40]=1;
    field.initialize(empty);warm.initialize(empty);
    field.advance(1,cells,N,1);warm.advance(1,cells,full,1);
    assert(field.taf()[40]>field.taf()[39] && field.taf()[39]>0);
    assert(field.perfused_volume()>0 && field.vessel_length()>0);
    assert(warm.perfused_volume()==0 && std::all_of(warm.taf().begin(),warm.taf().end(),[](double v){return v==0;}));
    for(double v:field.vessels()) assert(std::isfinite(v)&&v>=0&&v<=1);
    std::stringstream stream(std::ios::in|std::ios::out|std::ios::binary);field.save(stream);
    AngiogenesisField3D copy(law,{9,9,1},1,true);copy.load(stream);
    assert(copy.checksum()==field.checksum());

    auto nonlinear=law;nonlinear.tip_branching_per_hour=5;nonlinear.tip_anastomosis_per_hour=2;
    nonlinear.seed_tips_per_hour=10;nonlinear.maximum_tip_density=0.1;
    AngiogenesisField3D branches(nonlinear,{9,9,1},1,true);branches.initialize(empty);
    branches.advance(2,cells,N,1);
    for(double v:branches.tips()) assert(std::isfinite(v)&&v>=0&&v<=0.1);
    assert(branches.perfused_volume()>0);
    field.advance(1,cells,N,1);copy.advance(1,cells,N,1);
    assert(copy.checksum()==field.checksum());

    auto coupled=structured_pde::StructuredPdeConfig3D::load(std::filesystem::path(ATCG_SOURCE_DIR)/"ATCG3D_SharedRules/config/angiogenesis_v8.yaml");
    auto base=shared_rules::abm_config(coupled);
    Simulation3D source(base,std::make_unique<shared_rules::SharedResourceEnvironment3D>(coupled));source.initialize();
    structured_pde::StructuredPdeModel3D pde(coupled);pde.initialize_from_abm(source);
    for(int i=0;i<7;++i) assert(pde.step());
    const auto checkpoint=std::filesystem::temp_directory_path()/"atcg_vascular_field_pde.bin";std::filesystem::remove(checkpoint);
    pde.save_checkpoint(checkpoint);
    coupled.continuum.base.threads=4;
    structured_pde::StructuredPdeModel3D resumed(coupled);resumed.load_checkpoint(checkpoint);
    assert(pde.state_checksum()==resumed.state_checksum());
    while(pde.step()) { assert(resumed.step());assert(pde.state_checksum()==resumed.state_checksum()); }
    assert(pde.angiogenesis()->perfused_volume()>0);
    std::filesystem::remove(checkpoint);

    // The four-field continuum uses the same vascular solver and restart.
    ContinuumModel3D continuum(coupled.continuum);continuum.initialize_from_abm(source);
    for(int i=0;i<5;++i) assert(continuum.step());
    const auto other=std::filesystem::temp_directory_path()/"atcg_vascular_continuum.bin";std::filesystem::remove(other);
    continuum.save_checkpoint(other);ContinuumModel3D restored(coupled.continuum);restored.load_checkpoint(other);
    assert(continuum.state_checksum()==restored.state_checksum());
    while(continuum.step()) { assert(restored.step());assert(continuum.state_checksum()==restored.state_checksum()); }
    std::filesystem::remove(other);

    coupled.continuum.base.threads=1;
    auto resource=std::make_unique<shared_rules::SharedResourceEnvironment3D>(coupled);
    auto* original_field=resource.get();
    Simulation3D original(shared_rules::abm_config(coupled),std::move(resource));original.initialize();
    for(int i=0;i<100;++i) assert(original.step());
    const auto sidecar=std::filesystem::temp_directory_path()/"atcg_vascular_abm.resource.bin";std::filesystem::remove(sidecar);
    original_field->save_checkpoint(sidecar,original.state_checksum());
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
    const auto hdf5=std::filesystem::temp_directory_path()/"atcg_hypoxic_checkpoint.h5";std::filesystem::remove(hdf5);
    write_hdf5_checkpoint(hdf5,original);
    const auto loaded=read_hdf5_checkpoint(hdf5,shared_rules::abm_config(coupled));
    assert(loaded.state_checksum==original.state_checksum());
    std::filesystem::remove(hdf5);
#endif
    coupled.continuum.base.threads=4;
    auto reloaded=std::make_unique<shared_rules::SharedResourceEnvironment3D>(coupled);auto* restored_field=reloaded.get();
    restored_field->load_checkpoint(sidecar,original.state_checksum());
    Simulation3D restarted(shared_rules::abm_config(coupled),std::move(reloaded));
    restarted.restore(original.snapshot_cells(),original.next_uid(),original.clock(),original.stats(),original.lineage(),
        original.snapshot_vasculature(),original.cells().slot_count(),original.snapshot_cell_slots(),original.cells().free_slots());
    assert(original.state_checksum()==restarted.state_checksum());
    std::size_t continued=0;
    while(original.step()) {
        assert(restarted.step());
        assert(original.state_checksum()==restarted.state_checksum());
        assert(original_field->field_checksum()==restored_field->field_checksum());
        ++continued;
    }
    assert(!restarted.step() && continued>=10);
    std::filesystem::remove(sidecar);
}
