#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <optional>
#include <string>

#include "config/output_paths.hpp"
#include "engine/simulation.hpp"
#include "geometry/footprint.hpp"
#include "model/shared_resource_environment.hpp"
#include "model/structured_pde_model.hpp"
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
#include "io/checkpoint_hdf5.hpp"
#endif

namespace {
struct Summary {
    double r{}, K{}, active{}, time{};
    double vascular_length{}, vascular_path_length{}, perfused_volume{}, lesion_perfused_fraction{};
    double vascular_roots{}, vascular_anastomoses{}, tip_branching_rate{}, tip_anastomosis_rate{};
    std::vector<double> radial;
    std::uint64_t checksum{}, resource_checksum{};
    explicit Summary(std::size_t bins) : radial(bins, 0.0) {}
    void add(double radius, double mass) {
        const auto bin = std::min(radial.size() - 1, static_cast<std::size_t>(std::floor(radius)));
        radial[bin] += mass;
    }
    double quantile(double fraction) const {
        if (r + K <= 0.0) return 0.0;
        const double target = fraction * (r + K);
        double sum = 0.0;
        for (std::size_t bin = 0; bin < radial.size(); ++bin) {
            sum += radial[bin];
            if (sum >= target) return static_cast<double>(bin + 1);
        }
        return 0.0;
    }
    void write(std::ostream& out, const std::string& model, std::uint64_t seed) const {
        out << std::setprecision(17) << "{\"schema_version\":1,\"model\":\"" << model << "\",\"seed\":" << seed
            << ",\"time_hours\":" << time << ",\"total_mass\":" << r + K << ",\"r_mass\":" << r
            << ",\"K_mass\":" << K << ",\"r_K_ratio\":" << (K > 0.0 ? r / K : 0.0)
            << ",\"active_fraction\":" << (r > 0.0 ? active / r : 0.0)
            << ",\"r50\":" << quantile(0.5) << ",\"r90\":" << quantile(0.9) << ",\"r99\":" << quantile(0.99)
            << ",\"state_checksum\":" << checksum << ",\"resource_checksum\":" << resource_checksum
            << ",\"vascular_length\":" << vascular_length << ",\"vascular_path_length\":" << vascular_path_length
            << ",\"perfused_volume\":" << perfused_volume << ",\"lesion_perfused_fraction\":" << lesion_perfused_fraction
            << ",\"vascular_roots\":" << vascular_roots << ",\"vascular_anastomoses\":" << vascular_anastomoses
            << ",\"tip_branching_rate\":" << tip_branching_rate << ",\"tip_anastomosis_rate\":" << tip_anastomosis_rate
            << ",\"radial_mass\":[";
        for (std::size_t i = 0; i < radial.size(); ++i) { if (i != 0) out << ','; out << radial[i]; }
        out << "]}\n";
    }
};
void vascular_summary(Summary& summary, const std::vector<double>& vessels, const std::vector<double>& occupied,
                      const atcg3d::continuum::ContinuumModelConfig3D& config) {
    summary.perfused_volume=std::accumulate(vessels.begin(),vessels.end(),0.0);
    const double radius=config.angiogenesis.vessel_radius_voxels;
    const double cross_section=config.base.thin_layer ? 2*radius : std::acos(-1.0)*radius*radius;
    summary.vascular_length=summary.perfused_volume/cross_section;
    std::vector<std::uint8_t> mask(occupied.size(),0);
    if(config.base.thin_layer) {
        atcg3d::continuum::MovingTumorFrontWorkspace2D workspace;
        atcg3d::continuum::build_moving_tumor_front_mask_2d(occupied,config.grid.shape[0],config.grid.shape[1],
            config.nutrient.tumor_front_smoothing_radius_voxels,config.nutrient.tumor_front_density_threshold,workspace,mask);
    } else for(std::size_t i=0;i<mask.size();++i) mask[i]=occupied[i]>0;
    double area=0,perfused=0;
    for(std::size_t i=0;i<mask.size();++i) if(mask[i]) { area+=1;perfused+=vessels[i]; }
    summary.lesion_perfused_fraction=area>0 ? perfused/area : 0;
}
void vascular_abm(Summary& summary,const atcg3d::Simulation3D& simulation,
                  const atcg3d::shared_rules::SharedResourceEnvironment3D& resource,
                  const atcg3d::continuum::ContinuumModelConfig3D& config) {
    const auto geometry=config.shared_vascular_geometry();
    std::vector<double> occupied(resource.vessel_fraction().size(),0),vessels=resource.vessel_fraction();
    const auto index=[&](atcg3d::Vec3i site) {
        const auto x=static_cast<std::size_t>(std::floor(site.x-geometry.origin[0]));
        const auto y=static_cast<std::size_t>(std::floor(site.y-geometry.origin[1]));
        const auto z=geometry.thin_layer ? 0 : static_cast<std::size_t>(std::floor(site.z-geometry.origin[2]));
        return (z*geometry.shape[1]+y)*geometry.shape[0]+x;
    };
    for(const auto slot:simulation.cells().alive_slots()) {
        const auto add=[&](atcg3d::Vec3i site) { if(geometry.contains(site)) occupied[index(site)]+=1; };
        if(simulation.cells().stage(slot)==atcg3d::CellStage::large) for(auto site:atcg3d::large_footprint(simulation.cells().anchor(slot))) add(site);
        else add(simulation.cells().anchor(slot));
    }
    for(auto site:simulation.vessel_grid().occupied_sites()) if(geometry.contains(site)&&simulation.vessel_grid().perfused(site)) vessels[index(site)]=1;
    vascular_summary(summary,vessels,occupied,config);
    summary.vascular_roots = simulation.stats().angiogenesis_roots;
    summary.vascular_anastomoses = simulation.stats().vessel_anastomoses;
    for(const auto slot:simulation.vessel_nodes().alive_slots()) {
        const auto parent=simulation.vessel_nodes().parent_node_slot(slot);
        if(parent==atcg3d::kEmptyVesselNodeSlot) continue;
        const auto delta=simulation.vessel_nodes().position(slot)-simulation.vessel_nodes().position(parent);
        summary.vascular_path_length+=std::sqrt(static_cast<double>(atcg3d::squared_length(delta)));
    }
}
std::size_t radial_bin_count(const atcg3d::continuum::ContinuumModelConfig3D& config) {
    double squared=0;
    for(int edge:config.grid.shape) squared+=static_cast<double>(edge)*edge;
    return static_cast<std::size_t>(std::ceil(std::sqrt(squared)))+1;
}
}  // namespace

int main(int argc, char** argv) {
    try {
        std::filesystem::path yaml, report, output_root, checkpoint, resume;
        std::string model = "abm";
        std::optional<std::uint64_t> seed_override;
        int threads = 1;
        bool dry_run = false;
        for (int i = 1; i < argc; ++i) {
            const std::string argument = argv[i];
            if (argument == "--help") {
                std::cout << "Usage: atcg3d_shared_abm --config YAML [--model abm|pde] [--seed N] [--threads N]\n"
                          << "       [--report JSON] [--output-root PATH] [--checkpoint H5] [--resume-checkpoint H5] [--dry-run]\n";
                return 0;
            }
            if (argument == "--dry-run") { dry_run = true; continue; }
            if (i + 1 >= argc) throw std::invalid_argument("missing argument value");
            const std::string value = argv[++i];
            if (argument == "--config") yaml = value;
            else if (argument == "--model") model = value;
            else if (argument == "--seed") seed_override = std::stoull(value);
            else if (argument == "--threads") threads = std::stoi(value);
            else if (argument == "--report") report = value;
            else if (argument == "--output-root") output_root = value;
            else if (argument == "--checkpoint") checkpoint = value;
            else if (argument == "--resume-checkpoint") resume = value;
            else throw std::invalid_argument("unknown argument: " + argument);
        }
        if (yaml.empty() || (model != "abm" && model != "pde") || threads < 1) throw std::invalid_argument("invalid shared-model arguments");
        if (model == "pde" && (!checkpoint.empty() || !resume.empty())) throw std::invalid_argument("shared checkpoint options require --model abm");
        auto config = atcg3d::structured_pde::StructuredPdeConfig3D::load(yaml);
        if (seed_override) config.continuum.base.seed = *seed_override;
        const auto seed = config.continuum.base.seed;
        config.continuum.base.threads = threads;
        config.continuum.output.enabled = false;
        auto base = atcg3d::shared_rules::abm_config(config);
        if (dry_run) { std::cout << config.to_json() << '\n'; return 0; }
        auto environment = std::make_unique<atcg3d::shared_rules::SharedResourceEnvironment3D>(config);
        auto* field = environment.get();
        if (!resume.empty()) {
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
            base.run_mode = "resume";
            base.resume_checkpoint = resume;
            const auto saved = atcg3d::read_hdf5_checkpoint(resume, base);
            field->load_checkpoint(resume.string() + ".resource.bin", saved.state_checksum);
            atcg3d::Simulation3D restored(base, std::move(environment));
            restored.restore(saved.cells, saved.next_uid, saved.clock, saved.stats, saved.lineage,
                saved.vasculature, saved.cell_slot_count, saved.cell_slots, saved.cell_free_slots);
            restored.run();
            const auto bins = radial_bin_count(config.continuum);
            Summary result(bins);
            result.time = restored.clock().time_hours;
            for (const auto slot : restored.cells().alive_slots()) {
                const auto point = restored.cells().anchor(slot);
                (restored.cells().type(slot) == atcg3d::CellType::r ? result.r : result.K) += 1.0;
                if (restored.cells().type(slot) == atcg3d::CellType::r && (restored.cells().flags(slot) & atcg3d::kMigrationActive)) result.active += 1.0;
                result.add(std::sqrt((point.x + 0.5) * (point.x + 0.5) + (point.y + 0.5) * (point.y + 0.5) +
                    (base.thin_layer ? 0.0 : (point.z + 0.5) * (point.z + 0.5))), 1.0);
            }
            vascular_abm(result,restored,*field,config.continuum);
            result.checksum = restored.state_checksum(); result.resource_checksum = field->field_checksum();
            if (report.empty()) result.write(std::cout, model, seed);
            else {
                if (!report.parent_path().empty()) std::filesystem::create_directories(report.parent_path());
                std::ofstream out(report);
                if (!out) throw std::runtime_error("unable to write summary report");
                result.write(out, model, seed);
            }
            return 0;
#else
            throw std::runtime_error("shared ABM restart requires the HDF5 checkpoint build");
#endif
        }
        atcg3d::Simulation3D simulation(base, std::move(environment));
        simulation.initialize();
        const auto bins = radial_bin_count(config.continuum);
        Summary result(bins);
        if (model == "abm") {
            simulation.run();
            result.time = simulation.clock().time_hours;
            for (const auto slot : simulation.cells().alive_slots()) {
                const auto point = simulation.cells().anchor(slot);
                (simulation.cells().type(slot) == atcg3d::CellType::r ? result.r : result.K) += 1.0;
                if (simulation.cells().type(slot) == atcg3d::CellType::r && (simulation.cells().flags(slot) & atcg3d::kMigrationActive)) result.active += 1.0;
                result.add(std::sqrt((point.x + 0.5) * (point.x + 0.5) + (point.y + 0.5) * (point.y + 0.5) +
                    (base.thin_layer ? 0.0 : (point.z + 0.5) * (point.z + 0.5))), 1.0);
            }
            vascular_abm(result,simulation,*field,config.continuum);
            result.checksum = simulation.state_checksum(); result.resource_checksum = field->field_checksum();
            if (!checkpoint.empty()) {
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
                atcg3d::write_hdf5_checkpoint(checkpoint, simulation);
                field->save_checkpoint(checkpoint.string() + ".resource.bin", simulation.state_checksum());
#else
                throw std::runtime_error("shared ABM checkpoint requires the HDF5 build");
#endif
            }
        } else {
            atcg3d::structured_pde::StructuredPdeModel3D pde(config);
            pde.initialize_from_abm(simulation);
            while (pde.step()) {}
            std::vector<double> occupied(pde.voxel_count(),0);
            for(std::size_t i=0;i<pde.voxel_count();++i) occupied[i]=pde.occupied_fraction(i);
            vascular_summary(result,pde.vessel_fraction(),occupied,config.continuum);
            if (const auto* vascular = pde.angiogenesis()) {
                const auto& law = config.continuum.angiogenesis;
                for (std::size_t i = 0; i < pde.voxel_count(); ++i) {
                    const double tips = vascular->tips()[i];
                    result.tip_branching_rate += law.tip_branching_per_hour *
                        vascular->taf()[i] * tips * pde.voxel_measure();
                    result.tip_anastomosis_rate += law.tip_anastomosis_per_hour *
                        (pde.vessel_fraction()[i] + tips) * tips * pde.voxel_measure();
                }
            }
            const auto diagnostics = pde.diagnostics();
            result.r = diagnostics.r_total; result.K = diagnostics.K_total; result.active = diagnostics.r_active_total;
            result.time = pde.time_hours(); result.checksum = pde.state_checksum();
            for (std::size_t here = 0; here < pde.voxel_count(); ++here) {
                const auto point = pde.coordinate(here);
                double mass = 0.0;
                for (const auto stage : {atcg3d::structured_pde::StructuredStage3D::small, atcg3d::structured_pde::StructuredStage3D::large}) {
                    mass += pde.r_normal(stage, here) + pde.r_active(stage, here) + pde.K(stage, here);
                }
                result.add(std::sqrt(point[0] * point[0] + point[1] * point[1] + (base.thin_layer ? 0.0 : point[2] * point[2])), mass * pde.voxel_measure());
            }
        }
        if (report.empty()) {
            const auto directory = atcg3d::resolve_output_directory(config.continuum.output.directory.string() + "_" + model, output_root);
            std::filesystem::create_directories(directory);
            report = directory / "summary.json";
        } else if (!report.parent_path().empty()) std::filesystem::create_directories(report.parent_path());
        std::ofstream out(report);
        if (!out) throw std::runtime_error("unable to write summary report");
        result.write(out, model, seed);
        std::cout << "model=" << model << " time_hours=" << result.time << " mass=" << result.r + result.K << " checksum=" << result.checksum << '\n';
        return 0;
    } catch (const std::exception& error) { std::cerr << "shared rules: " << error.what() << '\n'; return 1; }
}
