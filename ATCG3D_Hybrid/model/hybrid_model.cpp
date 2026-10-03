#include "model/hybrid_model.hpp"

#include "geometry/footprint.hpp"
#include "rules/initial_rates.hpp"
#include <algorithm>
#include <bit>
#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>
#include <queue>
#include <unordered_set>
#include <yaml-cpp/yaml.h>
#ifdef ATCG3D_HAS_HDF5_CHECKPOINT
#include "io/checkpoint_hdf5.hpp"
#endif

namespace atcg3d::hybrid {
namespace {
std::uint64_t mix(std::uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}
template <class T> void put(std::ostream &out, const T &value) {
    out.write(reinterpret_cast<const char *>(&value), sizeof(T));
}
template <class T> void get(std::istream &in, T &value) {
    in.read(reinterpret_cast<char *>(&value), sizeof(T));
    if (!in)
        throw std::runtime_error("truncated hybrid checkpoint");
}
template <class T> void put_vector(std::ostream &out, const std::vector<T> &v) {
    put(out, std::uint64_t(v.size()));
    for (const auto &x : v)
        put(out, x);
}
template <class T> void get_vector(std::istream &in, std::vector<T> &v) {
    std::uint64_t n;
    get(in, n);
    if (n > 100000000)
        throw std::runtime_error("invalid hybrid checkpoint size");
    v.resize(n);
    for (auto &x : v)
        get(in, x);
}
void add_mixture(std::vector<double>& target, const std::vector<double>& source, double scale) {
    if (source.empty() || !(scale > 0.0)) return;
    if (target.empty()) target.resize(source.size());
    if (target.size() != source.size()) throw std::logic_error("hybrid distribution grid mismatch");
    for (std::size_t i = 0; i < source.size(); ++i) target[i] += source[i] * scale;
}
std::size_t sample_mixture(const std::vector<double>& values, std::uint64_t key) {
    const double mass = std::accumulate(values.begin(), values.end(), 0.0);
    if (!(mass > 0.0)) throw std::logic_error("hybrid cell lacks a clock distribution");
    double target = (mix(key) >> 11) * 0x1.0p-53 * mass;
    for (std::size_t i = 0; i < values.size(); ++i) {
        target -= values[i];
        if (target < 0.0) return i;
    }
    return values.size() - 1;
}
void debit_mixture(std::vector<double>& values, double fraction) {
    for (auto& value : values) value *= std::max(0.0, 1.0 - fraction);
}
} // namespace
HybridConfig3D HybridConfig3D::load(const std::filesystem::path &path) {
    auto y = YAML::LoadFile(path.string());
    HybridConfig3D c;
    if (y["schema"]["name"].as<std::string>() != "atcg3d.hybrid_config")
        throw std::invalid_argument("invalid hybrid schema");
    c.schema_version = y["schema"]["version"].as<int>();
    c.model = y["model"].as<std::string>();
    c.rules = structured_pde::StructuredPdeConfig3D::load(
        path.parent_path() / y["structured_config"].as<std::string>());
    if (y["output"] && y["output"]["directory"])
        c.output_directory = y["output"]["directory"].as<std::string>();
    auto h = y["hybrid"];
    c.mode = h["mode"].as<std::string>();
    c.exchange_every_hours = h["exchange_every_hours"].as<double>();
    c.smoothing_radius = h["smoothing_radius_voxels"].as<int>();
    c.core_on = h["core_on"].as<double>();
    c.core_off = h["core_off"].as<double>();
    if (c.schema_version >= 4) {
        c.front_band_voxels = h["front_band_voxels"].as<int>();
        c.front_hysteresis_voxels = h["front_hysteresis_voxels"].as<int>();
        c.active_guard_voxels = h["active_guard_voxels"].as<int>();
        c.gradient_on = h["nutrient_gradient_on"].as<double>();
        c.gradient_off = h["nutrient_gradient_off"].as<double>();
    }
    c.validate();
    return c;
}
void HybridConfig3D::validate() const {
    rules.validate();
    const bool distributions =
        (schema_version == 2 && model == "hybrid_distributions_v2") ||
        (schema_version == 3 && model == "hybrid_volume_coupling_v3") ||
        (schema_version == 4 && model == "hybrid_invasion_front_v4");
    if (schema_version >= 3 && rules.schema_version < 13)
        throw std::invalid_argument("hybrid v3 requires structured schema v13+");
    if (!distributions && rules.migration.activation_clock != "beta_mean_remaining_cycle_v1")
        throw std::invalid_argument("wrapper v1 does not carry activation duration distributions");
    if (!distributions && rules.division_clock_model != "mean_rate_v1")
        throw std::invalid_argument("hybrid v1 does not carry division work distributions");
    if (distributions && rules.division_clock_model != "transported_shifted_geometric_v1")
        throw std::invalid_argument("hybrid v2 requires transported division work");
    (void)shared_rules::abm_config(rules);
    if ((!distributions && (schema_version != 1 || model != "hybrid_shared_grid_v1")) ||
        (mode != "adaptive" && mode != "all_abm" && mode != "all_pde"))
        throw std::invalid_argument("unsupported hybrid model");
    if (schema_version >= 4 &&
        (rules.schema_version < 14 || front_band_voxels < 0 ||
         front_hysteresis_voxels < 0 || active_guard_voxels < 0 ||
         rules.continuum.base.ultrasmall_enabled ||
         !std::isfinite(gradient_on) || !std::isfinite(gradient_off) ||
         gradient_off < 0.0 || gradient_on <= gradient_off))
        throw std::invalid_argument("invalid hybrid invasion-front policy");
    const double ticks = exchange_every_hours / rules.continuum.time_step_hours;
    if (!std::isfinite(ticks) || ticks < 1 ||
        std::abs(ticks - std::round(ticks)) > 1e-10 || smoothing_radius < 0 ||
        core_off < 0 || core_on > 1 || core_off >= core_on ||
        rules.schema_version < 7)
        throw std::invalid_argument(
            "invalid hybrid exchange or core thresholds");
}
std::uint64_t HybridConfig3D::fingerprint() const {
    auto h = rules.dynamics_fingerprint();
    for (char c : model + mode)
        h = mix(h ^ std::uint8_t(c));
    for (double v :
         {exchange_every_hours, double(smoothing_radius), core_on, core_off})
        h = mix(h ^ std::bit_cast<std::uint64_t>(v));
    if (schema_version >= 4)
        for (double v : {double(front_band_voxels), double(front_hysteresis_voxels),
                         double(active_guard_voxels), gradient_on, gradient_off})
            h = mix(h ^ std::bit_cast<std::uint64_t>(v));
    return h;
}
HybridModel3D::HybridModel3D(HybridConfig3D config)
    : config_(std::move(config)) {
    config_.validate();
    auto density_rules = config_.rules;
    if (config_.schema_version >= 4 && config_.mode == "adaptive") {
        density_rules.operator_model = "shared_operator_switches_v1";
        density_rules.activation_operator_enabled = false;
    }
    pde_ = std::make_unique<structured_pde::StructuredPdeModel3D>(density_rules);
    auto env = std::make_unique<shared_rules::SharedResourceEnvironment3D>(
        config_.rules);
    environment_ = env.get();
    auto base = shared_rules::abm_config(config_.rules);
    if (config_.sample_every_hours > 0.0) {
        base.output_enabled = true;
        base.preview_every_hours = config_.sample_every_hours;
        base.full_every_hours = config_.sample_every_hours;
        base.checkpoint_every_hours = config_.sample_every_hours;
        base.preview_keyframe_every_hours = config_.sample_every_hours;
        base.full_keyframe_every_hours = std::max(base.full_keyframe_every_hours, config_.sample_every_hours);
        base.checkpoint_base_every_hours = std::max(base.checkpoint_base_every_hours, config_.sample_every_hours);
    }
    if (config_.mode != "all_abm")
        base.angiogenesis.enabled = false;
    if (config_.mode == "adaptive") {
        environment_->externally_driven_ = true;
        if (config_.schema_version >= 4 &&
            config_.rules.continuum.angiogenesis.model == "shared_vegf_lattice_v2")
            environment_->angiogenesis_.reset();
        environment_->next_refresh_ = 1.0e100;
        environment_->external_destination_ = [this](Vec3i p) {
            return available(p, config_.schema_version >= 4);
        };
        environment_->external_counts_ = [this](Vec3i p) {
            std::array<double, 2> counts{};
            const int edge =
                config_.rules.continuum.base.growth_density_window_edge;
            const int lower = (edge - 1) / 2, upper = edge - lower - 1;
            if (config_.schema_version >= 4) {
                const auto values = box_counts(p, lower, upper);
                return std::array<double, 2>{values[0], values[1]};
            }
            for (int z = config_.rules.continuum.base.thin_layer ? 0 : -lower;
                 z <= (config_.rules.continuum.base.thin_layer ? 0 : upper);
                 ++z)
                for (int y = -lower; y <= upper; ++y)
                    for (int x = -lower; x <= upper; ++x) {
                        const auto i = location(p + Vec3i{x, y, z});
                        if (i < own_r_.size()) {
                            counts[0] += own_r_[i];
                            counts[1] += own_K_[i];
                        }
                    }
            return counts;
        };
        environment_->external_activation_ = [this](Vec3i p, CellStage stage) {
            const auto &b = config_.rules.continuum.base;
            const int block = b.migration_activation_block_edge;
            auto center = [block](int v) {
                return int(std::floor(double(v) / block)) * block + block / 2;
            };
            p = {center(p.x), center(p.y), b.thin_layer ? 0 : center(p.z)};
            const int edge = b.migration_activation_window_edge,
                      lower = (edge - 1) / 2, upper = edge - lower - 1;
            double sum = 0;
            if (config_.schema_version >= 4) {
                const auto values = box_counts(p, lower, upper);
                sum = values[0] + values[1];
                const double volume = stage == CellStage::large
                    ? pde_->large_cell_volume_ : 1.0;
                return sum / (std::pow(double(edge), b.thin_layer ? 2 : 3) / volume);
            }
            for (int z = b.thin_layer ? 0 : -lower;
                 z <= (b.thin_layer ? 0 : upper); ++z)
                for (int y = -lower; y <= upper; ++y)
                    for (int x = -lower; x <= upper; ++x) {
                        const auto i = location(p + Vec3i{x, y, z});
                        if (i < own_r_.size())
                            sum += own_r_[i] + own_K_[i];
                    }
            const double volume = stage == CellStage::large
                ? pde_->large_cell_volume_ : 1.0;
            const double window = std::pow(double(edge), b.thin_layer ? 2 : 3);
            if (config_.schema_version >= 3)
                return sum / (window / volume);
            return sum / (window * volume);
        };
    }
    abm_ = std::make_unique<Simulation3D>(base, std::move(env));
    core_.resize(pde_->voxel_count_);
    own_r_.resize(core_.size());
    own_K_.resize(core_.size());
    own_occupied_.resize(core_.size());
}
std::size_t HybridModel3D::location(Vec3i point) const {
    int x, y, z;
    if (!pde_->grid_coordinate(point, x, y, z))
        return core_.size();
    return pde_->index(x, y, z);
}
Vec3i HybridModel3D::site(std::size_t i) const {
    auto c = pde_->coordinate(i);
    return {int(std::floor(c[0])), int(std::floor(c[1])),
            config_.rules.continuum.base.thin_layer ? 0
                                                    : int(std::floor(c[2]))};
}
std::vector<Vec3i> HybridModel3D::footprint(Vec3i anchor,
                                            CellStage stage) const {
    if (stage != CellStage::large)
        return {anchor};
    std::vector<Vec3i> out;
    for (auto p : large_footprint(anchor))
        if (!config_.rules.continuum.base.thin_layer || p.z == 0)
            out.push_back(p);
    return out;
}
bool HybridModel3D::available(Vec3i point, bool replaces_agent) const {
    auto resource_point = point;
    if (config_.schema_version >= 3 &&
        config_.rules.continuum.base.thin_layer && point.z == 1)
        resource_point.z = 0;
    const auto i = location(resource_point);
    if (config_.schema_version >= 3) {
        if (i >= core_.size() ||
            (config_.schema_version < 4 && core_[i]) || pde_->vessel_blocks_cells(i))
            return false;
        double pde_volume = 0.0;
        for (std::size_t stage = 0; stage < 2; ++stage) {
            const double volume = stage == 0 ? 1.0 : pde_->large_cell_volume_;
            pde_volume += volume * (pde_->r_normal_[stage][i] +
                pde_->active_total_[stage][i] + pde_->K_[stage][i]);
        }
        const double agent_volume = abm_->grid_.occupants(point).size();
        // The grid's native collision checks reject ordinary occupied targets,
        // retain a moving cell's overlapping footprint, and debit the partner
        // of an accepted swap. The external obstacle must allow replacement
        // of that existing unit, rather than counting it twice.
        const double incoming_volume = replaces_agent
            ? std::max(0.0, 1.0 - agent_volume) : 1.0;
        return pde_volume + agent_volume + incoming_volume <=
            config_.rules.continuum.reaction.maximum_occupied_fraction;
    }
    return i < core_.size() && !core_[i] &&
           own_occupied_[i] <= config_.rules.migration.minimum_density &&
           !pde_->vessel_blocks_cells(i);
}
void HybridModel3D::initialize() {
    if (initialized_)
        return;
    abm_->initialize();
    if (config_.mode == "all_pde")
        pde_->initialize_from_abm(*abm_);
    else if (config_.mode == "adaptive") {
        structured_pde::StructuredInitialFields3D f;
        for (int s = 0; s < 2; ++s) {
            f.r_normal[s].resize(core_.size());
            f.r_active[s].resize(core_.size());
            f.K[s].resize(core_.size());
            f.active_remaining_hours[s].resize(core_.size());
        }
        f.vessel_fraction = environment_->vessels_;
        pde_->initialize_from_arrays(std::move(f));
        abm_->grid_.attach_external_blocker(environment_);
        assemble();
        synchronize_environment(true);
        exchange();
        attach_coupling();
    }
    initialized_ = true;
    if (config_.schema_version >= 4)
        update_representation_coverage();
}
void HybridModel3D::attach_coupling() {
    abm_->grid_.attach_external_blocker(environment_);
    if (config_.schema_version >= 4 && config_.mode == "adaptive") {
        pde_->external_transport_destination_ = [this](std::size_t location) {
            return location < core_.size() && core_[location] != 0;
        };
    }
    pde_->external_after_vascular_advance_ = [this] {
        for (auto slot : abm_->cells_.alive_slots())
            for (auto point : footprint(abm_->cells_.anchor(slot),
                                        abm_->cells_.stage(slot))) {
                const auto i = location(point);
                if (i >= core_.size() || pde_->vessel_blocks_cells(i)) {
                    if (config_.schema_version >= 3) {
                        pde_->record_vascular_removal(
                            abm_->cells_.type(slot), abm_->cells_.stage(slot),
                            (abm_->cells_.flags(slot) & kMigrationActive) != 0,
                            1.0);
                    }
                    remove_agent(slot);
                    ++abm_->stats_.vascular_displacements;
                    break;
                }
            }
        assemble();
    };
}
void HybridModel3D::initialize_from_arrays(
    structured_pde::StructuredInitialFields3D fields) {
    if (initialized_ || config_.mode != "adaptive")
        throw std::logic_error(
            "hybrid array initialization requires fresh adaptive model");
    if (config_.schema_version >= 4)
        for (int stage = 0; stage < 2; ++stage)
            for (const auto* values : {&fields.r_normal[stage], &fields.r_active[stage]})
                if (std::any_of(values->begin(), values->end(),
                                [](double mass) { return mass != 0.0; }))
                    throw std::invalid_argument("invasion hybrid requires individual r initialization");
    abm_->restore({}, 1, {}, {}, {});
    pde_->initialize_from_arrays(std::move(fields));
    attach_coupling();
    initialized_ = true;
    assemble();
    canonicalize();
    synchronize_environment(false);
}
void HybridModel3D::assemble() {
    auto &p = *pde_;
    p.external_r_.assign(core_.size(), 0);
    p.external_K_.assign(core_.size(), 0);
    p.external_occupied_.assign(core_.size(), 0);
    if (config_.schema_version >= 3) {
        p.external_r_consumers_.assign(core_.size(), 0.0);
        p.external_K_consumers_.assign(core_.size(), 0.0);
    }
    for (auto slot : abm_->cells_.alive_slots()) {
        auto anchor = abm_->cells_.anchor(slot);
        const auto i = location(anchor);
        if (i >= core_.size())
            throw std::runtime_error("hybrid agent outside grid");
        (abm_->cells_.type(slot) == CellType::r ? p.external_r_[i]
                                                : p.external_K_[i]) += 1;
        const auto points = footprint(anchor, abm_->cells_.stage(slot));
        for (auto point : points) {
            auto j = location(point);
            if (j >= core_.size())
                throw std::runtime_error("hybrid footprint outside grid");
            p.external_occupied_[j] += 1;
            if (config_.schema_version >= 3) {
                auto& consumers = abm_->cells_.type(slot) == CellType::r
                    ? p.external_r_consumers_ : p.external_K_consumers_;
                consumers[j] += 1.0 / points.size();
                int x, y, z;
                p.grid_coordinate(point, x, y, z);
                p.include_population_location(x, y, z);
            }
        }
        int x, y, z;
        p.grid_coordinate(anchor, x, y, z);
        p.include_population_location(x, y, z);
    }
}
void HybridModel3D::synchronize_environment(bool refresh_rates) {
    auto &p = *pde_;
    for (std::size_t i = 0; i < core_.size(); ++i) {
        own_r_[i] = p.r_normal_[0][i] + p.r_normal_[1][i] +
                    p.r_active(structured_pde::StructuredStage3D::small, i) +
                    p.r_active(structured_pde::StructuredStage3D::large, i);
        own_K_[i] = p.K_[0][i] + p.K_[1][i];
        own_occupied_[i] =
            p.occupied_fraction(i) -
            (p.external_occupied_.empty() ? 0 : p.external_occupied_[i]);
    }
    if (config_.schema_version >= 4)
        rebuild_density_prefix();
    environment_->nutrient_ = p.nutrient_;
    environment_->vessels_ = p.vessel_;
    environment_->tumour_mask_ = p.tumour_mask_;
    environment_->rebuild_prefix();
    environment_->prepare_guidance(abm_->cells_, pde_->step_count_ + 1);
    environment_->last_refresh_ = abm_->clock_.time_hours;
    environment_->next_refresh_ = 1.0e100;
    environment_->refreshes_ = 1;
    if (environment_->angiogenesis_)
        environment_->angiogenesis_->initialize(p.vessel_);
    if (!refresh_rates)
        return;
    abm_->migration_proposal_cache_.clear();
    abm_->proposal_cache_horizon_ = -1;
    auto slots = abm_->cells_.alive_slots();
    std::sort(slots.begin(), slots.end(), [this](auto a, auto b) {
        return abm_->cells_.uid(a) < abm_->cells_.uid(b);
    });
    for (auto slot : slots) {
        auto f = refresh_growth_state(slot, abm_->clock_.time_hours,
                                      abm_->cells_, abm_->density_,
                                      abm_->config_, environment_, false);
        f.migration_activation_changed = abm_->apply_migration_activation_class(
            slot, abm_->migration_activation_class(
                      abm_->migration_activation_query_block(
                          abm_->cells_.anchor(slot))));
        abm_->apply_growth_refresh(slot, f);
    }
}
void HybridModel3D::remove_agent(Slot slot) {
    abm_->events_.cancel_cell(slot);
    remove_cell(slot, abm_->cells_, abm_->grid_, abm_->density_);
}
void HybridModel3D::classify_core() {
    synchronize_environment(false);
    if (config_.schema_version >= 4) {
        classify_invasion_core();
        return;
    }
    const int radius = config_.smoothing_radius;
    const bool thin = config_.rules.continuum.base.thin_layer;
    for (std::size_t i = 0; i < core_.size(); ++i) {
        double sum = 0;
        int count = 0;
        const auto origin = site(i);
        for (int z = thin ? 0 : -radius; z <= (thin ? 0 : radius); ++z)
            for (int y = -radius; y <= radius; ++y)
                for (int x = -radius; x <= radius; ++x) {
                    const auto j = location(origin + Vec3i{x, y, z});
                    if (j >= core_.size())
                        continue;
                    sum += own_occupied_[j] + pde_->external_occupied_[j];
                    ++count;
                }
        const double density = count ? sum / count : 0;
        core_[i] = density >= (core_[i] ? config_.core_off : config_.core_on);
    }
}
void HybridModel3D::canonicalize() {
    auto &p = *pde_;
    for (int s = 0; s < 2; ++s)
        for (std::size_t i = 0; i < core_.size(); ++i) {
            double total = 0;
            for (auto &bucket : p.active_direction_[s])
                total += bucket[i];
            p.active_total_[s][i] = float(total);
        }
    p.shrink_population_bounds();
    for (int s = 0; s < 2; ++s)
        p.shrink_active_bounds(s);
    p.rebuild_moving_tumour_front();
    p.build_activation_density();
}
void HybridModel3D::convert_to_density() {
    auto slots = abm_->cells_.alive_slots();
    std::sort(slots.begin(), slots.end(), [this](auto a, auto b) {
        return abm_->cells_.uid(a) < abm_->cells_.uid(b);
    });
    for (auto slot : slots) {
        const auto c = abm_->cells_.snapshot(slot);
        auto i = location(c.anchor);
        if (!core_[i])
            continue;
        if (config_.schema_version >= 4 && c.type == CellType::r)
            continue;
        // Ultrasmall colocations cannot be represented by this schema's two
        // stages.
        if (c.stage == CellStage::ultrasmall)
            continue;
        const auto footprint_sites = footprint(c.anchor, c.stage);
        bool valid = true;
        for (auto point : footprint_sites)
            if (location(point) >= core_.size() ||
                pde_->vessel_blocks_cells(location(point)) ||
                (config_.schema_version >= 4 && !core_[location(point)]))
                valid = false;
        if (!valid)
            continue;
        const int s = c.stage == CellStage::large ? 1 : 0;
        const double mass = 1.0 / footprint_sites.size();
        const bool distributions = config_.schema_version >= 2;
        const double work = std::max(0.0, double(c.division_work_remaining) -
            std::max(0.0, double(static_cast<float>(time_hours())) - c.last_update_time) *
            std::max(0.0, double(c.density_growth_rate)));
        // Large ABM footprints become a distributed density of one cell.
        for (auto point : footprint_sites) {
            const auto j = location(point);
            int x, y, z;
            pde_->grid_coordinate(point, x, y, z);
            pde_->include_population_location(x, y, z);
            if (distributions) {
                pde_->renewal_->add(j, (c.type == CellType::r ? 0 : 2) + s, mass, work);
                if (c.flags & kMigrationActive) {
                    if (pde_->duration_) pde_->duration_->add(j, s, mass,
                        std::max(0.0, c.migration_activation_end_time - time_hours()));
                    if (pde_->velocity_) pde_->velocity_->add(j, s, mass, c.normal_migration_rate);
                }
            }
            if (c.type == CellType::K)
                pde_->K_[s][j] += mass;
            else if (c.flags & kMigrationActive) {
                std::size_t bucket = 0;
                for (std::size_t d = 0; d < pde_->direction_ids_.size(); ++d)
                    if (pde_->direction_ids_[d] == c.last_direction)
                        bucket = d + 1;
                pde_->active_direction_[s][bucket][j] += float(mass);
                pde_->active_clock_[s][bucket][j] +=
                    float(mass * std::max(0.0, c.migration_activation_end_time -
                                                   time_hours()));
                pde_->include_active_location(s, x, y, z);
            } else {
                pde_->r_normal_[s][j] += mass;
                auto ref = environment_->refractory_.find(c.uid);
                if (ref != environment_->refractory_.end() &&
                    ref->second.until > time_hours()) {
                    pde_->r_refractory_[s][j] += mass;
                    pde_->refractory_clock_[s][j] +=
                        mass * (ref->second.until - time_hours());
                }
            }
        }
        environment_->refractory_.erase(c.uid);
        remove_agent(slot);
        ++to_pde_;
    }
}
void HybridModel3D::convert_to_agents() {
    // Fixed 16-voxel exchange blocks bound redistribution distance and keep
    // disconnected fronts from exchanging individual-cell mass globally.
    std::map<std::array<int, 3>, std::vector<std::size_t>> regions;
    for (std::size_t i = 0; i < core_.size(); ++i) {
        if (core_[i] || pde_->vessel_blocks_cells(i) || own_occupied_[i] <= 0)
            continue;
        const auto point = site(i);
        regions[{int(std::floor(double(point.x) / 16)),
                 int(std::floor(double(point.y) / 16)),
                 int(std::floor(double(point.z) / 16))}]
            .push_back(i);
    }
    for (const auto &[key, region] : regions)
        convert_region(region);
}
void HybridModel3D::convert_region(const std::vector<std::size_t> &region) {
    if (config_.schema_version >= 4) {
        convert_invasion_region(region);
        return;
    }
    // Dependent rounding debits exact units from regional phenotype/stage
    // budgets. Fractional remainders stay in PDE compartments. Clearing all
    // eligible cohorts together permits mixed r/K voxels to become agents.
    struct Pool {
        double mass{}, clock{};
        std::vector<double> work, duration, velocity, directions;
    };
    std::array<std::array<Pool, 4>, 2> pools{};
    std::vector<std::size_t> eligible = region;
    auto &p = *pde_;
    const bool distributions = config_.schema_version >= 2;
    if (distributions) {
        std::array<double, 8> budgets{};
        for (auto i : region) for (int s = 0; s < 2; ++s) {
            budgets[4 * s] += p.r_normal_[s][i] - p.r_refractory_[s][i];
            budgets[4 * s + 2] += p.K_[s][i];
            budgets[4 * s + 3] += p.r_refractory_[s][i];
            for (const auto& field : p.active_direction_[s]) budgets[4 * s + 1] += field[i];
        }
        if (*std::max_element(budgets.begin(), budgets.end()) < 1.0) return;
    }
    std::unordered_map<std::size_t, std::array<double, 8>> original_mass;
    for (auto i : eligible) {
        for (int stage = 0; stage < 2; ++stage) {
            pools[stage][0].mass +=
                p.r_normal_[stage][i] - p.r_refractory_[stage][i];
            pools[stage][3].mass += p.r_refractory_[stage][i];
            pools[stage][3].clock += p.refractory_clock_[stage][i];
            pools[stage][2].mass += p.K_[stage][i];
            double active_mass = 0.0;
            for (std::size_t bucket = 0;
                 bucket < p.active_direction_[stage].size(); ++bucket) {
                pools[stage][1].mass += p.active_direction_[stage][bucket][i];
                active_mass += p.active_direction_[stage][bucket][i];
                pools[stage][1].clock += p.active_clock_[stage][bucket][i];
                if (distributions) {
                    auto& directions = pools[stage][1].directions;
                    if (directions.empty()) directions.resize(p.active_direction_[stage].size());
                    directions[bucket] += p.active_direction_[stage][bucket][i];
                }
                p.active_direction_[stage][bucket][i] = 0;
                p.active_clock_[stage][bucket][i] = 0;
            }
            if (distributions) {
                const std::array<double, 4> mass{
                    p.r_normal_[stage][i] - p.r_refractory_[stage][i], active_mass,
                    p.K_[stage][i], p.r_refractory_[stage][i]};
                for (int kind = 0; kind < 4; ++kind) original_mass[i][4 * stage + kind] = mass[kind];
                for (int kind = 0; kind < 4; ++kind) {
                    const auto channel = (kind == 2 ? 2 : 0) + stage;
                    const double total = p.renewal_->mass(i, channel);
                    if (total > 0.0) add_mixture(pools[stage][kind].work,
                        p.renewal_->distribution(i, channel), mass[kind] / total);
                }
                for (auto [bank, target] : {std::pair{p.duration_.get(), &pools[stage][1].duration},
                                            std::pair{p.velocity_.get(), &pools[stage][1].velocity}})
                    if (bank && bank->mass(i, stage) > 0.0)
                        add_mixture(*target, bank->distribution(i, stage), active_mass / bank->mass(i, stage));
            }
            p.r_normal_[stage][i] = 0;
            p.r_refractory_[stage][i] = 0;
            p.refractory_clock_[stage][i] = 0;
            p.K_[stage][i] = 0;
            p.active_total_[stage][i] = 0;
        }
        if (distributions)
            for (auto* bank : {p.renewal_.get(), p.duration_.get(), p.velocity_.get()})
                if (bank) bank->erase(i);
    }
    synchronize_environment(false);
    const auto remaining_volume = [&] {
        double v = 0;
        for (int s = 0; s < 2; ++s)
            for (auto pool : pools[s])
                v += pool.mass * (s ? p.large_cell_volume_ : 1);
        return v;
    };
    for (int stage = 1; stage >= 0; --stage)
        for (int kind = 0; kind < 4; ++kind) {
            auto &pool = pools[stage][kind];
            if (pool.mass < 1)
                continue;
            const double mean_clock = pool.clock / pool.mass;
            std::sort(eligible.begin(), eligible.end(), [&](auto a, auto b) {
                return mix(a ^ config_.rules.continuum.base.seed ^
                           mix(exchange_count_ +
                               std::uint64_t(8 * stage + kind))) <
                       mix(b ^ config_.rules.continuum.base.seed ^
                           mix(exchange_count_ +
                               std::uint64_t(8 * stage + kind)));
            });
            const auto cell_stage = stage ? CellStage::large : CellStage::small;
            for (auto i : eligible) {
                if (pool.mass < 1)
                    break;
                auto anchor = site(i);
                auto points = footprint(anchor, cell_stage);
                bool free = true;
                for (auto point : points)
                    if (!available(point) ||
                        !abm_->grid_.occupants(point).empty())
                        free = false;
                if (!free)
                    continue;
                double capacity = 0;
                for (auto j : eligible)
                    if (abm_->grid_.occupants(site(j)).empty())
                        capacity += 1;
                std::size_t consumed = 0;
                for (auto point : points)
                    if (std::find(eligible.begin(), eligible.end(),
                                  location(point)) != eligible.end())
                        ++consumed;
                if (capacity - double(consumed) + 1e-10 <
                    remaining_volume() - (stage ? p.large_cell_volume_ : 1))
                    continue;
                CellInit cell;
                cell.anchor = anchor;
                cell.uid = abm_->next_uid_++;
                cell.clone_id = std::uint32_t(cell.uid);
                cell.type = kind == 2 ? CellType::K : CellType::r;
                cell.stage = cell_stage;
                cell.last_update_time = time_hours();
                auto rates = sample_initial_cell_rates(
                    abm_->config_, cell.type, abm_->config_.seed, cell.uid);
                cell.inherent_growth_rate = float(rates.inherent_growth_rate);
                cell.migration_rate = float(rates.migration_rate);
                cell.normal_migration_rate = sample_normal_migration_rate(
                    cell.type, cell.migration_rate, abm_->config_, cell.uid, 0);
                const auto sample = [&](const auto& values, std::uint64_t domain) {
                    return sample_mixture(values, config_.rules.continuum.base.seed ^
                        mix(cell.uid) ^ mix(exchange_count_) ^ domain);
                };
                if (distributions) {
                    cell.division_work_remaining = float(p.renewal_->bin_width() * sample(pool.work, 0x485942574f524bULL));
                    if (kind == 1 && p.velocity_)
                        cell.normal_migration_rate = float(p.velocity_->bin_width() * sample(pool.velocity, 0x48594252415445ULL));
                }
                if (cell.type == CellType::r &&
                    abm_->config_.activated_r_migration_rate_model ==
                        "normal_multiplier")
                    cell.migration_rate =
                        float(cell.normal_migration_rate *
                              abm_->config_.activated_r_normal_multiplier);
                if (kind == 1) {
                    cell.flags |= kMigrationActive;
                    cell.migration_activation_end_time =
                        time_hours() + (distributions && p.duration_ ?
                            p.duration_->bin_width() * sample(pool.duration, 0x48594254494d45ULL) : mean_clock);
                    if (distributions) {
                        const auto bucket = sample(pool.directions, 0x485942444952ULL);
                        cell.last_direction = bucket ? p.direction_ids_[bucket - 1] : kStayDirection;
                        if (!(cell.migration_activation_end_time > time_hours())) {
                            cell.flags &= std::uint8_t(~kMigrationActive);
                            cell.migration_activation_end_time = 0.0;
                            cell.last_direction = kStayDirection;
                        }
                    }
                }
                const double effective_rate = (cell.flags & kMigrationActive)
                                                  ? cell.migration_rate
                                                  : cell.normal_migration_rate;
                cell.next_migration_time =
                    effective_rate > 0 ? time_hours() + 1 / effective_rate : 0;
                const auto slot = abm_->cells_.create(cell);
                const bool placed =
                    stage ? abm_->grid_.place_large(anchor, slot)
                          : abm_->grid_.place_single(anchor, slot);
                if (!placed)
                    throw std::logic_error("hybrid rounded placement failed");
                abm_->density_.add(anchor, cell.type, slot);
                if (!distributions)
                    initialize_division_cycle(slot, time_hours(), abm_->cells_, abm_->config_);
                if (kind == 3)
                    environment_->refractory_[cell.uid] = {
                        time_hours() + mean_clock, false};
                if (distributions && kind == 1 && !(cell.flags & kMigrationActive))
                    environment_->refractory_[cell.uid] = {
                        time_hours() + config_.rules.migration.reactivation_cooldown_hours, false};
                abm_->schedule_cell(slot);
                if (distributions)
                    for (auto* values : {&pool.work, &pool.duration, &pool.velocity, &pool.directions})
                        debit_mixture(*values, 1 / pool.mass);
                pool.mass -= 1;
                pool.clock = pool.mass * mean_clock;
                ++to_abm_;
            }
        }
    // Residuals share vacant voxels up to the same biological volume limit.
    std::vector<double> used(core_.size(), 0);
    for (int stage = 0; stage < 2; ++stage)
        for (int kind = 0; kind < 4; ++kind) {
            auto &pool = pools[stage][kind];
            const double mean_clock =
                pool.mass > 0 ? pool.clock / pool.mass : 0;
            // Keep fractional fronts at their original sites. Repeatedly
            // moving sub-cell remainders to the first permuted site creates
            // an artificial block-scale diffusion, even without conversion.
            const auto& restore_order = distributions ? region : eligible;
            bool progressed;
            do {
                const double pass_mass = pool.mass;
                double weight_sum = 0.0;
                if (distributions) for (auto i : restore_order)
                    if (abm_->grid_.occupants(site(i)).empty() && used[i] < 1.0)
                        weight_sum += original_mass.at(i)[4 * stage + kind];
                const bool uniform = !(weight_sum > 0.0);
                if (distributions && uniform) for (auto i : restore_order)
                    if (abm_->grid_.occupants(site(i)).empty() && used[i] < 1.0) weight_sum += 1.0;
                progressed = false;
                for (auto i : restore_order) {
                    if (pool.mass <= 0)
                        break;
                    if (!abm_->grid_.occupants(site(i)).empty())
                        continue;
                    const double volume = stage ? p.large_cell_volume_ : 1;
                    const double desired = distributions ? (weight_sum > 0.0 ? pass_mass *
                        (uniform ? 1.0 : original_mass.at(i)[4 * stage + kind]) / weight_sum : 0.0) : pool.mass;
                    const double mass =
                        std::min({pool.mass, desired, std::max(0.0, 1 - used[i]) / volume});
                    if (distributions && !(mass > 0.0)) continue;
                    progressed = progressed || mass > 0.0;
                    double active_added = 0.0;
                    if (kind == 0 || kind == 3) {
                        p.r_normal_[stage][i] += mass;
                        if (kind == 3) {
                            p.r_refractory_[stage][i] += mass;
                            p.refractory_clock_[stage][i] += mass * mean_clock;
                        }
                    } else if (kind == 2)
                        p.K_[stage][i] += mass;
                    else {
                        if (distributions) {
                            const double total = std::accumulate(pool.directions.begin(), pool.directions.end(), 0.0);
                            for (std::size_t bucket = 0; bucket < pool.directions.size(); ++bucket) {
                                const float m = std::nextafter(float(mass * pool.directions[bucket] / total), 0.0F);
                                p.active_direction_[stage][bucket][i] += m;
                                p.active_clock_[stage][bucket][i] += float(m * mean_clock);
                                active_added += m;
                            }
                            p.active_total_[stage][i] += float(active_added);
                            p.r_normal_[stage][i] += mass - active_added;
                        } else {
                            const float m = std::nextafter(float(mass), 0.0F);
                            p.active_direction_[stage][0][i] += m;
                            p.active_clock_[stage][0][i] += float(m * mean_clock);
                            p.active_total_[stage][i] += m;
                            p.r_normal_[stage][i] += mass - double(m);
                        }
                    }
                    if (distributions) {
                        p.renewal_->add_distribution(i, (kind == 2 ? 2 : 0) + stage, pool.work, mass / pool.mass);
                        if (kind == 1) {
                            if (p.duration_) p.duration_->add_distribution(i, stage, pool.duration, active_added / pool.mass);
                            if (p.velocity_) p.velocity_->add_distribution(i, stage, pool.velocity, active_added / pool.mass);
                        }
                        for (auto* values : {&pool.work, &pool.duration, &pool.velocity, &pool.directions})
                            debit_mixture(*values, mass / pool.mass);
                    }
                    pool.mass -= mass;
                    used[i] += mass * volume;
                    int x, y, z;
                    p.grid_coordinate(site(i), x, y, z);
                    p.include_population_location(x, y, z);
                    if (kind == 1)
                        p.include_active_location(stage, x, y, z);
                }
            } while (distributions && progressed && pool.mass > 1e-12);
            if (pool.mass > 1e-9)
                throw std::logic_error(
                    "hybrid exchange exhausted residual capacity");
        }
}

void HybridModel3D::exchange() {
    if (config_.mode != "adaptive")
        return;
    assemble();
    classify_core();
    convert_to_density();
    assemble();
    canonicalize();
    convert_to_agents();
    assemble();
    canonicalize();
    synchronize_environment(true);
    ++exchange_count_;
    if (config_.schema_version >= 4 && initialized_)
        update_representation_coverage();
}
double HybridModel3D::time_hours() const noexcept {
    return config_.mode == "all_abm" ? abm_->clock().time_hours
                                     : pde_->time_hours();
}
bool HybridModel3D::step() {
    initialize();
    const auto &c = config_.rules.continuum;
    if (time_hours() >= c.end_time_hours)
        return false;
    if (config_.mode == "all_pde")
        return pde_->step();
    if (config_.mode == "all_abm") {
        abm_->run();
        return true;
    }
    const double target =
        std::min(c.end_time_hours, time_hours() + c.time_step_hours);
    if (config_.mode == "adaptive") {
        assemble();
        pde_->step();
        synchronize_environment(true);
    }
    const double end = abm_->config_.end_time_hours;
    abm_->config_.end_time_hours = target;
    abm_->run();
    abm_->config_.end_time_hours = end;
    if (config_.mode == "adaptive") {
        assemble();
        canonicalize();
        const auto ticks = std::uint64_t(
            std::llround(config_.exchange_every_hours / c.time_step_hours));
        if (pde_->step_count_ % ticks == 0)
            exchange();
    }
    if (config_.schema_version >= 4)
        update_representation_coverage();
    return true;
}
HybridDiagnostics3D HybridModel3D::diagnostics() const {
    HybridDiagnostics3D d;
    d.to_pde = to_pde_;
    d.to_abm = to_abm_;
    d.exchanges = exchange_count_;
    if (config_.mode != "all_pde") {
        d.abm_mass = abm_->cells().alive_count();
        for (auto slot : abm_->cells().alive_slots()) {
            (abm_->cells().type(slot) == CellType::r ? d.r_mass : d.K_mass) += 1;
            if (abm_->cells().flags(slot) & kMigrationActive)
                d.active_mass += 1;
        }
    }
    if (config_.mode != "all_abm") {
        auto p = pde_->diagnostics();
        d.pde_mass = p.r_total + p.K_total;
        d.r_mass += p.r_total;
        d.K_mass += p.K_total;
        d.active_mass += p.r_active_total;
        d.mean_nutrient = p.mean_nutrient;
    } else
        d.mean_nutrient = std::accumulate(environment_->nutrient().begin(),
                                          environment_->nutrient().end(), 0.0) /
                          environment_->nutrient().size();
    d.total_mass = d.abm_mass + d.pde_mass;
    if (config_.schema_version >= 4 && config_.mode == "adaptive")
        for (std::size_t i = 0; i < core_.size(); ++i)
            if (!core_[i])
                d.front_pde_mass += (pde_->K_[0][i] + pde_->K_[1][i]) * pde_->voxel_measure_;
    d.minimum_abm_fraction = minimum_abm_fraction_;
    d.minimum_active_fraction = minimum_active_fraction_;
    d.maximum_pde_fraction = maximum_pde_fraction_;
    return d;
}
std::vector<double> HybridModel3D::radial_mass() const {
    double squared = 0.0;
    for (int edge : config_.rules.continuum.grid.shape) squared += double(edge) * edge;
    std::vector<double> mass(std::size_t(std::ceil(std::sqrt(squared))) + 1);
    const auto add = [&](double x, double y, double z, double amount) {
        const double radius = std::sqrt(x * x + y * y +
            (config_.rules.continuum.base.thin_layer ? 0.0 : z * z));
        mass[std::min(mass.size() - 1, std::size_t(std::floor(radius)))] += amount;
    };
    if (config_.mode != "all_pde") for (auto slot : abm_->cells().alive_slots()) {
        const auto anchor = abm_->cells().anchor(slot);
        add(anchor.x + 0.5, anchor.y + 0.5, anchor.z + 0.5, 1);
    }
    if (config_.mode != "all_abm") for (std::size_t i = 0; i < core_.size(); ++i) {
        const auto center = pde_->coordinate(i);
        double amount = 0;
        for (int s = 0; s < 2; ++s) amount += pde_->r_normal_[s][i] + pde_->active_total_[s][i] + pde_->K_[s][i];
        add(center[0], center[1], center[2], amount * pde_->voxel_measure_);
    }
    return mass;
}
std::uint64_t HybridModel3D::state_checksum() const {
    if (config_.mode == "all_abm")
        return abm_->state_checksum();
    if (config_.mode == "all_pde")
        return pde_->state_checksum();
    auto h = mix(config_.fingerprint() ^ abm_->state_checksum());
    h = mix(h ^ pde_->state_checksum());
    h = mix(h ^ exchange_count_);
    h = mix(h ^ to_pde_);
    h = mix(h ^ to_abm_);
    for (auto x : core_)
        h = mix(h ^ x);
    if (config_.schema_version >= 4)
        for (double value : {minimum_abm_fraction_, minimum_active_fraction_, maximum_pde_fraction_})
            h = mix(h ^ std::bit_cast<std::uint64_t>(value));
    for (const auto &[uid, r] : environment_->refractory_) {
        h = mix(h ^ uid);
        h = mix(h ^ std::bit_cast<std::uint64_t>(r.until));
        h = mix(h ^ r.armed);
    }
    return h;
}
std::uint64_t HybridModel3D::checkpoint_magic() const noexcept {
    if (config_.schema_version >= 4)
        return 0x4154434748594234;
    if (config_.schema_version >= 3)
        return 0x4154434748594233;
    return config_.schema_version == 2 ? 0x4154434748594232 : 0x4154434748594231;
}
void HybridModel3D::save_checkpoint(const std::filesystem::path &path) const {
    if (!initialized_)
        throw std::logic_error("hybrid checkpoint before initialization");
    for (auto suffix : {"", ".pde.bin", ".resource.bin"})
        if (std::filesystem::exists(path.string() + suffix))
            throw std::runtime_error("refusing to overwrite hybrid checkpoint");
    if (!path.parent_path().empty())
        std::filesystem::create_directories(path.parent_path());
    if (config_.mode != "all_abm")
        pde_->save_checkpoint(path.string() + ".pde.bin");
    const auto temporary = path.string() + ".tmp";
    std::ofstream out(temporary, std::ios::binary);
    put(out, std::uint64_t(checkpoint_magic()));
    put(out, config_.fingerprint());
    put(out, state_checksum());
    put(out, abm_->state_checksum());
    put(out, exchange_count_);
    put(out, to_pde_);
    put(out, to_abm_);
    put_vector(out, core_);
    if (config_.schema_version >= 4) {
        put(out, minimum_abm_fraction_);
        put(out, minimum_active_fraction_);
        put(out, maximum_pde_fraction_);
    }
    put(out, abm_->clock_);
    put(out, abm_->stats_);
    put(out, abm_->next_uid_);
    put(out, std::uint64_t(abm_->cells_.slot_count()));
    auto slots = abm_->cells_.alive_slots();
    std::vector<CellInit> cells;
    for (auto slot : slots)
        cells.push_back(abm_->cells_.snapshot(slot));
    put_vector(out, slots);
    put_vector(out, cells);
    put_vector(out, abm_->cells_.free_slots());
    put_vector(out, abm_->lineage_);
    const auto vascular = abm_->snapshot_vasculature();
    put(out, vascular.process);
    put(out, vascular.next_vessel_id);
    put(out, vascular.next_node_uid);
    put(out, vascular.next_tip_uid);
    put(out, vascular.lesions.next_lesion_id);
    put(out, vascular.lesions.last_refresh_time_hours);
    put(out, vascular.lesions.next_refresh_time_hours);
    put(out, vascular.lesions.refresh_schedule_generation);
    put_vector(out, vascular.lesions.core_identity);
    put_vector(out, vascular.lesions.dirty_blocks);
    put_vector(out, vascular.lesions.processes);
    put_vector(out, vascular.lesions.source_ownership);
    put_vector(out, vascular.nodes);
    put_vector(out, vascular.tips);
    put_vector(out, vascular.perfused_vessels);
    if (!out)
        throw std::runtime_error("hybrid checkpoint write failed");
    out.close();
    environment_->save_checkpoint(path.string() + ".resource.bin",
                                  abm_->state_checksum());
    std::filesystem::rename(temporary, path);
}
void HybridModel3D::load_checkpoint(const std::filesystem::path &path) {
    std::ifstream in(path, std::ios::binary);
    std::uint64_t magic, fingerprint, checksum, abm_checksum;
    get(in, magic);
    get(in, fingerprint);
    get(in, checksum);
    get(in, abm_checksum);
    if (magic != checkpoint_magic() ||
        fingerprint != config_.fingerprint())
        throw std::runtime_error("hybrid checkpoint configuration mismatch");
    get(in, exchange_count_);
    get(in, to_pde_);
    get(in, to_abm_);
    get_vector(in, core_);
    if (config_.schema_version >= 4) {
        get(in, minimum_abm_fraction_);
        get(in, minimum_active_fraction_);
        get(in, maximum_pde_fraction_);
    }
    SimulationClock3D clock;
    SimulationStats3D stats;
    CellUid uid;
    std::uint64_t count;
    get(in, clock);
    get(in, stats);
    get(in, uid);
    get(in, count);
    std::vector<Slot> slots, free;
    std::vector<CellInit> cells;
    std::vector<LineageEdge> lineage;
    get_vector(in, slots);
    get_vector(in, cells);
    get_vector(in, free);
    get_vector(in, lineage);
    VasculatureState3D vascular;
    get(in, vascular.process);
    get(in, vascular.next_vessel_id);
    get(in, vascular.next_node_uid);
    get(in, vascular.next_tip_uid);
    get(in, vascular.lesions.next_lesion_id);
    get(in, vascular.lesions.last_refresh_time_hours);
    get(in, vascular.lesions.next_refresh_time_hours);
    get(in, vascular.lesions.refresh_schedule_generation);
    get_vector(in, vascular.lesions.core_identity);
    get_vector(in, vascular.lesions.dirty_blocks);
    get_vector(in, vascular.lesions.processes);
    get_vector(in, vascular.lesions.source_ownership);
    get_vector(in, vascular.nodes);
    get_vector(in, vascular.tips);
    get_vector(in, vascular.perfused_vessels);
    if (in.peek() != std::char_traits<char>::eof())
        throw std::runtime_error("hybrid checkpoint has trailing data");
    if (config_.mode != "all_abm")
        pde_->load_checkpoint(path.string() + ".pde.bin");
    environment_->load_checkpoint(path.string() + ".resource.bin",
                                  abm_checksum);
    abm_->restore(cells, uid, clock, stats, lineage, vascular, count, slots,
                  free);
    initialized_ = true;
    if (config_.mode == "adaptive") {
        abm_->grid_.attach_external_blocker(environment_);
        assemble();
        canonicalize();
        synchronize_environment(false);
        attach_coupling();
    }
    if (checksum != state_checksum())
        throw std::runtime_error("hybrid checkpoint checksum mismatch");
}
} // namespace atcg3d::hybrid
