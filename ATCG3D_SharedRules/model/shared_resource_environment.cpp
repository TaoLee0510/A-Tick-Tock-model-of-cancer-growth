#include "model/shared_resource_environment.hpp"

#include <bit>
#include <cmath>
#include <fstream>
#include <set>
#include <stdexcept>

#include "core/cell_store.hpp"
#include "engine/parallelism.hpp"
#include "geometry/footprint.hpp"
#include "model/shared_resource.hpp"
#include "vasculature/vessel_grid.hpp"

namespace atcg3d::shared_rules {
namespace {
std::uint64_t mix(std::uint64_t state, std::uint64_t value) noexcept {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    return state ^ ((value ^ (value >> 31U)) + (state << 6U) + (state >> 2U));
}
template <class T> void write(std::ostream& out, T value) {
    out.write(reinterpret_cast<const char*>(&value), sizeof(T));
    if (!out) throw std::runtime_error("unable to write shared resource checkpoint");
}
template <class T> T read(std::istream& in) {
    T value{};
    in.read(reinterpret_cast<char*>(&value), sizeof(T));
    if (!in) throw std::runtime_error("truncated shared resource checkpoint");
    return value;
}
}  // namespace

Model3DConfig abm_config(const structured_pde::StructuredPdeConfig3D& config) {
    if (config.schema_version < 5) throw std::invalid_argument("shared ABM requires structured schema v5+");
    auto base = config.continuum.base;
    base.direction_guidance_model = "nutrient_gradient_shared_resource_v3";
    base.static_vasculature = config.continuum.shared_vascular_geometry();
    const double mapping_scale = base.legacy_mapping.density_count_scale_2d_to_3d;
    base.legacy_mapping.source_r_limit = base.legacy_mapping.source_K_limit = config.continuum.nutrient.common_density_limit;
    base.legacy_mapping.source_carrying_capacity_r =
        base.legacy_mapping.source_carrying_capacity_K = config.continuum.nutrient.common_carrying_capacity;
    base.r_limit = base.K_limit = base.legacy_mapping.source_r_limit * mapping_scale;
    base.carrying_capacity_r = base.carrying_capacity_K = base.legacy_mapping.source_carrying_capacity_r * mapping_scale;
    base.migration_swap_enabled = config.migration.crowding_exchange != "none";
    base.end_time_hours = config.continuum.end_time_hours;
    base.output_enabled = false;
    base.control_enabled = false;
    base.run_mode = "new";
    base.resume_checkpoint.clear();
    if (base.static_vasculature.spacing_voxels != 1.0) {
        throw std::invalid_argument("shared ABM requires unit lattice spacing");
    }
    base.validate();
    return base;
}

SharedResourceEnvironment3D::SharedResourceEnvironment3D(structured_pde::StructuredPdeConfig3D config)
    : config_(std::move(config)), fingerprint_(config_.dynamics_fingerprint()),
      geometry_(config_.continuum.shared_vascular_geometry()) {
    config_.validate();
    (void)abm_config(config_);
    const auto& grid = config_.continuum.grid;
    const std::size_t size = static_cast<std::size_t>(grid.shape[0]) * grid.shape[1] * grid.shape[2];
    const int dimensions = geometry_.thin_layer ? 2 : 3;
    if (config_.continuum.nutrient.diffusion_voxels2_per_hour *
        config_.continuum.time_step_hours * 2.0 * dimensions > 1.0 + 1.0e-12) {
        throw std::invalid_argument("shared resource diffusion violates CFL bound");
    }
    nutrient_.assign(size, config_.continuum.nutrient.initial_value);
    for (auto* field : {&next_, &consumers_, &occupied_, &vessels_}) field->assign(size, 0.0);
    tumour_mask_.assign(size, 0U);
    if(config_.continuum.angiogenesis.model != "disabled") angiogenesis_=std::make_unique<continuum::AngiogenesisField3D>(config_.continuum.angiogenesis,grid.shape,grid.spacing_voxels,geometry_.thin_layer);
    const int edge = config_.migration.direction_nutrient_window_edge;
    const int lower = (edge - 1) / 2;
    const int upper = edge - lower - 1;
    const double cosine_limit = std::cos(config_.continuum.base.direction_density_half_angle_degrees * std::acos(-1.0) / 180.0);
    for (DirectionId direction = 1; direction <= 26; ++direction) {
        const Vec3i forward = direction_vector(direction);
        if (geometry_.thin_layer && forward.z != 0) continue;
        const double length = std::sqrt(static_cast<double>(squared_length(forward)));
        for (int z = geometry_.thin_layer ? 0 : -lower; z <= (geometry_.thin_layer ? 0 : upper); ++z) {
            for (int y = -lower; y <= upper; ++y) {
                int start = -lower;
                bool in_span = false;
                for (int x = -lower; x <= upper + 1; ++x) {
                    const Vec3i offset{x, y, z};
                    const bool included = x <= upper && squared_length(offset) > 0 &&
                        static_cast<double>(dot(offset, forward)) /
                            (length * std::sqrt(static_cast<double>(squared_length(offset)))) + 1.0e-12 >= cosine_limit;
                    if (included && !geometry_.thin_layer) sector_offsets_[direction].push_back(offset);
                    if (geometry_.thin_layer) {
                        if (included && !in_span) { start = x; in_span = true; }
                        if (!included && in_span) { sector_rows_[direction].push_back({y, start, x - 1}); in_span = false; }
                    }
                }
            }
        }
    }
    rebuild_prefix();
}

std::size_t SharedResourceEnvironment3D::index(int x, int y, int z) const noexcept {
    return (static_cast<std::size_t>(z) * geometry_.shape[1] + y) * geometry_.shape[0] + x;
}
std::size_t SharedResourceEnvironment3D::location(Vec3i site) const noexcept {
    return index(static_cast<int>(std::floor(site.x - geometry_.origin[0])),
        static_cast<int>(std::floor(site.y - geometry_.origin[1])),
        geometry_.thin_layer ? 0 : static_cast<int>(std::floor(site.z - geometry_.origin[2])));
}
bool SharedResourceEnvironment3D::contains_resource_site(Vec3i site) const noexcept {
    return (!geometry_.thin_layer || site.z == 0) && geometry_.contains(site);
}
double SharedResourceEnvironment3D::normalized_resource(Vec3i site) const noexcept {
    return contains_resource_site(site) ? nutrient_[location(site)] / config_.continuum.nutrient.vessel_value : 0.0;
}
double SharedResourceEnvironment3D::shared_density_limit() const noexcept { return config_.continuum.nutrient.common_density_limit; }
double SharedResourceEnvironment3D::shared_carrying_capacity() const noexcept { return config_.continuum.nutrient.common_carrying_capacity; }
double SharedResourceEnvironment3D::growth_resource_scale(Vec3i site) const noexcept {
    const double value = contains_resource_site(site) ? nutrient_[location(site)] : 0.0;
    return value / (config_.continuum.nutrient.growth_half_saturation + value);
}

void SharedResourceEnvironment3D::assemble(const CellStore3D& cells, const SparseVesselGrid3D& vessels) {
    std::fill(consumers_.begin(), consumers_.end(), 0.0);
    std::fill(occupied_.begin(), occupied_.end(), 0.0);
    std::fill(vessels_.begin(), vessels_.end(), 0.0);
    std::set<CellUid> alive;
    const double large_volume = geometry_.thin_layer ? 4.0 : 8.0;
    for (const auto slot : cells.alive_slots()) {
        alive.insert(cells.uid(slot));
        const bool large = cells.stage(slot) == CellStage::large;
        const auto add = [&](Vec3i site) {
            if (!contains_resource_site(site)) return;
            consumers_[location(site)] += large ? 1.0 / large_volume : 1.0;
            occupied_[location(site)] += 1.0;
        };
        if (large) { for (const auto site : large_footprint(cells.anchor(slot))) add(site); }
        else add(cells.anchor(slot));
    }
    for (auto it = refractory_.begin(); it != refractory_.end();) {
        if (!alive.contains(it->first)) it = refractory_.erase(it);
        else ++it;
    }
    if (config_.schema_version>=8 || geometry_.source_mode == "abm_perfusion" || geometry_.source_mode == "abm_plus_synthetic_line") {
        for (const auto site : vessels.occupied_sites()) {
            if (vessels.perfused(site) && contains_resource_site(site)) vessels_[location(site)] = 1.0;
        }
    }
    for (int z = 0; z < geometry_.shape[2]; ++z) {
        for (int y = 0; y < geometry_.shape[1]; ++y) {
            for (int x = 0; x < geometry_.shape[0]; ++x) {
                if (geometry_.source_voxel({x, y, z})) vessels_[index(x, y, z)] = 1.0;
            }
        }
    }
    rebuild_sources();
}

void SharedResourceEnvironment3D::rebuild_sources() {
    const auto& n = config_.continuum.nutrient;
    if (n.boundary_mode == "moving_tumor_front_dirichlet_v2" ||
        n.boundary_mode == "moving_tumor_front_and_vessels_dirichlet_v2") {
        continuum::build_moving_tumor_front_mask_2d(occupied_, geometry_.shape[0], geometry_.shape[1],
            n.tumor_front_smoothing_radius_voxels, n.tumor_front_density_threshold,
            front_workspace_, tumour_mask_);
    }
}
bool SharedResourceEnvironment3D::source(int x, int y, int z, std::size_t here) const noexcept {
    const auto& mode = config_.continuum.nutrient.boundary_mode;
    if (mode == "moving_tumor_front_dirichlet_v2" || mode == "moving_tumor_front_and_vessels_dirichlet_v2") {
        return tumour_mask_[here] == 0U || (mode == "moving_tumor_front_and_vessels_dirichlet_v2" && vessels_[here] > 0.0);
    }
    const bool edge = x == 0 || y == 0 || x + 1 == geometry_.shape[0] || y + 1 == geometry_.shape[1] ||
        (!geometry_.thin_layer && (z == 0 || z + 1 == geometry_.shape[2]));
    return (mode != "vessels_dirichlet_v1" && edge) || (mode != "planar_edges_dirichlet_v1" && vessels_[here] > 0.0);
}

EnvironmentInitializationResult3D SharedResourceEnvironment3D::initialize(double now, const CellStore3D& cells,
                                                                          const SparseVesselGrid3D& vessels) {
    if (externally_driven_) return {false};
    if (restored_) {
        if (now < last_refresh_ - 1.0e-10 || now > next_refresh_ + 1.0e-10) throw std::runtime_error("shared resource restore time mismatch");
        restored_ = false;
        return {false};
    }
    assemble(cells, vessels);
    if(angiogenesis_) angiogenesis_->initialize(vessels_);
    for (int z = 0; z < geometry_.shape[2]; ++z) {
        for (int y = 0; y < geometry_.shape[1]; ++y) {
            for (int x = 0; x < geometry_.shape[0]; ++x) {
                if (source(x, y, z, index(x, y, z))) nutrient_[index(x, y, z)] = config_.continuum.nutrient.vessel_value;
            }
        }
    }
    last_refresh_ = now;
    next_refresh_ = now + config_.continuum.time_step_hours;
    refreshes_ = 1;
    rebuild_prefix();
    return {true};
}

void SharedResourceEnvironment3D::refresh(double now, const CellStore3D& cells, const SparseVesselGrid3D& vessels) {
    if (std::abs(now - next_refresh_) > 1.0e-9) throw std::logic_error("shared resource refresh outside schedule");
    const double dt = now - last_refresh_;
    assemble(cells, vessels);
    const auto& n = config_.continuum.nutrient;
    if(angiogenesis_) { angiogenesis_->initialize(vessels_); angiogenesis_->advance(dt,consumers_,nutrient_,n.vessel_value,false); }
    const double mu = n.diffusion_voxels2_per_hour * dt;
    const double decay = std::exp(-n.decay_per_hour * dt);
    const int dimensions = geometry_.thin_layer ? 2 : 3;
    const int nx = geometry_.shape[0], ny = geometry_.shape[1], nz = geometry_.shape[2];
    deterministic_parallel_for(nutrient_.size(), config_.continuum.base.threads, [&](std::size_t here) {
        const int x = static_cast<int>(here % nx);
        const auto yz = here / nx;
        const int y = static_cast<int>(yz % ny), z = static_cast<int>(yz / ny);
        if (source(x, y, z, here)) { next_[here] = n.vessel_value; return; }
        const double old = nutrient_[here];
        double neighbours = (x > 0 ? nutrient_[index(x - 1, y, z)] : old) +
            (x + 1 < nx ? nutrient_[index(x + 1, y, z)] : old) +
            (y > 0 ? nutrient_[index(x, y - 1, z)] : old) +
            (y + 1 < ny ? nutrient_[index(x, y + 1, z)] : old);
        if (!geometry_.thin_layer) neighbours += (z > 0 ? nutrient_[index(x, y, z - 1)] : old) +
            (z + 1 < nz ? nutrient_[index(x, y, z + 1)] : old);
        const double diffused = std::clamp(old + mu * (neighbours - 2.0 * dimensions * old), 0.0, n.vessel_value) * decay;
        next_[here] = continuum::resource_after_uptake(diffused, consumers_[here], n.K_consumption_rate_per_hour,
            n.K_consumption_half_saturation, dt, n.vessel_value);
        if(angiogenesis_) next_[here]=n.vessel_value-(n.vessel_value-next_[here])*
            std::exp(-dt*config_.continuum.angiogenesis.perfusion_exchange_per_hour*vessels_[here]);
    });
    nutrient_.swap(next_);
    last_refresh_ = now;
    next_refresh_ = now + config_.continuum.time_step_hours;
    ++generation_;
    ++refreshes_;
    rebuild_prefix();
}

void SharedResourceEnvironment3D::rebuild_prefix() {
    if (!geometry_.thin_layer) return;
    const int nx = geometry_.shape[0], ny = geometry_.shape[1];
    row_prefix_.assign(static_cast<std::size_t>(nx + 1) * ny, 0.0);
    for (int y = 0; y < ny; ++y) {
        double sum = 0.0;
        for (int x = 0; x < nx; ++x) {
            sum += nutrient_[index(x, y, 0)];
            row_prefix_[static_cast<std::size_t>(y) * (nx + 1) + x + 1] = sum;
        }
    }
}

double SharedResourceEnvironment3D::nutrient_direction_weight(Vec3i site, DirectionId direction) const {
    if (direction == 0 || direction > 26 || !contains_resource_site(site)) return 0.0;
    const auto here = location(site);
    const int nx = geometry_.shape[0], ny = geometry_.shape[1];
    const int x = static_cast<int>(here % nx), y = static_cast<int>((here / nx) % ny);
    double sum = 0.0;
    std::size_t count = 0;
    if (geometry_.thin_layer) {
        for (const auto span : sector_rows_[direction]) {
            const int row = y + span.dy;
            if (row < 0 || row >= ny) continue;
            const int left = std::max(0, x + span.dx0), right = std::min(nx, x + span.dx1 + 1);
            if (right <= left) continue;
            const auto start = static_cast<std::size_t>(row) * (nx + 1);
            sum += row_prefix_[start + right] - row_prefix_[start + left];
            count += right - left;
        }
    } else {
        for (const auto offset : sector_offsets_[direction]) {
            if (contains_resource_site(site + offset)) { sum += nutrient_[location(site + offset)]; ++count; }
        }
    }
    if (count == 0) return 0.0;
    double gradient = sum / count - nutrient_[here];
    if (std::abs(gradient) <= config_.migration.zero_gradient_tolerance) gradient = 0.0;
    return std::exp(std::clamp(config_.migration.chemotaxis_strength * gradient, -40.0, 40.0));
}

bool SharedResourceEnvironment3D::activation_ready(CellUid uid, double now, double density) {
    auto& state = refractory_[uid];
    if (!state.armed && now >= state.until && density <= config_.migration.reactivation_density_threshold) state.armed = true;
    return state.armed;
}
void SharedResourceEnvironment3D::activation_expired(CellUid uid, double now) {
    refractory_[uid] = {now + config_.migration.reactivation_cooldown_hours, false};
}
std::size_t SharedResourceEnvironment3D::allocated_bytes() const noexcept {
    std::size_t result = (nutrient_.capacity() + next_.capacity() + consumers_.capacity() + occupied_.capacity() +
        vessels_.capacity() + row_prefix_.capacity()) * sizeof(double) + tumour_mask_.capacity();
    for (const auto& offsets : sector_offsets_) result += offsets.capacity() * sizeof(Vec3i);
    if(angiogenesis_) result+=5*nutrient_.size()*sizeof(double);
    return result + refractory_.size() * sizeof(std::pair<CellUid, Refractory>);
}
std::uint64_t SharedResourceEnvironment3D::field_checksum() const noexcept {
    auto state = mix(fingerprint_, std::bit_cast<std::uint64_t>(last_refresh_));
    state = mix(state, std::bit_cast<std::uint64_t>(next_refresh_));
    state = mix(mix(state, generation_), refreshes_);
    for (const auto* field : {&nutrient_, &vessels_}) for (const double value : *field) state = mix(state, std::bit_cast<std::uint64_t>(value));
    for (const auto value : tumour_mask_) state = mix(state, value);
    for (const auto& [uid, refractory] : refractory_) {
        state = mix(mix(mix(state, uid), std::bit_cast<std::uint64_t>(refractory.until)), refractory.armed);
    }
    if(angiogenesis_) state=mix(state,angiogenesis_->checksum());
    return state;
}

void SharedResourceEnvironment3D::save_checkpoint(const std::filesystem::path& path, std::uint64_t abm_checksum) const {
    if (std::filesystem::exists(path)) throw std::runtime_error("refusing to overwrite shared resource checkpoint");
    const auto temporary = path.string() + ".tmp";
    std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
    write(out, std::uint64_t{0x4154434753525631ULL});
    write(out, fingerprint_); write(out, abm_checksum);
    write(out, last_refresh_); write(out, next_refresh_); write(out, generation_); write(out, refreshes_);
    for (const auto* field : {&nutrient_, &vessels_}) {
        write(out, static_cast<std::uint64_t>(field->size()));
        for (const double value : *field) write(out, value);
    }
    for (const auto value : tumour_mask_) write(out, value);
    write(out, static_cast<std::uint64_t>(refractory_.size()));
    for (const auto& [uid, refractory] : refractory_) {
        write(out, uid); write(out, refractory.until); write(out, static_cast<std::uint8_t>(refractory.armed));
    }
    if(angiogenesis_) angiogenesis_->save(out);
    write(out, field_checksum()); out.close();
    std::filesystem::rename(temporary, path);
}
void SharedResourceEnvironment3D::load_checkpoint(const std::filesystem::path& path, std::uint64_t abm_checksum) {
    std::ifstream in(path, std::ios::binary);
    if (read<std::uint64_t>(in) != 0x4154434753525631ULL || read<std::uint64_t>(in) != fingerprint_ ||
        read<std::uint64_t>(in) != abm_checksum) throw std::runtime_error("shared resource checkpoint identity mismatch");
    last_refresh_ = read<double>(in); next_refresh_ = read<double>(in);
    generation_ = read<std::uint32_t>(in); refreshes_ = read<std::uint64_t>(in);
    if (!std::isfinite(last_refresh_) || !std::isfinite(next_refresh_) || last_refresh_ < 0.0 ||
        next_refresh_ <= last_refresh_ || generation_ == 0 || refreshes_ == 0) {
        throw std::runtime_error("invalid shared resource checkpoint clock");
    }
    for (auto* field : {&nutrient_, &vessels_}) {
        const double maximum = field == &nutrient_ ? config_.continuum.nutrient.vessel_value : 1.0;
        if (read<std::uint64_t>(in) != field->size()) throw std::runtime_error("shared resource checkpoint size mismatch");
        for (double& value : *field) {
            value = read<double>(in);
            if (!std::isfinite(value) || value < 0.0 || value > maximum) {
                throw std::runtime_error("invalid shared resource checkpoint field");
            }
        }
    }
    for (auto& value : tumour_mask_) {
        value = read<std::uint8_t>(in);
        if (value > 1U) throw std::runtime_error("invalid shared tumour mask");
    }
    const auto count = read<std::uint64_t>(in);
    if (count > 100000000ULL) throw std::runtime_error("shared refractory checkpoint too large");
    refractory_.clear();
    for (std::uint64_t i = 0; i < count; ++i) {
        const auto uid = read<CellUid>(in); const double until = read<double>(in); const auto armed = read<std::uint8_t>(in);
        if (!std::isfinite(until) || until < 0.0 || armed > 1 || !refractory_.emplace(uid, Refractory{until, armed != 0}).second) {
            throw std::runtime_error("invalid shared refractory checkpoint");
        }
    }
    if(angiogenesis_) { angiogenesis_->load(in); if(angiogenesis_->vessels()!=vessels_) throw std::runtime_error("shared vascular checkpoint mismatch"); }
    if (read<std::uint64_t>(in) != field_checksum() || in.peek() != std::char_traits<char>::eof()) {
        throw std::runtime_error("shared resource checkpoint checksum mismatch");
    }
    restored_ = true;
    rebuild_prefix();
}

}  // namespace atcg3d::shared_rules
