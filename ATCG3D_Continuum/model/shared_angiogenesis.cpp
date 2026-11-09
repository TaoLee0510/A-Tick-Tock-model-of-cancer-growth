#include "model/shared_angiogenesis.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <istream>
#include <limits>
#include <numeric>
#include <ostream>
#include <stdexcept>

#include "engine/parallelism.hpp"

namespace atcg3d::continuum {
namespace {
constexpr double kActivityCutoff = 1.0e-14;
constexpr std::size_t kMaximumIndividualTips = 2000000;

template <class T> void write_value(std::ostream& out, T value) {
    out.write(reinterpret_cast<const char*>(&value), sizeof(T));
    if (!out) throw std::runtime_error("unable to save shared angiogenesis");
}

template <class T> T read_value(std::istream& in) {
    T value{};
    in.read(reinterpret_cast<char*>(&value), sizeof(T));
    if (!in) throw std::runtime_error("truncated shared angiogenesis checkpoint");
    return value;
}

std::uint64_t mix(std::uint64_t hash, std::uint64_t value) noexcept {
    hash ^= value;
    return hash * 1099511628211ULL;
}
}  // namespace

SharedAngiogenesis3D::SharedAngiogenesis3D(
    AngiogenesisFieldConfig3D config, std::array<int, 3> shape,
    double spacing, bool thin, int threads)
    : config_(std::move(config)), shape_(shape), spacing_(spacing),
      measure_(std::pow(spacing, thin ? 2 : 3)),
      cross_section_(thin ? 2.0 * config_.vessel_radius_voxels :
          std::acos(-1.0) * config_.vessel_radius_voxels * config_.vessel_radius_voxels),
      dimensions_(thin ? 2 : 3), threads_(std::max(1, threads)) {
    config_.validate();
    if (config_.model != "shared_vegf_lattice_v2" || !std::isfinite(spacing) || spacing <= 0.0 ||
        shape[0] <= 0 || shape[1] <= 0 || shape[2] <= 0 || (thin && shape[2] != 1) || threads < 1) {
        throw std::invalid_argument("invalid shared vascular grid or representation");
    }
    const auto size = static_cast<std::size_t>(shape[0]) * shape[1] * shape[2];
    for (auto* field : {&taf_, &tips_, &vessels_, &work_taf_, &work_tips_,
                       &hypoxia_, &surface_, &rates_, &moved_, &branch_counts_,
                       &connection_counts_, &discarded_taf_, &discarded_tips_}) {
        field->assign(size, 0.0);
    }
    for (int axis = 0; axis < dimensions_; ++axis) {
        edges_[axis].assign(size, 0.0);
        work_edges_[axis].assign(size, 0.0);
    }
}

void SharedAngiogenesis3D::use_individual_tips(std::uint64_t seed) {
    if (support_.valid || random_counter_ != 0 || !individuals_.empty()) {
        throw std::logic_error("cannot change an initialized vascular representation");
    }
    individual_ = true;
    seed_ = seed;
}

void SharedAngiogenesis3D::initialize(std::vector<double> vessels) {
    if (vessels.size() != vessels_.size()) {
        throw std::invalid_argument("shared vascular field shape mismatch");
    }
    for (const double value : vessels) {
        if (!std::isfinite(value) || value < 0.0 || value > 1.0) {
            throw std::invalid_argument("invalid shared vessel fraction");
        }
    }
    vessels_ = std::move(vessels);
    diagnostics_.initial_perfused_volume = measure_ *
        std::accumulate(vessels_.begin(), vessels_.end(), 0.0);
}

void SharedAngiogenesis3D::initialize_tip_density(const std::vector<double>& tips) {
    if (tips.size() != tips_.size() || diagnostics_.seeded_tips != 0.0 || !individuals_.empty()) {
        throw std::invalid_argument("tip initialization requires an unseeded vascular field");
    }
    double total = 0.0;
    for (std::size_t here = 0; here < tips.size(); ++here) {
        const double mass = tips[here] * measure_;
        if (!std::isfinite(mass) || mass < 0.0 ||
            (individual_ && mass != std::floor(mass))) {
            throw std::invalid_argument("individual tip initialization requires integer site counts");
        }
        total += mass;
        if (total > kMaximumIndividualTips) throw std::invalid_argument("initial vascular tip budget exceeded");
    }
    tips_ = tips;
    for (std::size_t here = 0; here < tips.size(); ++here) {
        if (tips[here] > 0.0) include(support_, here);
        if (individual_) {
            const auto count = static_cast<std::size_t>(tips[here] * measure_);
            for (std::size_t tip = 0; tip < count; ++tip) individuals_.push_back({next_uid_++, here});
        }
    }
    diagnostics_.seeded_tips = total;
}

std::size_t SharedAngiogenesis3D::index(int x, int y, int z) const noexcept {
    return (static_cast<std::size_t>(z) * shape_[1] + y) * shape_[0] + x;
}

std::array<int, 3> SharedAngiogenesis3D::coordinate(std::size_t location) const noexcept {
    const int x = static_cast<int>(location % shape_[0]);
    const auto yz = location / shape_[0];
    return {x, static_cast<int>(yz % shape_[1]), static_cast<int>(yz / shape_[1])};
}

bool SharedAngiogenesis3D::neighbor(std::size_t here, int face,
                                   std::size_t& other) const noexcept {
    auto point = coordinate(here);
    const int axis = face / 2;
    point[axis] += face % 2 == 0 ? -1 : 1;
    if (point[axis] < 0 || point[axis] >= shape_[axis]) return false;
    other = index(point[0], point[1], point[2]);
    return true;
}

void SharedAngiogenesis3D::include(Bounds& bounds, std::size_t location) const noexcept {
    const auto point = coordinate(location);
    if (!bounds.valid) {
        bounds.lower = bounds.upper = point;
        bounds.valid = true;
        return;
    }
    for (int axis = 0; axis < 3; ++axis) {
        bounds.lower[axis] = std::min(bounds.lower[axis], point[axis]);
        bounds.upper[axis] = std::max(bounds.upper[axis], point[axis]);
    }
}

SharedAngiogenesis3D::Bounds SharedAngiogenesis3D::expanded(Bounds bounds) const noexcept {
    if (!bounds.valid) return bounds;
    for (int axis = 0; axis < dimensions_; ++axis) {
        bounds.lower[axis] = std::max(0, bounds.lower[axis] - 1);
        bounds.upper[axis] = std::min(shape_[axis] - 1, bounds.upper[axis] + 1);
    }
    return bounds;
}

std::size_t SharedAngiogenesis3D::count(const Bounds& bounds) const noexcept {
    if (!bounds.valid) return 0;
    return static_cast<std::size_t>(bounds.upper[0] - bounds.lower[0] + 1) *
        (bounds.upper[1] - bounds.lower[1] + 1) *
        (bounds.upper[2] - bounds.lower[2] + 1);
}

std::size_t SharedAngiogenesis3D::bounded_location(
    const Bounds& bounds, std::size_t offset) const noexcept {
    const int nx = bounds.upper[0] - bounds.lower[0] + 1;
    const int ny = bounds.upper[1] - bounds.lower[1] + 1;
    return index(bounds.lower[0] + static_cast<int>(offset % nx),
        bounds.lower[1] + static_cast<int>((offset / nx) % ny),
        bounds.lower[2] + static_cast<int>(offset / (static_cast<std::size_t>(nx) * ny)));
}

SharedAngiogenesis3D::Bounds SharedAngiogenesis3D::prepare_sources(
    const std::vector<double>& cells, const std::vector<double>& nutrient, double maximum,
    const Bounds* consumer_bounds) {
    Bounds sources;
    double total = 0.0, hypoxic = 0.0;
    Bounds input{{0, 0, 0}, {shape_[0] - 1, shape_[1] - 1, shape_[2] - 1}, true};
    if (consumer_bounds) input = *consumer_bounds;
    for (int axis = 0; axis < 3; ++axis) {
        if (input.valid && (input.lower[axis] < 0 || input.upper[axis] >= shape_[axis] ||
                            input.lower[axis] > input.upper[axis])) {
            throw std::invalid_argument("invalid vascular consumer range");
        }
    }
    for (std::size_t offset = 0; offset < count(previous_sources_); ++offset) {
        const auto here = bounded_location(previous_sources_, offset);
        hypoxia_[here] = surface_[here] = 0.0;
    }
    for (std::size_t offset = 0; offset < count(input); ++offset) {
        const auto here = bounded_location(input, offset);
        if (!std::isfinite(cells[here]) || cells[here] < 0.0 ||
            !std::isfinite(nutrient[here]) || nutrient[here] < 0.0) {
            throw std::invalid_argument("invalid shared vascular consumer/resource field");
        }
        hypoxia_[here] = std::clamp(1.0 - nutrient[here] /
            (maximum * config_.hypoxia_threshold), 0.0, 1.0) * cells[here];
        total += cells[here];
        if (nutrient[here] < maximum * config_.hypoxia_threshold) hypoxic += cells[here];
        if (cells[here] > 0.0) include(sources, here);
    }
    seed_rate_ = total > 0.0 ? config_.seed_tips_per_hour * hypoxic / total : 0.0;
    surface_sum_ = 0.0;
    seed_sites_.clear();
    seed_prefix_.clear();
    for (std::size_t offset = 0; offset < count(sources); ++offset) {
        const auto here = bounded_location(sources, offset);
        bool boundary = false;
        for (int face = 0; face < 2 * dimensions_; ++face) {
            std::size_t other{};
            if (!neighbor(here, face, other) || cells[other] < 0.1 * cells[here]) {
                boundary = true;
            }
        }
        surface_[here] = boundary ? hypoxia_[here] : 0.0;
        surface_sum_ += surface_[here];
    }
    const bool use_all = surface_sum_ == 0.0;
    surface_sum_ = 0.0;
    for (std::size_t offset = 0; offset < count(sources); ++offset) {
        const auto here = bounded_location(sources, offset);
        if (use_all) surface_[here] = hypoxia_[here];
        if (surface_[here] > 0.0) {
            surface_sum_ += surface_[here];
            seed_sites_.push_back(here);
            seed_prefix_.push_back(surface_sum_);
        }
    }
    if (surface_sum_ == 0.0) seed_rate_ = 0.0;
    previous_sources_ = sources;
    return sources;
}

double SharedAngiogenesis3D::jump_rate(std::size_t here, std::size_t other) const noexcept {
    return (config_.tip_diffusion_voxels2_per_hour + config_.tip_chemotaxis *
        std::max(0.0, taf_[other] - taf_[here])) / (spacing_ * spacing_);
}

double SharedAngiogenesis3D::outgoing_rate(std::size_t here) const noexcept {
    double result = 0.0;
    for (int face = 0; face < 2 * dimensions_; ++face) {
        std::size_t other{};
        if (neighbor(here, face, other)) result += jump_rate(here, other);
    }
    return result;
}

double SharedAngiogenesis3D::substep(double remaining, const Bounds& bounds) {
    deterministic_parallel_for(count(bounds), threads_, [&](std::size_t offset) {
        const auto here = bounded_location(bounds, offset);
        rates_[here] = outgoing_rate(here);
    });
    double largest = 2.0 * dimensions_ * config_.taf_diffusion_voxels2_per_hour /
        (spacing_ * spacing_);
    for (std::size_t offset = 0; offset < count(bounds); ++offset) {
        largest = std::max(largest, rates_[bounded_location(bounds, offset)]);
    }
    return std::min(remaining, largest > 0.0 ? 0.45 / largest : remaining);
}

void SharedAngiogenesis3D::advance_taf(double dt, const Bounds& bounds) {
    const double decay = std::exp(-dt * config_.taf_decay_per_hour);
    deterministic_parallel_for(count(bounds), threads_, [&](std::size_t offset) {
        const auto here = bounded_location(bounds, offset);
        double laplacian = 0.0;
        for (int face = 0; face < 2 * dimensions_; ++face) {
            std::size_t other{};
            if (neighbor(here, face, other)) laplacian += taf_[other] - taf_[here];
        }
        const double value = std::max(0.0, taf_[here] + dt *
            config_.taf_diffusion_voxels2_per_hour * laplacian / (spacing_ * spacing_)) *
            decay + dt * config_.taf_production_per_cell_hour * hypoxia_[here];
        discarded_taf_[here] = value < kActivityCutoff ? value * measure_ : 0.0;
        work_taf_[here] = value < kActivityCutoff ? 0.0 : value;
    });
}

void SharedAngiogenesis3D::advance_density_tips(double dt, const Bounds& bounds) {
    deterministic_parallel_for(count(bounds), threads_, [&](std::size_t offset) {
        const auto here = bounded_location(bounds, offset);
        double incoming = 0.0;
        for (int face = 0; face < 2 * dimensions_; ++face) {
            std::size_t other{};
            if (neighbor(here, face, other)) incoming += tips_[other] * jump_rate(other, here);
        }
        const double moved = std::max(0.0, tips_[here] * (1.0 - dt * rates_[here]) + dt * incoming);
        moved_[here] = moved;
        const double branched = moved * std::exp(dt * config_.tip_branching_per_hour * taf_[here]);
        const double surviving = branched * std::exp(-dt * config_.tip_anastomosis_per_hour *
            (vessels_[here] + moved));
        const double seeded = surface_sum_ > 0.0 ?
            dt * seed_rate_ * surface_[here] / (surface_sum_ * measure_) : 0.0;
        const double value = surviving + seeded;
        branch_counts_[here] = (branched - moved) * measure_;
        connection_counts_[here] = (branched - surviving) * measure_;
        discarded_tips_[here] = value < kActivityCutoff ? value * measure_ : 0.0;
        work_tips_[here] = value < kActivityCutoff ? 0.0 : value;
    });
    diagnostics_.seeded_tips += dt * seed_rate_;
}

double SharedAngiogenesis3D::uniform() {
    if (random_counter_ == std::numeric_limits<std::uint64_t>::max()) {
        throw std::runtime_error("vascular random counter exhausted");
    }
    std::uint64_t value = seed_ + (++random_counter_) * 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
    value ^= value >> 31U;
    return (static_cast<double>(value >> 12U) + 0.5) * 0x1.0p-52;
}

std::uint64_t SharedAngiogenesis3D::poisson(double mean) {
    if (!(mean > 0.0)) return 0;
    if (!std::isfinite(mean) || mean > 100.0) {
        throw std::runtime_error("vascular seeding requires a smaller substep");
    }
    const double limit = std::exp(-mean);
    std::uint64_t count = 0;
    double product = 1.0;
    do {
        ++count;
        product *= uniform();
    } while (product > limit);
    return count - 1;
}

std::size_t SharedAngiogenesis3D::seed_location() {
    const double target = uniform() * surface_sum_;
    const auto found = std::lower_bound(seed_prefix_.begin(), seed_prefix_.end(), target);
    const auto offset = std::min(seed_sites_.size() - 1,
        static_cast<std::size_t>(found - seed_prefix_.begin()));
    return seed_sites_[offset];
}

void SharedAngiogenesis3D::advance_individual_tips(double dt, const Bounds& bounds) {
    for (std::size_t offset = 0; offset < count(bounds); ++offset) {
        const auto here = bounded_location(bounds, offset);
        moved_[here] = work_tips_[here] = branch_counts_[here] = connection_counts_[here] = 0.0;
        discarded_tips_[here] = 0.0;
        for (int axis = 0; axis < dimensions_; ++axis) work_edges_[axis][here] = 0.0;
    }
    next_individuals_.clear();
    for (const auto tip : individuals_) {
        auto moved_tip = tip;
        const double choice = uniform();
        double cumulative = 0.0;
        for (int face = 0; face < 2 * dimensions_; ++face) {
            std::size_t other{};
            if (!neighbor(tip.location, face, other)) continue;
            cumulative += dt * jump_rate(tip.location, other);
            if (choice < cumulative) {
                moved_tip.location = other;
                work_edges_[face / 2][std::min(tip.location, other)] = 1.0;
                break;
            }
        }
        next_individuals_.push_back(moved_tip);
        moved_[moved_tip.location] += 1.0 / measure_;
    }
    individuals_.clear();
    for (const auto tip : next_individuals_) {
        const auto here = tip.location;
        const double probability = std::exp(-dt * config_.tip_branching_per_hour * taf_[here]);
        double offspring = 0.0;
        if (probability < 1.0) offspring = std::floor(std::log(uniform()) / std::log1p(-probability));
        if (!std::isfinite(offspring) || offspring > kMaximumIndividualTips) {
            throw std::runtime_error("shared vascular offspring budget exceeded");
        }
        const auto children = static_cast<std::size_t>(offspring);
        branch_counts_[here] += static_cast<double>(children);
        const double survival = std::exp(-dt * config_.tip_anastomosis_per_hour *
            (vessels_[here] + moved_[here]));
        for (std::size_t child = 0; child <= children; ++child) {
            const Tip candidate{child == 0 ? tip.uid : next_uid_++, here};
            if (uniform() < survival) {
                if (individuals_.size() >= kMaximumIndividualTips) {
                    throw std::runtime_error("shared vascular individual budget exceeded");
                }
                individuals_.push_back(candidate);
                work_tips_[here] += 1.0 / measure_;
            } else {
                connection_counts_[here] += 1.0;
            }
        }
    }
    const auto seeds = surface_sum_ > 0.0 ? poisson(dt * seed_rate_) : 0;
    if (seeds > kMaximumIndividualTips - individuals_.size()) {
        throw std::runtime_error("shared vascular seeding budget exceeded");
    }
    for (std::uint64_t seed = 0; seed < seeds; ++seed) {
        const auto here = seed_location();
        individuals_.push_back({next_uid_++, here});
        work_tips_[here] += 1.0 / measure_;
    }
    std::sort(individuals_.begin(), individuals_.end(), [](const Tip& first, const Tip& second) {
        return first.uid < second.uid;
    });
    diagnostics_.seeded_tips += static_cast<double>(seeds);
}

void SharedAngiogenesis3D::deposit_centerlines(double dt, const Bounds& bounds) {
    deterministic_parallel_for(count(bounds), threads_, [&](std::size_t offset) {
        const auto here = bounded_location(bounds, offset);
        for (int axis = 0; axis < dimensions_; ++axis) {
            std::size_t other{};
            if (!neighbor(here, 2 * axis + 1, other)) {
                work_edges_[axis][here] = 0.0;
                continue;
            }
            const double old = edges_[axis][here];
            const double traversals = dt * measure_ * (tips_[here] * jump_rate(here, other) +
                tips_[other] * jump_rate(other, here));
            const double next = individual_ ? (work_edges_[axis][here] > 0.0 ? 1.0 : old) :
                1.0 - (1.0 - old) * std::exp(-traversals);
            work_edges_[axis][here] = next - old;
            edges_[axis][here] = next;
        }
    });
    deterministic_parallel_for(count(bounds), threads_, [&](std::size_t offset) {
        const auto here = bounded_location(bounds, offset);
        double new_edges = 0.0;
        for (int axis = 0; axis < dimensions_; ++axis) {
            new_edges += work_edges_[axis][here];
            std::size_t other{};
            if (neighbor(here, 2 * axis, other)) new_edges += work_edges_[axis][other];
        }
        const double dose = 0.5 * cross_section_ * spacing_ * new_edges / measure_;
        vessels_[here] = 1.0 - (1.0 - vessels_[here]) * std::exp(-dose);
    });
}

void SharedAngiogenesis3D::reduce_step(const Bounds& bounds) {
    support_ = {};
    for (std::size_t offset = 0; offset < count(bounds); ++offset) {
        const auto here = bounded_location(bounds, offset);
        diagnostics_.branches += branch_counts_[here];
        diagnostics_.anastomoses += connection_counts_[here];
        diagnostics_.discarded_taf_mass += discarded_taf_[here];
        diagnostics_.discarded_tip_mass += discarded_tips_[here];
        for (int axis = 0; axis < dimensions_; ++axis) {
            diagnostics_.centerline_growth += spacing_ * work_edges_[axis][here];
        }
        if (work_taf_[here] > 0.0 || work_tips_[here] > 0.0) include(support_, here);
    }
    diagnostics_.active_voxels = count(support_);
}

void SharedAngiogenesis3D::clear_previous_workspace(const Bounds& bounds) {
    deterministic_parallel_for(count(workspace_bounds_), threads_, [&](std::size_t offset) {
        const auto here = bounded_location(workspace_bounds_, offset);
        const auto point = coordinate(here);
        bool outside = !bounds.valid;
        for (int axis = 0; axis < 3; ++axis) {
            outside = outside || point[axis] < bounds.lower[axis] || point[axis] > bounds.upper[axis];
        }
        if (!outside) return;
        work_taf_[here] = work_tips_[here] = 0.0;
        for (int axis = 0; axis < dimensions_; ++axis) work_edges_[axis][here] = 0.0;
    });
}

void SharedAngiogenesis3D::advance(double dt, const std::vector<double>& cells,
                                  const std::vector<double>& nutrient,
                                  double maximum, bool grow,
                                  const VascularConsumerBounds3D* consumer_bounds) {
    if (dt < 0.0 || !std::isfinite(dt) || cells.size() != taf_.size() ||
        nutrient.size() != taf_.size() || !std::isfinite(maximum) || maximum <= 0.0) {
        throw std::invalid_argument("invalid shared vascular update");
    }
    const auto sources = prepare_sources(cells, nutrient, maximum, consumer_bounds);
    double elapsed = 0.0;
    while (elapsed < dt) {
        Bounds bounds = support_;
        if (sources.valid) {
            include(bounds, index(sources.lower[0], sources.lower[1], sources.lower[2]));
            include(bounds, index(sources.upper[0], sources.upper[1], sources.upper[2]));
        }
        bounds = expanded(bounds);
        clear_previous_workspace(bounds);
        if (!bounds.valid) return;
        double step = substep(dt - elapsed, bounds);
        if (seed_rate_ > 0.0) step = std::min(step, 100.0 / seed_rate_);
        if (!(step > 0.0) || elapsed + step == elapsed) {
            throw std::runtime_error("shared vascular substep underflow");
        }
        advance_taf(step, bounds);
        for (std::size_t offset = 0; offset < count(bounds); ++offset) {
            if (!std::isfinite(work_taf_[bounded_location(bounds, offset)])) {
                throw std::runtime_error("nonfinite shared VEGF field");
            }
        }
        if (grow) {
            if (individual_) advance_individual_tips(step, bounds);
            else advance_density_tips(step, bounds);
            for (std::size_t offset = 0; offset < count(bounds); ++offset) {
                const auto here = bounded_location(bounds, offset);
                if (!std::isfinite(work_tips_[here]) || !std::isfinite(branch_counts_[here]) ||
                    !std::isfinite(connection_counts_[here])) {
                    throw std::runtime_error("nonfinite shared tip density");
                }
            }
            deposit_centerlines(step, bounds);
        } else {
            for (std::size_t offset = 0; offset < count(bounds); ++offset) {
                const auto here = bounded_location(bounds, offset);
                work_tips_[here] = branch_counts_[here] = connection_counts_[here] = discarded_tips_[here] = 0.0;
                for (int axis = 0; axis < dimensions_; ++axis) work_edges_[axis][here] = 0.0;
            }
        }
        reduce_step(bounds);
        taf_.swap(work_taf_);
        tips_.swap(work_tips_);
        workspace_bounds_ = bounds;
        elapsed += step;
    }
}

std::size_t SharedAngiogenesis3D::allocated_bytes() const noexcept {
    std::size_t result = (individuals_.capacity() + next_individuals_.capacity()) * sizeof(Tip);
    for (const auto* field : {&taf_, &tips_, &vessels_, &work_taf_, &work_tips_,
                              &hypoxia_, &surface_, &rates_, &moved_, &branch_counts_,
                              &connection_counts_, &discarded_taf_, &discarded_tips_}) {
        result += field->capacity() * sizeof(double);
    }
    for (int axis = 0; axis < dimensions_; ++axis) {
        result += (edges_[axis].capacity() + work_edges_[axis].capacity()) * sizeof(double);
    }
    return result + seed_sites_.capacity() * sizeof(std::size_t) + seed_prefix_.capacity() * sizeof(double);
}

std::uint64_t SharedAngiogenesis3D::checksum() const {
    auto hash = mix(config_.fingerprint(), individual_);
    hash = mix(mix(mix(hash, seed_), random_counter_), next_uid_);
    for (const auto* field : {&taf_, &tips_, &vessels_}) {
        for (const auto value : *field) hash = mix(hash, std::bit_cast<std::uint64_t>(value));
    }
    for (int axis = 0; axis < dimensions_; ++axis) {
        for (const auto value : edges_[axis]) hash = mix(hash, std::bit_cast<std::uint64_t>(value));
    }
    for (const auto value : {diagnostics_.centerline_growth, diagnostics_.initial_perfused_volume,
                             diagnostics_.seeded_tips, diagnostics_.branches, diagnostics_.anastomoses,
                             diagnostics_.discarded_taf_mass, diagnostics_.discarded_tip_mass}) {
        hash = mix(hash, std::bit_cast<std::uint64_t>(value));
    }
    for (const auto tip : individuals_) hash = mix(mix(hash, tip.uid), tip.location);
    hash = mix(hash, support_.valid);
    for (int axis = 0; axis < 3; ++axis) {
        hash = mix(mix(hash, support_.lower[axis]), support_.upper[axis]);
    }
    return hash;
}

void SharedAngiogenesis3D::save(std::ostream& out) const {
    write_value(out, std::uint64_t{0x4154434756415332ULL});
    write_value(out, config_.fingerprint());
    for (const auto value : shape_) write_value(out, value);
    write_value(out, spacing_);
    write_value(out, static_cast<std::uint8_t>(individual_));
    write_value(out, seed_);
    write_value(out, random_counter_);
    write_value(out, next_uid_);
    write_value(out, static_cast<std::uint8_t>(support_.valid));
    for (int axis = 0; axis < 3; ++axis) {
        write_value(out, support_.lower[axis]);
        write_value(out, support_.upper[axis]);
    }
    for (const auto value : {diagnostics_.centerline_growth, diagnostics_.initial_perfused_volume,
                             diagnostics_.seeded_tips, diagnostics_.branches, diagnostics_.anastomoses,
                             diagnostics_.discarded_taf_mass, diagnostics_.discarded_tip_mass}) {
        write_value(out, value);
    }
    for (const auto* field : {&taf_, &tips_, &vessels_}) {
        out.write(reinterpret_cast<const char*>(field->data()), field->size() * sizeof(double));
    }
    for (int axis = 0; axis < dimensions_; ++axis) {
        out.write(reinterpret_cast<const char*>(edges_[axis].data()), edges_[axis].size() * sizeof(double));
    }
    write_value(out, static_cast<std::uint64_t>(individuals_.size()));
    for (const auto tip : individuals_) {
        write_value(out, tip.uid);
        write_value(out, static_cast<std::uint64_t>(tip.location));
    }
    if (!out) throw std::runtime_error("unable to save shared vascular fields");
}

void SharedAngiogenesis3D::load(std::istream& in) {
    if (read_value<std::uint64_t>(in) != 0x4154434756415332ULL ||
        read_value<std::uint64_t>(in) != config_.fingerprint()) {
        throw std::runtime_error("shared vascular checkpoint configuration mismatch");
    }
    for (const auto value : shape_) {
        if (read_value<int>(in) != value) {
            throw std::runtime_error("shared vascular checkpoint shape mismatch");
        }
    }
    if (read_value<double>(in) != spacing_ ||
        read_value<std::uint8_t>(in) != static_cast<std::uint8_t>(individual_)) {
        throw std::runtime_error("shared vascular checkpoint representation mismatch");
    }
    seed_ = read_value<std::uint64_t>(in);
    random_counter_ = read_value<std::uint64_t>(in);
    next_uid_ = read_value<std::uint64_t>(in);
    support_.valid = read_value<std::uint8_t>(in) != 0;
    for (int axis = 0; axis < 3; ++axis) {
        support_.lower[axis] = read_value<int>(in);
        support_.upper[axis] = read_value<int>(in);
        if (support_.valid && (support_.lower[axis] < 0 || support_.upper[axis] >= shape_[axis] ||
                              support_.lower[axis] > support_.upper[axis])) {
            throw std::runtime_error("invalid shared vascular checkpoint bounds");
        }
    }
    for (auto* value : {&diagnostics_.centerline_growth, &diagnostics_.initial_perfused_volume,
                       &diagnostics_.seeded_tips, &diagnostics_.branches, &diagnostics_.anastomoses,
                       &diagnostics_.discarded_taf_mass, &diagnostics_.discarded_tip_mass}) {
        *value = read_value<double>(in);
        if (!std::isfinite(*value) || *value < 0.0) {
            throw std::runtime_error("invalid shared vascular checkpoint diagnostics");
        }
    }
    for (auto* field : {&taf_, &tips_, &vessels_}) {
        in.read(reinterpret_cast<char*>(field->data()), field->size() * sizeof(double));
        for (const auto value : *field) {
            if (!std::isfinite(value) || value < 0.0 || (field == &vessels_ && value > 1.0)) {
                throw std::runtime_error("invalid shared vascular checkpoint field");
            }
        }
    }
    for (int axis = 0; axis < dimensions_; ++axis) {
        auto& field = edges_[axis];
        in.read(reinterpret_cast<char*>(field.data()), field.size() * sizeof(double));
        for (const auto value : field) {
            if (!std::isfinite(value) || value < 0.0 || value > 1.0) {
                throw std::runtime_error("invalid vascular centerline checkpoint");
            }
        }
    }
    const auto size = read_value<std::uint64_t>(in);
    if (size > kMaximumIndividualTips || (!individual_ && size != 0)) {
        throw std::runtime_error("invalid vascular individual checkpoint size");
    }
    individuals_.clear();
    for (std::uint64_t tip = 0; tip < size; ++tip) {
        const auto uid = read_value<std::uint64_t>(in);
        const auto location = read_value<std::uint64_t>(in);
        if (location >= taf_.size() || uid == 0 || uid >= next_uid_ ||
            (!individuals_.empty() && uid <= individuals_.back().uid)) {
            throw std::runtime_error("invalid vascular individual checkpoint");
        }
        individuals_.push_back({uid, static_cast<std::size_t>(location)});
    }
    diagnostics_.active_voxels = count(support_);
    workspace_bounds_ = {};
    previous_sources_ = {};
    std::fill(hypoxia_.begin(), hypoxia_.end(), 0.0);
    std::fill(surface_.begin(), surface_.end(), 0.0);
    std::fill(work_taf_.begin(), work_taf_.end(), 0.0);
    std::fill(work_tips_.begin(), work_tips_.end(), 0.0);
    for (int axis = 0; axis < dimensions_; ++axis) {
        std::fill(work_edges_[axis].begin(), work_edges_[axis].end(), 0.0);
    }
    if (!in) throw std::runtime_error("truncated shared vascular checkpoint fields");
}

}  // namespace atcg3d::continuum
