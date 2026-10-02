#include "model/division_renewal.hpp"

#include <algorithm>
#include <bit>
#include <cmath>
#include <istream>
#include <numeric>
#include <ostream>
#include <stdexcept>

namespace atcg3d::structured_pde {
namespace {
template<class T> void write(std::ostream& out, T value) {
    out.write(reinterpret_cast<const char*>(&value), sizeof(T));
    if (!out) throw std::runtime_error("division renewal checkpoint write failed");
}
template<class T> T read(std::istream& in) {
    T value{};
    in.read(reinterpret_cast<char*>(&value), sizeof(T));
    if (!in) throw std::runtime_error("truncated division renewal checkpoint");
    return value;
}
}

DivisionRenewal3D::DivisionRenewal3D(DivisionTimingConfig timing, double width,
                                   double maximum, std::array<double, 2> inherent)
    : width_(width) {
    if (!(width > 0.0) || !(maximum > width) || !std::isfinite(width) ||
        !std::isfinite(maximum) || maximum / width > 4096 ||
        !(timing.base_cycle_hours > 0.0) || !(timing.stochastic_time_quantum_hours > 0.0)) {
        throw std::invalid_argument("invalid division renewal work grid/timing");
    }
    bins_ = std::size_t(std::ceil(maximum / width)) + 1;
    for (std::size_t type = 0; type < 2; ++type) {
        kernel_[type].resize(bins_);
        const double rate = inherent[type];
        if (!(rate > 0.0) || !std::isfinite(rate))
            throw std::invalid_argument("invalid division renewal inherent rate");
        const double quantum = timing.stochastic_time_quantum_hours * rate;
        const double minimum = timing.minimum_fraction * timing.base_cycle_hours;
        const double p = std::clamp(quantum /
            std::max(timing.stochastic_tail_fraction * timing.base_cycle_hours, quantum),
            1.0e-12, 1.0);
        double remainder = 1.0;
        for (std::size_t k = 1; remainder > 1.0e-13; ++k) {
            const double work = minimum + k * quantum;
            if (work > maximum) {
                if (remainder > 1.0e-9) {
                    throw std::invalid_argument("division maximum_work truncates geometric tail");
                }
                deposit(kernel_[type], remainder, maximum);
                break;
            }
            const double probability = remainder * p;
            deposit(kernel_[type], probability, work);
            remainder -= probability;
        }
        const double sum = std::accumulate(kernel_[type].begin(), kernel_[type].end(), 0.0);
        for (auto& value : kernel_[type]) value /= sum;
    }
}

DivisionRenewal3D::Node& DivisionRenewal3D::node(
    std::map<std::size_t, Node>& state, std::size_t location) {
    auto& result = state[location];
    for (auto& field : result) if (field.empty()) field.resize(bins_);
    return result;
}

void DivisionRenewal3D::deposit(std::vector<double>& values, double mass, double work) const {
    if (!(mass >= 0.0) || !std::isfinite(mass) || !(work >= 0.0) ||
        !std::isfinite(work) || work > (bins_ - 1) * width_) {
        throw std::invalid_argument("invalid division renewal mass/work");
    }
    const double coordinate = work / width_;
    const auto lower = std::min(bins_ - 1, std::size_t(std::floor(coordinate)));
    const auto upper = std::min(bins_ - 1, lower + 1);
    const double fraction = coordinate - lower;
    values[lower] += mass * (1.0 - fraction);
    values[upper] += mass * fraction;
}

void DivisionRenewal3D::add(std::size_t location, std::size_t channel, double mass, double work) {
    deposit(node(state_, location).at(channel), mass, work);
}
void DivisionRenewal3D::add_fresh(std::size_t location, std::size_t channel, double mass) {
    if (!(mass >= 0.0) || !std::isfinite(mass)) throw std::invalid_argument("invalid renewal birth mass");
    auto& values = node(state_, location).at(channel);
    const auto& kernel = kernel_[channel / 2];
    for (std::size_t i = 0; i < bins_; ++i) values[i] += mass * kernel[i];
}
void DivisionRenewal3D::erase(std::size_t location) { state_.erase(location); }
double DivisionRenewal3D::mass(std::size_t location, std::size_t channel) const {
    const auto found = state_.find(location);
    if (found == state_.end()) return 0.0;
    const auto& values = found->second.at(channel);
    return std::accumulate(values.begin(), values.end(), 0.0);
}
double DivisionRenewal3D::mean_work(std::size_t location, std::size_t channel) const {
    const auto found = state_.find(location);
    if (found == state_.end()) return 0.0;
    double numerator = 0.0;
    for (std::size_t i = 0; i < bins_; ++i) numerator += found->second.at(channel)[i] * i * width_;
    const double total = mass(location, channel);
    return total > 0.0 ? numerator / total : 0.0;
}
void DivisionRenewal3D::begin_transport() {
    work_ = state_;
    totals_.clear();
    for (const auto& [location, fields] : state_) {
        auto& totals = totals_[location];
        for (std::size_t c = 0; c < 4; ++c)
            totals[c] = std::accumulate(fields[c].begin(), fields[c].end(), 0.0);
    }
}
void DivisionRenewal3D::transfer(std::size_t source, std::size_t target,
                               std::size_t channel, double amount) {
    if (!(amount > 0.0) || source == target) return;
    const auto found = state_.find(source);
    if (found == state_.end() || !(totals_.at(source)[channel] > 0.0)) {
        throw std::logic_error("population flux lacks division work");
    }
    const auto& from = found->second[channel];
    auto& outgoing = node(work_, source)[channel];
    auto& incoming = node(work_, target)[channel];
    const double factor = amount / totals_.at(source)[channel];
    for (std::size_t b = 0; b < bins_; ++b) {
        const double moved = from[b] * factor;
        outgoing[b] -= moved;
        incoming[b] += moved;
    }
}
void DivisionRenewal3D::finish_transport(
    const std::function<double(std::size_t, std::size_t)>& total) {
    state_.swap(work_);
    work_.clear();
    for (auto it = state_.begin(); it != state_.end();) {
        bool populated = false;
        for (std::size_t c = 0; c < 4; ++c) {
            auto& values = it->second[c];
            for (auto& value : values) {
                if (value < -1.0e-9) throw std::runtime_error("negative transported division mass");
                value = std::max(0.0, value);
            }
            reconcile(it->first, c, total(it->first, c));
            populated = populated || total(it->first, c) > 0.0;
        }
        if (!populated) it = state_.erase(it); else ++it;
    }
}
double DivisionRenewal3D::advance(std::size_t location, std::size_t channel, double completed_work) {
    auto& values = node(state_, location).at(channel);
    std::vector<double> next(bins_);
    for (std::size_t b = 0; b < bins_; ++b)
        deposit(next, values[b], std::max(0.0, b * width_ - completed_work));
    const double completed = next[0];
    next[0] = 0.0;
    values.swap(next);
    return completed;
}
void DivisionRenewal3D::reconcile(std::size_t location, std::size_t channel, double total) {
    const double present = mass(location, channel);
    if (total > present) add_fresh(location, channel, total - present);
    else {
        auto& values = node(state_, location).at(channel);
        const double factor = present > 0.0 ? std::max(0.0, total) / present : 0.0;
        for (auto& value : values) value *= factor;
    }
}
std::uint64_t DivisionRenewal3D::checksum() const {
    std::uint64_t hash = 1469598103934665603ULL;
    const auto mix = [&](std::uint64_t value) { hash = (hash ^ value) * 1099511628211ULL; };
    for (const auto& [location, fields] : state_) {
        mix(location);
        for (const auto& field : fields) for (double value : field) mix(std::bit_cast<std::uint64_t>(value));
    }
    return hash;
}
void DivisionRenewal3D::save(std::ostream& out) const {
    write(out, std::uint64_t(bins_));
    write(out, std::uint64_t(state_.size()));
    for (const auto& [location, fields] : state_) {
        write(out, std::uint64_t(location));
        for (const auto& field : fields) for (double value : field) write(out, value);
    }
}
void DivisionRenewal3D::load(std::istream& in, std::size_t voxels) {
    if (read<std::uint64_t>(in) != bins_) throw std::runtime_error("division work grid mismatch");
    const auto count = read<std::uint64_t>(in);
    if (count > voxels) throw std::runtime_error("invalid division renewal node count");
    state_.clear();
    for (std::size_t i = 0; i < count; ++i) {
        const auto location = read<std::uint64_t>(in);
        if (location >= voxels || state_.contains(location)) throw std::runtime_error("invalid division renewal node");
        auto& fields = node(state_, location);
        for (auto& field : fields) for (auto& value : field) {
            value = read<double>(in);
            if (!(value >= 0.0) || !std::isfinite(value)) throw std::runtime_error("invalid division renewal density");
        }
    }
}
} // namespace atcg3d::structured_pde
