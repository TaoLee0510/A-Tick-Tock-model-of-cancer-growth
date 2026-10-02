#pragma once

#include <array>
#include <cstdint>
#include <functional>
#include <iosfwd>
#include <map>
#include <unordered_map>
#include <vector>

#include "config/model_config.hpp"

namespace atcg3d::structured_pde {

// Cell-number density resolved by remaining division work. Channels are
// r-small, r-large, K-small and K-large. Activity is a separate marginal;
// conservative native fluxes carry the local work mixture with their mass.
class DivisionRenewal3D {
public:
    DivisionRenewal3D(DivisionTimingConfig timing, double width, double maximum,
                      std::array<double, 2> inherent);
    DivisionRenewal3D(double width, double maximum,
                      std::array<std::vector<double>, 2> kernel);
    void add(std::size_t location, std::size_t channel, double mass, double work);
    void add_fresh(std::size_t location, std::size_t channel, double mass);
    void erase(std::size_t location);
    void begin_transport();
    void begin_weighted_transport(std::array<std::vector<double>, 4> weights);
    double weighted_transport_mass(std::size_t location, std::size_t channel) const;
    void transfer(std::size_t source, std::size_t target, std::size_t channel,
                  double mass);
    void finish_transport(const std::function<double(std::size_t, std::size_t)>& total);
    double advance(std::size_t location, std::size_t channel, double work);
    void reconcile(std::size_t location, std::size_t channel, double total);
    double mass(std::size_t location, std::size_t channel) const;
    double mean_work(std::size_t location, std::size_t channel) const;
    std::uint64_t checksum() const;
    void save(std::ostream& out) const;
    void load(std::istream& in, std::size_t voxels);
private:
    using Node = std::array<std::vector<double>, 4>;
    Node& node(std::map<std::size_t, Node>& state, std::size_t location);
    void deposit(std::vector<double>& values, double mass, double work) const;
    double width_;
    std::size_t bins_;
    std::array<std::vector<double>, 2> kernel_;
    std::map<std::size_t, Node> state_, work_;
    std::map<std::size_t, std::array<double, 4>> totals_;
    bool compact_channels_{};
    std::unordered_map<std::size_t, const Node*> source_cache_;
    std::unordered_map<std::size_t, Node*> transport_cache_;
    std::unordered_map<std::size_t, std::array<double, 4>> compact_totals_;
    std::array<std::vector<double>, 4> transport_weights_;
};

} // namespace atcg3d::structured_pde
