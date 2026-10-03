#include "model/hybrid_model.hpp"

#include "rules/initial_rates.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <unordered_map>

namespace atcg3d::hybrid {
namespace {
std::uint64_t conversion_key(std::uint64_t value) {
    value += 0x9e3779b97f4a7c15ULL;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
    return value ^ (value >> 31);
}

std::size_t sample_work(const std::vector<double>& work, std::uint64_t key) {
    const double total = std::accumulate(work.begin(), work.end(), 0.0);
    if (!(total > 0.0))
        throw std::logic_error("hybrid K cohort lacks a work distribution");
    double target = (conversion_key(key) >> 11) * 0x1.0p-53 * total;
    for (std::size_t bin = 0; bin < work.size(); ++bin) {
        target -= work[bin];
        if (target < 0.0)
            return bin;
    }
    return work.size() - 1;
}

struct CapacityWeight {
    std::size_t index{};
    double capacity{}, weight{};
};

std::vector<double> water_fill(std::vector<CapacityWeight> candidates,
                               std::size_t size, double mass) {
    std::vector<double> allocation(size);
    double weights = 0.0, capacity = 0.0;
    for (const auto& candidate : candidates) {
        weights += candidate.weight;
        capacity += candidate.capacity;
    }
    const double target = std::min(mass, capacity);
    if (!(weights > 0.0) || !(target > 0.0))
        return allocation;
    std::sort(candidates.begin(), candidates.end(), [](const auto& a, const auto& b) {
        const double left = a.capacity / a.weight;
        const double right = b.capacity / b.weight;
        return left != right ? left < right : a.index < b.index;
    });
    double saturated = 0.0;
    double level = target / weights;
    for (const auto& candidate : candidates) {
        if (level <= candidate.capacity / candidate.weight)
            break;
        saturated += candidate.capacity;
        weights -= candidate.weight;
        level = weights > 0.0 ? std::max(0.0, target - saturated) / weights : 0.0;
    }
    for (const auto& candidate : candidates)
        allocation[candidate.index] = std::min(candidate.capacity, level * candidate.weight);
    return allocation;
}
}

void HybridModel3D::convert_invasion_region(const std::vector<std::size_t>& region) {
    struct KPool {
        double mass{};
        std::vector<double> work;
    };
    auto& p = *pde_;
    std::array<KPool, 2> pools;
    for (const auto i : region)
        for (int stage = 0; stage < 2; ++stage) {
            if (p.r_normal_[stage][i] != 0.0 || p.active_total_[stage][i] != 0.0)
                throw std::logic_error("invasion hybrid contains PDE r mass");
            pools[stage].mass += p.K_[stage][i];
        }
    if (std::max(pools[0].mass, pools[1].mass) < 1.0)
        return;
    std::vector<std::array<double, 2>> original(region.size());
    std::unordered_map<std::size_t, std::size_t> membership;
    membership.reserve(region.size());
    for (std::size_t local = 0; local < region.size(); ++local) {
        const auto i = region[local];
        membership.emplace(i, local);
        for (int stage = 0; stage < 2; ++stage) {
            const double mass = p.K_[stage][i];
            original[local][stage] = mass;
            const auto values = p.renewal_->distribution(i, 2 + stage);
            if (pools[stage].work.empty())
                pools[stage].work.resize(values.size());
            for (std::size_t bin = 0; bin < values.size(); ++bin)
                pools[stage].work[bin] += values[bin];
            p.K_[stage][i] = 0.0;
        }
        p.renewal_->erase(i);
    }
    // Placement reads the live PDE fields directly. Rebuilding the entire
    // coupling grid per region would make exchange quadratic in grid size.
    // The exchange driver rebuilds it once after all regions are committed.
    const double maximum = config_.rules.continuum.reaction.maximum_occupied_fraction;
    double vacant_volume = 0.0;
    for (const auto i : region)
        if (abm_->grid_.occupants(site(i)).empty())
            vacant_volume += maximum;
    double remaining_volume = pools[0].mass + p.large_cell_volume_ * pools[1].mass;
    auto eligible = region;
    for (int stage = 1; stage >= 0; --stage) {
        auto& pool = pools[stage];
        const double volume = stage ? p.large_cell_volume_ : 1.0;
        std::sort(eligible.begin(), eligible.end(), [&](auto a, auto b) {
            const auto key = config_.rules.continuum.base.seed ^
                conversion_key(exchange_count_ + std::uint64_t(8 * stage + 2));
            const auto left = conversion_key(a ^ key), right = conversion_key(b ^ key);
            return left != right ? left < right : a < b;
        });
        for (const auto i : eligible) {
            if (pool.mass < 1.0)
                break;
            const auto anchor = site(i);
            const auto cell_stage = stage ? CellStage::large : CellStage::small;
            const auto points = footprint(anchor, cell_stage);
            bool free = true;
            double consumed_volume = 0.0;
            for (const auto point : points) {
                if (!available(point) || !abm_->grid_.occupants(point).empty())
                    free = false;
                if (membership.contains(location(point)))
                    consumed_volume += maximum;
            }
            if (!free || vacant_volume - consumed_volume + 1e-10 < remaining_volume - volume)
                continue;
            CellInit cell;
            cell.anchor = anchor;
            cell.uid = abm_->next_uid_++;
            cell.clone_id = std::uint32_t(cell.uid);
            cell.type = CellType::K;
            cell.stage = cell_stage;
            cell.last_update_time = time_hours();
            const auto rates = sample_initial_cell_rates(
                abm_->config_, cell.type, abm_->config_.seed, cell.uid);
            cell.inherent_growth_rate = float(rates.inherent_growth_rate);
            cell.migration_rate = float(rates.migration_rate);
            cell.normal_migration_rate = sample_normal_migration_rate(
                cell.type, cell.migration_rate, abm_->config_, cell.uid, 0);
            const auto work_key = config_.rules.continuum.base.seed ^ conversion_key(cell.uid) ^
                conversion_key(exchange_count_) ^ 0x485942574f524bULL;
            cell.division_work_remaining = float(p.renewal_->bin_width() * sample_work(pool.work, work_key));
            cell.next_migration_time = cell.normal_migration_rate > 0.0
                ? time_hours() + 1.0 / cell.normal_migration_rate : 0.0;
            const auto slot = abm_->cells_.create(cell);
            const bool placed = stage ? abm_->grid_.place_large(anchor, slot)
                                      : abm_->grid_.place_single(anchor, slot);
            if (!placed)
                throw std::logic_error("invasion hybrid rounded placement failed");
            abm_->density_.add(anchor, cell.type, slot);
            abm_->schedule_cell(slot);
            for (auto& value : pool.work)
                value *= std::max(0.0, 1.0 - 1.0 / pool.mass);
            pool.mass -= 1.0;
            remaining_volume -= volume;
            vacant_volume -= consumed_volume;
            ++to_abm_;
        }
    }
    std::vector<double> used(region.size());
    for (int stage = 0; stage < 2; ++stage) {
        auto& pool = pools[stage];
        const double volume = stage ? p.large_cell_volume_ : 1.0;
        for (bool uniform : {false, true}) {
            if (!(pool.mass > 0.0))
                break;
            std::vector<CapacityWeight> candidates;
            for (std::size_t local = 0; local < region.size(); ++local) {
                if (!abm_->grid_.occupants(site(region[local])).empty())
                    continue;
                const double capacity = std::max(0.0, maximum - used[local]) / volume;
                const double weight = uniform ? 1.0 : original[local][stage];
                if (capacity > 0.0 && weight > 0.0)
                    candidates.push_back({local, capacity, weight});
            }
            const auto allocation = water_fill(candidates, region.size(), pool.mass);
            for (std::size_t local = 0; local < region.size(); ++local) {
                const double mass = std::min(pool.mass, allocation[local]);
                if (!(mass > 0.0))
                    continue;
                const auto i = region[local];
                p.K_[stage][i] += mass;
                p.renewal_->add_distribution(i, 2 + stage, pool.work, mass / pool.mass);
                for (auto& value : pool.work)
                    value *= std::max(0.0, 1.0 - mass / pool.mass);
                pool.mass -= mass;
                used[local] += volume * mass;
                int x, y, z;
                p.grid_coordinate(site(i), x, y, z);
                p.include_population_location(x, y, z);
            }
        }
        // Preserve the final rounding residue without a third global scan.
        if (pool.mass > 0.0 && pool.mass <= 1e-9)
            for (std::size_t local = 0; local < region.size(); ++local) {
                if (!abm_->grid_.occupants(site(region[local])).empty() ||
                    maximum - used[local] < volume * pool.mass)
                    continue;
                const auto i = region[local];
                p.K_[stage][i] += pool.mass;
                p.renewal_->add_distribution(i, 2 + stage, pool.work, 1.0);
                used[local] += volume * pool.mass;
                pool.mass = 0.0;
                break;
            }
        if (pool.mass > 0.0)
            throw std::logic_error("invasion hybrid exchange exhausted residual capacity");
    }
}
} // namespace atcg3d::hybrid
