#include "model/structured_pde_model.hpp"
#include "model/beta_duration.hpp"
#include "model/feasible_jump.hpp"
#include "model/shared_resource.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <utility>

#include "common/density_growth_rule.hpp"
#include "engine/parallelism.hpp"
#include "engine/simulation.hpp"
#include "geometry/directions.hpp"
#include "geometry/footprint.hpp"
#include "rules/initial_rates.hpp"

#include "model/operator_math.hpp"

namespace atcg3d::structured_pde {

void StructuredPdeModel3D::advance_angiogenesis(double dt) {
    if(!angiogenesis_) return;
    if (config_.continuum.angiogenesis.model == "shared_vegf_lattice_v2") {
        vascular_consumers_work_.resize(voxel_count_);
        continuum::VascularConsumerBounds3D bounds;
        const auto nx = config_.continuum.grid.shape[0], ny = config_.continuum.grid.shape[1];
        for (std::size_t here = 0; here < voxel_count_; ++here) {
            double consumers = 0.0;
            for (std::size_t stage = 0; stage < 2; ++stage) {
                consumers += r_normal_[stage][here] + active_total_[stage][here] + K_[stage][here];
            }
            if (!external_r_.empty()) consumers += external_consumers(here);
            vascular_consumers_work_[here] = consumers;
            if (!(consumers > 0.0)) continue;
            const std::array<int, 3> point{static_cast<int>(here % nx),
                static_cast<int>((here / nx) % ny), static_cast<int>(here / (static_cast<std::size_t>(nx) * ny))};
            if (!bounds.valid) {
                bounds.lower = bounds.upper = point;
                bounds.valid = true;
            } else {
                for (int axis = 0; axis < 3; ++axis) {
                    bounds.lower[axis] = std::min(bounds.lower[axis], point[axis]);
                    bounds.upper[axis] = std::max(bounds.upper[axis], point[axis]);
                }
            }
        }
        angiogenesis_->advance(dt, vascular_consumers_work_, nutrient_,
            config_.continuum.nutrient.vessel_value, true, &bounds);
        vessel_ = angiogenesis_->vessels();
        clear_cells_from_vessels();
        return;
    }
    std::vector<double> consumers(voxel_count_,0.0);
    for(std::size_t here=0;here<voxel_count_;++here) {
        for(std::size_t stage=0;stage<2;++stage) consumers[here]+=r_normal_[stage][here]+active_total_[stage][here]+K_[stage][here];
        if (!external_r_.empty()) consumers[here] += external_consumers(here);
    }
    angiogenesis_->advance(dt,consumers,nutrient_,config_.continuum.nutrient.vessel_value);
    vessel_=angiogenesis_->vessels();
    clear_cells_from_vessels();
}

void StructuredPdeModel3D::add_synthetic_vessel() {
    if (operators_.vascular.shared_geometry) {
        const auto geometry = config_.continuum.shared_vascular_geometry();
        const int nx = geometry.shape[0];
        const int ny = geometry.shape[1];
        for (std::size_t location = 0; location < voxel_count_; ++location) {
            const int x = static_cast<int>(location % nx);
            const auto yz = location / nx;
            const int y = static_cast<int>(yz % ny);
            const int z = static_cast<int>(yz / ny);
            if (geometry.source_voxel({x, y, z})) vessel_[location] = 1.0;
        }
        return;
    }
    const auto& vascular = config_.continuum.vascular;
    const int axis = vascular.synthetic_axis == "x" ? 0
        : (vascular.synthetic_axis == "y" ? 1 : 2);
    const double radius_squared = vascular.synthetic_radius_voxels *
        vascular.synthetic_radius_voxels;
    for (std::size_t location = 0; location < voxel_count_; ++location) {
        const auto point = coordinate(location);
        double distance_squared = 0.0;
        for (int dimension = 0; dimension < 3; ++dimension) {
            if (dimension == axis ||
                (config_.continuum.base.thin_layer && dimension == 2)) continue;
            const double offset =
                point[dimension] - vascular.synthetic_center[dimension];
            distance_squared += offset * offset;
        }
        if (distance_squared <= radius_squared) vessel_[location] = 1.0;
    }
}

bool StructuredPdeModel3D::vessel_blocks_cells(
    std::size_t location) const noexcept {
    return config_.migration.vessel_exclusion && (!operators_.vascular.density_exclusion
        ? vessel_[location] > 0.0 : vessel_[location] >= config_.continuum.angiogenesis.exclusion_fraction);
}

bool StructuredPdeModel3D::transport_destination_blocks_cells(
    std::size_t location, StructuredStage3D stage) const {
    if (vessel_blocks_cells(location)) return true;
    if (!external_transport_destination_) return false;
    if (!external_transport_destination_(location)) return true;
    if (stage == StructuredStage3D::small) return false;
    const auto& shape = config_.continuum.grid.shape;
    const int x = static_cast<int>(location % shape[0]);
    const int y = static_cast<int>((location / shape[0]) % shape[1]);
    const int z = static_cast<int>(location /
        (static_cast<std::size_t>(shape[0]) * shape[1]));
    const int depth = config_.continuum.base.thin_layer ? 1 : 2;
    for (int dz = 0; dz < depth; ++dz) {
        for (int dy = 0; dy < 2; ++dy) {
            for (int dx = 0; dx < 2; ++dx) {
                if (x + dx >= shape[0] || y + dy >= shape[1] ||
                    z + dz >= shape[2] ||
                    !external_transport_destination_(index(x + dx, y + dy, z + dz)))
                    return true;
            }
        }
    }
    return false;
}

void StructuredPdeModel3D::clear_cells_from_vessels() {
    if (!config_.migration.vessel_exclusion) return;
    for (std::size_t location = 0; location < voxel_count_; ++location) {
        if (!vessel_blocks_cells(location)) continue;
        if(operators_.vascular.skip_empty_cells&&occupied_fraction(location)==0)continue;
        if (renewal_) renewal_->erase(location);
        if (duration_) duration_->erase(location);
        if (velocity_) velocity_->erase(location);
        if (operators_.vascular.shared_geometry && !initialized_ && occupied_fraction(location) > 1.0e-10) {
            throw std::invalid_argument("structured initial population overlaps an excluded vessel");
        }
        for (std::size_t stage = 0; stage < kStructuredStageCount3D; ++stage) {
            if (operators_.vascular.removal_diagnostics) {
                const auto cell_stage = stage == 0
                    ? CellStage::small : CellStage::large;
                record_vascular_removal(CellType::r, cell_stage, false,
                    r_normal_[stage][location] * voxel_measure_);
                record_vascular_removal(CellType::r, cell_stage, true,
                    active_total_[stage][location] * voxel_measure_);
                record_vascular_removal(CellType::K, cell_stage, false,
                    K_[stage][location] * voxel_measure_);
            }
            r_normal_[stage][location] = 0.0;
            if (operators_.vascular.shared_geometry) {
                r_refractory_[stage][location] = 0.0;
                refractory_clock_[stage][location] = 0.0;
            }
            K_[stage][location] = 0.0;
            active_total_[stage][location] = 0.0F;
            activation_cooldown_[stage][location] = 0.0F;
            activation_armed_[stage][location] = 0U;
            for (std::size_t bucket = 0;
                 bucket < active_direction_[stage].size(); ++bucket) {
                active_direction_[stage][bucket][location] = 0.0F;
                active_clock_[stage][bucket][location] = 0.0F;
            }
        }
    }
    for (std::size_t stage = 0; stage < kStructuredStageCount3D; ++stage) {
        shrink_active_bounds(stage);
    }
    shrink_population_bounds();
}

void StructuredPdeModel3D::record_vascular_removal(
    CellType type, CellStage stage, bool active, double mass) {
    const std::size_t phase = stage == CellStage::large ? 1 : 0;
    if (type == CellType::K)
        vascular_removed_mass_.K[phase] += mass;
    else if (active)
        vascular_removed_mass_.r_active[phase] += mass;
    else
        vascular_removed_mass_.r_normal[phase] += mass;
}

}  // namespace atcg3d::structured_pde
