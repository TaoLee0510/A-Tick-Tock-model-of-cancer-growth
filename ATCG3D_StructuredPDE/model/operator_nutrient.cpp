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

double StructuredPdeModel3D::external_consumers(
    std::size_t location) const noexcept {
    if (!external_r_consumers_.empty())
        return external_r_consumers_[location] + external_K_consumers_[location];
    return external_r_.empty()
        ? 0.0 : external_r_[location] + external_K_[location];
}

double StructuredPdeModel3D::capacity_multiplier(double value) const noexcept {
    if (operators_.nutrient.transient_resources) return 1.0;
    const auto& nutrient = config_.continuum.nutrient;
    const double local = std::clamp(value, 0.0, nutrient.vessel_value);
    const double raw = local / (nutrient.capacity_half_saturation + local);
    const double at_vessel = nutrient.vessel_value /
        (nutrient.capacity_half_saturation + nutrient.vessel_value);
    const double saturation = at_vessel > 0.0
        ? std::clamp(raw / at_vessel, 0.0, 1.0) : 0.0;
    return 1.0 + (nutrient.maximum_capacity_multiplier - 1.0) * saturation;
}

void StructuredPdeModel3D::rebuild_moving_tumour_front() {
    if (!operators_.nutrient.moving_front) return;
    const auto& nutrient = config_.continuum.nutrient;
    const auto& mode = nutrient.boundary_mode;
    if (mode != "moving_tumor_front_dirichlet_v2" &&
        mode != "moving_tumor_front_and_vessels_dirichlet_v2") return;

    const int nx = config_.continuum.grid.shape[0];
    const int ny = config_.continuum.grid.shape[1];
    const StructuredActiveBounds3D previous_bounds = tumour_mask_bounds_;
    if (tumour_mask_bounds_.valid) {
        for (int y = tumour_mask_bounds_.y0; y < tumour_mask_bounds_.y1; ++y) {
            for (int x = tumour_mask_bounds_.x0; x < tumour_mask_bounds_.x1;
                 ++x) {
                tumour_mask_[index(x, y, 0)] = 0U;
            }
        }
    }
    tumour_mask_bounds_ = {};
    tumour_voxel_count_ = 0U;
    tumour_front_voxel_count_ = 0U;
    if (!population_bounds_.valid) {
        nutrient_update_bounds_ = previous_bounds;
        return;
    }

    const int margin = nutrient.tumor_front_smoothing_radius_voxels + 2;
    const int x0 = std::max(0, population_bounds_.x0 - margin);
    const int y0 = std::max(0, population_bounds_.y0 - margin);
    const int x1 = std::min(nx, population_bounds_.x1 + margin);
    const int y1 = std::min(ny, population_bounds_.y1 + margin);
    const int width = x1 - x0;
    const int height = y1 - y0;
    if (width <= 0 || height <= 0) {
        nutrient_update_bounds_ = previous_bounds;
        return;
    }

    tumour_occupancy_work_.resize(
        static_cast<std::size_t>(width) * height);
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            tumour_occupancy_work_[
                static_cast<std::size_t>(y) * width + x] =
                occupied_fraction(index(x0 + x, y0 + y, 0));
        }
    }
    const auto summary = continuum::build_moving_tumor_front_mask_2d(
        tumour_occupancy_work_, width, height,
        nutrient.tumor_front_smoothing_radius_voxels,
        nutrient.tumor_front_density_threshold,
        tumour_front_workspace_, tumour_local_mask_work_);
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            if (tumour_local_mask_work_[
                    static_cast<std::size_t>(y) * width + x] == 0U) continue;
            tumour_mask_[index(x0 + x, y0 + y, 0)] = 1U;
        }
    }
    tumour_mask_bounds_ = {x0, y0, 0, x1, y1, 1, true};
    nutrient_update_bounds_ = union_bounds(previous_bounds, tumour_mask_bounds_);
    tumour_voxel_count_ = summary.tumour_voxels;
    tumour_front_voxel_count_ = summary.front_voxels;
}

bool StructuredPdeModel3D::nutrient_source(
    int x, int y, int z, std::size_t location) const noexcept {
    if (!operators_.nutrient.transient_resources) return false;
    const auto& mode = config_.continuum.nutrient.boundary_mode;
    if (operators_.nutrient.moving_front &&
        (mode == "moving_tumor_front_dirichlet_v2" ||
         mode == "moving_tumor_front_and_vessels_dirichlet_v2")) {
        // The in-grid host exterior is maintained at Nmax. When tumour reaches
        // the computational box, tumour voxels on the box edge are not sources
        // and the existing ghost-cell rule therefore remains zero-flux.
        const bool use_vessels = mode ==
            "moving_tumor_front_and_vessels_dirichlet_v2";
        return tumour_mask_[location] == 0U ||
            (use_vessels && (!operators_.vascular.density_exclusion ? vessel_[location] > 0.0 : vessel_[location] >= config_.continuum.angiogenesis.exclusion_fraction));
    }
    const auto& grid = config_.continuum.grid;
    const bool planar_edge = x == 0 || x + 1 == grid.shape[0] ||
        y == 0 || y + 1 == grid.shape[1] ||
        (!config_.continuum.base.thin_layer &&
         (z == 0 || z + 1 == grid.shape[2]));
    const bool use_edges = mode != "vessels_dirichlet_v1";
    const bool use_vessels = mode != "planar_edges_dirichlet_v1";
    return (use_edges && planar_edge) ||
        (use_vessels && (!operators_.vascular.density_exclusion ? vessel_[location] > 0.0 : vessel_[location] >= config_.continuum.angiogenesis.exclusion_fraction));
}

void StructuredPdeModel3D::advance_transient_nutrient(double dt) {
    const auto& continuum = config_.continuum;
    const auto& nutrient = continuum.nutrient;
    const int nx = continuum.grid.shape[0];
    const int ny = continuum.grid.shape[1];
    const int nz = continuum.grid.shape[2];
    const int dimensions = continuum.base.thin_layer ? 2 : 3;
    const double inverse_h2 = 1.0 /
        (continuum.grid.spacing_voxels * continuum.grid.spacing_voxels);
    const double mu = nutrient.diffusion_voxels2_per_hour * dt * inverse_h2;
    if (mu > 1.0 / (2.0 * dimensions) + 1.0e-12) {
        throw std::runtime_error(
            "transient nutrient diffusion violates the explicit CFL limit");
    }
    const double source_value = nutrient.vessel_value;
    const double decay = std::exp(-nutrient.decay_per_hour * dt);
    const int workers = std::max(
        1, std::min(continuum.base.threads, available_worker_threads()));
    const bool moving_front = operators_.nutrient.moving_front &&
        (nutrient.boundary_mode == "moving_tumor_front_dirichlet_v2" ||
         nutrient.boundary_mode ==
             "moving_tumor_front_and_vessels_dirichlet_v2");
    const StructuredActiveBounds3D update_bounds = moving_front
        ? nutrient_update_bounds_
        : StructuredActiveBounds3D{0, 0, 0, nx, ny, nz, true};
    if (!update_bounds.valid) {
        ++nutrient_solve_count_;
        return;
    }
    const std::size_t bx = static_cast<std::size_t>(
        update_bounds.x1 - update_bounds.x0);
    const std::size_t by = static_cast<std::size_t>(
        update_bounds.y1 - update_bounds.y0);
    const std::size_t bz = static_cast<std::size_t>(
        update_bounds.z1 - update_bounds.z0);
    const std::size_t update_count = bx * by * bz;
    const auto coordinates = [&](std::size_t offset, int& x, int& y, int& z) {
        x = update_bounds.x0 + static_cast<int>(offset % bx);
        const std::size_t local_yz = offset / bx;
        y = update_bounds.y0 + static_cast<int>(local_yz % by);
        z = update_bounds.z0 + static_cast<int>(local_yz / by);
    };
    deterministic_parallel_for(update_count, workers, [&](std::size_t offset) {
        int x{}, y{}, z{};
        coordinates(offset, x, y, z);
        const std::size_t here = index(x, y, z);
        if (nutrient_source(x, y, z, here)) {
            nutrient_next_[here] = source_value;
            return;
        }

        const double old = nutrient_[here];
        // A boundary that is not selected as a fixed nutrient source is
        // reflecting. The centre-valued ghost cell gives zero normal flux and
        // avoids any domain-exterior access in vessel-only controls.
        double neighbor_sum =
            (x > 0 ? nutrient_[index(x - 1, y, z)] : old) +
            (x + 1 < nx ? nutrient_[index(x + 1, y, z)] : old) +
            (y > 0 ? nutrient_[index(x, y - 1, z)] : old) +
            (y + 1 < ny ? nutrient_[index(x, y + 1, z)] : old);
        if (!continuum.base.thin_layer) {
            neighbor_sum +=
                (z > 0 ? nutrient_[index(x, y, z - 1)] : old) +
                (z + 1 < nz ? nutrient_[index(x, y, z + 1)] : old);
        }
        const double diffused = std::clamp(
            old + mu * (neighbor_sum - 2.0 * dimensions * old),
            0.0, source_value) * decay;

        // Schema v3 consumption is per biological cell, independent of
        // phenotype and footprint. Solve the local Michaelis-Menten sink
        // implicitly so consumption cannot drive the field negative.
        double consumers = r_normal_[0][here] +
            r_active(StructuredStage3D::small, here) +
            r_normal_[1][here] +
            r_active(StructuredStage3D::large, here) +
            K_[0][here] + K_[1][here];
        if (!external_r_.empty()) consumers += external_consumers(here);
        const double demand =
            nutrient.K_consumption_rate_per_hour * consumers;
        const double half = nutrient.K_consumption_half_saturation;
        const double b = half + dt * demand - diffused;
        const double discriminant = std::max(
            0.0, b * b + 4.0 * half * diffused);
        const double consumed = 0.5 * (-b + std::sqrt(discriminant));
        nutrient_next_[here] = operators_.nutrient.footprint_consumers
            ? continuum::resource_after_uptake(diffused, consumers, nutrient.K_consumption_rate_per_hour,
                half, dt, source_value)
            : std::clamp(consumed, 0.0, source_value);
        if(angiogenesis_) nutrient_next_[here]=source_value-(source_value-nutrient_next_[here])*
            std::exp(-dt*continuum.angiogenesis.perfusion_exchange_per_hour*vessel_[here]);
    });
    nutrient_.swap(nutrient_next_);
    if (moving_front) {
        // The explicit stencil above must read the previous value of a voxel
        // that has just become exterior. After the swap, synchronize those
        // Dirichlet sites in the inactive buffer so the solve box may shrink
        // again on the next step without leaving stale nutrient behind.
        deterministic_parallel_for(
            update_count, workers, [&](std::size_t offset) {
                int x{}, y{}, z{};
                coordinates(offset, x, y, z);
                const std::size_t here = index(x, y, z);
                if (nutrient_source(x, y, z, here)) {
                    nutrient_next_[here] = source_value;
                }
            });
    }
    ++nutrient_solve_count_;
    validate_resources();
}

void StructuredPdeModel3D::solve_nutrient() {
    const auto& continuum = config_.continuum;
    if (operators_.nutrient.transient_resources) {
        const int nx = continuum.grid.shape[0];
        const int ny = continuum.grid.shape[1];
        for (std::size_t here = 0; here < voxel_count_; ++here) {
            const int x = static_cast<int>(
                here % static_cast<std::size_t>(nx));
            const std::size_t yz = here / static_cast<std::size_t>(nx);
            const int y = static_cast<int>(
                yz % static_cast<std::size_t>(ny));
            const int z = static_cast<int>(
                yz / static_cast<std::size_t>(ny));
            if (nutrient_source(x, y, z, here)) {
                nutrient_[here] = continuum.nutrient.vessel_value;
            }
        }
        ++nutrient_solve_count_;
        validate_resources();
        if (operators_.nutrient.moving_front) nutrient_next_ = nutrient_;
        return;
    }
    const int nx = continuum.grid.shape[0];
    const int ny = continuum.grid.shape[1];
    const int nz = continuum.grid.shape[2];
    const int dimensions = continuum.base.thin_layer ? 2 : 3;
    const double diffusion = continuum.nutrient.diffusion_voxels2_per_hour /
        (continuum.grid.spacing_voxels * continuum.grid.spacing_voxels);
    const double laplacian_diagonal = 2.0 * dimensions * diffusion;
    const int workers = std::max(
        1, std::min(continuum.base.threads, available_worker_threads()));
    for (int iteration = 0; iteration < continuum.nutrient.solver_iterations;
         ++iteration) {
        deterministic_parallel_for(voxel_count_, workers, [&](std::size_t here) {
            const int x = static_cast<int>(here % static_cast<std::size_t>(nx));
            const std::size_t yz = here / static_cast<std::size_t>(nx);
            const int y = static_cast<int>(yz % static_cast<std::size_t>(ny));
            const int z = static_cast<int>(yz / static_cast<std::size_t>(ny));
            double neighbor_sum = 0.0;
            if (x > 0) neighbor_sum += nutrient_[index(x - 1, y, z)];
            if (x + 1 < nx) neighbor_sum += nutrient_[index(x + 1, y, z)];
            if (y > 0) neighbor_sum += nutrient_[index(x, y - 1, z)];
            if (y + 1 < ny) neighbor_sum += nutrient_[index(x, y + 1, z)];
            if (!continuum.base.thin_layer) {
                if (z > 0) neighbor_sum += nutrient_[index(x, y, z - 1)];
                if (z + 1 < nz) neighbor_sum += nutrient_[index(x, y, z + 1)];
            }
            const double old = nutrient_[here];
            const bool per_cell = continuum.nutrient.consumption_model ==
                "per_cell_ratio_v2";
            const double large_consumption_weight =
                per_cell ? 1.0 : large_cell_volume_;
            const double r_consumers = r_normal_[0][here] +
                r_active(StructuredStage3D::small, here) +
                large_consumption_weight * (r_normal_[1][here] +
                    r_active(StructuredStage3D::large, here));
            const double K_consumers =
                K_[0][here] + large_consumption_weight * K_[1][here];
            const double r_sink =
                continuum.nutrient.r_consumption_rate_per_hour *
                r_consumers /
                (continuum.nutrient.r_consumption_half_saturation + old);
            const double K_sink =
                continuum.nutrient.K_consumption_rate_per_hour *
                K_consumers /
                (continuum.nutrient.K_consumption_half_saturation + old);
            const double exchange = continuum.nutrient.vessel_exchange_per_hour *
                std::clamp(vessel_[here], 0.0, 1.0);
            const double denominator = laplacian_diagonal +
                continuum.nutrient.decay_per_hour + exchange + r_sink + K_sink;
            const double candidate = denominator > 0.0
                ? (diffusion * neighbor_sum +
                   exchange * continuum.nutrient.vessel_value) / denominator
                : 0.0;
            nutrient_next_[here] = std::clamp(
                old + continuum.nutrient.relaxation * (candidate - old),
                0.0, continuum.nutrient.vessel_value);
        });
        nutrient_.swap(nutrient_next_);
    }
    ++nutrient_solve_count_;
    validate_resources();
}

}  // namespace atcg3d::structured_pde
