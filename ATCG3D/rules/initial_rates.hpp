#pragma once

#include <cstdint>

#include "config/model_config.hpp"
#include "core/types.hpp"

namespace atcg3d {

struct InitialCellRates3D {
    double inherent_growth_rate{};
    double migration_rate{};
};

// These functions are pure and stateless. Their output is keyed only by the
// effective configuration, seed, immutable cell UID, and cell type; they never
// read or advance the simulation event RNG sequence.
double sample_initial_growth_rate(const Model3DConfig& config,
                                  CellType type,
                                  std::uint64_t seed,
                                  CellUid uid);

double sample_initial_migration_rate(const Model3DConfig& config,
                                     CellType type,
                                     std::uint64_t seed,
                                     CellUid uid);

InitialCellRates3D sample_initial_cell_rates(const Model3DConfig& config,
                                             CellType type,
                                             std::uint64_t seed,
                                             CellUid uid);

// Pure event-keyed Beta draw used by finite migration-state transitions. The
// result is in [0,1] and does not advance any mutable RNG engine.
double sample_stateless_beta_fraction(double alpha,
                                      double beta,
                                      std::uint64_t seed,
                                      CellUid uid,
                                      std::uint64_t event_domain,
                                      std::uint64_t event_sequence);

}  // namespace atcg3d
