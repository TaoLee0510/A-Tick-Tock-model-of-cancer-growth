#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <vector>

#include "config/model_config.hpp"
#include "engine/simulation.hpp"

namespace atcg3d {

struct CheckpointData3D {
    std::vector<CellInit> cells;
    std::vector<Slot> cell_slots;
    std::size_t cell_slot_count{};
    std::vector<Slot> cell_free_slots;
    CellUid next_uid{};
    SimulationClock3D clock;
    SimulationStats3D stats;
    std::vector<LineageEdge> lineage;
    VasculatureState3D vasculature;
    std::uint64_t state_checksum{};
};

void write_hdf5_checkpoint(const std::filesystem::path& path,
                           const Simulation3D& simulation);
CheckpointData3D read_hdf5_checkpoint(const std::filesystem::path& path,
                                      const Model3DConfig& expected_config);

}  // namespace atcg3d
