#include <cassert>
#include <filesystem>
#include <fstream>
#include <stdexcept>

#include <H5Cpp.h>

#include "config/model_config.hpp"
#include "engine/simulation.hpp"
#include "io/checkpoint_hdf5.hpp"

int main() {
    using namespace atcg3d;
    H5::Exception::dontPrint();
    Model3DConfig short_config;
    short_config.output_enabled = false;
    short_config.initial_r_cells = 4;
    short_config.initial_K_cells = 4;
    short_config.initial_radius = 7;
    short_config.end_time_hours = 5.0;
    short_config.max_events = 10000;
    short_config.density_block_edge = 2;

    Simulation3D original(short_config);
    original.run();
    const auto expected_checksum = original.state_checksum();
    const std::filesystem::path path = std::filesystem::temp_directory_path() / "atcg3d_checkpoint_test.h5";
    std::filesystem::remove(path);
    write_hdf5_checkpoint(path, original);
    const CheckpointData3D data = read_hdf5_checkpoint(path, short_config);
    assert(data.state_checksum == expected_checksum);

    Simulation3D restored(short_config);
    restored.restore(data.cells, data.next_uid, data.clock, data.stats, data.lineage);
    assert(restored.state_checksum() == expected_checksum);
    assert(restored.clock().time_hours == original.clock().time_hours);
    assert(restored.next_uid() == original.next_uid());

    // Stopping at a checkpoint and resuming must reproduce an uninterrupted
    // run. End time, event cap, thread count, and output policy may change;
    // biological and numerical parameters may not.
    Model3DConfig long_config = short_config;
    long_config.end_time_hours = 10.0;
    long_config.threads = 2;
    Simulation3D continuous(long_config);
    continuous.run();

    const CheckpointData3D resume_data = read_hdf5_checkpoint(path, long_config);
    Simulation3D resumed(long_config);
    resumed.restore(resume_data.cells, resume_data.next_uid, resume_data.clock,
                    resume_data.stats, resume_data.lineage);
    resumed.run();
    assert(resumed.state_checksum() == continuous.state_checksum());
    assert(resumed.clock().completed_events == continuous.clock().completed_events);
    assert(resumed.stats().migration_commits == continuous.stats().migration_commits);

    Model3DConfig incompatible = long_config;
    incompatible.continue_probability = 0.5;
    bool mismatch_rejected = false;
    try {
        (void)read_hdf5_checkpoint(path, incompatible);
    } catch (const std::runtime_error&) {
        mismatch_rejected = true;
    }
    assert(mismatch_rejected);

    const std::filesystem::path corrupt =
        std::filesystem::temp_directory_path() / "atcg3d_checkpoint_corrupt_test.h5";
    std::filesystem::remove(corrupt);
    {
        std::ofstream stream(corrupt, std::ios::binary);
        stream << "not an HDF5 checkpoint";
    }
    bool corrupt_rejected = false;
    try {
        (void)read_hdf5_checkpoint(corrupt, long_config);
    } catch (const std::runtime_error&) {
        corrupt_rejected = true;
    }
    assert(corrupt_rejected);
    std::filesystem::remove(path);
    std::filesystem::remove(corrupt);
}
