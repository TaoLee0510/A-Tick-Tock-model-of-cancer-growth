#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include <H5Cpp.h>
#include <hdf5.h>

#include "config/model_config.hpp"
#include "engine/simulation.hpp"
#include "io/checkpoint_hdf5.hpp"

namespace {

template <class Callback>
void expect_rejected(Callback&& callback) {
    bool rejected = false;
    try {
        callback();
    } catch (const std::runtime_error&) {
        rejected = true;
    }
    assert(rejected);
}

atcg3d::Model3DConfig cell_only_config(double end_time) {
    atcg3d::Model3DConfig config;
    config.output_enabled = false;
    config.initial_r_cells = 4;
    config.initial_K_cells = 4;
    config.initial_radius = 7;
    config.end_time_hours = end_time;
    config.max_events = 10000;
    config.density_block_edge = 2;
    config.validate();
    return config;
}

atcg3d::Model3DConfig vascular_config(double end_time) {
    atcg3d::Model3DConfig config;
    config.output_enabled = false;
    config.initial_r_cells = 20;
    config.initial_K_cells = 20;
    config.initial_radius = 8;
    config.initial_shell_inner_radius = 6;
    config.initial_shell_thickness = 2;
    config.initial_inner_small_radius = 6;
    config.initial_K_migration_rate = 0.0;
    config.end_time_hours = end_time;
    config.max_events = 100000;

    auto& vessels = config.angiogenesis;
    vessels.enabled = true;
    vessels.trigger_activation_volume_voxels3 = 1.0;
    vessels.trigger_deactivation_volume_voxels3 = 0.0;
    vessels.seed_rate_sites_per_30_days = 720000.0;
    vessels.seed_rate_sites_per_hour = vessels.seed_rate_sites_per_30_days / 720.0;
    vessels.max_total_roots = 1;
    vessels.max_active_tips = 2;
    vessels.surface_min_separation_voxels = 0;
    vessels.diameter_voxels = 1.0;
    vessels.inward_speed_voxels_per_hour = 10.0;
    vessels.outward_speed_voxels_per_hour = 10.0;
    vessels.inward_max_length_voxels = 12;
    vessels.outward_max_length_voxels = 12;
    vessels.outward_external_connection_distance_voxels = 1.0;
    vessels.influence_activation = "immediate";
    vessels.influence_cutoff_radius_voxels = 3.0;
    vessels.influence_max_relief_fraction = 0.5;
    vessels.lesion_block_edge = 8;
    vessels.lesion_core_activation_occupied_fraction = 1.0 / 512.0;
    vessels.lesion_core_deactivation_occupied_fraction = 0.0;
    vessels.lesion_minimum_cells_per_core_block = 1;
    vessels.lesion_minimum_biological_volume_per_core_block = 0.0;
    vessels.trigger_minimum_core_blocks = 1;
    config.validate();
    return config;
}

atcg3d::CellInit lesion_cell(atcg3d::Vec3i anchor,
                             atcg3d::CellUid uid) {
    atcg3d::CellInit cell;
    cell.anchor = anchor;
    cell.uid = uid;
    cell.clone_id = static_cast<std::uint32_t>(uid);
    cell.type = atcg3d::CellType::r;
    cell.stage = atcg3d::CellStage::small;
    cell.inherent_growth_rate = 1.0F;
    cell.density_growth_rate = 1.0F;
    cell.migration_rate = 1.0F;
    cell.normal_migration_rate = 0.25F;
    cell.division_work_remaining = 100.0F;
    cell.next_migration_time = 100.0;
    cell.next_division_time = 100.0;
    cell.event_sequence = 11;
    return cell;
}

std::vector<atcg3d::CellInit> lesion_cluster(
    atcg3d::Vec3i origin, atcg3d::CellUid& next_uid) {
    std::vector<atcg3d::CellInit> cells;
    cells.reserve(8);
    for (int z = 0; z < 2; ++z) {
        for (int y = 0; y < 2; ++y) {
            for (int x = 0; x < 2; ++x) {
                cells.push_back(lesion_cell(
                    origin + atcg3d::Vec3i{x, y, z}, next_uid++));
            }
        }
    }
    return cells;
}

atcg3d::Model3DConfig two_lesion_config() {
    atcg3d::Model3DConfig config = vascular_config(2.0);
    config.initial_r_cells = 0;
    config.initial_K_cells = 0;
    config.max_events = 2;
    auto& vessels = config.angiogenesis;
    vessels.lesion_block_edge = 1;
    vessels.lesion_core_activation_occupied_fraction = 1.0;
    vessels.lesion_core_deactivation_occupied_fraction = 1.0;
    vessels.lesion_halo_blocks = 0;
    vessels.trigger_activation_volume_voxels3 = 4.0;
    vessels.seed_rate_sites_per_30_days = 720000000.0;
    vessels.seed_rate_sites_per_hour =
        vessels.seed_rate_sites_per_30_days / 720.0;
    vessels.max_total_roots = 2;
    vessels.max_roots_per_lesion = 1;
    vessels.max_active_tips = 4;
    vessels.max_active_tips_per_lesion = 2;
    config.validate();
    return config;
}

void corrupt_schema_version(const std::filesystem::path& source,
                            const std::filesystem::path& destination,
                            std::uint32_t version) {
    std::filesystem::copy_file(source, destination,
                               std::filesystem::copy_options::overwrite_existing);
    H5::H5File file(destination.string(), H5F_ACC_RDWR);
    H5::Group meta = file.openGroup("/meta");
    H5::Attribute attribute = meta.openAttribute("schema_version");
    attribute.write(H5::PredType::NATIVE_UINT32, &version);
    file.flush(H5F_SCOPE_GLOBAL);
}

void corrupt_cell_column_length(const std::filesystem::path& source,
                                const std::filesystem::path& destination) {
    std::filesystem::copy_file(source, destination,
                               std::filesystem::copy_options::overwrite_existing);
    H5::H5File file(destination.string(), H5F_ACC_RDWR);
    H5::Group cells = file.openGroup("/cells");
    const std::size_t count = cells.openDataSet("uid").getSpace().getSimpleExtentNpoints();
    assert(count > 1);
    if (H5Ldelete(cells.getId(), "x", H5P_DEFAULT) < 0) {
        throw std::runtime_error("failed to corrupt checkpoint test dataset");
    }
    const hsize_t dimensions[1] = {static_cast<hsize_t>(count - 1)};
    H5::DataSpace space(1, dimensions);
    H5::DataSet x = cells.createDataSet("x", H5::PredType::STD_I32LE, space);
    const std::vector<std::int32_t> values(count - 1, 0);
    x.write(values.data(), H5::PredType::NATIVE_INT32);
    file.flush(H5F_SCOPE_GLOBAL);
}

void corrupt_division_work(const std::filesystem::path& source,
                           const std::filesystem::path& destination) {
    std::filesystem::copy_file(source, destination,
                               std::filesystem::copy_options::overwrite_existing);
    H5::H5File file(destination.string(), H5F_ACC_RDWR);
    H5::DataSet dataset = file.openDataSet("/cells/division_work_remaining");
    const std::size_t count = dataset.getSpace().getSimpleExtentNpoints();
    assert(count > 0);
    std::vector<float> values(count);
    dataset.read(values.data(), H5::PredType::NATIVE_FLOAT);
    // Keep the value finite, nonnegative, and otherwise schema-valid. The
    // reader must reject it because reconstructed state no longer matches the
    // stored checksum, rather than because a scalar validator catches it.
    values.front() += 1.0F;
    dataset.write(values.data(), H5::PredType::NATIVE_FLOAT);
    file.flush(H5F_SCOPE_GLOBAL);
}

void corrupt_migration_activation_pair(const std::filesystem::path& source,
                                       const std::filesystem::path& destination) {
    std::filesystem::copy_file(source, destination,
                               std::filesystem::copy_options::overwrite_existing);
    H5::H5File file(destination.string(), H5F_ACC_RDWR);
    H5::DataSet flags_data = file.openDataSet("/cells/flags");
    H5::DataSet end_data =
        file.openDataSet("/cells/migration_activation_end_time");
    const std::size_t count = flags_data.getSpace().getSimpleExtentNpoints();
    assert(count > 0 &&
           static_cast<std::size_t>(
               end_data.getSpace().getSimpleExtentNpoints()) == count);
    std::vector<std::uint8_t> flags(count);
    std::vector<double> ends(count);
    flags_data.read(flags.data(), H5::PredType::NATIVE_UINT8);
    end_data.read(ends.data(), H5::PredType::NATIVE_DOUBLE);
    flags.front() &= static_cast<std::uint8_t>(~atcg3d::kMigrationActive);
    ends.front() = 1.0;
    flags_data.write(flags.data(), H5::PredType::NATIVE_UINT8);
    end_data.write(ends.data(), H5::PredType::NATIVE_DOUBLE);
    file.flush(H5F_SCOPE_GLOBAL);
}

void corrupt_lesion_process_aggregate(
    const std::filesystem::path& source,
    const std::filesystem::path& destination) {
    std::filesystem::copy_file(source, destination,
                               std::filesystem::copy_options::overwrite_existing);
    H5::H5File file(destination.string(), H5F_ACC_RDWR);
    for (const char* name : {"attempted_events", "rejected_events"}) {
        H5::DataSet dataset = file.openDataSet(
            std::string("/vasculature/lesions/processes/") + name);
        const std::size_t count = dataset.getSpace().getSimpleExtentNpoints();
        assert(count > 0);
        std::vector<std::uint64_t> values(count);
        dataset.read(values.data(), H5::PredType::NATIVE_UINT64);
        ++values.front();
        dataset.write(values.data(), H5::PredType::NATIVE_UINT64);
    }
    file.flush(H5F_SCOPE_GLOBAL);
}

void corrupt_child_node_source_lesion(
    const std::filesystem::path& source,
    const std::filesystem::path& destination) {
    std::filesystem::copy_file(source, destination,
                               std::filesystem::copy_options::overwrite_existing);
    H5::H5File file(destination.string(), H5F_ACC_RDWR);
    H5::DataSet parent_dataset =
        file.openDataSet("/vasculature/nodes/parent_node_slot");
    H5::DataSet source_dataset =
        file.openDataSet("/vasculature/nodes/source_lesion_id");
    const std::size_t count = parent_dataset.getSpace().getSimpleExtentNpoints();
    assert(count > 1 &&
           source_dataset.getSpace().getSimpleExtentNpoints() == count);
    std::vector<std::uint32_t> parents(count);
    std::vector<std::uint64_t> sources(count);
    parent_dataset.read(parents.data(), H5::PredType::NATIVE_UINT32);
    source_dataset.read(sources.data(), H5::PredType::NATIVE_UINT64);
    const auto child = std::find_if(
        parents.begin(), parents.end(), [](std::uint32_t parent) {
            return parent != atcg3d::kEmptyVesselNodeSlot;
        });
    assert(child != parents.end());
    const std::size_t index = static_cast<std::size_t>(child - parents.begin());
    ++sources[index];
    source_dataset.write(sources.data(), H5::PredType::NATIVE_UINT64);
    file.flush(H5F_SCOPE_GLOBAL);
}

void corrupt_tip_source_lesion(const std::filesystem::path& source,
                               const std::filesystem::path& destination) {
    std::filesystem::copy_file(source, destination,
                               std::filesystem::copy_options::overwrite_existing);
    H5::H5File file(destination.string(), H5F_ACC_RDWR);
    H5::DataSet dataset =
        file.openDataSet("/vasculature/tips/source_lesion_id");
    const std::size_t count = dataset.getSpace().getSimpleExtentNpoints();
    assert(count > 0);
    std::vector<std::uint64_t> sources(count);
    dataset.read(sources.data(), H5::PredType::NATIVE_UINT64);
    ++sources.front();
    dataset.write(sources.data(), H5::PredType::NATIVE_UINT64);
    file.flush(H5F_SCOPE_GLOBAL);
}

void corrupt_next_lesion_refresh_time(
    const std::filesystem::path& source,
    const std::filesystem::path& destination) {
    std::filesystem::copy_file(source, destination,
                               std::filesystem::copy_options::overwrite_existing);
    H5::H5File file(destination.string(), H5F_ACC_RDWR);
    H5::Group lesions = file.openGroup("/vasculature/lesions");
    H5::Attribute attribute =
        lesions.openAttribute("next_refresh_time_hours");
    const double invalid = -1.0;
    attribute.write(H5::PredType::NATIVE_DOUBLE, &invalid);
    file.flush(H5F_SCOPE_GLOBAL);
}

void corrupt_lesion_dirty_block_exists(
    const std::filesystem::path& source,
    const std::filesystem::path& destination) {
    std::filesystem::copy_file(source, destination,
                               std::filesystem::copy_options::overwrite_existing);
    H5::H5File file(destination.string(), H5F_ACC_RDWR);
    H5::DataSet dataset =
        file.openDataSet("/vasculature/lesions/dirty_blocks/exists");
    const std::size_t count = dataset.getSpace().getSimpleExtentNpoints();
    assert(count > 0);
    std::vector<std::uint8_t> values(count);
    dataset.read(values.data(), H5::PredType::NATIVE_UINT8);
    values.front() = 2;
    dataset.write(values.data(), H5::PredType::NATIVE_UINT8);
    file.flush(H5F_SCOPE_GLOBAL);
}

void remove_effective_config_provenance(
    const std::filesystem::path& source,
    const std::filesystem::path& destination) {
    std::filesystem::copy_file(source, destination,
                               std::filesystem::copy_options::overwrite_existing);
    H5::H5File file(destination.string(), H5F_ACC_RDWR);
    H5::Group meta = file.openGroup("/meta");
    if (H5Adelete(meta.getId(), "effective_config_json") < 0) {
        throw std::runtime_error(
            "failed to remove checkpoint effective-config provenance");
    }
    file.flush(H5F_SCOPE_GLOBAL);
}

void corrupt_effective_config_provenance(
    const std::filesystem::path& source,
    const std::filesystem::path& destination) {
    std::filesystem::copy_file(source, destination,
                               std::filesystem::copy_options::overwrite_existing);
    H5::H5File file(destination.string(), H5F_ACC_RDWR);
    H5::Group meta = file.openGroup("/meta");
    if (H5Adelete(meta.getId(), "effective_config_json") < 0) {
        throw std::runtime_error(
            "failed to replace checkpoint effective-config provenance");
    }
    H5::StrType type(H5::PredType::C_S1, H5T_VARIABLE);
    type.setCset(H5T_CSET_UTF8);
    H5::DataSpace scalar(H5S_SCALAR);
    H5::Attribute attribute =
        meta.createAttribute("effective_config_json", type, scalar);
    const std::string invalid = "{not-json";
    attribute.write(type, invalid);
    file.flush(H5F_SCOPE_GLOBAL);
}

void remove_source_ownership_table(
    const std::filesystem::path& source,
    const std::filesystem::path& destination) {
    std::filesystem::copy_file(source, destination,
                               std::filesystem::copy_options::overwrite_existing);
    H5::H5File file(destination.string(), H5F_ACC_RDWR);
    H5::Group lesions = file.openGroup("/vasculature/lesions");
    if (H5Ldelete(lesions.getId(), "source_ownership", H5P_DEFAULT) < 0) {
        throw std::runtime_error(
            "failed to remove checkpoint lesion-source ownership table");
    }
    file.flush(H5F_SCOPE_GLOBAL);
}

void corrupt_source_ownership_target(
    const std::filesystem::path& source,
    const std::filesystem::path& destination) {
    std::filesystem::copy_file(source, destination,
                               std::filesystem::copy_options::overwrite_existing);
    H5::H5File file(destination.string(), H5F_ACC_RDWR);
    H5::DataSet dataset = file.openDataSet(
        "/vasculature/lesions/source_ownership/current_lesion_id");
    const std::size_t count = dataset.getSpace().getSimpleExtentNpoints();
    assert(count > 0);
    std::vector<std::uint64_t> targets(count);
    dataset.read(targets.data(), H5::PredType::NATIVE_UINT64);
    targets.front() = std::numeric_limits<std::uint64_t>::max();
    dataset.write(targets.data(), H5::PredType::NATIVE_UINT64);
    file.flush(H5F_SCOPE_GLOBAL);
}

}  // namespace

int main() {
    using namespace atcg3d;
    H5::Exception::dontPrint();
    const std::filesystem::path directory =
        std::filesystem::temp_directory_path() / "atcg3d_checkpoint_v3_test";
    std::filesystem::create_directories(directory);
    const std::filesystem::path cell_path = directory / "cells.h5";
    const std::filesystem::path active_path = directory / "active_migration.h5";
    const std::filesystem::path non_float_clock_path =
        directory / "non_float_clock.h5";
    const std::filesystem::path holes_path = directory / "holes.h5";
    const std::filesystem::path vascular_path = directory / "vascular.h5";
    const std::filesystem::path two_lesion_path =
        directory / "two_lesions.h5";
    const std::filesystem::path historical_source_path =
        directory / "historical_lesion_source.h5";
    const std::filesystem::path corrupt_path = directory / "corrupt.h5";
    const std::filesystem::path version_path = directory / "wrong_version.h5";
    const std::filesystem::path columns_path = directory / "wrong_columns.h5";
    const std::filesystem::path work_path = directory / "wrong_work.h5";
    const std::filesystem::path activation_path =
        directory / "wrong_migration_activation.h5";
    const std::filesystem::path aggregate_path =
        directory / "wrong_lesion_aggregate.h5";
    const std::filesystem::path node_source_path =
        directory / "wrong_node_source_lesion.h5";
    const std::filesystem::path tip_source_path =
        directory / "wrong_tip_source_lesion.h5";
    const std::filesystem::path refresh_time_path =
        directory / "wrong_lesion_refresh_time.h5";
    const std::filesystem::path dirty_block_path =
        directory / "wrong_lesion_dirty_block.h5";
    const std::filesystem::path missing_effective_config_path =
        directory / "missing_effective_config.h5";
    const std::filesystem::path corrupt_effective_config_path =
        directory / "corrupt_effective_config.h5";
    const std::filesystem::path missing_source_ownership_path =
        directory / "missing_source_ownership.h5";
    const std::filesystem::path corrupt_source_ownership_path =
        directory / "corrupt_source_ownership.h5";
    std::filesystem::remove(cell_path);
    std::filesystem::remove(active_path);
    std::filesystem::remove(non_float_clock_path);
    std::filesystem::remove(holes_path);
    std::filesystem::remove(vascular_path);
    std::filesystem::remove(two_lesion_path);
    std::filesystem::remove(historical_source_path);
    std::filesystem::remove(corrupt_path);
    std::filesystem::remove(version_path);
    std::filesystem::remove(columns_path);
    std::filesystem::remove(work_path);
    std::filesystem::remove(activation_path);
    std::filesystem::remove(aggregate_path);
    std::filesystem::remove(node_source_path);
    std::filesystem::remove(tip_source_path);
    std::filesystem::remove(refresh_time_path);
    std::filesystem::remove(dirty_block_path);
    std::filesystem::remove(missing_effective_config_path);
    std::filesystem::remove(corrupt_effective_config_path);
    std::filesystem::remove(missing_source_ownership_path);
    std::filesystem::remove(corrupt_source_ownership_path);

    // Schema v3 stores the three event-kind generations independently, even
    // when angiogenesis is disabled and the vascular tables are empty.
    Model3DConfig short_config = cell_only_config(5.0);
    Simulation3D original(short_config);
    original.run();
    const auto expected_checksum = original.state_checksum();
    write_hdf5_checkpoint(cell_path, original);
    assert(original.state_checksum() == expected_checksum);
    assert(!std::filesystem::exists(cell_path.string() + ".tmp"));
    const CheckpointData3D data = read_hdf5_checkpoint(cell_path, short_config);
    assert(data.state_checksum == expected_checksum);
    assert(data.vasculature.nodes.empty());
    assert(data.vasculature.tips.empty());
    assert(data.vasculature.next_vessel_id == 1);
    assert(data.vasculature.lesions.next_refresh_time_hours == 0.0);
    assert(data.vasculature.lesions.refresh_schedule_generation == 0);
    assert(!data.cells.empty());
    assert(data.cells.front().migration_schedule_generation != 0);
    assert(data.cells.front().division_schedule_generation != 0);
    assert(data.cells.front().death_schedule_generation != 0);
    assert(data.cells.front().division_work_remaining > 0.0F);
    assert(data.cells.front().normal_migration_rate > 0.0F);
    {
        H5::H5File file(cell_path.string(), H5F_ACC_RDONLY);
        H5::Group meta = file.openGroup("/meta");
        std::uint32_t schema{};
        meta.openAttribute("schema_version")
            .read(H5::PredType::NATIVE_UINT32, &schema);
        assert(schema == 3);
        assert(static_cast<std::size_t>(
                   file.openDataSet("/cells/normal_migration_rate")
                       .getSpace().getSimpleExtentNpoints()) == data.cells.size());
        assert(static_cast<std::size_t>(
                   file.openDataSet("/cells/migration_activation_end_time")
                       .getSpace().getSimpleExtentNpoints()) == data.cells.size());
    }

    Simulation3D restored(short_config);
    restored.restore(data.cells, data.next_uid, data.clock, data.stats, data.lineage,
                     data.vasculature, data.cell_slot_count, data.cell_slots,
                     data.cell_free_slots);
    assert(restored.state_checksum() == expected_checksum);
    assert(restored.clock().time_hours == original.clock().time_hours);
    assert(restored.next_uid() == original.next_uid());
    const Slot checksum_slot = restored.cells().alive_slots().front();
    restored.cells().set_division_work_remaining(
        checksum_slot,
        restored.cells().division_work_remaining(checksum_slot) + 1.0F);
    assert(restored.state_checksum() != expected_checksum);
    Simulation3D migration_state_changed(short_config);
    migration_state_changed.restore(
        data.cells, data.next_uid, data.clock, data.stats, data.lineage,
        data.vasculature, data.cell_slot_count, data.cell_slots,
        data.cell_free_slots);
    const Slot migration_slot =
        migration_state_changed.cells().alive_slots().front();
    migration_state_changed.cells().set_normal_migration_rate(
        migration_slot,
        migration_state_changed.cells().normal_migration_rate(migration_slot) +
            0.01F);
    assert(migration_state_changed.state_checksum() != expected_checksum);

    // A checkpoint may be taken in the middle of a finite active interval.
    // Both rates, the end event time, and the shared migration generation must
    // round-trip exactly.
    CellInit active_cell;
    active_cell.uid = 9000;
    active_cell.flags = kDirtyDensity | kMigrationActive;
    active_cell.migration_rate = 1.5F;
    active_cell.normal_migration_rate = 0.2F;
    active_cell.migration_activation_end_time = 4.0;
    active_cell.density_growth_rate = 1.0F;
    active_cell.division_work_remaining = 9.0F;
    active_cell.next_migration_time = 2.0;
    active_cell.next_division_time = 10.0;
    active_cell.last_update_time = 1.0;
    active_cell.migration_schedule_generation = 7;
    active_cell.division_schedule_generation = 8;
    active_cell.death_schedule_generation = 9;
    SimulationClock3D active_clock;
    active_clock.time_hours = 1.0;
    Simulation3D active_simulation(short_config);
    active_simulation.restore({active_cell}, 9001, active_clock, {}, {});
    const std::uint64_t active_checksum = active_simulation.state_checksum();
    write_hdf5_checkpoint(active_path, active_simulation);
    const CheckpointData3D active_data =
        read_hdf5_checkpoint(active_path, short_config);
    assert(active_data.cells.size() == 1);
    assert(active_data.cells.front().migration_rate == 1.5F);
    assert(active_data.cells.front().normal_migration_rate == 0.2F);
    assert(active_data.cells.front().migration_activation_end_time == 4.0);
    Simulation3D active_restored(short_config);
    active_restored.restore(
        active_data.cells, active_data.next_uid, active_data.clock,
        active_data.stats, active_data.lineage, active_data.vasculature,
        active_data.cell_slot_count, active_data.cell_slots,
        active_data.cell_free_slots);
    assert(active_restored.state_checksum() == active_checksum);

    // A non-float-exact global clock can occur when a vessel event refreshes a
    // cell. Nearest-rounded last_update may then be slightly after the clock;
    // restore accepts that bounded representation error, while future event
    // times round upward and cannot be silently dropped as past events.
    const float lower_clock = 1.0F;
    const float upper_clock = std::nextafter(
        lower_clock, std::numeric_limits<float>::infinity());
    const double non_float_clock = static_cast<double>(lower_clock) +
        0.75 * (static_cast<double>(upper_clock) - lower_clock);
    CellInit non_float_cell;
    non_float_cell.uid = 9100;
    non_float_cell.density_growth_rate = 1.0F;
    non_float_cell.division_work_remaining = 20.0F;
    non_float_cell.last_update_time = non_float_clock;
    non_float_cell.next_division_time = non_float_clock + 20.0;
    non_float_cell.migration_schedule_generation = 1;
    non_float_cell.division_schedule_generation = 1;
    non_float_cell.death_schedule_generation = 1;
    SimulationClock3D non_float_state_clock;
    non_float_state_clock.time_hours = non_float_clock;
    Simulation3D non_float_simulation(short_config);
    non_float_simulation.restore(
        {non_float_cell}, 9101, non_float_state_clock, {}, {});
    assert(non_float_simulation.cells().last_update_time(0) ==
           static_cast<double>(upper_clock));
    assert(non_float_simulation.cells().last_update_time(0) > non_float_clock);
    assert(non_float_simulation.cells().next_division_time(0) >=
           non_float_cell.next_division_time);
    const std::uint64_t non_float_checksum =
        non_float_simulation.state_checksum();
    write_hdf5_checkpoint(non_float_clock_path, non_float_simulation);
    const CheckpointData3D non_float_data = read_hdf5_checkpoint(
        non_float_clock_path, short_config);
    Simulation3D non_float_restored(short_config);
    non_float_restored.restore(
        non_float_data.cells, non_float_data.next_uid, non_float_data.clock,
        non_float_data.stats, non_float_data.lineage,
        non_float_data.vasculature, non_float_data.cell_slot_count,
        non_float_data.cell_slots, non_float_data.cell_free_slots);
    assert(non_float_restored.state_checksum() == non_float_checksum);

    const auto checksum_for = [&](SimulationStats3D stats,
                                  std::vector<LineageEdge> lineage) {
        Simulation3D candidate(short_config);
        candidate.restore(data.cells, data.next_uid, data.clock, stats,
                          std::move(lineage), data.vasculature,
                          data.cell_slot_count, data.cell_slots,
                          data.cell_free_slots);
        return candidate.state_checksum();
    };
    using StatCounter = std::uint64_t SimulationStats3D::*;
    constexpr std::array<StatCounter, 12> stat_counters{
        &SimulationStats3D::migration_attempts,
        &SimulationStats3D::migration_commits,
        &SimulationStats3D::divisions,
        &SimulationStats3D::deaths,
        &SimulationStats3D::conflict_rejections,
        &SimulationStats3D::angiogenesis_seed_attempts,
        &SimulationStats3D::angiogenesis_roots,
        &SimulationStats3D::angiogenesis_seed_rejections,
        &SimulationStats3D::vessel_growth_attempts,
        &SimulationStats3D::vessel_growth_commits,
        &SimulationStats3D::vessel_anastomoses,
        &SimulationStats3D::vascular_displacements,
    };
    for (const StatCounter counter : stat_counters) {
        SimulationStats3D changed = data.stats;
        ++(changed.*counter);
        assert(checksum_for(changed, data.lineage) != expected_checksum);
    }
    std::vector<LineageEdge> changed_lineage = data.lineage;
    changed_lineage.push_back({0.0, data.cells.front().uid, 0, 17,
                               data.cells.front().type});
    assert(checksum_for(data.stats, std::move(changed_lineage)) != expected_checksum);

    Model3DConfig long_config = cell_only_config(10.0);
    long_config.threads = 2;
    Simulation3D continuous(long_config);
    continuous.run();
    const CheckpointData3D resume_data = read_hdf5_checkpoint(cell_path, long_config);
    Simulation3D resumed(long_config);
    resumed.restore(resume_data.cells, resume_data.next_uid, resume_data.clock,
                    resume_data.stats, resume_data.lineage, resume_data.vasculature,
                    resume_data.cell_slot_count, resume_data.cell_slots,
                    resume_data.cell_free_slots);
    resumed.run();
    assert(resumed.state_checksum() == continuous.state_checksum());
    assert(resumed.clock().completed_events == continuous.clock().completed_events);

    // A checkpoint taken while bidirectional vessel tips have pending events
    // must be bit-for-bit continuous after resume, including reconstructed
    // vessel occupancy and perfusion influence.
    Model3DConfig short_vascular = vascular_config(0.45);
    Simulation3D vascular(short_vascular);
    vascular.run();
    assert(vascular.vessel_nodes().alive_count() > 1);
    assert(vascular.active_vessel_tip_count() > 0);
    assert(vascular.stats().angiogenesis_roots == 1);
    assert(vascular.stats().vessel_growth_commits > 0);
    const VasculatureState3D vascular_snapshot =
        vascular.snapshot_vasculature();
    assert(!vascular_snapshot.lesions.dirty_blocks.empty());
    const std::uint64_t vascular_checkpoint_checksum = vascular.state_checksum();
    write_hdf5_checkpoint(vascular_path, vascular);
    const CheckpointData3D vascular_data =
        read_hdf5_checkpoint(vascular_path, short_vascular);
    assert(!vascular_data.vasculature.nodes.empty());
    assert(vascular_data.vasculature.tips.size() == 2);
    assert(vascular_data.vasculature.process.committed_roots == 1);
    assert(vascular_data.vasculature.lesions.next_lesion_id > 1);
    assert(vascular_data.vasculature.lesions.next_refresh_time_hours ==
           vascular_snapshot.lesions.next_refresh_time_hours);
    assert(vascular_data.vasculature.lesions.refresh_schedule_generation ==
           vascular_snapshot.lesions.refresh_schedule_generation);
    assert(vascular_data.vasculature.lesions.dirty_blocks ==
           vascular_snapshot.lesions.dirty_blocks);
    assert(!vascular_data.vasculature.lesions.core_identity.empty());
    assert(vascular_data.vasculature.lesions.processes.size() == 1);
    const LesionId source_lesion_id =
        vascular_data.vasculature.lesions.processes.front().lesion_id;
    assert(source_lesion_id != kNoLesionId);
    assert(vascular_data.vasculature.lesions.processes.front()
               .process.committed_roots == 1);
    for (const VesselNodeInit3D& node : vascular_data.vasculature.nodes) {
        assert(node.source_lesion_id == source_lesion_id);
    }
    for (const VesselTipInit3D& tip : vascular_data.vasculature.tips) {
        assert(tip.source_lesion_id == source_lesion_id);
    }
    assert(vascular_data.vasculature.next_vessel_id == 2);
    assert(!vascular_data.vasculature.perfused_vessels.empty());
    {
        H5::H5File file(vascular_path.string(), H5F_ACC_RDONLY);
        assert(file.openDataSet("/vasculature/lesions/core_identity/lesion_id")
                   .getSpace().getSimpleExtentNpoints() > 0);
        assert(file.openDataSet("/vasculature/lesions/processes/lesion_id")
                   .getSpace().getSimpleExtentNpoints() == 1);
        assert(file.openDataSet("/vasculature/lesions/dirty_blocks/block_x")
                   .getSpace().getSimpleExtentNpoints() ==
               vascular_snapshot.lesions.dirty_blocks.size());
        assert(file.openDataSet("/vasculature/nodes/source_lesion_id")
                   .getSpace().getSimpleExtentNpoints() ==
               vascular_data.vasculature.nodes.size());
        assert(file.openDataSet("/vasculature/tips/source_lesion_id")
                   .getSpace().getSimpleExtentNpoints() ==
               vascular_data.vasculature.tips.size());
    }

    Simulation3D vascular_round_trip(short_vascular);
    vascular_round_trip.restore(
        vascular_data.cells, vascular_data.next_uid, vascular_data.clock,
        vascular_data.stats, vascular_data.lineage, vascular_data.vasculature,
        vascular_data.cell_slot_count, vascular_data.cell_slots,
        vascular_data.cell_free_slots);
    assert(vascular_round_trip.state_checksum() == vascular_checkpoint_checksum);
    assert(vascular_round_trip.vessel_grid().occupied_voxel_count() ==
           vascular.vessel_grid().occupied_voxel_count());
    assert(vascular_round_trip.vascular_influence().influenced_voxel_count() ==
           vascular.vascular_influence().influenced_voxel_count());

    // A merge transfers predecessor counters into the retained process but
    // intentionally leaves vessel provenance unchanged. Persist the sparse
    // historical-source -> current-owner alias and verify that per-current-
    // lesion tip limits include tips born under the retired source ID.
    VasculatureState3D historical_source = vascular_data.vasculature;
    const LesionId retired_source = historical_source.lesions.next_lesion_id++;
    historical_source.lesions.source_ownership.push_back(
        {retired_source, source_lesion_id});
    for (VesselNodeInit3D& node : historical_source.nodes) {
        node.source_lesion_id = retired_source;
    }
    for (VesselTipInit3D& tip : historical_source.tips) {
        tip.source_lesion_id = retired_source;
    }
    Simulation3D historical_source_simulation(short_vascular);
    historical_source_simulation.restore(
        vascular_data.cells, vascular_data.next_uid, vascular_data.clock,
        vascular_data.stats, vascular_data.lineage, historical_source,
        vascular_data.cell_slot_count, vascular_data.cell_slots,
        vascular_data.cell_free_slots);
    assert(historical_source_simulation.active_vessel_tip_count(
               source_lesion_id) ==
           historical_source_simulation.active_vessel_tip_count());
    write_hdf5_checkpoint(historical_source_path,
                          historical_source_simulation);
    const CheckpointData3D historical_data =
        read_hdf5_checkpoint(historical_source_path, short_vascular);
    assert(std::all_of(
        historical_data.vasculature.nodes.begin(),
        historical_data.vasculature.nodes.end(),
        [retired_source](const VesselNodeInit3D& node) {
            return node.source_lesion_id == retired_source;
        }));
    assert(std::all_of(
        historical_data.vasculature.tips.begin(),
        historical_data.vasculature.tips.end(),
        [retired_source](const VesselTipInit3D& tip) {
            return tip.source_lesion_id == retired_source;
        }));
    const std::vector<LesionSourceOwnership3D> expected_source_ownership{
        {retired_source, source_lesion_id}};
    assert(historical_data.vasculature.lesions.source_ownership ==
           expected_source_ownership);
    assert(std::none_of(
        historical_data.vasculature.lesions.processes.begin(),
        historical_data.vasculature.lesions.processes.end(),
        [retired_source](const LesionAngiogenesisState3D& process) {
            return process.lesion_id == retired_source;
        }));
    {
        H5::H5File file(historical_source_path.string(), H5F_ACC_RDONLY);
        H5::DataSet sources = file.openDataSet(
            "/vasculature/lesions/source_ownership/source_lesion_id");
        H5::DataSet owners = file.openDataSet(
            "/vasculature/lesions/source_ownership/current_lesion_id");
        assert(sources.getSpace().getSimpleExtentNpoints() == 1);
        assert(owners.getSpace().getSimpleExtentNpoints() == 1);
        std::uint64_t persisted_source{};
        std::uint64_t persisted_owner{};
        sources.read(&persisted_source, H5::PredType::NATIVE_UINT64);
        owners.read(&persisted_owner, H5::PredType::NATIVE_UINT64);
        assert(persisted_source == retired_source);
        assert(persisted_owner == source_lesion_id);
    }
    Simulation3D historical_source_round_trip(short_vascular);
    historical_source_round_trip.restore(
        historical_data.cells, historical_data.next_uid,
        historical_data.clock, historical_data.stats,
        historical_data.lineage, historical_data.vasculature,
        historical_data.cell_slot_count, historical_data.cell_slots,
        historical_data.cell_free_slots);
    assert(historical_source_round_trip.active_vessel_tip_count(
               source_lesion_id) ==
           historical_source_round_trip.active_vessel_tip_count());

    VasculatureState3D no_current_owner = historical_data.vasculature;
    no_current_owner.lesions.source_ownership.front().current_lesion_id =
        kNoLesionId;
    Simulation3D no_current_owner_simulation(short_vascular);
    no_current_owner_simulation.restore(
        historical_data.cells, historical_data.next_uid,
        historical_data.clock, historical_data.stats,
        historical_data.lineage, no_current_owner,
        historical_data.cell_slot_count, historical_data.cell_slots,
        historical_data.cell_free_slots);
    assert(no_current_owner_simulation.active_vessel_tip_count(
               source_lesion_id) == 0);
    assert(no_current_owner_simulation.state_checksum() !=
           historical_source_round_trip.state_checksum());

    VasculatureState3D missing_historical_owner =
        historical_data.vasculature;
    missing_historical_owner.lesions.source_ownership.clear();
    expect_rejected([&] {
        Simulation3D missing_historical_owner_simulation(short_vascular);
        missing_historical_owner_simulation.restore(
            historical_data.cells, historical_data.next_uid,
            historical_data.clock, historical_data.stats,
            historical_data.lineage, missing_historical_owner,
            historical_data.cell_slot_count, historical_data.cell_slots,
            historical_data.cell_free_slots);
    });

    // Disconnected primary/metastatic lesions keep independent process rows.
    // Their aggregate counters and historical source IDs must survive one
    // checkpoint without collapsing back to the old global-process model.
    const Model3DConfig two_config = two_lesion_config();
    CellUid next_lesion_uid = 2000;
    std::vector<CellInit> two_lesion_cells =
        lesion_cluster({0, 0, 0}, next_lesion_uid);
    std::vector<CellInit> distant_cells =
        lesion_cluster({30, 0, 0}, next_lesion_uid);
    two_lesion_cells.insert(two_lesion_cells.end(), distant_cells.begin(),
                            distant_cells.end());
    Simulation3D two_lesions(two_config);
    two_lesions.restore(two_lesion_cells, next_lesion_uid, {}, {}, {});
    two_lesions.run();
    assert(two_lesions.stats().angiogenesis_roots == 2);
    const std::uint64_t two_lesion_checksum = two_lesions.state_checksum();
    write_hdf5_checkpoint(two_lesion_path, two_lesions);
    const CheckpointData3D two_data =
        read_hdf5_checkpoint(two_lesion_path, two_config);
    // Root displacement may trigger stage recovery and deterministic lesion
    // splits before this checkpoint, so additional child-process rows are
    // valid; the two original source lesions must not collapse to one.
    assert(two_data.vasculature.lesions.processes.size() >= 2U);
    assert(two_data.vasculature.process.committed_roots == 2);
    std::set<LesionId> source_lesions;
    for (const VesselNodeInit3D& node : two_data.vasculature.nodes) {
        if (node.role == VesselBranchRole::root) {
            source_lesions.insert(node.source_lesion_id);
        }
    }
    assert(source_lesions.size() == 2);
    Simulation3D two_lesions_round_trip(two_config);
    two_lesions_round_trip.restore(
        two_data.cells, two_data.next_uid, two_data.clock, two_data.stats,
        two_data.lineage, two_data.vasculature, two_data.cell_slot_count,
        two_data.cell_slots, two_data.cell_free_slots);
    assert(two_lesions_round_trip.state_checksum() == two_lesion_checksum);

    Model3DConfig long_vascular = vascular_config(2.0);
    long_vascular.threads = 4;
    Simulation3D vascular_continuous(long_vascular);
    vascular_continuous.run();
    const CheckpointData3D vascular_resume_data =
        read_hdf5_checkpoint(vascular_path, long_vascular);
    Simulation3D vascular_resumed(long_vascular);
    vascular_resumed.restore(
        vascular_resume_data.cells, vascular_resume_data.next_uid,
        vascular_resume_data.clock, vascular_resume_data.stats,
        vascular_resume_data.lineage, vascular_resume_data.vasculature,
        vascular_resume_data.cell_slot_count, vascular_resume_data.cell_slots,
        vascular_resume_data.cell_free_slots);
    vascular_resumed.run();
    assert(vascular_resumed.state_checksum() == vascular_continuous.state_checksum());
    assert(vascular_resumed.stats().vessel_growth_commits ==
           vascular_continuous.stats().vessel_growth_commits);
    assert(vascular_resumed.stats().vascular_displacements ==
           vascular_continuous.stats().vascular_displacements);

    // Slot topology is part of the deterministic continuation state. Preserve
    // holes and the free-list LIFO order, rather than densely recreating live
    // cells and changing the slot chosen by the next birth.
    std::vector<CellInit> hole_cells(3);
    for (std::size_t index = 0; index < hole_cells.size(); ++index) {
        CellInit& cell = hole_cells[index];
        cell.anchor = {static_cast<int>(index) * 3, 30, 0};
        cell.uid = static_cast<CellUid>(100 + index);
        cell.clone_id = static_cast<std::uint32_t>(index + 1);
        cell.density_growth_rate = 1.0F;
        cell.migration_rate = 0.0F;
        cell.division_work_remaining = 10.0F;
        cell.next_migration_time = 0.0;
        cell.next_division_time = 10.0;
        cell.death_deadline = 0.0;
        cell.migration_schedule_generation = 1;
        cell.division_schedule_generation = 1;
        cell.death_schedule_generation = 1;
    }
    const std::vector<Slot> hole_alive_slots{0, 2, 4};
    const std::vector<Slot> hole_free_slots{1, 3};
    Simulation3D holes(short_config);
    holes.restore(hole_cells, 200, {}, {}, {}, {}, 5, hole_alive_slots,
                  hole_free_slots);
    const std::uint64_t holes_checksum = holes.state_checksum();
    write_hdf5_checkpoint(holes_path, holes);
    const CheckpointData3D holes_data =
        read_hdf5_checkpoint(holes_path, short_config);
    assert(holes_data.cell_slot_count == 5);
    assert(holes_data.cell_slots == hole_alive_slots);
    assert(holes_data.cell_free_slots == hole_free_slots);
    Simulation3D holes_restored(short_config);
    holes_restored.restore(
        holes_data.cells, holes_data.next_uid, holes_data.clock,
        holes_data.stats, holes_data.lineage, holes_data.vasculature,
        holes_data.cell_slot_count, holes_data.cell_slots,
        holes_data.cell_free_slots);
    assert(holes_restored.state_checksum() == holes_checksum);
    CellInit next_cell = hole_cells.front();
    next_cell.uid = 200;
    next_cell.anchor = {20, 30, 0};
    assert(holes_restored.cells().create(next_cell) == 3);

    Model3DConfig incompatible = long_config;
    incompatible.continue_probability = 0.5;
    expect_rejected([&] { (void)read_hdf5_checkpoint(cell_path, incompatible); });
    Model3DConfig incompatible_activated_r_scale = long_config;
    incompatible_activated_r_scale.activated_r_migration_beta.scale = 2.0;
    expect_rejected([&] {
        (void)read_hdf5_checkpoint(cell_path,
                                   incompatible_activated_r_scale);
    });

    {
        std::ofstream stream(corrupt_path, std::ios::binary);
        stream << "not an HDF5 checkpoint";
    }
    expect_rejected([&] { (void)read_hdf5_checkpoint(corrupt_path, long_config); });

    corrupt_schema_version(cell_path, version_path, 99);
    expect_rejected([&] { (void)read_hdf5_checkpoint(version_path, long_config); });
    corrupt_schema_version(cell_path, version_path, 1);
    expect_rejected([&] { (void)read_hdf5_checkpoint(version_path, long_config); });
    corrupt_schema_version(cell_path, version_path, 2);
    expect_rejected([&] { (void)read_hdf5_checkpoint(version_path, long_config); });

    corrupt_cell_column_length(cell_path, columns_path);
    expect_rejected([&] { (void)read_hdf5_checkpoint(columns_path, long_config); });

    corrupt_division_work(cell_path, work_path);
    expect_rejected([&] { (void)read_hdf5_checkpoint(work_path, long_config); });

    corrupt_migration_activation_pair(cell_path, activation_path);
    expect_rejected([&] {
        (void)read_hdf5_checkpoint(activation_path, long_config);
    });

    corrupt_lesion_process_aggregate(vascular_path, aggregate_path);
    expect_rejected([&] {
        (void)read_hdf5_checkpoint(aggregate_path, short_vascular);
    });
    corrupt_child_node_source_lesion(vascular_path, node_source_path);
    expect_rejected([&] {
        (void)read_hdf5_checkpoint(node_source_path, short_vascular);
    });
    corrupt_tip_source_lesion(vascular_path, tip_source_path);
    expect_rejected([&] {
        (void)read_hdf5_checkpoint(tip_source_path, short_vascular);
    });
    corrupt_next_lesion_refresh_time(vascular_path, refresh_time_path);
    expect_rejected([&] {
        (void)read_hdf5_checkpoint(refresh_time_path, short_vascular);
    });
    corrupt_lesion_dirty_block_exists(vascular_path, dirty_block_path);
    expect_rejected([&] {
        (void)read_hdf5_checkpoint(dirty_block_path, short_vascular);
    });
    remove_effective_config_provenance(cell_path,
                                       missing_effective_config_path);
    expect_rejected([&] {
        (void)read_hdf5_checkpoint(missing_effective_config_path,
                                   long_config);
    });
    corrupt_effective_config_provenance(cell_path,
                                        corrupt_effective_config_path);
    expect_rejected([&] {
        (void)read_hdf5_checkpoint(corrupt_effective_config_path,
                                   long_config);
    });
    remove_source_ownership_table(historical_source_path,
                                  missing_source_ownership_path);
    expect_rejected([&] {
        (void)read_hdf5_checkpoint(missing_source_ownership_path,
                                   short_vascular);
    });
    corrupt_source_ownership_target(historical_source_path,
                                    corrupt_source_ownership_path);
    expect_rejected([&] {
        (void)read_hdf5_checkpoint(corrupt_source_ownership_path,
                                   short_vascular);
    });

    std::filesystem::remove_all(directory);
    return 0;
}
