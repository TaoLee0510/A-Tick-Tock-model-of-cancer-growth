#include "io/checkpoint_hdf5.hpp"

#include <H5Cpp.h>

#include <algorithm>
#include <array>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace atcg3d {
namespace {

template <class T>
const H5::PredType& hdf_type();

template <> const H5::PredType& hdf_type<std::uint8_t>() { return H5::PredType::NATIVE_UINT8; }
template <> const H5::PredType& hdf_type<std::uint32_t>() { return H5::PredType::NATIVE_UINT32; }
template <> const H5::PredType& hdf_type<std::uint64_t>() { return H5::PredType::NATIVE_UINT64; }
template <> const H5::PredType& hdf_type<std::int32_t>() { return H5::PredType::NATIVE_INT32; }
template <> const H5::PredType& hdf_type<float>() { return H5::PredType::NATIVE_FLOAT; }
template <> const H5::PredType& hdf_type<double>() { return H5::PredType::NATIVE_DOUBLE; }

template <class T>
void write_vector(H5::Group& group, const std::string& name, const std::vector<T>& values) {
    const hsize_t dimensions[1] = {static_cast<hsize_t>(values.size())};
    H5::DataSpace space(1, dimensions);
    H5::DataSet dataset = group.createDataSet(name, hdf_type<T>(), space);
    if (!values.empty()) {
        dataset.write(values.data(), hdf_type<T>());
    }
}

template <class T>
std::vector<T> read_vector(H5::Group& group, const std::string& name) {
    H5::DataSet dataset = group.openDataSet(name);
    H5::DataSpace space = dataset.getSpace();
    if (space.getSimpleExtentNdims() != 1) {
        throw std::runtime_error("checkpoint dataset is not one-dimensional: " + name);
    }
    hsize_t dimensions[1]{};
    space.getSimpleExtentDims(dimensions);
    std::vector<T> values(static_cast<std::size_t>(dimensions[0]));
    if (!values.empty()) {
        dataset.read(values.data(), hdf_type<T>());
    }
    return values;
}

template <class T>
void write_scalar_attribute(H5::H5Object& object, const std::string& name, const T& value) {
    H5::DataSpace scalar(H5S_SCALAR);
    H5::Attribute attribute = object.createAttribute(name, hdf_type<T>(), scalar);
    attribute.write(hdf_type<T>(), &value);
}

template <class T>
T read_scalar_attribute(H5::H5Object& object, const std::string& name) {
    T value{};
    H5::Attribute attribute = object.openAttribute(name);
    attribute.read(hdf_type<T>(), &value);
    return value;
}

void write_string_attribute(H5::H5Object& object, const std::string& name, const std::string& value) {
    H5::StrType string_type(H5::PredType::C_S1, H5T_VARIABLE);
    H5::DataSpace scalar(H5S_SCALAR);
    H5::Attribute attribute = object.createAttribute(name, string_type, scalar);
    attribute.write(string_type, value);
}

std::string read_string_attribute(H5::H5Object& object, const std::string& name) {
    H5::Attribute attribute = object.openAttribute(name);
    H5::StrType string_type = attribute.getStrType();
    std::string value;
    attribute.read(string_type, value);
    return value;
}

template <class Member, class Row, class Getter>
std::vector<Member> column(const std::vector<Row>& cells, Getter getter) {
    std::vector<Member> values;
    values.reserve(cells.size());
    for (const Row& cell : cells) {
        values.push_back(getter(cell));
    }
    return values;
}

void require_equal_sizes(std::size_t expected,
                         const std::vector<std::pair<std::string, std::size_t>>& sizes) {
    for (const auto& [name, size] : sizes) {
        if (size != expected) {
            throw std::runtime_error("checkpoint cell column length mismatch: " + name);
        }
    }
}

}  // namespace

void write_hdf5_checkpoint(const std::filesystem::path& path,
                           const Simulation3D& simulation) {
    std::filesystem::create_directories(path.parent_path());
    const std::filesystem::path temporary = path.string() + ".tmp";
    if (std::filesystem::exists(path)) {
        throw std::runtime_error("refusing to overwrite checkpoint: " + path.string());
    }
    std::filesystem::remove(temporary);
    try {
        H5::H5File file(temporary.string(), H5F_ACC_TRUNC);
        H5::Group meta = file.createGroup("/meta");
        H5::Group cells_group = file.createGroup("/cells");
        H5::Group lineage_group = file.createGroup("/lineage");

        write_scalar_attribute(meta, "schema_version", std::uint32_t{1});
        write_scalar_attribute(meta, "dimension", std::uint32_t{3});
        write_scalar_attribute(meta, "completed_events", simulation.clock().completed_events);
        write_scalar_attribute(meta, "time_hours", simulation.clock().time_hours);
        write_scalar_attribute(meta, "next_uid", simulation.next_uid());
        write_scalar_attribute(meta, "state_checksum", simulation.state_checksum());
        write_scalar_attribute(meta, "migration_attempts", simulation.stats().migration_attempts);
        write_scalar_attribute(meta, "migration_commits", simulation.stats().migration_commits);
        write_scalar_attribute(meta, "divisions", simulation.stats().divisions);
        write_scalar_attribute(meta, "deaths", simulation.stats().deaths);
        write_scalar_attribute(meta, "conflict_rejections", simulation.stats().conflict_rejections);
        write_string_attribute(meta, "effective_config_json", simulation.config().to_json());
        write_string_attribute(meta, "dynamics_config_json", simulation.config().dynamics_json());

        const std::vector<CellInit> cells = simulation.snapshot_cells();
        write_vector(cells_group, "x", column<std::int32_t>(cells, [](const CellInit& c) { return c.anchor.x; }));
        write_vector(cells_group, "y", column<std::int32_t>(cells, [](const CellInit& c) { return c.anchor.y; }));
        write_vector(cells_group, "z", column<std::int32_t>(cells, [](const CellInit& c) { return c.anchor.z; }));
        write_vector(cells_group, "uid", column<std::uint64_t>(cells, [](const CellInit& c) { return c.uid; }));
        write_vector(cells_group, "parent_uid", column<std::uint64_t>(cells, [](const CellInit& c) { return c.parent_uid; }));
        write_vector(cells_group, "clone_id", column<std::uint32_t>(cells, [](const CellInit& c) { return c.clone_id; }));
        write_vector(cells_group, "type", column<std::uint8_t>(cells, [](const CellInit& c) { return static_cast<std::uint8_t>(c.type); }));
        write_vector(cells_group, "stage", column<std::uint8_t>(cells, [](const CellInit& c) { return static_cast<std::uint8_t>(c.stage); }));
        write_vector(cells_group, "viability", column<std::uint8_t>(cells, [](const CellInit& c) { return c.viability; }));
        write_vector(cells_group, "flags", column<std::uint8_t>(cells, [](const CellInit& c) { return c.flags; }));
        write_vector(cells_group, "last_direction", column<std::uint8_t>(cells, [](const CellInit& c) { return c.last_direction; }));
        write_vector(cells_group, "inherent_growth_rate", column<float>(cells, [](const CellInit& c) { return c.inherent_growth_rate; }));
        write_vector(cells_group, "density_growth_rate", column<float>(cells, [](const CellInit& c) { return c.density_growth_rate; }));
        write_vector(cells_group, "migration_rate", column<float>(cells, [](const CellInit& c) { return c.migration_rate; }));
        write_vector(cells_group, "next_migration_time", column<double>(cells, [](const CellInit& c) { return c.next_migration_time; }));
        write_vector(cells_group, "next_division_time", column<double>(cells, [](const CellInit& c) { return c.next_division_time; }));
        write_vector(cells_group, "death_deadline", column<float>(cells, [](const CellInit& c) { return static_cast<float>(c.death_deadline); }));
        write_vector(cells_group, "last_update_time", column<double>(cells, [](const CellInit& c) { return c.last_update_time; }));
        write_vector(cells_group, "event_sequence", column<std::uint64_t>(cells, [](const CellInit& c) { return c.event_sequence; }));
        write_vector(cells_group, "schedule_generation", column<std::uint32_t>(cells, [](const CellInit& c) { return c.schedule_generation; }));

        const auto& lineage = simulation.lineage();
        write_vector(lineage_group, "birth_time", column<double>(lineage, [](const LineageEdge& e) { return e.birth_time; }));
        write_vector(lineage_group, "child_uid", column<std::uint64_t>(lineage, [](const LineageEdge& e) { return e.child_uid; }));
        write_vector(lineage_group, "parent_uid", column<std::uint64_t>(lineage, [](const LineageEdge& e) { return e.parent_uid; }));
        write_vector(lineage_group, "clone_id", column<std::uint32_t>(lineage, [](const LineageEdge& e) { return e.clone_id; }));
        write_vector(lineage_group, "type", column<std::uint8_t>(lineage, [](const LineageEdge& e) { return static_cast<std::uint8_t>(e.type); }));
        file.flush(H5F_SCOPE_GLOBAL);
        file.close();
        std::filesystem::rename(temporary, path);
    } catch (...) {
        std::filesystem::remove(temporary);
        throw;
    }
}

CheckpointData3D read_hdf5_checkpoint(const std::filesystem::path& path,
                                      const Model3DConfig& expected_config) {
    try {
        H5::H5File file(path.string(), H5F_ACC_RDONLY);
        H5::Group meta = file.openGroup("/meta");
        if (read_scalar_attribute<std::uint32_t>(meta, "schema_version") != 1 ||
            read_scalar_attribute<std::uint32_t>(meta, "dimension") != 3) {
            throw std::runtime_error("unsupported checkpoint schema or dimension");
        }
        if (read_string_attribute(meta, "dynamics_config_json") != expected_config.dynamics_json()) {
            throw std::runtime_error("checkpoint biological/numerical configuration does not match");
        }

        H5::Group cells_group = file.openGroup("/cells");
        const auto x = read_vector<std::int32_t>(cells_group, "x");
        const auto y = read_vector<std::int32_t>(cells_group, "y");
        const auto z = read_vector<std::int32_t>(cells_group, "z");
        const auto uid = read_vector<std::uint64_t>(cells_group, "uid");
        const auto parent_uid = read_vector<std::uint64_t>(cells_group, "parent_uid");
        const auto clone_id = read_vector<std::uint32_t>(cells_group, "clone_id");
        const auto type = read_vector<std::uint8_t>(cells_group, "type");
        const auto stage = read_vector<std::uint8_t>(cells_group, "stage");
        const auto viability = read_vector<std::uint8_t>(cells_group, "viability");
        const auto flags = read_vector<std::uint8_t>(cells_group, "flags");
        const auto last_direction = read_vector<std::uint8_t>(cells_group, "last_direction");
        const auto inherent_growth_rate = read_vector<float>(cells_group, "inherent_growth_rate");
        const auto density_growth_rate = read_vector<float>(cells_group, "density_growth_rate");
        const auto migration_rate = read_vector<float>(cells_group, "migration_rate");
        const auto next_migration_time = read_vector<double>(cells_group, "next_migration_time");
        const auto next_division_time = read_vector<double>(cells_group, "next_division_time");
        const auto death_deadline = read_vector<float>(cells_group, "death_deadline");
        const auto last_update_time = read_vector<double>(cells_group, "last_update_time");
        const auto event_sequence = read_vector<std::uint64_t>(cells_group, "event_sequence");
        const auto schedule_generation = read_vector<std::uint32_t>(cells_group, "schedule_generation");
        const std::size_t count = uid.size();
        require_equal_sizes(count, {{"x", x.size()}, {"y", y.size()}, {"z", z.size()},
                                    {"parent_uid", parent_uid.size()}, {"clone_id", clone_id.size()},
                                    {"type", type.size()}, {"stage", stage.size()},
                                    {"viability", viability.size()}, {"flags", flags.size()},
                                    {"last_direction", last_direction.size()},
                                    {"inherent_growth_rate", inherent_growth_rate.size()},
                                    {"density_growth_rate", density_growth_rate.size()},
                                    {"migration_rate", migration_rate.size()},
                                    {"next_migration_time", next_migration_time.size()},
                                    {"next_division_time", next_division_time.size()},
                                    {"death_deadline", death_deadline.size()},
                                    {"last_update_time", last_update_time.size()},
                                    {"event_sequence", event_sequence.size()},
                                    {"schedule_generation", schedule_generation.size()}});
        std::vector<CellInit> cells(count);
        for (std::size_t index = 0; index < count; ++index) {
            cells[index] = {{x[index], y[index], z[index]}, uid[index], parent_uid[index], clone_id[index],
                            static_cast<CellType>(type[index]), static_cast<CellStage>(stage[index]),
                            viability[index], flags[index], last_direction[index], inherent_growth_rate[index],
                            density_growth_rate[index], migration_rate[index], next_migration_time[index],
                            next_division_time[index], death_deadline[index], last_update_time[index],
                            event_sequence[index], schedule_generation[index]};
        }
        std::sort(cells.begin(), cells.end(), [](const CellInit& lhs, const CellInit& rhs) { return lhs.uid < rhs.uid; });
        for (std::size_t index = 1; index < cells.size(); ++index) {
            if (cells[index - 1].uid == cells[index].uid || cells[index].uid == 0) {
                throw std::runtime_error("checkpoint contains duplicate or zero cell uid");
            }
        }

        H5::Group lineage_group = file.openGroup("/lineage");
        const auto birth_time = read_vector<double>(lineage_group, "birth_time");
        const auto child_uid = read_vector<std::uint64_t>(lineage_group, "child_uid");
        const auto lineage_parent_uid = read_vector<std::uint64_t>(lineage_group, "parent_uid");
        const auto lineage_clone_id = read_vector<std::uint32_t>(lineage_group, "clone_id");
        const auto lineage_type = read_vector<std::uint8_t>(lineage_group, "type");
        require_equal_sizes(birth_time.size(), {{"child_uid", child_uid.size()},
                                                {"parent_uid", lineage_parent_uid.size()},
                                                {"clone_id", lineage_clone_id.size()},
                                                {"type", lineage_type.size()}});
        std::vector<LineageEdge> lineage(birth_time.size());
        for (std::size_t index = 0; index < lineage.size(); ++index) {
            lineage[index] = {birth_time[index], child_uid[index], lineage_parent_uid[index],
                              lineage_clone_id[index], static_cast<CellType>(lineage_type[index])};
        }

        CheckpointData3D result;
        result.cells = std::move(cells);
        result.next_uid = read_scalar_attribute<std::uint64_t>(meta, "next_uid");
        result.clock.completed_events = read_scalar_attribute<std::uint64_t>(meta, "completed_events");
        result.clock.time_hours = read_scalar_attribute<double>(meta, "time_hours");
        result.state_checksum = read_scalar_attribute<std::uint64_t>(meta, "state_checksum");
        result.stats.migration_attempts = read_scalar_attribute<std::uint64_t>(meta, "migration_attempts");
        result.stats.migration_commits = read_scalar_attribute<std::uint64_t>(meta, "migration_commits");
        result.stats.divisions = read_scalar_attribute<std::uint64_t>(meta, "divisions");
        result.stats.deaths = read_scalar_attribute<std::uint64_t>(meta, "deaths");
        result.stats.conflict_rejections = read_scalar_attribute<std::uint64_t>(meta, "conflict_rejections");
        result.lineage = std::move(lineage);
        return result;
    } catch (const H5::Exception& error) {
        throw std::runtime_error("invalid or unreadable HDF5 checkpoint: " + error.getDetailMsg());
    }
}

}  // namespace atcg3d
