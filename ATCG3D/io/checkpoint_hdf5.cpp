#include "io/checkpoint_hdf5.hpp"

#include <H5Cpp.h>
#include <hdf5.h>
#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iterator>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace atcg3d {
namespace {

constexpr std::uint32_t kLegacyCheckpointSchemaVersion = 4;
constexpr std::uint32_t kLegacyCheckpointDeltaSchemaVersion = 5;
constexpr std::uint32_t kCheckpointSchemaVersion = 6;
constexpr std::uint32_t kCheckpointDeltaSchemaVersion = 7;
constexpr std::uint32_t kCheckpointJournalDeltaSchemaVersion = 8;
constexpr std::uint32_t kCheckpointDimension = 3;

struct DatasetWriteOptions {
    int compression_level{};
    hsize_t chunk_elements{262144};
};

thread_local DatasetWriteOptions g_dataset_write_options{};

class ScopedDatasetWriteOptions {
public:
    explicit ScopedDatasetWriteOptions(const Model3DConfig& config)
        : previous_(g_dataset_write_options) {
        g_dataset_write_options.compression_level =
            config.hdf5_compression_level;
        g_dataset_write_options.chunk_elements = static_cast<hsize_t>(
            config.hdf5_chunk_elements);
    }
    ~ScopedDatasetWriteOptions() { g_dataset_write_options = previous_; }

private:
    DatasetWriteOptions previous_;
};

class CheckpointReadGuard {
public:
    explicit CheckpointReadGuard(const std::filesystem::path& path)
        : key_(std::filesystem::absolute(path).lexically_normal().string()) {
        if (stack_.size() >= 4096 ||
            std::find(stack_.begin(), stack_.end(), key_) != stack_.end()) {
            throw std::runtime_error(
                "checkpoint delta chain is cyclic or exceeds the hard limit");
        }
        stack_.push_back(key_);
        outermost_ = stack_.size() == 1;
    }
    ~CheckpointReadGuard() { stack_.pop_back(); }
    bool outermost() const noexcept { return outermost_; }

private:
    std::string key_;
    bool outermost_{};
    static thread_local std::vector<std::string> stack_;
};

thread_local std::vector<std::string> CheckpointReadGuard::stack_;

double float_storage_time_tolerance(double lhs, double rhs) noexcept {
    return 2.0 * static_cast<double>(std::numeric_limits<float>::epsilon()) *
           std::max({1.0, std::abs(lhs), std::abs(rhs)});
}

template <class T>
const H5::PredType& file_hdf_type();

template <class T>
const H5::PredType& native_hdf_type();

template <> const H5::PredType& file_hdf_type<std::uint8_t>() {
    return H5::PredType::STD_U8LE;
}
template <> const H5::PredType& file_hdf_type<std::uint32_t>() {
    return H5::PredType::STD_U32LE;
}
template <> const H5::PredType& file_hdf_type<std::uint64_t>() {
    return H5::PredType::STD_U64LE;
}
template <> const H5::PredType& file_hdf_type<std::int32_t>() {
    return H5::PredType::STD_I32LE;
}
template <> const H5::PredType& file_hdf_type<std::int64_t>() {
    return H5::PredType::STD_I64LE;
}
template <> const H5::PredType& file_hdf_type<float>() {
    return H5::PredType::IEEE_F32LE;
}
template <> const H5::PredType& file_hdf_type<double>() {
    return H5::PredType::IEEE_F64LE;
}

template <> const H5::PredType& native_hdf_type<std::uint8_t>() {
    return H5::PredType::NATIVE_UINT8;
}
template <> const H5::PredType& native_hdf_type<std::uint32_t>() {
    return H5::PredType::NATIVE_UINT32;
}
template <> const H5::PredType& native_hdf_type<std::uint64_t>() {
    return H5::PredType::NATIVE_UINT64;
}
template <> const H5::PredType& native_hdf_type<std::int32_t>() {
    return H5::PredType::NATIVE_INT32;
}
template <> const H5::PredType& native_hdf_type<std::int64_t>() {
    return H5::PredType::NATIVE_INT64;
}
template <> const H5::PredType& native_hdf_type<float>() {
    return H5::PredType::NATIVE_FLOAT;
}
template <> const H5::PredType& native_hdf_type<double>() {
    return H5::PredType::NATIVE_DOUBLE;
}

template <class T>
void write_vector(H5::Group& group, const std::string& name,
                  std::span<const T> values) {
    const hsize_t dimensions[1] = {static_cast<hsize_t>(values.size())};
    H5::DataSpace space(1, dimensions);
    H5::DSetCreatPropList properties;
    if (!values.empty()) {
        const hsize_t chunk[1] = {
            std::min<hsize_t>(dimensions[0],
                              g_dataset_write_options.chunk_elements)};
        properties.setChunk(1, chunk);
        if (g_dataset_write_options.compression_level > 0) {
            properties.setShuffle();
            properties.setDeflate(
                g_dataset_write_options.compression_level);
        }
        properties.setFletcher32();
    }
    H5::DataSet dataset = values.empty()
        ? group.createDataSet(name, file_hdf_type<T>(), space)
        : group.createDataSet(name, file_hdf_type<T>(), space, properties);
    if (!values.empty()) dataset.write(values.data(), native_hdf_type<T>());
}

template <class T>
void write_vector(H5::Group& group, const std::string& name,
                  const std::vector<T>& values) {
    write_vector(group, name, std::span<const T>(values));
}

template <class T>
std::vector<T> read_vector(H5::Group& group, const std::string& name) {
    H5::DataSet dataset = group.openDataSet(name);
    if (H5Tequal(dataset.getDataType().getId(), file_hdf_type<T>().getId()) <= 0) {
        throw std::runtime_error("checkpoint dataset has the wrong fixed-width type: " + name);
    }
    H5::DataSpace space = dataset.getSpace();
    if (space.getSimpleExtentNdims() != 1) {
        throw std::runtime_error("checkpoint dataset is not one-dimensional: " + name);
    }
    hsize_t dimensions[1]{};
    space.getSimpleExtentDims(dimensions);
    if (dimensions[0] > static_cast<hsize_t>(std::numeric_limits<std::size_t>::max())) {
        throw std::runtime_error("checkpoint dataset is too large: " + name);
    }
    std::vector<T> values(static_cast<std::size_t>(dimensions[0]));
    if (!values.empty()) dataset.read(values.data(), native_hdf_type<T>());
    return values;
}

template <class T>
void write_scalar_attribute(H5::H5Object& object, const std::string& name,
                            const T& value) {
    H5::DataSpace scalar(H5S_SCALAR);
    H5::Attribute attribute = object.createAttribute(name, file_hdf_type<T>(), scalar);
    attribute.write(native_hdf_type<T>(), &value);
}

template <class T>
T read_scalar_attribute(H5::H5Object& object, const std::string& name) {
    H5::Attribute attribute = object.openAttribute(name);
    if (H5Tequal(attribute.getDataType().getId(), file_hdf_type<T>().getId()) <= 0) {
        throw std::runtime_error("checkpoint attribute has the wrong fixed-width type: " + name);
    }
    if (attribute.getSpace().getSimpleExtentType() != H5S_SCALAR) {
        throw std::runtime_error("checkpoint attribute is not scalar: " + name);
    }
    T value{};
    attribute.read(native_hdf_type<T>(), &value);
    return value;
}

void write_bool_attribute(H5::H5Object& object, const std::string& name, bool value) {
    write_scalar_attribute(object, name, static_cast<std::uint8_t>(value ? 1U : 0U));
}

bool read_bool_attribute(H5::H5Object& object, const std::string& name) {
    const std::uint8_t value = read_scalar_attribute<std::uint8_t>(object, name);
    if (value > 1U) throw std::runtime_error("checkpoint boolean attribute is invalid: " + name);
    return value != 0;
}

void write_string_attribute(H5::H5Object& object, const std::string& name,
                            const std::string& value) {
    H5::StrType string_type(H5::PredType::C_S1, H5T_VARIABLE);
    string_type.setCset(H5T_CSET_UTF8);
    H5::DataSpace scalar(H5S_SCALAR);
    H5::Attribute attribute = object.createAttribute(name, string_type, scalar);
    attribute.write(string_type, value);
}

std::string read_string_attribute(H5::H5Object& object, const std::string& name) {
    H5::Attribute attribute = object.openAttribute(name);
    if (attribute.getTypeClass() != H5T_STRING ||
        attribute.getSpace().getSimpleExtentType() != H5S_SCALAR) {
        throw std::runtime_error("checkpoint string attribute is invalid: " + name);
    }
    H5::StrType string_type = attribute.getStrType();
    std::string value;
    attribute.read(string_type, value);
    return value;
}

void validate_effective_config_provenance(
    const std::string& text,
    const Model3DConfig& expected_config) {
    if (text.empty()) {
        throw std::runtime_error(
            "checkpoint effective configuration provenance is empty");
    }
    try {
        const YAML::Node root = YAML::Load(text);
        if (!root || !root.IsMap()) {
            throw std::runtime_error(
                "checkpoint effective configuration provenance is not an object");
        }
        const auto required_scalar = [&root](const char* key) -> YAML::Node {
            const YAML::Node value = root[key];
            if (!value || !value.IsScalar()) {
                throw std::runtime_error(
                    std::string("checkpoint effective configuration is missing scalar ") +
                    key);
            }
            return value;
        };
        const std::string schema_name =
            required_scalar("schema_name").as<std::string>();
        const std::uint32_t schema_version =
            required_scalar("schema_version").as<std::uint32_t>();
        const std::string profile =
            required_scalar("profile").as<std::string>();
        if (schema_name != expected_config.schema_name ||
            schema_version != expected_config.schema_version ||
            profile != expected_config.profile) {
            throw std::runtime_error(
                "checkpoint effective configuration provenance identity does not match");
        }
    } catch (const YAML::Exception& error) {
        throw std::runtime_error(
            "checkpoint effective configuration provenance is malformed: " +
            std::string(error.what()));
    }
}

template <class Member, class Rows, class Getter>
std::vector<Member> column(const Rows& rows, Getter getter) {
    std::vector<Member> values;
    values.reserve(rows.size());
    for (const auto& row : rows) values.push_back(getter(row));
    return values;
}

void require_equal_sizes(
    const std::string& table, std::size_t expected,
    const std::vector<std::pair<std::string, std::size_t>>& sizes) {
    for (const auto& [name, size] : sizes) {
        if (size != expected) {
            throw std::runtime_error("checkpoint " + table +
                                     " column length mismatch: " + name);
        }
    }
}

void require_finite_nonnegative(double value, const std::string& name) {
    if (!std::isfinite(value) || value < 0.0) {
        throw std::runtime_error("checkpoint value must be finite and nonnegative: " + name);
    }
}

void require_finite(float value, const std::string& name) {
    if (!std::isfinite(value)) {
        throw std::runtime_error("checkpoint value must be finite: " + name);
    }
}

bool valid_cell_type(std::uint8_t value) {
    return value == static_cast<std::uint8_t>(CellType::r) ||
           value == static_cast<std::uint8_t>(CellType::K);
}

bool valid_cell_stage(std::uint8_t value) {
    return value <= static_cast<std::uint8_t>(CellStage::ultrasmall);
}

bool valid_vessel_role(std::uint8_t value) {
    return value <= static_cast<std::uint8_t>(VesselBranchRole::outward);
}

bool valid_tip_status(std::uint8_t value) {
    return value <= static_cast<std::uint8_t>(VesselTipStatus::transiting);
}

bool same_cell(const CellInit& lhs, const CellInit& rhs) noexcept {
    return lhs.anchor == rhs.anchor && lhs.uid == rhs.uid &&
           lhs.parent_uid == rhs.parent_uid && lhs.clone_id == rhs.clone_id &&
           lhs.type == rhs.type && lhs.stage == rhs.stage &&
           lhs.viability == rhs.viability && lhs.flags == rhs.flags &&
           lhs.last_direction == rhs.last_direction &&
           lhs.inherent_growth_rate == rhs.inherent_growth_rate &&
           lhs.density_growth_rate == rhs.density_growth_rate &&
           lhs.migration_rate == rhs.migration_rate &&
           lhs.normal_migration_rate == rhs.normal_migration_rate &&
           lhs.migration_activation_end_time ==
               rhs.migration_activation_end_time &&
           lhs.division_work_remaining == rhs.division_work_remaining &&
           lhs.next_migration_time == rhs.next_migration_time &&
           lhs.next_division_time == rhs.next_division_time &&
           lhs.death_deadline == rhs.death_deadline &&
           lhs.last_update_time == rhs.last_update_time &&
           lhs.swap_ready_time == rhs.swap_ready_time &&
           lhs.swap_wait_state == rhs.swap_wait_state &&
           lhs.pending_swap_direction == rhs.pending_swap_direction &&
           lhs.event_sequence == rhs.event_sequence &&
           lhs.migration_schedule_generation ==
               rhs.migration_schedule_generation &&
           lhs.division_schedule_generation ==
               rhs.division_schedule_generation &&
           lhs.death_schedule_generation == rhs.death_schedule_generation;
}

constexpr std::uint64_t field_bit(CheckpointCellField3D field) noexcept {
    return checkpoint_cell_field_mask(field);
}

constexpr std::uint64_t kKnownCellFieldMask =
    (1ULL << 25U) - 1ULL;

std::uint64_t changed_cell_field_mask(const CellInit& current,
                                      const CellInit& parent) noexcept {
    std::uint64_t mask = 0;
    const auto changed = [&](CheckpointCellField3D field, bool value) {
        if (value) mask |= field_bit(field);
    };
    changed(CheckpointCellField3D::anchor, current.anchor != parent.anchor);
    changed(CheckpointCellField3D::parent_uid,
            current.parent_uid != parent.parent_uid);
    changed(CheckpointCellField3D::clone_id,
            current.clone_id != parent.clone_id);
    changed(CheckpointCellField3D::type, current.type != parent.type);
    changed(CheckpointCellField3D::stage, current.stage != parent.stage);
    changed(CheckpointCellField3D::viability,
            current.viability != parent.viability);
    changed(CheckpointCellField3D::flags, current.flags != parent.flags);
    changed(CheckpointCellField3D::last_direction,
            current.last_direction != parent.last_direction);
    changed(CheckpointCellField3D::inherent_growth_rate,
            current.inherent_growth_rate != parent.inherent_growth_rate);
    changed(CheckpointCellField3D::density_growth_rate,
            current.density_growth_rate != parent.density_growth_rate);
    changed(CheckpointCellField3D::migration_rate,
            current.migration_rate != parent.migration_rate);
    changed(CheckpointCellField3D::normal_migration_rate,
            current.normal_migration_rate != parent.normal_migration_rate);
    changed(CheckpointCellField3D::migration_activation_end_time,
            current.migration_activation_end_time !=
                parent.migration_activation_end_time);
    changed(CheckpointCellField3D::division_work_remaining,
            current.division_work_remaining != parent.division_work_remaining);
    changed(CheckpointCellField3D::next_migration_time,
            current.next_migration_time != parent.next_migration_time);
    changed(CheckpointCellField3D::next_division_time,
            current.next_division_time != parent.next_division_time);
    changed(CheckpointCellField3D::death_deadline,
            current.death_deadline != parent.death_deadline);
    changed(CheckpointCellField3D::last_update_time,
            current.last_update_time != parent.last_update_time);
    changed(CheckpointCellField3D::swap_ready_time,
            current.swap_ready_time != parent.swap_ready_time);
    changed(CheckpointCellField3D::swap_wait_state,
            current.swap_wait_state != parent.swap_wait_state);
    changed(CheckpointCellField3D::pending_swap_direction,
            current.pending_swap_direction != parent.pending_swap_direction);
    changed(CheckpointCellField3D::event_sequence,
            current.event_sequence != parent.event_sequence);
    changed(CheckpointCellField3D::migration_schedule_generation,
            current.migration_schedule_generation !=
                parent.migration_schedule_generation);
    changed(CheckpointCellField3D::division_schedule_generation,
            current.division_schedule_generation !=
                parent.division_schedule_generation);
    changed(CheckpointCellField3D::death_schedule_generation,
            current.death_schedule_generation !=
                parent.death_schedule_generation);
    return mask;
}

struct CellFieldDelta3D {
    std::vector<CellInit> births;
    std::vector<Slot> birth_slots;
    std::vector<CellUid> removed_uids;
    std::vector<CellUid> update_uids;
    std::vector<Slot> update_slots;
    std::vector<std::uint64_t> field_masks;

    std::vector<std::int32_t> x, y, z;
    std::vector<std::uint64_t> parent_uid;
    std::vector<std::uint32_t> clone_id;
    std::vector<std::uint8_t> type, stage, viability, flags, last_direction;
    std::vector<float> inherent_growth_rate, density_growth_rate,
        migration_rate, normal_migration_rate, division_work_remaining,
        death_deadline;
    std::vector<double> migration_activation_end_time, next_migration_time,
        next_division_time, last_update_time, swap_ready_time;
    std::vector<std::uint8_t> swap_wait_state, pending_swap_direction;
    std::vector<std::uint64_t> event_sequence;
    std::vector<std::uint32_t> migration_schedule_generation,
        division_schedule_generation, death_schedule_generation;
};

bool mask_has(std::uint64_t mask, CheckpointCellField3D field) noexcept {
    return (mask & field_bit(field)) != 0;
}

void append_cell_field_values(CellFieldDelta3D& delta,
                              const CellInit& cell,
                              std::uint64_t mask) {
    if (mask_has(mask, CheckpointCellField3D::anchor)) {
        delta.x.push_back(cell.anchor.x);
        delta.y.push_back(cell.anchor.y);
        delta.z.push_back(cell.anchor.z);
    }
    if (mask_has(mask, CheckpointCellField3D::parent_uid))
        delta.parent_uid.push_back(cell.parent_uid);
    if (mask_has(mask, CheckpointCellField3D::clone_id))
        delta.clone_id.push_back(cell.clone_id);
    if (mask_has(mask, CheckpointCellField3D::type))
        delta.type.push_back(static_cast<std::uint8_t>(cell.type));
    if (mask_has(mask, CheckpointCellField3D::stage))
        delta.stage.push_back(static_cast<std::uint8_t>(cell.stage));
    if (mask_has(mask, CheckpointCellField3D::viability))
        delta.viability.push_back(cell.viability);
    if (mask_has(mask, CheckpointCellField3D::flags))
        delta.flags.push_back(cell.flags);
    if (mask_has(mask, CheckpointCellField3D::last_direction))
        delta.last_direction.push_back(cell.last_direction);
    if (mask_has(mask, CheckpointCellField3D::inherent_growth_rate))
        delta.inherent_growth_rate.push_back(cell.inherent_growth_rate);
    if (mask_has(mask, CheckpointCellField3D::density_growth_rate))
        delta.density_growth_rate.push_back(cell.density_growth_rate);
    if (mask_has(mask, CheckpointCellField3D::migration_rate))
        delta.migration_rate.push_back(cell.migration_rate);
    if (mask_has(mask, CheckpointCellField3D::normal_migration_rate))
        delta.normal_migration_rate.push_back(cell.normal_migration_rate);
    if (mask_has(mask, CheckpointCellField3D::migration_activation_end_time))
        delta.migration_activation_end_time.push_back(
            cell.migration_activation_end_time);
    if (mask_has(mask, CheckpointCellField3D::division_work_remaining))
        delta.division_work_remaining.push_back(cell.division_work_remaining);
    if (mask_has(mask, CheckpointCellField3D::next_migration_time))
        delta.next_migration_time.push_back(cell.next_migration_time);
    if (mask_has(mask, CheckpointCellField3D::next_division_time))
        delta.next_division_time.push_back(cell.next_division_time);
    if (mask_has(mask, CheckpointCellField3D::death_deadline))
        delta.death_deadline.push_back(
            static_cast<float>(cell.death_deadline));
    if (mask_has(mask, CheckpointCellField3D::last_update_time))
        delta.last_update_time.push_back(cell.last_update_time);
    if (mask_has(mask, CheckpointCellField3D::swap_ready_time))
        delta.swap_ready_time.push_back(cell.swap_ready_time);
    if (mask_has(mask, CheckpointCellField3D::swap_wait_state))
        delta.swap_wait_state.push_back(cell.swap_wait_state);
    if (mask_has(mask, CheckpointCellField3D::pending_swap_direction))
        delta.pending_swap_direction.push_back(cell.pending_swap_direction);
    if (mask_has(mask, CheckpointCellField3D::event_sequence))
        delta.event_sequence.push_back(cell.event_sequence);
    if (mask_has(mask, CheckpointCellField3D::migration_schedule_generation))
        delta.migration_schedule_generation.push_back(
            cell.migration_schedule_generation);
    if (mask_has(mask, CheckpointCellField3D::division_schedule_generation))
        delta.division_schedule_generation.push_back(
            cell.division_schedule_generation);
    if (mask_has(mask, CheckpointCellField3D::death_schedule_generation))
        delta.death_schedule_generation.push_back(
            cell.death_schedule_generation);
}

bool same_lineage(const LineageEdge& lhs, const LineageEdge& rhs) noexcept {
    return lhs.birth_time == rhs.birth_time && lhs.child_uid == rhs.child_uid &&
           lhs.parent_uid == rhs.parent_uid && lhs.clone_id == rhs.clone_id &&
           lhs.type == rhs.type;
}

void validate_snapshot_order(const CheckpointSnapshotView3D& snapshot,
                             const char* name) {
    if (snapshot.cells.size() != snapshot.cell_slots.size()) {
        throw std::invalid_argument(std::string(name) +
                                    " cell rows and slots are misaligned");
    }
    CellUid previous = 0;
    for (std::size_t index = 0; index < snapshot.cells.size(); ++index) {
        const CellUid uid = snapshot.cells[index].uid;
        if (uid == 0 || (index != 0 && uid <= previous)) {
            throw std::invalid_argument(std::string(name) +
                                        " cells must be sorted by unique uid");
        }
        previous = uid;
    }
}

template <class T>
void require_unique_nonzero_ids(const std::vector<T>& ids, const std::string& name) {
    std::unordered_set<T> unique;
    unique.reserve(ids.size());
    for (const T id : ids) {
        if (id == 0 || !unique.insert(id).second) {
            throw std::runtime_error("checkpoint contains duplicate or zero " + name);
        }
    }
}

void validate_cell(const CellInit& cell) {
    if (!valid_cell_type(static_cast<std::uint8_t>(cell.type)) ||
        !valid_cell_stage(static_cast<std::uint8_t>(cell.stage)) ||
        cell.viability > 1U ||
        (cell.flags & static_cast<std::uint8_t>(~(kMigrationActive | kDirtyDensity))) != 0 ||
        cell.last_direction > 26 || cell.pending_swap_direction > 26) {
        throw std::runtime_error("checkpoint contains an invalid cell enum, flag, or direction");
    }
    require_finite(cell.inherent_growth_rate, "cells.inherent_growth_rate");
    require_finite(cell.density_growth_rate, "cells.density_growth_rate");
    require_finite(cell.migration_rate, "cells.migration_rate");
    require_finite(cell.normal_migration_rate, "cells.normal_migration_rate");
    require_finite(cell.division_work_remaining,
                   "cells.division_work_remaining");
    if (!(cell.inherent_growth_rate > 0.0F) || cell.migration_rate < 0.0F ||
        cell.normal_migration_rate < 0.0F ||
        cell.division_work_remaining < 0.0F) {
        throw std::runtime_error("checkpoint contains an invalid cell rate");
    }
    require_finite_nonnegative(cell.next_migration_time, "cells.next_migration_time");
    require_finite_nonnegative(cell.migration_activation_end_time,
                               "cells.migration_activation_end_time");
    require_finite_nonnegative(cell.next_division_time, "cells.next_division_time");
    require_finite_nonnegative(cell.death_deadline, "cells.death_deadline");
    require_finite_nonnegative(cell.last_update_time, "cells.last_update_time");
    require_finite_nonnegative(cell.swap_ready_time, "cells.swap_ready_time");
}

void validate_process(const AngiogenesisProcessState3D& state,
                      double checkpoint_time) {
    require_finite_nonnegative(state.next_seed_time_hours,
                               "angiogenesis.process.next_seed_time_hours");
    require_finite_nonnegative(state.eligibility_started_hours,
                               "angiogenesis.process.eligibility_started_hours");
    require_finite_nonnegative(state.accumulated_eligible_hours,
                               "angiogenesis.process.accumulated_eligible_hours");
    require_finite_nonnegative(state.remaining_hazard,
                               "angiogenesis.process.remaining_hazard");
    require_finite_nonnegative(state.hazard_last_update_hours,
                               "angiogenesis.process.hazard_last_update_hours");
    require_finite_nonnegative(state.hazard_not_before_hours,
                               "angiogenesis.process.hazard_not_before_hours");
    require_finite_nonnegative(state.current_rate_sites_per_30_days,
                               "angiogenesis.process.current_rate_sites_per_30_days");
    if (!std::isfinite(state.current_density_stress) ||
        state.current_density_stress < 0.0 || state.current_density_stress > 1.0) {
        throw std::runtime_error(
            "checkpoint angiogenesis density stress is outside [0,1]");
    }
    if (!state.eligible && state.next_seed_time_hours != 0.0) {
        throw std::runtime_error("ineligible checkpoint angiogenesis process has a pending event");
    }
    if (state.eligible &&
        (!(state.remaining_hazard > 0.0) ||
         state.eligibility_started_hours > checkpoint_time ||
         state.hazard_last_update_hours > checkpoint_time ||
         state.hazard_not_before_hours < state.eligibility_started_hours ||
         ((state.current_rate_sites_per_30_days > 0.0) &&
          !(state.next_seed_time_hours > checkpoint_time)) ||
         ((state.current_rate_sites_per_30_days == 0.0) &&
          state.next_seed_time_hours != 0.0))) {
        throw std::runtime_error("eligible checkpoint angiogenesis timing is invalid");
    }
    if (state.committed_roots > state.attempted_events ||
        state.rejected_events != state.attempted_events - state.committed_roots) {
        throw std::runtime_error("checkpoint angiogenesis process counters are inconsistent");
    }
}

void write_stats(H5::Group& group, const SimulationStats3D& stats) {
    write_scalar_attribute(group, "migration_attempts", stats.migration_attempts);
    write_scalar_attribute(group, "migration_commits", stats.migration_commits);
    write_scalar_attribute(group, "migration_swap_waits",
                           stats.migration_swap_waits);
    write_scalar_attribute(group, "migration_swap_attempts",
                           stats.migration_swap_attempts);
    write_scalar_attribute(group, "migration_swap_commits",
                           stats.migration_swap_commits);
    write_scalar_attribute(group, "migration_swap_rejections",
                           stats.migration_swap_rejections);
    write_scalar_attribute(group, "divisions", stats.divisions);
    write_scalar_attribute(group, "deaths", stats.deaths);
    write_scalar_attribute(group, "conflict_rejections", stats.conflict_rejections);
    write_scalar_attribute(group, "angiogenesis_seed_attempts", stats.angiogenesis_seed_attempts);
    write_scalar_attribute(group, "angiogenesis_roots", stats.angiogenesis_roots);
    write_scalar_attribute(group, "angiogenesis_seed_rejections",
                           stats.angiogenesis_seed_rejections);
    write_scalar_attribute(group, "vessel_growth_attempts", stats.vessel_growth_attempts);
    write_scalar_attribute(group, "vessel_growth_commits", stats.vessel_growth_commits);
    write_scalar_attribute(group, "vessel_anastomoses", stats.vessel_anastomoses);
    write_scalar_attribute(group, "vascular_displacements", stats.vascular_displacements);
}

SimulationStats3D read_stats(H5::Group& group, bool has_swap_stats) {
    SimulationStats3D stats;
    stats.migration_attempts = read_scalar_attribute<std::uint64_t>(group, "migration_attempts");
    stats.migration_commits = read_scalar_attribute<std::uint64_t>(group, "migration_commits");
    if (has_swap_stats) {
        stats.migration_swap_waits =
            read_scalar_attribute<std::uint64_t>(
                group, "migration_swap_waits");
        stats.migration_swap_attempts =
            read_scalar_attribute<std::uint64_t>(
                group, "migration_swap_attempts");
        stats.migration_swap_commits =
            read_scalar_attribute<std::uint64_t>(
                group, "migration_swap_commits");
        stats.migration_swap_rejections =
            read_scalar_attribute<std::uint64_t>(
                group, "migration_swap_rejections");
    }
    stats.divisions = read_scalar_attribute<std::uint64_t>(group, "divisions");
    stats.deaths = read_scalar_attribute<std::uint64_t>(group, "deaths");
    stats.conflict_rejections =
        read_scalar_attribute<std::uint64_t>(group, "conflict_rejections");
    stats.angiogenesis_seed_attempts =
        read_scalar_attribute<std::uint64_t>(group, "angiogenesis_seed_attempts");
    stats.angiogenesis_roots =
        read_scalar_attribute<std::uint64_t>(group, "angiogenesis_roots");
    stats.angiogenesis_seed_rejections =
        read_scalar_attribute<std::uint64_t>(group, "angiogenesis_seed_rejections");
    stats.vessel_growth_attempts =
        read_scalar_attribute<std::uint64_t>(group, "vessel_growth_attempts");
    stats.vessel_growth_commits =
        read_scalar_attribute<std::uint64_t>(group, "vessel_growth_commits");
    stats.vessel_anastomoses =
        read_scalar_attribute<std::uint64_t>(group, "vessel_anastomoses");
    stats.vascular_displacements =
        read_scalar_attribute<std::uint64_t>(group, "vascular_displacements");
    return stats;
}

void write_process(H5::Group& group, const AngiogenesisProcessState3D& state) {
    write_bool_attribute(group, "eligible", state.eligible);
    write_scalar_attribute(group, "next_seed_time_hours", state.next_seed_time_hours);
    write_scalar_attribute(group, "eligibility_started_hours", state.eligibility_started_hours);
    write_scalar_attribute(group, "accumulated_eligible_hours", state.accumulated_eligible_hours);
    write_scalar_attribute(group, "remaining_hazard", state.remaining_hazard);
    write_scalar_attribute(group, "hazard_last_update_hours", state.hazard_last_update_hours);
    write_scalar_attribute(group, "hazard_not_before_hours", state.hazard_not_before_hours);
    write_scalar_attribute(group, "current_rate_sites_per_30_days",
                           state.current_rate_sites_per_30_days);
    write_scalar_attribute(group, "current_density_stress",
                           state.current_density_stress);
    write_scalar_attribute(group, "event_sequence", state.event_sequence);
    write_scalar_attribute(group, "schedule_generation", state.schedule_generation);
    write_scalar_attribute(group, "attempted_events", state.attempted_events);
    write_scalar_attribute(group, "committed_roots", state.committed_roots);
    write_scalar_attribute(group, "rejected_events", state.rejected_events);
}

AngiogenesisProcessState3D read_process(H5::Group& group) {
    AngiogenesisProcessState3D state;
    state.eligible = read_bool_attribute(group, "eligible");
    state.next_seed_time_hours =
        read_scalar_attribute<double>(group, "next_seed_time_hours");
    state.eligibility_started_hours =
        read_scalar_attribute<double>(group, "eligibility_started_hours");
    state.accumulated_eligible_hours =
        read_scalar_attribute<double>(group, "accumulated_eligible_hours");
    state.remaining_hazard =
        read_scalar_attribute<double>(group, "remaining_hazard");
    state.hazard_last_update_hours =
        read_scalar_attribute<double>(group, "hazard_last_update_hours");
    state.hazard_not_before_hours =
        read_scalar_attribute<double>(group, "hazard_not_before_hours");
    state.current_rate_sites_per_30_days =
        read_scalar_attribute<double>(group, "current_rate_sites_per_30_days");
    state.current_density_stress =
        read_scalar_attribute<double>(group, "current_density_stress");
    state.event_sequence = read_scalar_attribute<std::uint64_t>(group, "event_sequence");
    state.schedule_generation =
        read_scalar_attribute<std::uint32_t>(group, "schedule_generation");
    state.attempted_events = read_scalar_attribute<std::uint64_t>(group, "attempted_events");
    state.committed_roots = read_scalar_attribute<std::uint64_t>(group, "committed_roots");
    state.rejected_events = read_scalar_attribute<std::uint64_t>(group, "rejected_events");
    return state;
}

void write_lesion_core_identity(
    H5::Group& group,
    const std::vector<LesionCoreIdentity3D>& identity) {
    write_vector(group, "block_x",
                 column<std::int32_t>(identity,
                                      [](const auto& entry) {
                                          return entry.block.x;
                                      }));
    write_vector(group, "block_y",
                 column<std::int32_t>(identity,
                                      [](const auto& entry) {
                                          return entry.block.y;
                                      }));
    write_vector(group, "block_z",
                 column<std::int32_t>(identity,
                                      [](const auto& entry) {
                                          return entry.block.z;
                                      }));
    write_vector(group, "lesion_id",
                 column<std::uint64_t>(identity,
                                       [](const auto& entry) {
                                           return entry.lesion_id;
                                       }));
}

std::vector<LesionCoreIdentity3D> read_lesion_core_identity(
    H5::Group& group, LesionId next_lesion_id) {
    const auto x = read_vector<std::int32_t>(group, "block_x");
    const auto y = read_vector<std::int32_t>(group, "block_y");
    const auto z = read_vector<std::int32_t>(group, "block_z");
    const auto lesion_id = read_vector<std::uint64_t>(group, "lesion_id");
    require_equal_sizes("lesion core identity", lesion_id.size(),
                        {{"block_x", x.size()}, {"block_y", y.size()},
                         {"block_z", z.size()}});

    std::vector<LesionCoreIdentity3D> identity(lesion_id.size());
    Vec3i previous{};
    for (std::size_t index = 0; index < identity.size(); ++index) {
        const Vec3i block{x[index], y[index], z[index]};
        if (lesion_id[index] == kNoLesionId ||
            lesion_id[index] >= next_lesion_id) {
            throw std::runtime_error(
                "checkpoint lesion core identity contains an invalid lesion ID");
        }
        if (index != 0 && !(previous < block)) {
            throw std::runtime_error(
                "checkpoint lesion core identity must be sorted by unique block");
        }
        identity[index] = {block, lesion_id[index]};
        previous = block;
    }
    return identity;
}

void write_lesion_dirty_blocks(
    H5::Group& group,
    const std::vector<LesionDirtyBlockState3D>& blocks) {
    write_vector(group, "block_x",
                 column<std::int32_t>(blocks, [](const auto& block) {
                     return block.block.x;
                 }));
    write_vector(group, "block_y",
                 column<std::int32_t>(blocks, [](const auto& block) {
                     return block.block.y;
                 }));
    write_vector(group, "block_z",
                 column<std::int32_t>(blocks, [](const auto& block) {
                     return block.block.z;
                 }));
    write_vector(group, "exists",
                 column<std::uint8_t>(blocks, [](const auto& block) {
                     return static_cast<std::uint8_t>(block.exists ? 1U : 0U);
                 }));
    write_vector(group, "cell_count",
                 column<std::uint64_t>(blocks, [](const auto& block) {
                     return block.cell_count;
                 }));
    write_vector(group, "occupied_voxel_count",
                 column<std::uint64_t>(blocks, [](const auto& block) {
                     return block.occupied_voxel_count;
                 }));
    write_vector(group, "biological_volume",
                 column<double>(blocks, [](const auto& block) {
                     return block.biological_volume;
                 }));
    write_vector(group, "cell_coordinate_sum_x",
                 column<std::int64_t>(blocks, [](const auto& block) {
                     return block.cell_coordinate_sum_x;
                 }));
    write_vector(group, "cell_coordinate_sum_y",
                 column<std::int64_t>(blocks, [](const auto& block) {
                     return block.cell_coordinate_sum_y;
                 }));
    write_vector(group, "cell_coordinate_sum_z",
                 column<std::int64_t>(blocks, [](const auto& block) {
                     return block.cell_coordinate_sum_z;
                 }));
    write_vector(group, "occupied_coordinate_sum_x",
                 column<std::int64_t>(blocks, [](const auto& block) {
                     return block.occupied_coordinate_sum_x;
                 }));
    write_vector(group, "occupied_coordinate_sum_y",
                 column<std::int64_t>(blocks, [](const auto& block) {
                     return block.occupied_coordinate_sum_y;
                 }));
    write_vector(group, "occupied_coordinate_sum_z",
                 column<std::int64_t>(blocks, [](const auto& block) {
                     return block.occupied_coordinate_sum_z;
                 }));
}

std::int32_t clamp_lattice_coordinate(std::int64_t value) noexcept {
    return static_cast<std::int32_t>(std::clamp(
        value,
        static_cast<std::int64_t>(std::numeric_limits<std::int32_t>::min()),
        static_cast<std::int64_t>(std::numeric_limits<std::int32_t>::max())));
}

void validate_coordinate_sum(std::int64_t sum, std::uint64_t count,
                             std::int32_t minimum, std::int32_t maximum,
                             const char* name) {
    if (count == 0) {
        if (sum != 0) {
            throw std::runtime_error(std::string("checkpoint ") + name +
                                     " must be zero for an empty count");
        }
        return;
    }
    const std::int64_t signed_count = static_cast<std::int64_t>(count);
    const std::int64_t lower =
        static_cast<std::int64_t>(minimum) * signed_count;
    const std::int64_t upper =
        static_cast<std::int64_t>(maximum) * signed_count;
    if (sum < lower || sum > upper) {
        throw std::runtime_error(std::string("checkpoint ") + name +
                                 " lies outside its lesion block");
    }
}

std::vector<LesionDirtyBlockState3D> read_lesion_dirty_blocks(
    H5::Group& group, const Model3DConfig& config) {
    const auto x = read_vector<std::int32_t>(group, "block_x");
    const auto y = read_vector<std::int32_t>(group, "block_y");
    const auto z = read_vector<std::int32_t>(group, "block_z");
    const auto exists = read_vector<std::uint8_t>(group, "exists");
    const auto cell_count = read_vector<std::uint64_t>(group, "cell_count");
    const auto occupied_count =
        read_vector<std::uint64_t>(group, "occupied_voxel_count");
    const auto biological_volume =
        read_vector<double>(group, "biological_volume");
    const auto cell_sum_x =
        read_vector<std::int64_t>(group, "cell_coordinate_sum_x");
    const auto cell_sum_y =
        read_vector<std::int64_t>(group, "cell_coordinate_sum_y");
    const auto cell_sum_z =
        read_vector<std::int64_t>(group, "cell_coordinate_sum_z");
    const auto occupied_sum_x =
        read_vector<std::int64_t>(group, "occupied_coordinate_sum_x");
    const auto occupied_sum_y =
        read_vector<std::int64_t>(group, "occupied_coordinate_sum_y");
    const auto occupied_sum_z =
        read_vector<std::int64_t>(group, "occupied_coordinate_sum_z");
    const std::size_t count = x.size();
    require_equal_sizes(
        "lesion dirty block", count,
        {{"block_y", y.size()}, {"block_z", z.size()},
         {"exists", exists.size()}, {"cell_count", cell_count.size()},
         {"occupied_voxel_count", occupied_count.size()},
         {"biological_volume", biological_volume.size()},
         {"cell_coordinate_sum_x", cell_sum_x.size()},
         {"cell_coordinate_sum_y", cell_sum_y.size()},
         {"cell_coordinate_sum_z", cell_sum_z.size()},
         {"occupied_coordinate_sum_x", occupied_sum_x.size()},
         {"occupied_coordinate_sum_y", occupied_sum_y.size()},
         {"occupied_coordinate_sum_z", occupied_sum_z.size()}});

    const std::uint64_t edge = static_cast<std::uint64_t>(
        config.angiogenesis.lesion_block_edge);
    if (edge == 0 ||
        edge > std::numeric_limits<std::uint64_t>::max() / edge ||
        edge * edge > std::numeric_limits<std::uint64_t>::max() / edge) {
        throw std::runtime_error(
            "checkpoint lesion block capacity is invalid");
    }
    const std::uint64_t voxel_capacity = edge * edge * edge;
    std::vector<LesionDirtyBlockState3D> blocks(count);
    Vec3i previous{};
    for (std::size_t index = 0; index < count; ++index) {
        const Vec3i coordinate{x[index], y[index], z[index]};
        if (index != 0 && !(previous < coordinate)) {
            throw std::runtime_error(
                "checkpoint lesion dirty blocks must be sorted and unique");
        }
        if (exists[index] > 1U ||
            !std::isfinite(biological_volume[index]) ||
            biological_volume[index] < 0.0 ||
            cell_count[index] > std::numeric_limits<Slot>::max() ||
            occupied_count[index] > voxel_capacity) {
            throw std::runtime_error(
                "checkpoint lesion dirty block contains an invalid value");
        }
        const bool payload_is_zero =
            cell_count[index] == 0 && occupied_count[index] == 0 &&
            biological_volume[index] == 0.0 && cell_sum_x[index] == 0 &&
            cell_sum_y[index] == 0 && cell_sum_z[index] == 0 &&
            occupied_sum_x[index] == 0 && occupied_sum_y[index] == 0 &&
            occupied_sum_z[index] == 0;
        if ((exists[index] == 0 && !payload_is_zero) ||
            (exists[index] != 0 && cell_count[index] == 0 &&
             occupied_count[index] == 0) ||
            (cell_count[index] == 0 && biological_volume[index] != 0.0) ||
            (cell_count[index] != 0 && biological_volume[index] <= 0.0)) {
            throw std::runtime_error(
                "checkpoint lesion dirty block payload is inconsistent");
        }

        const std::int64_t coordinate_x = coordinate.x;
        const std::int64_t coordinate_y = coordinate.y;
        const std::int64_t coordinate_z = coordinate.z;
        const std::int64_t signed_edge =
            config.angiogenesis.lesion_block_edge;
        const Vec3i minimum{
            clamp_lattice_coordinate(coordinate_x * signed_edge),
            clamp_lattice_coordinate(coordinate_y * signed_edge),
            clamp_lattice_coordinate(coordinate_z * signed_edge)};
        const Vec3i maximum{
            clamp_lattice_coordinate(coordinate_x * signed_edge +
                                     signed_edge - 1),
            clamp_lattice_coordinate(coordinate_y * signed_edge +
                                     signed_edge - 1),
            clamp_lattice_coordinate(coordinate_z * signed_edge +
                                     signed_edge - 1)};
        validate_coordinate_sum(cell_sum_x[index], cell_count[index],
                                minimum.x, maximum.x,
                                "lesion dirty cell x sum");
        validate_coordinate_sum(cell_sum_y[index], cell_count[index],
                                minimum.y, maximum.y,
                                "lesion dirty cell y sum");
        validate_coordinate_sum(cell_sum_z[index], cell_count[index],
                                minimum.z, maximum.z,
                                "lesion dirty cell z sum");
        validate_coordinate_sum(occupied_sum_x[index], occupied_count[index],
                                minimum.x, maximum.x,
                                "lesion dirty occupied x sum");
        validate_coordinate_sum(occupied_sum_y[index], occupied_count[index],
                                minimum.y, maximum.y,
                                "lesion dirty occupied y sum");
        validate_coordinate_sum(occupied_sum_z[index], occupied_count[index],
                                minimum.z, maximum.z,
                                "lesion dirty occupied z sum");

        blocks[index] = {
            .block = coordinate,
            .exists = exists[index] != 0,
            .cell_count = cell_count[index],
            .occupied_voxel_count = occupied_count[index],
            .biological_volume = biological_volume[index],
            .cell_coordinate_sum_x = cell_sum_x[index],
            .cell_coordinate_sum_y = cell_sum_y[index],
            .cell_coordinate_sum_z = cell_sum_z[index],
            .occupied_coordinate_sum_x = occupied_sum_x[index],
            .occupied_coordinate_sum_y = occupied_sum_y[index],
            .occupied_coordinate_sum_z = occupied_sum_z[index],
        };
        previous = coordinate;
    }
    return blocks;
}

void write_lesion_processes(
    H5::Group& group,
    const std::vector<LesionAngiogenesisState3D>& processes) {
    write_vector(group, "lesion_id",
                 column<std::uint64_t>(processes,
                                       [](const auto& entry) {
                                           return entry.lesion_id;
                                       }));
    write_vector(group, "eligible",
                 column<std::uint8_t>(processes,
                                      [](const auto& entry) {
                                          return static_cast<std::uint8_t>(
                                              entry.process.eligible ? 1U : 0U);
                                      }));
    write_vector(group, "next_seed_time_hours",
                 column<double>(processes,
                                [](const auto& entry) {
                                    return entry.process.next_seed_time_hours;
                                }));
    write_vector(group, "eligibility_started_hours",
                 column<double>(processes,
                                [](const auto& entry) {
                                    return entry.process.eligibility_started_hours;
                                }));
    write_vector(group, "accumulated_eligible_hours",
                 column<double>(processes,
                                [](const auto& entry) {
                                    return entry.process.accumulated_eligible_hours;
                                }));
    write_vector(group, "remaining_hazard",
                 column<double>(processes, [](const auto& entry) {
                     return entry.process.remaining_hazard;
                 }));
    write_vector(group, "hazard_last_update_hours",
                 column<double>(processes, [](const auto& entry) {
                     return entry.process.hazard_last_update_hours;
                 }));
    write_vector(group, "hazard_not_before_hours",
                 column<double>(processes, [](const auto& entry) {
                     return entry.process.hazard_not_before_hours;
                 }));
    write_vector(group, "current_rate_sites_per_30_days",
                 column<double>(processes, [](const auto& entry) {
                     return entry.process.current_rate_sites_per_30_days;
                 }));
    write_vector(group, "current_density_stress",
                 column<double>(processes, [](const auto& entry) {
                     return entry.process.current_density_stress;
                 }));
    write_vector(group, "event_sequence",
                 column<std::uint64_t>(processes,
                                       [](const auto& entry) {
                                           return entry.process.event_sequence;
                                       }));
    write_vector(group, "schedule_generation",
                 column<std::uint32_t>(processes,
                                       [](const auto& entry) {
                                           return entry.process.schedule_generation;
                                       }));
    write_vector(group, "attempted_events",
                 column<std::uint64_t>(processes,
                                       [](const auto& entry) {
                                           return entry.process.attempted_events;
                                       }));
    write_vector(group, "committed_roots",
                 column<std::uint64_t>(processes,
                                       [](const auto& entry) {
                                           return entry.process.committed_roots;
                                       }));
    write_vector(group, "rejected_events",
                 column<std::uint64_t>(processes,
                                       [](const auto& entry) {
                                           return entry.process.rejected_events;
                                       }));
}

std::vector<LesionAngiogenesisState3D> read_lesion_processes(
    H5::Group& group, LesionId next_lesion_id, double checkpoint_time) {
    const auto lesion_id = read_vector<std::uint64_t>(group, "lesion_id");
    const auto eligible = read_vector<std::uint8_t>(group, "eligible");
    const auto next_seed_time =
        read_vector<double>(group, "next_seed_time_hours");
    const auto eligibility_started =
        read_vector<double>(group, "eligibility_started_hours");
    const auto accumulated =
        read_vector<double>(group, "accumulated_eligible_hours");
    const auto remaining_hazard =
        read_vector<double>(group, "remaining_hazard");
    const auto hazard_last_update =
        read_vector<double>(group, "hazard_last_update_hours");
    const auto hazard_not_before =
        read_vector<double>(group, "hazard_not_before_hours");
    const auto current_rate =
        read_vector<double>(group, "current_rate_sites_per_30_days");
    const auto current_stress =
        read_vector<double>(group, "current_density_stress");
    const auto event_sequence =
        read_vector<std::uint64_t>(group, "event_sequence");
    const auto generation =
        read_vector<std::uint32_t>(group, "schedule_generation");
    const auto attempted =
        read_vector<std::uint64_t>(group, "attempted_events");
    const auto committed =
        read_vector<std::uint64_t>(group, "committed_roots");
    const auto rejected =
        read_vector<std::uint64_t>(group, "rejected_events");
    const std::size_t count = lesion_id.size();
    require_equal_sizes(
        "lesion angiogenesis process", count,
        {{"eligible", eligible.size()},
         {"next_seed_time_hours", next_seed_time.size()},
         {"eligibility_started_hours", eligibility_started.size()},
         {"accumulated_eligible_hours", accumulated.size()},
         {"remaining_hazard", remaining_hazard.size()},
         {"hazard_last_update_hours", hazard_last_update.size()},
         {"hazard_not_before_hours", hazard_not_before.size()},
         {"current_rate_sites_per_30_days", current_rate.size()},
         {"current_density_stress", current_stress.size()},
         {"event_sequence", event_sequence.size()},
         {"schedule_generation", generation.size()},
         {"attempted_events", attempted.size()},
         {"committed_roots", committed.size()},
         {"rejected_events", rejected.size()}});

    std::vector<LesionAngiogenesisState3D> processes(count);
    LesionId previous = kNoLesionId;
    for (std::size_t index = 0; index < count; ++index) {
        if (lesion_id[index] == kNoLesionId ||
            lesion_id[index] >= next_lesion_id ||
            (index != 0 && lesion_id[index] <= previous) ||
            eligible[index] > 1U) {
            throw std::runtime_error(
                "checkpoint lesion process IDs/boolean must be sorted and valid");
        }
        AngiogenesisProcessState3D process;
        process.eligible = eligible[index] != 0;
        process.next_seed_time_hours = next_seed_time[index];
        process.eligibility_started_hours = eligibility_started[index];
        process.accumulated_eligible_hours = accumulated[index];
        process.remaining_hazard = remaining_hazard[index];
        process.hazard_last_update_hours = hazard_last_update[index];
        process.hazard_not_before_hours = hazard_not_before[index];
        process.current_rate_sites_per_30_days = current_rate[index];
        process.current_density_stress = current_stress[index];
        process.event_sequence = event_sequence[index];
        process.schedule_generation = generation[index];
        process.attempted_events = attempted[index];
        process.committed_roots = committed[index];
        process.rejected_events = rejected[index];
        validate_process(process, checkpoint_time);
        processes[index] = {lesion_id[index], process};
        previous = lesion_id[index];
    }
    return processes;
}

void write_lesion_source_ownership(
    H5::Group& group,
    const std::vector<LesionSourceOwnership3D>& ownership) {
    write_vector(group, "source_lesion_id",
                 column<std::uint64_t>(ownership, [](const auto& entry) {
                     return entry.source_lesion_id;
                 }));
    write_vector(group, "current_lesion_id",
                 column<std::uint64_t>(ownership, [](const auto& entry) {
                     return entry.current_lesion_id;
                 }));
}

std::vector<LesionSourceOwnership3D> read_lesion_source_ownership(
    H5::Group& group, LesionId next_lesion_id) {
    const auto source =
        read_vector<std::uint64_t>(group, "source_lesion_id");
    const auto current =
        read_vector<std::uint64_t>(group, "current_lesion_id");
    require_equal_sizes("lesion source ownership", source.size(),
                        {{"current_lesion_id", current.size()}});
    std::vector<LesionSourceOwnership3D> result(source.size());
    LesionId previous = kNoLesionId;
    for (std::size_t index = 0; index < source.size(); ++index) {
        if (source[index] == kNoLesionId || source[index] >= next_lesion_id ||
            (index != 0 && source[index] <= previous) ||
            current[index] >= next_lesion_id ||
            source[index] == current[index]) {
            throw std::runtime_error(
                "checkpoint lesion source ownership IDs must be sorted and valid");
        }
        result[index] = {source[index], current[index]};
        previous = source[index];
    }
    return result;
}

void write_cells(H5::Group& group,
                 std::span<const CellInit> cells,
                 std::span<const Slot> slots,
                 std::size_t slot_count,
                 std::span<const Slot> free_slots) {
    if (cells.size() != slots.size()) {
        throw std::logic_error("cell checkpoint rows and slots are misaligned");
    }
    write_scalar_attribute(
        group, "slot_count", static_cast<std::uint64_t>(slot_count));
    write_vector<Slot>(group, "slot", slots);
    write_vector<Slot>(group, "free_slots", free_slots);
    write_vector(group, "x", column<std::int32_t>(cells, [](const auto& c) { return c.anchor.x; }));
    write_vector(group, "y", column<std::int32_t>(cells, [](const auto& c) { return c.anchor.y; }));
    write_vector(group, "z", column<std::int32_t>(cells, [](const auto& c) { return c.anchor.z; }));
    write_vector(group, "uid", column<std::uint64_t>(cells, [](const auto& c) { return c.uid; }));
    write_vector(group, "parent_uid", column<std::uint64_t>(cells, [](const auto& c) { return c.parent_uid; }));
    write_vector(group, "clone_id", column<std::uint32_t>(cells, [](const auto& c) { return c.clone_id; }));
    write_vector(group, "type", column<std::uint8_t>(cells, [](const auto& c) { return static_cast<std::uint8_t>(c.type); }));
    write_vector(group, "stage", column<std::uint8_t>(cells, [](const auto& c) { return static_cast<std::uint8_t>(c.stage); }));
    write_vector(group, "viability", column<std::uint8_t>(cells, [](const auto& c) { return c.viability; }));
    write_vector(group, "flags", column<std::uint8_t>(cells, [](const auto& c) { return c.flags; }));
    write_vector(group, "last_direction", column<std::uint8_t>(cells, [](const auto& c) { return c.last_direction; }));
    write_vector(group, "inherent_growth_rate", column<float>(cells, [](const auto& c) { return c.inherent_growth_rate; }));
    write_vector(group, "density_growth_rate", column<float>(cells, [](const auto& c) { return c.density_growth_rate; }));
    write_vector(group, "migration_rate", column<float>(cells, [](const auto& c) { return c.migration_rate; }));
    write_vector(group, "normal_migration_rate", column<float>(cells, [](const auto& c) { return c.normal_migration_rate; }));
    write_vector(group, "migration_activation_end_time", column<double>(cells, [](const auto& c) { return c.migration_activation_end_time; }));
    write_vector(group, "division_work_remaining", column<float>(cells, [](const auto& c) { return c.division_work_remaining; }));
    write_vector(group, "next_migration_time", column<double>(cells, [](const auto& c) { return c.next_migration_time; }));
    write_vector(group, "next_division_time", column<double>(cells, [](const auto& c) { return c.next_division_time; }));
    write_vector(group, "death_deadline", column<float>(cells, [](const auto& c) { return static_cast<float>(c.death_deadline); }));
    write_vector(group, "last_update_time", column<double>(cells, [](const auto& c) { return c.last_update_time; }));
    write_vector(group, "swap_ready_time", column<double>(cells, [](const auto& c) { return c.swap_ready_time; }));
    write_vector(group, "swap_wait_state", column<std::uint8_t>(cells, [](const auto& c) { return c.swap_wait_state; }));
    write_vector(group, "pending_swap_direction", column<std::uint8_t>(cells, [](const auto& c) { return c.pending_swap_direction; }));
    write_vector(group, "event_sequence", column<std::uint64_t>(cells, [](const auto& c) { return c.event_sequence; }));
    write_vector(group, "migration_schedule_generation", column<std::uint32_t>(cells, [](const auto& c) { return c.migration_schedule_generation; }));
    write_vector(group, "division_schedule_generation", column<std::uint32_t>(cells, [](const auto& c) { return c.division_schedule_generation; }));
    write_vector(group, "death_schedule_generation", column<std::uint32_t>(cells, [](const auto& c) { return c.death_schedule_generation; }));
}

std::vector<CellInit> read_cells(H5::Group& group, CellUid next_uid,
                                 const Model3DConfig& config,
                                 std::vector<Slot>& slots,
                                 std::size_t& slot_count,
                                 std::vector<Slot>& free_slots,
                                 bool has_swap_fields) {
    const std::uint64_t stored_slot_count =
        read_scalar_attribute<std::uint64_t>(group, "slot_count");
    if (stored_slot_count > static_cast<std::uint64_t>(kEmptySlot) ||
        stored_slot_count >
            static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
        throw std::runtime_error("checkpoint cell slot count is too large");
    }
    slot_count = static_cast<std::size_t>(stored_slot_count);
    slots = read_vector<Slot>(group, "slot");
    free_slots = read_vector<Slot>(group, "free_slots");
    const auto x = read_vector<std::int32_t>(group, "x");
    const auto y = read_vector<std::int32_t>(group, "y");
    const auto z = read_vector<std::int32_t>(group, "z");
    const auto uid = read_vector<std::uint64_t>(group, "uid");
    const auto parent_uid = read_vector<std::uint64_t>(group, "parent_uid");
    const auto clone_id = read_vector<std::uint32_t>(group, "clone_id");
    const auto type = read_vector<std::uint8_t>(group, "type");
    const auto stage = read_vector<std::uint8_t>(group, "stage");
    const auto viability = read_vector<std::uint8_t>(group, "viability");
    const auto flags = read_vector<std::uint8_t>(group, "flags");
    const auto last_direction = read_vector<std::uint8_t>(group, "last_direction");
    const auto inherent_growth_rate = read_vector<float>(group, "inherent_growth_rate");
    const auto density_growth_rate = read_vector<float>(group, "density_growth_rate");
    const auto migration_rate = read_vector<float>(group, "migration_rate");
    const auto normal_migration_rate = read_vector<float>(group, "normal_migration_rate");
    const auto migration_activation_end_time =
        read_vector<double>(group, "migration_activation_end_time");
    const auto division_work = read_vector<float>(group, "division_work_remaining");
    const auto next_migration_time = read_vector<double>(group, "next_migration_time");
    const auto next_division_time = read_vector<double>(group, "next_division_time");
    const auto death_deadline = read_vector<float>(group, "death_deadline");
    const auto last_update_time = read_vector<double>(group, "last_update_time");
    const auto swap_ready_time = has_swap_fields
        ? read_vector<double>(group, "swap_ready_time")
        : std::vector<double>(uid.size(), 0.0);
    const auto swap_wait_state = has_swap_fields
        ? read_vector<std::uint8_t>(group, "swap_wait_state")
        : std::vector<std::uint8_t>(uid.size(), 0U);
    const auto pending_swap_direction = has_swap_fields
        ? read_vector<std::uint8_t>(group, "pending_swap_direction")
        : std::vector<std::uint8_t>(uid.size(), kStayDirection);
    const auto event_sequence = read_vector<std::uint64_t>(group, "event_sequence");
    const auto migration_generation =
        read_vector<std::uint32_t>(group, "migration_schedule_generation");
    const auto division_generation =
        read_vector<std::uint32_t>(group, "division_schedule_generation");
    const auto death_generation =
        read_vector<std::uint32_t>(group, "death_schedule_generation");
    const std::size_t count = uid.size();
    require_equal_sizes("cell", count,
                        {{"slot", slots.size()},
                         {"x", x.size()}, {"y", y.size()}, {"z", z.size()},
                         {"parent_uid", parent_uid.size()}, {"clone_id", clone_id.size()},
                         {"type", type.size()}, {"stage", stage.size()},
                         {"viability", viability.size()}, {"flags", flags.size()},
                         {"last_direction", last_direction.size()},
                         {"inherent_growth_rate", inherent_growth_rate.size()},
                         {"density_growth_rate", density_growth_rate.size()},
                         {"migration_rate", migration_rate.size()},
                         {"normal_migration_rate", normal_migration_rate.size()},
                         {"migration_activation_end_time", migration_activation_end_time.size()},
                         {"division_work_remaining", division_work.size()},
                         {"next_migration_time", next_migration_time.size()},
                         {"next_division_time", next_division_time.size()},
                         {"death_deadline", death_deadline.size()},
                         {"last_update_time", last_update_time.size()},
                         {"swap_ready_time", swap_ready_time.size()},
                         {"swap_wait_state", swap_wait_state.size()},
                         {"pending_swap_direction", pending_swap_direction.size()},
                         {"event_sequence", event_sequence.size()},
                         {"migration_schedule_generation", migration_generation.size()},
                         {"division_schedule_generation", division_generation.size()},
                         {"death_schedule_generation", death_generation.size()}});
    require_unique_nonzero_ids(uid, "cell uid");
    if (!std::is_sorted(uid.begin(), uid.end())) {
        throw std::runtime_error("checkpoint cell rows are not sorted by uid");
    }

    std::vector<CellInit> cells(count);
    for (std::size_t index = 0; index < count; ++index) {
        CellInit& cell = cells[index];
        cell.anchor = {x[index], y[index], z[index]};
        cell.uid = uid[index];
        cell.parent_uid = parent_uid[index];
        cell.clone_id = clone_id[index];
        cell.type = static_cast<CellType>(type[index]);
        cell.stage = static_cast<CellStage>(stage[index]);
        cell.viability = viability[index];
        cell.flags = flags[index];
        cell.last_direction = last_direction[index];
        cell.inherent_growth_rate = inherent_growth_rate[index];
        cell.density_growth_rate = density_growth_rate[index];
        cell.migration_rate = migration_rate[index];
        cell.normal_migration_rate = normal_migration_rate[index];
        cell.migration_activation_end_time = migration_activation_end_time[index];
        cell.division_work_remaining = division_work[index];
        cell.next_migration_time = next_migration_time[index];
        cell.next_division_time = next_division_time[index];
        cell.death_deadline = death_deadline[index];
        cell.last_update_time = last_update_time[index];
        cell.swap_ready_time = swap_ready_time[index];
        cell.swap_wait_state = swap_wait_state[index];
        cell.pending_swap_direction = pending_swap_direction[index];
        cell.event_sequence = event_sequence[index];
        cell.migration_schedule_generation = migration_generation[index];
        cell.division_schedule_generation = division_generation[index];
        cell.death_schedule_generation = death_generation[index];
        cell.schedule_generation = std::max({migration_generation[index],
                                             division_generation[index],
                                             death_generation[index]});
        validate_cell(cell);
        const bool active =
            (cell.flags & static_cast<std::uint8_t>(kMigrationActive)) != 0;
        if ((!config.migration_activation_enabled && active) ||
            active != (cell.migration_activation_end_time > 0.0)) {
            throw std::runtime_error(
                "checkpoint migration activation flag/end time are inconsistent");
        }
        if (cell.uid >= next_uid) {
            throw std::runtime_error("checkpoint next cell uid does not exceed all cell uids");
        }
    }
    return cells;
}

std::size_t field_value_count(std::span<const std::uint64_t> masks,
                              CheckpointCellField3D field) {
    return static_cast<std::size_t>(std::count_if(
        masks.begin(), masks.end(),
        [field](std::uint64_t mask) { return mask_has(mask, field); }));
}

void write_cell_field_delta(H5::Group& cells_group,
                            const CellFieldDelta3D& delta,
                            std::size_t slot_count,
                            std::span<const Slot> free_slots) {
    write_scalar_attribute(
        cells_group, "slot_count", static_cast<std::uint64_t>(slot_count));
    write_vector<Slot>(cells_group, "free_slots", free_slots);
    write_vector(cells_group, "removed_uid", delta.removed_uids);

    H5::Group births_group = cells_group.createGroup("births");
    write_cells(births_group, delta.births, delta.birth_slots, slot_count, {});

    H5::Group updates_group = cells_group.createGroup("updates");
    write_vector(updates_group, "uid", delta.update_uids);
    write_vector(updates_group, "slot", delta.update_slots);
    write_vector(updates_group, "field_mask", delta.field_masks);
    write_vector(updates_group, "x", delta.x);
    write_vector(updates_group, "y", delta.y);
    write_vector(updates_group, "z", delta.z);
    write_vector(updates_group, "parent_uid", delta.parent_uid);
    write_vector(updates_group, "clone_id", delta.clone_id);
    write_vector(updates_group, "type", delta.type);
    write_vector(updates_group, "stage", delta.stage);
    write_vector(updates_group, "viability", delta.viability);
    write_vector(updates_group, "flags", delta.flags);
    write_vector(updates_group, "last_direction", delta.last_direction);
    write_vector(updates_group, "inherent_growth_rate",
                 delta.inherent_growth_rate);
    write_vector(updates_group, "density_growth_rate",
                 delta.density_growth_rate);
    write_vector(updates_group, "migration_rate", delta.migration_rate);
    write_vector(updates_group, "normal_migration_rate",
                 delta.normal_migration_rate);
    write_vector(updates_group, "migration_activation_end_time",
                 delta.migration_activation_end_time);
    write_vector(updates_group, "division_work_remaining",
                 delta.division_work_remaining);
    write_vector(updates_group, "next_migration_time",
                 delta.next_migration_time);
    write_vector(updates_group, "next_division_time",
                 delta.next_division_time);
    write_vector(updates_group, "death_deadline", delta.death_deadline);
    write_vector(updates_group, "last_update_time", delta.last_update_time);
    write_vector(updates_group, "swap_ready_time", delta.swap_ready_time);
    write_vector(updates_group, "swap_wait_state", delta.swap_wait_state);
    write_vector(updates_group, "pending_swap_direction",
                 delta.pending_swap_direction);
    write_vector(updates_group, "event_sequence", delta.event_sequence);
    write_vector(updates_group, "migration_schedule_generation",
                 delta.migration_schedule_generation);
    write_vector(updates_group, "division_schedule_generation",
                 delta.division_schedule_generation);
    write_vector(updates_group, "death_schedule_generation",
                 delta.death_schedule_generation);
}

CellFieldDelta3D read_cell_field_delta(H5::Group& cells_group,
                                      CellUid next_uid,
                                      const Model3DConfig& config,
                                      std::size_t& slot_count,
                                      std::vector<Slot>& free_slots) {
    const std::uint64_t stored_slot_count =
        read_scalar_attribute<std::uint64_t>(cells_group, "slot_count");
    if (stored_slot_count > static_cast<std::uint64_t>(kEmptySlot) ||
        stored_slot_count >
            static_cast<std::uint64_t>(
                std::numeric_limits<std::size_t>::max())) {
        throw std::runtime_error("checkpoint cell slot count is too large");
    }
    slot_count = static_cast<std::size_t>(stored_slot_count);
    free_slots = read_vector<Slot>(cells_group, "free_slots");

    CellFieldDelta3D delta;
    delta.removed_uids =
        read_vector<CellUid>(cells_group, "removed_uid");
    require_unique_nonzero_ids(delta.removed_uids, "removed cell uid");
    if (!std::is_sorted(delta.removed_uids.begin(),
                        delta.removed_uids.end())) {
        throw std::runtime_error(
            "incremental checkpoint removed cell uids are not sorted");
    }

    H5::Group births_group = cells_group.openGroup("births");
    std::size_t ignored_slot_count{};
    std::vector<Slot> ignored_free_slots;
    delta.births = read_cells(
        births_group, next_uid, config, delta.birth_slots,
        ignored_slot_count, ignored_free_slots, true);
    if (!ignored_free_slots.empty() || ignored_slot_count != slot_count) {
        throw std::runtime_error(
            "field delta birth table has invalid slot metadata");
    }

    H5::Group updates_group = cells_group.openGroup("updates");
    delta.update_uids = read_vector<CellUid>(updates_group, "uid");
    delta.update_slots = read_vector<Slot>(updates_group, "slot");
    delta.field_masks =
        read_vector<std::uint64_t>(updates_group, "field_mask");
    require_equal_sizes(
        "cell field update", delta.update_uids.size(),
        {{"slot", delta.update_slots.size()},
         {"field_mask", delta.field_masks.size()}});
    require_unique_nonzero_ids(delta.update_uids, "updated cell uid");
    if (!std::is_sorted(delta.update_uids.begin(), delta.update_uids.end())) {
        throw std::runtime_error(
            "incremental checkpoint updated cell uids are not sorted");
    }
    for (const std::uint64_t mask : delta.field_masks) {
        if (mask == 0 || (mask & ~kKnownCellFieldMask) != 0) {
            throw std::runtime_error(
                "incremental checkpoint contains an invalid cell field mask");
        }
    }

    delta.x = read_vector<std::int32_t>(updates_group, "x");
    delta.y = read_vector<std::int32_t>(updates_group, "y");
    delta.z = read_vector<std::int32_t>(updates_group, "z");
    delta.parent_uid =
        read_vector<std::uint64_t>(updates_group, "parent_uid");
    delta.clone_id =
        read_vector<std::uint32_t>(updates_group, "clone_id");
    delta.type = read_vector<std::uint8_t>(updates_group, "type");
    delta.stage = read_vector<std::uint8_t>(updates_group, "stage");
    delta.viability =
        read_vector<std::uint8_t>(updates_group, "viability");
    delta.flags = read_vector<std::uint8_t>(updates_group, "flags");
    delta.last_direction =
        read_vector<std::uint8_t>(updates_group, "last_direction");
    delta.inherent_growth_rate =
        read_vector<float>(updates_group, "inherent_growth_rate");
    delta.density_growth_rate =
        read_vector<float>(updates_group, "density_growth_rate");
    delta.migration_rate =
        read_vector<float>(updates_group, "migration_rate");
    delta.normal_migration_rate =
        read_vector<float>(updates_group, "normal_migration_rate");
    delta.migration_activation_end_time =
        read_vector<double>(updates_group,
                            "migration_activation_end_time");
    delta.division_work_remaining =
        read_vector<float>(updates_group, "division_work_remaining");
    delta.next_migration_time =
        read_vector<double>(updates_group, "next_migration_time");
    delta.next_division_time =
        read_vector<double>(updates_group, "next_division_time");
    delta.death_deadline =
        read_vector<float>(updates_group, "death_deadline");
    delta.last_update_time =
        read_vector<double>(updates_group, "last_update_time");
    delta.swap_ready_time =
        read_vector<double>(updates_group, "swap_ready_time");
    delta.swap_wait_state =
        read_vector<std::uint8_t>(updates_group, "swap_wait_state");
    delta.pending_swap_direction =
        read_vector<std::uint8_t>(updates_group,
                                  "pending_swap_direction");
    delta.event_sequence =
        read_vector<std::uint64_t>(updates_group, "event_sequence");
    delta.migration_schedule_generation =
        read_vector<std::uint32_t>(
            updates_group, "migration_schedule_generation");
    delta.division_schedule_generation =
        read_vector<std::uint32_t>(
            updates_group, "division_schedule_generation");
    delta.death_schedule_generation =
        read_vector<std::uint32_t>(
            updates_group, "death_schedule_generation");

    const auto expected = [&](CheckpointCellField3D field) {
        return field_value_count(delta.field_masks, field);
    };
    require_equal_sizes(
        "cell field update anchor",
        expected(CheckpointCellField3D::anchor),
        {{"x", delta.x.size()}, {"y", delta.y.size()},
         {"z", delta.z.size()}});
#define ATCG3D_REQUIRE_FIELD_SIZE(member, field_name)                         \
    require_equal_sizes(                                                     \
        "cell field update", expected(CheckpointCellField3D::field_name),    \
        {{#member, delta.member.size()}})
    ATCG3D_REQUIRE_FIELD_SIZE(parent_uid, parent_uid);
    ATCG3D_REQUIRE_FIELD_SIZE(clone_id, clone_id);
    ATCG3D_REQUIRE_FIELD_SIZE(type, type);
    ATCG3D_REQUIRE_FIELD_SIZE(stage, stage);
    ATCG3D_REQUIRE_FIELD_SIZE(viability, viability);
    ATCG3D_REQUIRE_FIELD_SIZE(flags, flags);
    ATCG3D_REQUIRE_FIELD_SIZE(last_direction, last_direction);
    ATCG3D_REQUIRE_FIELD_SIZE(inherent_growth_rate, inherent_growth_rate);
    ATCG3D_REQUIRE_FIELD_SIZE(density_growth_rate, density_growth_rate);
    ATCG3D_REQUIRE_FIELD_SIZE(migration_rate, migration_rate);
    ATCG3D_REQUIRE_FIELD_SIZE(normal_migration_rate, normal_migration_rate);
    ATCG3D_REQUIRE_FIELD_SIZE(migration_activation_end_time,
                              migration_activation_end_time);
    ATCG3D_REQUIRE_FIELD_SIZE(division_work_remaining,
                              division_work_remaining);
    ATCG3D_REQUIRE_FIELD_SIZE(next_migration_time, next_migration_time);
    ATCG3D_REQUIRE_FIELD_SIZE(next_division_time, next_division_time);
    ATCG3D_REQUIRE_FIELD_SIZE(death_deadline, death_deadline);
    ATCG3D_REQUIRE_FIELD_SIZE(last_update_time, last_update_time);
    ATCG3D_REQUIRE_FIELD_SIZE(swap_ready_time, swap_ready_time);
    ATCG3D_REQUIRE_FIELD_SIZE(swap_wait_state, swap_wait_state);
    ATCG3D_REQUIRE_FIELD_SIZE(pending_swap_direction,
                              pending_swap_direction);
    ATCG3D_REQUIRE_FIELD_SIZE(event_sequence, event_sequence);
    ATCG3D_REQUIRE_FIELD_SIZE(migration_schedule_generation,
                              migration_schedule_generation);
    ATCG3D_REQUIRE_FIELD_SIZE(division_schedule_generation,
                              division_schedule_generation);
    ATCG3D_REQUIRE_FIELD_SIZE(death_schedule_generation,
                              death_schedule_generation);
#undef ATCG3D_REQUIRE_FIELD_SIZE
    return delta;
}

void apply_cell_field_updates(std::vector<CellInit>& cells,
                              const std::vector<Slot>& slots,
                              const CellFieldDelta3D& delta) {
    std::size_t parent_index = 0;
    std::size_t x_index = 0, parent_uid_index = 0, clone_id_index = 0;
    std::size_t type_index = 0, stage_index = 0, viability_index = 0;
    std::size_t flags_index = 0, direction_index = 0;
    std::size_t inherent_index = 0, density_index = 0, migration_index = 0;
    std::size_t normal_migration_index = 0, activation_index = 0;
    std::size_t division_work_index = 0, next_migration_index = 0;
    std::size_t next_division_index = 0, death_index = 0;
    std::size_t last_update_index = 0, swap_ready_index = 0;
    std::size_t swap_state_index = 0, swap_direction_index = 0;
    std::size_t event_sequence_index = 0, migration_generation_index = 0;
    std::size_t division_generation_index = 0, death_generation_index = 0;

    for (std::size_t update = 0; update < delta.update_uids.size();
         ++update) {
        const CellUid uid = delta.update_uids[update];
        while (parent_index < cells.size() &&
               cells[parent_index].uid < uid) {
            ++parent_index;
        }
        if (parent_index == cells.size() ||
            cells[parent_index].uid != uid) {
            throw std::runtime_error(
                "incremental checkpoint updates a missing cell uid");
        }
        if (slots[parent_index] != delta.update_slots[update]) {
            throw std::runtime_error(
                "incremental checkpoint changes a stable cell slot");
        }
        CellInit& cell = cells[parent_index];
        const std::uint64_t mask = delta.field_masks[update];
        if (mask_has(mask, CheckpointCellField3D::anchor)) {
            cell.anchor = {delta.x[x_index], delta.y[x_index],
                           delta.z[x_index]};
            ++x_index;
        }
#define ATCG3D_APPLY_FIELD(field_name, member, values, index_name)            \
        if (mask_has(mask, CheckpointCellField3D::field_name)) {             \
            cell.member = delta.values[index_name++];                        \
        }
        ATCG3D_APPLY_FIELD(parent_uid, parent_uid, parent_uid,
                           parent_uid_index);
        ATCG3D_APPLY_FIELD(clone_id, clone_id, clone_id, clone_id_index);
        if (mask_has(mask, CheckpointCellField3D::type))
            cell.type = static_cast<CellType>(delta.type[type_index++]);
        if (mask_has(mask, CheckpointCellField3D::stage))
            cell.stage = static_cast<CellStage>(delta.stage[stage_index++]);
        ATCG3D_APPLY_FIELD(viability, viability, viability, viability_index);
        ATCG3D_APPLY_FIELD(flags, flags, flags, flags_index);
        ATCG3D_APPLY_FIELD(last_direction, last_direction, last_direction,
                           direction_index);
        ATCG3D_APPLY_FIELD(inherent_growth_rate, inherent_growth_rate,
                           inherent_growth_rate, inherent_index);
        ATCG3D_APPLY_FIELD(density_growth_rate, density_growth_rate,
                           density_growth_rate, density_index);
        ATCG3D_APPLY_FIELD(migration_rate, migration_rate, migration_rate,
                           migration_index);
        ATCG3D_APPLY_FIELD(normal_migration_rate, normal_migration_rate,
                           normal_migration_rate, normal_migration_index);
        ATCG3D_APPLY_FIELD(migration_activation_end_time,
                           migration_activation_end_time,
                           migration_activation_end_time, activation_index);
        ATCG3D_APPLY_FIELD(division_work_remaining,
                           division_work_remaining,
                           division_work_remaining, division_work_index);
        ATCG3D_APPLY_FIELD(next_migration_time, next_migration_time,
                           next_migration_time, next_migration_index);
        ATCG3D_APPLY_FIELD(next_division_time, next_division_time,
                           next_division_time, next_division_index);
        ATCG3D_APPLY_FIELD(death_deadline, death_deadline, death_deadline,
                           death_index);
        ATCG3D_APPLY_FIELD(last_update_time, last_update_time,
                           last_update_time, last_update_index);
        ATCG3D_APPLY_FIELD(swap_ready_time, swap_ready_time,
                           swap_ready_time, swap_ready_index);
        ATCG3D_APPLY_FIELD(swap_wait_state, swap_wait_state,
                           swap_wait_state, swap_state_index);
        ATCG3D_APPLY_FIELD(pending_swap_direction, pending_swap_direction,
                           pending_swap_direction, swap_direction_index);
        ATCG3D_APPLY_FIELD(event_sequence, event_sequence, event_sequence,
                           event_sequence_index);
        ATCG3D_APPLY_FIELD(migration_schedule_generation,
                           migration_schedule_generation,
                           migration_schedule_generation,
                           migration_generation_index);
        ATCG3D_APPLY_FIELD(division_schedule_generation,
                           division_schedule_generation,
                           division_schedule_generation,
                           division_generation_index);
        ATCG3D_APPLY_FIELD(death_schedule_generation,
                           death_schedule_generation,
                           death_schedule_generation,
                           death_generation_index);
#undef ATCG3D_APPLY_FIELD
        cell.schedule_generation =
            std::max({cell.migration_schedule_generation,
                      cell.division_schedule_generation,
                      cell.death_schedule_generation});
        validate_cell(cell);
    }
}

void write_lineage(H5::Group& group, std::span<const LineageEdge> lineage) {
    write_vector(group, "birth_time", column<double>(lineage, [](const auto& e) { return e.birth_time; }));
    write_vector(group, "child_uid", column<std::uint64_t>(lineage, [](const auto& e) { return e.child_uid; }));
    write_vector(group, "parent_uid", column<std::uint64_t>(lineage, [](const auto& e) { return e.parent_uid; }));
    write_vector(group, "clone_id", column<std::uint32_t>(lineage, [](const auto& e) { return e.clone_id; }));
    write_vector(group, "type", column<std::uint8_t>(lineage, [](const auto& e) { return static_cast<std::uint8_t>(e.type); }));
}

std::vector<LineageEdge> read_lineage(H5::Group& group) {
    const auto birth_time = read_vector<double>(group, "birth_time");
    const auto child_uid = read_vector<std::uint64_t>(group, "child_uid");
    const auto parent_uid = read_vector<std::uint64_t>(group, "parent_uid");
    const auto clone_id = read_vector<std::uint32_t>(group, "clone_id");
    const auto type = read_vector<std::uint8_t>(group, "type");
    require_equal_sizes("lineage", birth_time.size(),
                        {{"child_uid", child_uid.size()}, {"parent_uid", parent_uid.size()},
                         {"clone_id", clone_id.size()}, {"type", type.size()}});
    require_unique_nonzero_ids(child_uid, "lineage child uid");
    std::vector<LineageEdge> lineage(birth_time.size());
    for (std::size_t index = 0; index < lineage.size(); ++index) {
        require_finite_nonnegative(birth_time[index], "lineage.birth_time");
        if (!valid_cell_type(type[index])) {
            throw std::runtime_error("checkpoint contains an invalid lineage cell type");
        }
        lineage[index] = {birth_time[index], child_uid[index], parent_uid[index],
                          clone_id[index], static_cast<CellType>(type[index])};
    }
    return lineage;
}

void write_nodes(H5::Group& group, const std::vector<VesselNodeInit3D>& nodes) {
    write_vector(group, "x", column<std::int32_t>(nodes, [](const auto& n) { return n.position.x; }));
    write_vector(group, "y", column<std::int32_t>(nodes, [](const auto& n) { return n.position.y; }));
    write_vector(group, "z", column<std::int32_t>(nodes, [](const auto& n) { return n.position.z; }));
    write_vector(group, "uid", column<std::uint64_t>(nodes, [](const auto& n) { return n.uid; }));
    write_vector(group, "parent_uid", column<std::uint64_t>(nodes, [](const auto& n) { return n.parent_uid; }));
    write_vector(group, "parent_node_slot", column<std::uint32_t>(nodes, [](const auto& n) { return n.parent_node_slot; }));
    write_vector(group, "vessel_id", column<std::uint64_t>(nodes, [](const auto& n) { return n.vessel_id; }));
    write_vector(group, "source_lesion_id", column<std::uint64_t>(nodes, [](const auto& n) { return n.source_lesion_id; }));
    write_vector(group, "role", column<std::uint8_t>(nodes, [](const auto& n) { return static_cast<std::uint8_t>(n.role); }));
    write_vector(group, "perfused", column<std::uint8_t>(nodes, [](const auto& n) { return static_cast<std::uint8_t>(n.perfused ? 1U : 0U); }));
    write_vector(group, "diameter_voxels", column<float>(nodes, [](const auto& n) { return n.diameter_voxels; }));
    write_vector(group, "created_time_hours", column<double>(nodes, [](const auto& n) { return n.created_time_hours; }));
}

std::vector<VesselNodeInit3D> read_nodes(H5::Group& group,
                                         VesselNodeUid next_node_uid,
                                         VesselId next_vessel_id,
                                         double checkpoint_time) {
    const auto x = read_vector<std::int32_t>(group, "x");
    const auto y = read_vector<std::int32_t>(group, "y");
    const auto z = read_vector<std::int32_t>(group, "z");
    const auto uid = read_vector<std::uint64_t>(group, "uid");
    const auto parent_uid = read_vector<std::uint64_t>(group, "parent_uid");
    const auto parent_slot = read_vector<std::uint32_t>(group, "parent_node_slot");
    const auto vessel_id = read_vector<std::uint64_t>(group, "vessel_id");
    const auto source_lesion_id = read_vector<std::uint64_t>(group, "source_lesion_id");
    const auto role = read_vector<std::uint8_t>(group, "role");
    const auto perfused = read_vector<std::uint8_t>(group, "perfused");
    const auto diameter = read_vector<float>(group, "diameter_voxels");
    const auto created_time = read_vector<double>(group, "created_time_hours");
    const std::size_t count = uid.size();
    require_equal_sizes("vessel node", count,
                        {{"x", x.size()}, {"y", y.size()}, {"z", z.size()},
                         {"parent_uid", parent_uid.size()},
                         {"parent_node_slot", parent_slot.size()},
                         {"vessel_id", vessel_id.size()},
                         {"source_lesion_id", source_lesion_id.size()},
                         {"role", role.size()},
                         {"perfused", perfused.size()}, {"diameter_voxels", diameter.size()},
                         {"created_time_hours", created_time.size()}});
    require_unique_nonzero_ids(uid, "vessel node uid");
    std::vector<VesselNodeInit3D> nodes(count);
    for (std::size_t index = 0; index < count; ++index) {
        if (!valid_vessel_role(role[index]) || perfused[index] > 1U ||
            !(diameter[index] > 0.0F) || !std::isfinite(diameter[index])) {
            throw std::runtime_error("checkpoint contains an invalid vessel node value");
        }
        require_finite_nonnegative(created_time[index], "vessel_nodes.created_time_hours");
        if (uid[index] >= next_node_uid || vessel_id[index] == 0 ||
            vessel_id[index] >= next_vessel_id || created_time[index] > checkpoint_time) {
            throw std::runtime_error("checkpoint vessel node id or creation time is invalid");
        }
        nodes[index] = {{x[index], y[index], z[index]}, uid[index], parent_uid[index],
                        parent_slot[index], vessel_id[index], source_lesion_id[index],
                        static_cast<VesselBranchRole>(role[index]), perfused[index] != 0,
                        diameter[index], created_time[index]};
        if (nodes[index].role == VesselBranchRole::root) {
            if (nodes[index].parent_uid != 0 ||
                nodes[index].parent_node_slot != kEmptyVesselNodeSlot) {
                throw std::runtime_error("checkpoint root vessel node has a parent");
            }
        } else if (nodes[index].parent_uid == 0 ||
                   nodes[index].parent_node_slot >= index ||
                   nodes[nodes[index].parent_node_slot].uid != nodes[index].parent_uid ||
                   nodes[nodes[index].parent_node_slot].vessel_id != nodes[index].vessel_id ||
                   nodes[nodes[index].parent_node_slot].source_lesion_id !=
                       nodes[index].source_lesion_id) {
            throw std::runtime_error("checkpoint vessel node parent reference is invalid");
        }
    }
    return nodes;
}

void write_tips(H5::Group& group, const std::vector<VesselTipInit3D>& tips) {
    write_vector(group, "x", column<std::int32_t>(tips, [](const auto& t) { return t.position.x; }));
    write_vector(group, "y", column<std::int32_t>(tips, [](const auto& t) { return t.position.y; }));
    write_vector(group, "z", column<std::int32_t>(tips, [](const auto& t) { return t.position.z; }));
    write_vector(group, "bias_x", column<std::int32_t>(tips, [](const auto& t) { return t.bias_axis.x; }));
    write_vector(group, "bias_y", column<std::int32_t>(tips, [](const auto& t) { return t.bias_axis.y; }));
    write_vector(group, "bias_z", column<std::int32_t>(tips, [](const auto& t) { return t.bias_axis.z; }));
    write_vector(group, "target_x", column<std::int32_t>(tips, [](const auto& t) { return t.target.x; }));
    write_vector(group, "target_y", column<std::int32_t>(tips, [](const auto& t) { return t.target.y; }));
    write_vector(group, "target_z", column<std::int32_t>(tips, [](const auto& t) { return t.target.z; }));
    write_vector(group, "uid", column<std::uint64_t>(tips, [](const auto& t) { return t.uid; }));
    write_vector(group, "vessel_id", column<std::uint64_t>(tips, [](const auto& t) { return t.vessel_id; }));
    write_vector(group, "source_lesion_id", column<std::uint64_t>(tips, [](const auto& t) { return t.source_lesion_id; }));
    write_vector(group, "current_node_uid", column<std::uint64_t>(tips, [](const auto& t) { return t.current_node_uid; }));
    write_vector(group, "current_node_slot", column<std::uint32_t>(tips, [](const auto& t) { return t.current_node_slot; }));
    write_vector(group, "role", column<std::uint8_t>(tips, [](const auto& t) { return static_cast<std::uint8_t>(t.role); }));
    write_vector(group, "status", column<std::uint8_t>(tips, [](const auto& t) { return static_cast<std::uint8_t>(t.status); }));
    write_vector(group, "perfused", column<std::uint8_t>(tips, [](const auto& t) { return static_cast<std::uint8_t>(t.perfused ? 1U : 0U); }));
    write_vector(group, "last_direction", column<std::uint8_t>(tips, [](const auto& t) { return t.last_direction; }));
    write_vector(group, "pending_direction", column<std::uint8_t>(tips, [](const auto& t) { return t.pending_direction; }));
    write_vector(group, "diameter_voxels", column<float>(tips, [](const auto& t) { return t.diameter_voxels; }));
    write_vector(group, "speed_voxels_per_hour", column<float>(tips, [](const auto& t) { return t.speed_voxels_per_hour; }));
    write_vector(group, "max_length_voxels", column<float>(tips, [](const auto& t) { return t.max_length_voxels; }));
    write_vector(group, "grown_length_voxels", column<float>(tips, [](const auto& t) { return t.grown_length_voxels; }));
    write_vector(group, "next_growth_time", column<double>(tips, [](const auto& t) { return t.next_growth_time; }));
    write_vector(group, "event_sequence", column<std::uint64_t>(tips, [](const auto& t) { return t.event_sequence; }));
    write_vector(group, "schedule_generation", column<std::uint32_t>(tips, [](const auto& t) { return t.schedule_generation; }));
}

std::vector<VesselTipInit3D> read_tips(H5::Group& group,
                                       const std::vector<VesselNodeInit3D>& nodes,
                                       VesselTipUid next_tip_uid,
                                       VesselId next_vessel_id,
                                       double checkpoint_time) {
    const auto x = read_vector<std::int32_t>(group, "x");
    const auto y = read_vector<std::int32_t>(group, "y");
    const auto z = read_vector<std::int32_t>(group, "z");
    const auto bias_x = read_vector<std::int32_t>(group, "bias_x");
    const auto bias_y = read_vector<std::int32_t>(group, "bias_y");
    const auto bias_z = read_vector<std::int32_t>(group, "bias_z");
    const auto target_x = read_vector<std::int32_t>(group, "target_x");
    const auto target_y = read_vector<std::int32_t>(group, "target_y");
    const auto target_z = read_vector<std::int32_t>(group, "target_z");
    const auto uid = read_vector<std::uint64_t>(group, "uid");
    const auto vessel_id = read_vector<std::uint64_t>(group, "vessel_id");
    const auto source_lesion_id = read_vector<std::uint64_t>(group, "source_lesion_id");
    const auto current_node_uid = read_vector<std::uint64_t>(group, "current_node_uid");
    const auto current_node_slot = read_vector<std::uint32_t>(group, "current_node_slot");
    const auto role = read_vector<std::uint8_t>(group, "role");
    const auto status = read_vector<std::uint8_t>(group, "status");
    const auto perfused = read_vector<std::uint8_t>(group, "perfused");
    const auto last_direction = read_vector<std::uint8_t>(group, "last_direction");
    const auto pending_direction = read_vector<std::uint8_t>(group, "pending_direction");
    const auto diameter = read_vector<float>(group, "diameter_voxels");
    const auto speed = read_vector<float>(group, "speed_voxels_per_hour");
    const auto max_length = read_vector<float>(group, "max_length_voxels");
    const auto grown_length = read_vector<float>(group, "grown_length_voxels");
    const auto next_growth_time = read_vector<double>(group, "next_growth_time");
    const auto event_sequence = read_vector<std::uint64_t>(group, "event_sequence");
    const auto generation = read_vector<std::uint32_t>(group, "schedule_generation");
    const std::size_t count = uid.size();
    require_equal_sizes("vessel tip", count,
                        {{"x", x.size()}, {"y", y.size()}, {"z", z.size()},
                         {"bias_x", bias_x.size()}, {"bias_y", bias_y.size()},
                         {"bias_z", bias_z.size()}, {"target_x", target_x.size()},
                         {"target_y", target_y.size()}, {"target_z", target_z.size()},
                         {"vessel_id", vessel_id.size()},
                         {"source_lesion_id", source_lesion_id.size()},
                         {"current_node_uid", current_node_uid.size()},
                         {"current_node_slot", current_node_slot.size()},
                         {"role", role.size()}, {"status", status.size()},
                         {"perfused", perfused.size()},
                         {"last_direction", last_direction.size()},
                         {"pending_direction", pending_direction.size()},
                         {"diameter_voxels", diameter.size()},
                         {"speed_voxels_per_hour", speed.size()},
                         {"max_length_voxels", max_length.size()},
                         {"grown_length_voxels", grown_length.size()},
                         {"next_growth_time", next_growth_time.size()},
                         {"event_sequence", event_sequence.size()},
                         {"schedule_generation", generation.size()}});
    require_unique_nonzero_ids(uid, "vessel tip uid");
    std::vector<VesselTipInit3D> tips(count);
    for (std::size_t index = 0; index < count; ++index) {
        if (uid[index] >= next_tip_uid || vessel_id[index] == 0 ||
            vessel_id[index] >= next_vessel_id ||
            current_node_slot[index] >= nodes.size() ||
            nodes[current_node_slot[index]].uid != current_node_uid[index] ||
            nodes[current_node_slot[index]].vessel_id != vessel_id[index] ||
            nodes[current_node_slot[index]].source_lesion_id != source_lesion_id[index] ||
            !valid_vessel_role(role[index]) || !valid_tip_status(status[index]) ||
            perfused[index] > 1U || last_direction[index] > 26 ||
            pending_direction[index] > 26 || !(diameter[index] > 0.0F) ||
            !std::isfinite(diameter[index]) || speed[index] < 0.0F ||
            !std::isfinite(speed[index]) || max_length[index] < 0.0F ||
            !std::isfinite(max_length[index]) || grown_length[index] < 0.0F ||
            !std::isfinite(grown_length[index]) ||
            grown_length[index] > max_length[index]) {
            throw std::runtime_error("checkpoint contains an invalid vessel tip value");
        }
        require_finite_nonnegative(next_growth_time[index],
                                   "vessel_tips.next_growth_time");
        const auto tip_status = static_cast<VesselTipStatus>(status[index]);
        if (vessel_tip_growing(tip_status) &&
            (!(speed[index] > 0.0F) || !(max_length[index] > grown_length[index]) ||
             !(next_growth_time[index] > checkpoint_time))) {
            throw std::runtime_error("checkpoint active vessel tip scheduling state is invalid");
        }
        // Terminal tips may retain the just-consumed proposal fields. Their
        // status and bumped generation make those fields inert, but they are
        // still persisted because the checksum deliberately covers them.
        VesselTipInit3D& tip = tips[index];
        tip.position = {x[index], y[index], z[index]};
        tip.bias_axis = {bias_x[index], bias_y[index], bias_z[index]};
        tip.target = {target_x[index], target_y[index], target_z[index]};
        tip.uid = uid[index];
        tip.vessel_id = vessel_id[index];
        tip.source_lesion_id = source_lesion_id[index];
        tip.current_node_uid = current_node_uid[index];
        tip.current_node_slot = current_node_slot[index];
        tip.role = static_cast<VesselBranchRole>(role[index]);
        tip.status = tip_status;
        tip.perfused = perfused[index] != 0;
        tip.last_direction = last_direction[index];
        tip.pending_direction = pending_direction[index];
        tip.diameter_voxels = diameter[index];
        tip.speed_voxels_per_hour = speed[index];
        tip.max_length_voxels = max_length[index];
        tip.grown_length_voxels = grown_length[index];
        tip.next_growth_time = next_growth_time[index];
        tip.event_sequence = event_sequence[index];
        tip.schedule_generation = generation[index];
    }
    return tips;
}

void checked_add(std::uint64_t& total, std::uint64_t value,
                 const char* field) {
    if (value > std::numeric_limits<std::uint64_t>::max() - total) {
        throw std::runtime_error(std::string("checkpoint angiogenesis ") +
                                 field + " counter overflows");
    }
    total += value;
}

bool same_process_state(const AngiogenesisProcessState3D& lhs,
                        const AngiogenesisProcessState3D& rhs,
                        std::size_t process_count) noexcept {
    const double sum_tolerance =
        static_cast<double>(std::max<std::size_t>(process_count, 1U)) *
        std::numeric_limits<double>::epsilon() *
        std::max({1.0, std::abs(lhs.accumulated_eligible_hours),
                  std::abs(rhs.accumulated_eligible_hours)});
    return lhs.eligible == rhs.eligible &&
           lhs.next_seed_time_hours == rhs.next_seed_time_hours &&
           lhs.eligibility_started_hours == rhs.eligibility_started_hours &&
           std::abs(lhs.accumulated_eligible_hours -
                    rhs.accumulated_eligible_hours) <= sum_tolerance &&
           lhs.remaining_hazard == rhs.remaining_hazard &&
           lhs.hazard_last_update_hours == rhs.hazard_last_update_hours &&
           lhs.hazard_not_before_hours == rhs.hazard_not_before_hours &&
           lhs.current_rate_sites_per_30_days ==
               rhs.current_rate_sites_per_30_days &&
           lhs.current_density_stress == rhs.current_density_stress &&
           lhs.event_sequence == rhs.event_sequence &&
           lhs.schedule_generation == rhs.schedule_generation &&
           lhs.attempted_events == rhs.attempted_events &&
           lhs.committed_roots == rhs.committed_roots &&
           lhs.rejected_events == rhs.rejected_events;
}

void validate_vasculature(VasculatureState3D& state,
                          const Model3DConfig& config,
                          const SimulationStats3D& stats,
                          double checkpoint_time) {
    if (state.next_vessel_id == 0 || state.next_node_uid == 0 || state.next_tip_uid == 0) {
        throw std::runtime_error("checkpoint next vasculature id is zero");
    }
    if (state.lesions.next_lesion_id == kNoLesionId) {
        throw std::runtime_error("checkpoint next lesion id is zero");
    }
    require_finite_nonnegative(state.lesions.last_refresh_time_hours,
                               "vasculature.lesions.last_refresh_time_hours");
    require_finite_nonnegative(state.lesions.next_refresh_time_hours,
                               "vasculature.lesions.next_refresh_time_hours");
    if (state.lesions.last_refresh_time_hours > checkpoint_time) {
        throw std::runtime_error(
            "checkpoint lesion refresh time exceeds the checkpoint clock");
    }
    const double refresh_tolerance =
        1e-10 * std::max({1.0, checkpoint_time,
                          state.lesions.next_refresh_time_hours});
    if (state.lesions.next_refresh_time_hours > 0.0 &&
        state.lesions.next_refresh_time_hours + refresh_tolerance <
            checkpoint_time) {
        throw std::runtime_error(
            "checkpoint next lesion refresh event is in the past");
    }
    if (state.lesions.dirty_blocks.empty() !=
        (state.lesions.next_refresh_time_hours == 0.0)) {
        throw std::runtime_error(
            "checkpoint lesion dirty blocks and refresh event disagree");
    }
    if (state.lesions.next_refresh_time_hours > 0.0) {
        const double expected_refresh =
            state.lesions.last_refresh_time_hours +
            config.angiogenesis.lesion_refresh_interval_hours;
        const double interval_tolerance =
            1e-10 * std::max({1.0, expected_refresh,
                              state.lesions.next_refresh_time_hours});
        if (std::abs(state.lesions.next_refresh_time_hours -
                     expected_refresh) > interval_tolerance) {
            throw std::runtime_error(
                "checkpoint lesion refresh event does not match the configured interval");
        }
    }

    validate_process(state.process, checkpoint_time);
    const AngiogenesisProcessState3D aggregate =
        aggregate_angiogenesis_process_states(
            state.lesions.processes, checkpoint_time);
    if (!same_process_state(state.process, aggregate,
                            state.lesions.processes.size())) {
        throw std::runtime_error(
            "checkpoint aggregate angiogenesis process does not equal the "
            "per-lesion processes");
    }
    if (aggregate.attempted_events != stats.angiogenesis_seed_attempts ||
        aggregate.committed_roots != stats.angiogenesis_roots ||
        aggregate.rejected_events != stats.angiogenesis_seed_rejections) {
        throw std::runtime_error("checkpoint angiogenesis process and simulation stats disagree");
    }

    std::unordered_set<LesionId> process_ids;
    process_ids.reserve(state.lesions.processes.size());
    for (const LesionAngiogenesisState3D& entry :
         state.lesions.processes) {
        process_ids.insert(entry.lesion_id);
    }
    std::unordered_set<LesionId> current_lesion_ids;
    current_lesion_ids.reserve(state.lesions.core_identity.size());
    for (const LesionCoreIdentity3D& entry :
         state.lesions.core_identity) {
        current_lesion_ids.insert(entry.lesion_id);
        if (!process_ids.contains(entry.lesion_id)) {
            throw std::runtime_error(
                "checkpoint current lesion has no angiogenesis process");
        }
    }
    for (const LesionAngiogenesisState3D& entry :
         state.lesions.processes) {
        if (entry.process.eligible &&
            !current_lesion_ids.contains(entry.lesion_id)) {
            throw std::runtime_error(
                "checkpoint eligible process references a retired lesion");
        }
    }
    for (const LesionSourceOwnership3D& ownership :
         state.lesions.source_ownership) {
        if (current_lesion_ids.contains(ownership.source_lesion_id) ||
            (ownership.current_lesion_id != kNoLesionId &&
             !current_lesion_ids.contains(ownership.current_lesion_id))) {
            throw std::runtime_error(
                "checkpoint lesion source ownership is not a direct historical-to-current mapping");
        }
    }
    std::unordered_set<LesionId> known_vessel_sources = current_lesion_ids;
    known_vessel_sources.reserve(current_lesion_ids.size() +
                                 state.lesions.source_ownership.size());
    for (const LesionSourceOwnership3D& ownership :
         state.lesions.source_ownership) {
        known_vessel_sources.insert(ownership.source_lesion_id);
    }

    if (!std::is_sorted(state.perfused_vessels.begin(), state.perfused_vessels.end()) ||
        std::adjacent_find(state.perfused_vessels.begin(), state.perfused_vessels.end()) !=
            state.perfused_vessels.end()) {
        throw std::runtime_error("checkpoint perfused vessel ids must be sorted and unique");
    }
    std::unordered_set<VesselId> vessels;
    std::unordered_map<VesselId, LesionId> vessel_sources;
    std::uint64_t root_count = 0;
    for (const auto& node : state.nodes) {
        vessels.insert(node.vessel_id);
        if (node.source_lesion_id == kNoLesionId ||
            node.source_lesion_id >= state.lesions.next_lesion_id ||
            !known_vessel_sources.contains(node.source_lesion_id)) {
            throw std::runtime_error(
                "checkpoint vessel node source is neither current nor historically owned");
        }
        const auto [source, inserted] = vessel_sources.try_emplace(
            node.vessel_id, node.source_lesion_id);
        if (!inserted && source->second != node.source_lesion_id) {
            throw std::runtime_error(
                "checkpoint vessel contains nodes from different source lesions");
        }
        if (node.role == VesselBranchRole::root) {
            checked_add(root_count, 1, "root-node");
        }
    }
    const std::unordered_set<VesselId> perfused(state.perfused_vessels.begin(),
                                                state.perfused_vessels.end());
    for (const VesselId id : state.perfused_vessels) {
        if (id == 0 || id >= state.next_vessel_id || !vessels.contains(id)) {
            throw std::runtime_error("checkpoint perfused vessel id is invalid");
        }
    }
    for (const auto& node : state.nodes) {
        if (node.perfused != perfused.contains(node.vessel_id)) {
            throw std::runtime_error("checkpoint vessel node perfusion state is inconsistent");
        }
    }
    for (const auto& tip : state.tips) {
        const auto source = vessel_sources.find(tip.vessel_id);
        if (tip.source_lesion_id == kNoLesionId ||
            tip.source_lesion_id >= state.lesions.next_lesion_id ||
            source == vessel_sources.end() ||
            source->second != tip.source_lesion_id) {
            throw std::runtime_error(
                "checkpoint vessel tip source lesion is inconsistent");
        }
        if (tip.perfused != perfused.contains(tip.vessel_id)) {
            throw std::runtime_error("checkpoint vessel tip perfusion state is inconsistent");
        }
    }
    if (root_count != aggregate.committed_roots) {
        throw std::runtime_error("checkpoint root-node and angiogenesis counts disagree");
    }
    if (!config.angiogenesis.enabled &&
        (!state.nodes.empty() || !state.tips.empty() || !state.perfused_vessels.empty() ||
         !state.lesions.core_identity.empty() ||
         !state.lesions.dirty_blocks.empty() ||
         !state.lesions.processes.empty() ||
         !state.lesions.source_ownership.empty() ||
         state.lesions.next_lesion_id != 1 ||
         state.lesions.last_refresh_time_hours != 0.0 ||
         state.lesions.next_refresh_time_hours != 0.0 ||
         state.lesions.refresh_schedule_generation != 0 ||
         state.process.eligible || state.process.event_sequence != 0 ||
         state.process.attempted_events != 0 || state.next_vessel_id != 1 ||
         state.next_node_uid != 1 || state.next_tip_uid != 1)) {
        throw std::runtime_error(
            "checkpoint contains vasculature but angiogenesis is disabled");
    }
}

void read_vasculature_payload(H5::H5File& file,
                              CheckpointData3D& result,
                              const Model3DConfig& expected_config) {
    H5::Group vasculature_group = file.openGroup("/vasculature");
    result.vasculature.next_vessel_id =
        read_scalar_attribute<std::uint64_t>(vasculature_group,
                                             "next_vessel_id");
    result.vasculature.next_node_uid =
        read_scalar_attribute<std::uint64_t>(vasculature_group,
                                             "next_node_uid");
    result.vasculature.next_tip_uid =
        read_scalar_attribute<std::uint64_t>(vasculature_group,
                                             "next_tip_uid");
    result.vasculature.perfused_vessels =
        read_vector<std::uint64_t>(vasculature_group,
                                   "perfused_vessel_ids");

    H5::Group process_group = file.openGroup("/vasculature/process");
    result.vasculature.process = read_process(process_group);
    H5::Group lesions_group = file.openGroup("/vasculature/lesions");
    result.vasculature.lesions.next_lesion_id =
        read_scalar_attribute<std::uint64_t>(lesions_group,
                                             "next_lesion_id");
    result.vasculature.lesions.last_refresh_time_hours =
        read_scalar_attribute<double>(lesions_group,
                                      "last_refresh_time_hours");
    result.vasculature.lesions.next_refresh_time_hours =
        read_scalar_attribute<double>(lesions_group,
                                      "next_refresh_time_hours");
    result.vasculature.lesions.refresh_schedule_generation =
        read_scalar_attribute<std::uint32_t>(
            lesions_group, "refresh_schedule_generation");
    H5::Group lesion_core_group =
        file.openGroup("/vasculature/lesions/core_identity");
    result.vasculature.lesions.core_identity =
        read_lesion_core_identity(
            lesion_core_group,
            result.vasculature.lesions.next_lesion_id);
    H5::Group lesion_dirty_group =
        file.openGroup("/vasculature/lesions/dirty_blocks");
    result.vasculature.lesions.dirty_blocks =
        read_lesion_dirty_blocks(lesion_dirty_group, expected_config);
    H5::Group lesion_processes_group =
        file.openGroup("/vasculature/lesions/processes");
    result.vasculature.lesions.processes =
        read_lesion_processes(
            lesion_processes_group,
            result.vasculature.lesions.next_lesion_id,
            result.clock.time_hours);
    H5::Group lesion_source_ownership_group =
        file.openGroup("/vasculature/lesions/source_ownership");
    result.vasculature.lesions.source_ownership =
        read_lesion_source_ownership(
            lesion_source_ownership_group,
            result.vasculature.lesions.next_lesion_id);
    H5::Group nodes_group = file.openGroup("/vasculature/nodes");
    result.vasculature.nodes =
        read_nodes(nodes_group, result.vasculature.next_node_uid,
                   result.vasculature.next_vessel_id,
                   result.clock.time_hours);
    H5::Group tips_group = file.openGroup("/vasculature/tips");
    result.vasculature.tips =
        read_tips(tips_group, result.vasculature.nodes,
                  result.vasculature.next_tip_uid,
                  result.vasculature.next_vessel_id,
                  result.clock.time_hours);
}

void validate_complete_checkpoint(CheckpointData3D& result,
                                  const Model3DConfig& expected_config,
                                  bool verify_state_checksum = true) {
    if (result.next_uid == 0) {
        throw std::runtime_error("checkpoint next cell uid is zero");
    }
    require_finite_nonnegative(result.clock.time_hours, "meta.time_hours");
    for (const CellInit& cell : result.cells) {
        const double tolerance = float_storage_time_tolerance(
            result.clock.time_hours, cell.last_update_time);
        if (cell.last_update_time > result.clock.time_hours + tolerance) {
            throw std::runtime_error(
                "checkpoint cell was updated after the checkpoint clock");
        }
        for (const double event_time : {cell.next_migration_time,
                                        cell.migration_activation_end_time,
                                        cell.next_division_time,
                                        cell.death_deadline}) {
            if (event_time > 0.0 && event_time < result.clock.time_hours) {
                throw std::runtime_error(
                    "checkpoint contains a living cell event in the past");
            }
        }
    }
    for (const LineageEdge& edge : result.lineage) {
        if (edge.birth_time > result.clock.time_hours ||
            edge.child_uid >= result.next_uid) {
            throw std::runtime_error(
                "checkpoint lineage time or child uid is invalid");
        }
    }
    validate_vasculature(result.vasculature, expected_config, result.stats,
                         result.clock.time_hours);

    Simulation3D reconstructed(expected_config);
    reconstructed.restore(result.cells, result.next_uid, result.clock,
                          result.stats, result.lineage, result.vasculature,
                          result.cell_slot_count, result.cell_slots,
                          result.cell_free_slots);
    if (verify_state_checksum &&
        reconstructed.state_checksum() != result.state_checksum) {
        throw std::runtime_error(
            "checkpoint reconstructed state checksum does not match metadata");
    }
}

}  // namespace

CheckpointDeltaSummary3D summarize_checkpoint_delta(
    const CheckpointSnapshotView3D& current,
    const CheckpointSnapshotView3D& parent) {
    validate_snapshot_order(current, "current checkpoint");
    validate_snapshot_order(parent, "parent checkpoint");
    if (current.clock.time_hours < parent.clock.time_hours ||
        current.next_uid < parent.next_uid ||
        current.lineage.size() < parent.lineage.size()) {
        throw std::invalid_argument(
            "checkpoint delta state precedes its parent");
    }
    for (std::size_t index = 0; index < parent.lineage.size(); ++index) {
        if (!same_lineage(current.lineage[index], parent.lineage[index])) {
            throw std::invalid_argument(
                "checkpoint lineage is not an append-only parent prefix");
        }
    }

    CheckpointDeltaSummary3D summary;
    std::size_t current_index = 0;
    std::size_t parent_index = 0;
    while (current_index < current.cells.size() &&
           parent_index < parent.cells.size()) {
        const CellUid current_uid = current.cells[current_index].uid;
        const CellUid parent_uid = parent.cells[parent_index].uid;
        if (current_uid < parent_uid) {
            ++summary.changed_cells;
            ++current_index;
        } else if (parent_uid < current_uid) {
            ++summary.removed_cells;
            ++parent_index;
        } else {
            if (current.cell_slots[current_index] !=
                    parent.cell_slots[parent_index] ||
                !same_cell(current.cells[current_index],
                           parent.cells[parent_index])) {
                ++summary.changed_cells;
            }
            ++current_index;
            ++parent_index;
        }
    }
    summary.changed_cells += current.cells.size() - current_index;
    summary.removed_cells += parent.cells.size() - parent_index;
    summary.appended_lineage_edges =
        current.lineage.size() - parent.lineage.size();
    const std::size_t denominator =
        std::max(current.cells.size(), parent.cells.size());
    summary.changed_fraction = denominator == 0
        ? 0.0
        : static_cast<double>(summary.changed_cells + summary.removed_cells) /
              static_cast<double>(denominator);
    return summary;
}

void write_hdf5_checkpoint(const std::filesystem::path& path,
                           const Simulation3D& simulation) {
    const std::vector<Slot> cell_slots = simulation.snapshot_cell_slots();
    std::vector<CellInit> cells;
    cells.reserve(cell_slots.size());
    for (const Slot slot : cell_slots) {
        cells.push_back(simulation.cells().snapshot(slot));
    }
    const VasculatureState3D vasculature = simulation.snapshot_vasculature();
    const CheckpointSnapshotView3D snapshot{
        cells,
        cell_slots,
        simulation.cells().slot_count(),
        simulation.cells().free_slots(),
        simulation.next_uid(),
        simulation.clock(),
        simulation.stats(),
        simulation.lineage(),
        vasculature,
        simulation.state_checksum()};
    write_hdf5_checkpoint(path, snapshot, simulation.config());
}

void write_hdf5_checkpoint(const std::filesystem::path& path,
                           const CheckpointSnapshotView3D& snapshot,
                           const Model3DConfig& config) {
    if (path.empty()) throw std::invalid_argument("checkpoint path must not be empty");
    if (!path.parent_path().empty()) std::filesystem::create_directories(path.parent_path());
    const std::filesystem::path temporary = path.string() + ".tmp";
    if (std::filesystem::exists(path)) {
        throw std::runtime_error("refusing to overwrite checkpoint: " + path.string());
    }
    std::filesystem::remove(temporary);
    try {
        const ScopedDatasetWriteOptions dataset_options(config);
        H5::H5File file(temporary.string(), H5F_ACC_TRUNC);
        H5::Group meta = file.createGroup("/meta");
        H5::Group stats_group = file.createGroup("/stats");
        H5::Group cells_group = file.createGroup("/cells");
        H5::Group lineage_group = file.createGroup("/lineage");
        H5::Group vasculature_group = file.createGroup("/vasculature");
        H5::Group process_group = file.createGroup("/vasculature/process");
        H5::Group lesions_group = file.createGroup("/vasculature/lesions");
        H5::Group lesion_core_group =
            file.createGroup("/vasculature/lesions/core_identity");
        H5::Group lesion_dirty_group =
            file.createGroup("/vasculature/lesions/dirty_blocks");
        H5::Group lesion_processes_group =
            file.createGroup("/vasculature/lesions/processes");
        H5::Group lesion_source_ownership_group =
            file.createGroup("/vasculature/lesions/source_ownership");
        H5::Group nodes_group = file.createGroup("/vasculature/nodes");
        H5::Group tips_group = file.createGroup("/vasculature/tips");

        write_scalar_attribute(meta, "schema_version", kCheckpointSchemaVersion);
        write_scalar_attribute(meta, "dimension", kCheckpointDimension);
        write_scalar_attribute(meta, "completed_events", snapshot.clock.completed_events);
        write_scalar_attribute(meta, "time_hours", snapshot.clock.time_hours);
        write_scalar_attribute(meta, "next_uid", snapshot.next_uid);
        write_scalar_attribute(meta, "state_checksum", snapshot.state_checksum);
        write_string_attribute(meta, "effective_config_json", config.to_json());
        write_string_attribute(meta, "dynamics_config_json", config.dynamics_json());

        write_stats(stats_group, snapshot.stats);
        write_cells(cells_group, snapshot.cells, snapshot.cell_slots,
                    snapshot.cell_slot_count, snapshot.cell_free_slots);
        write_lineage(lineage_group, snapshot.lineage);

        const VasculatureState3D& vasculature = snapshot.vasculature;
        write_scalar_attribute(vasculature_group, "next_vessel_id",
                               vasculature.next_vessel_id);
        write_scalar_attribute(vasculature_group, "next_node_uid", vasculature.next_node_uid);
        write_scalar_attribute(vasculature_group, "next_tip_uid", vasculature.next_tip_uid);
        write_process(process_group, vasculature.process);
        write_scalar_attribute(lesions_group, "next_lesion_id",
                               vasculature.lesions.next_lesion_id);
        write_scalar_attribute(lesions_group, "last_refresh_time_hours",
                               vasculature.lesions.last_refresh_time_hours);
        write_scalar_attribute(lesions_group, "next_refresh_time_hours",
                               vasculature.lesions.next_refresh_time_hours);
        write_scalar_attribute(
            lesions_group, "refresh_schedule_generation",
            vasculature.lesions.refresh_schedule_generation);
        write_lesion_core_identity(lesion_core_group,
                                   vasculature.lesions.core_identity);
        write_lesion_dirty_blocks(lesion_dirty_group,
                                  vasculature.lesions.dirty_blocks);
        write_lesion_processes(lesion_processes_group,
                               vasculature.lesions.processes);
        write_lesion_source_ownership(
            lesion_source_ownership_group,
            vasculature.lesions.source_ownership);
        write_nodes(nodes_group, vasculature.nodes);
        write_tips(tips_group, vasculature.tips);
        write_vector(vasculature_group, "perfused_vessel_ids",
                     vasculature.perfused_vessels);

        file.flush(H5F_SCOPE_GLOBAL);
        file.close();

        {
            H5::H5File verified(temporary.string(), H5F_ACC_RDONLY);
            H5::Group verified_meta = verified.openGroup("/meta");
            H5::Group verified_cells = verified.openGroup("/cells");
            if (read_scalar_attribute<std::uint32_t>(
                    verified_meta, "schema_version") !=
                    kCheckpointSchemaVersion ||
                read_vector<CellUid>(verified_cells, "uid").size() !=
                    snapshot.cells.size() ||
                read_vector<double>(
                    verified_cells, "swap_ready_time").size() !=
                    snapshot.cells.size()) {
                throw std::runtime_error(
                    "checkpoint local schema or dataset verification failed");
            }
        }
        std::filesystem::rename(temporary, path);
    } catch (...) {
        std::filesystem::remove(temporary);
        throw;
    }
}

void write_hdf5_delta_checkpoint(
    const std::filesystem::path& path,
    const CheckpointSnapshotView3D& snapshot,
    const CheckpointSnapshotView3D& parent_snapshot,
    const std::filesystem::path& parent_path,
    std::uint64_t chain_length,
    const Model3DConfig& config) {
    const CheckpointDeltaSummary3D summary =
        summarize_checkpoint_delta(snapshot, parent_snapshot);
    if (path.empty() || parent_path.empty() || chain_length == 0 ||
        chain_length > config.checkpoint_max_delta_chain) {
        throw std::invalid_argument("incremental checkpoint path or chain is invalid");
    }
    if (!std::filesystem::is_regular_file(parent_path)) {
        throw std::invalid_argument("incremental checkpoint parent does not exist");
    }
    const std::filesystem::path output_directory =
        std::filesystem::absolute(path.parent_path().empty()
                                      ? std::filesystem::path{"."}
                                      : path.parent_path())
            .lexically_normal();
    const std::filesystem::path parent_directory =
        std::filesystem::absolute(parent_path.parent_path().empty()
                                      ? std::filesystem::path{"."}
                                      : parent_path.parent_path())
            .lexically_normal();
    if (output_directory != parent_directory ||
        parent_path.filename().empty()) {
        throw std::invalid_argument(
            "incremental checkpoint parent must be in the same directory");
    }

    CellFieldDelta3D delta;
    delta.births.reserve(summary.changed_cells);
    delta.birth_slots.reserve(summary.changed_cells);
    delta.update_uids.reserve(summary.changed_cells);
    delta.update_slots.reserve(summary.changed_cells);
    delta.field_masks.reserve(summary.changed_cells);
    delta.removed_uids.reserve(summary.removed_cells);
    std::size_t current_index = 0;
    std::size_t parent_index = 0;
    while (current_index < snapshot.cells.size() &&
           parent_index < parent_snapshot.cells.size()) {
        const CellUid current_uid = snapshot.cells[current_index].uid;
        const CellUid parent_uid = parent_snapshot.cells[parent_index].uid;
        if (current_uid < parent_uid) {
            delta.births.push_back(snapshot.cells[current_index]);
            delta.birth_slots.push_back(snapshot.cell_slots[current_index]);
            ++current_index;
        } else if (parent_uid < current_uid) {
            delta.removed_uids.push_back(parent_uid);
            ++parent_index;
        } else {
            if (snapshot.cell_slots[current_index] !=
                parent_snapshot.cell_slots[parent_index]) {
                throw std::invalid_argument(
                    "incremental checkpoint changes a stable cell slot");
            }
            const std::uint64_t mask = changed_cell_field_mask(
                snapshot.cells[current_index],
                parent_snapshot.cells[parent_index]);
            if (mask != 0) {
                delta.update_uids.push_back(current_uid);
                delta.update_slots.push_back(
                    snapshot.cell_slots[current_index]);
                delta.field_masks.push_back(mask);
                append_cell_field_values(
                    delta, snapshot.cells[current_index], mask);
            }
            ++current_index;
            ++parent_index;
        }
    }
    for (; current_index < snapshot.cells.size(); ++current_index) {
        delta.births.push_back(snapshot.cells[current_index]);
        delta.birth_slots.push_back(snapshot.cell_slots[current_index]);
    }
    for (; parent_index < parent_snapshot.cells.size(); ++parent_index) {
        delta.removed_uids.push_back(
            parent_snapshot.cells[parent_index].uid);
    }

    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    const std::filesystem::path temporary = path.string() + ".tmp";
    if (std::filesystem::exists(path)) {
        throw std::runtime_error("refusing to overwrite checkpoint: " +
                                 path.string());
    }
    std::filesystem::remove(temporary);
    try {
        const ScopedDatasetWriteOptions dataset_options(config);
        H5::H5File file(temporary.string(), H5F_ACC_TRUNC);
        H5::Group meta = file.createGroup("/meta");
        H5::Group stats_group = file.createGroup("/stats");
        H5::Group cells_group = file.createGroup("/cells");
        H5::Group lineage_group = file.createGroup("/lineage");
        H5::Group vasculature_group = file.createGroup("/vasculature");
        H5::Group process_group = file.createGroup("/vasculature/process");
        H5::Group lesions_group = file.createGroup("/vasculature/lesions");
        H5::Group lesion_core_group =
            file.createGroup("/vasculature/lesions/core_identity");
        H5::Group lesion_dirty_group =
            file.createGroup("/vasculature/lesions/dirty_blocks");
        H5::Group lesion_processes_group =
            file.createGroup("/vasculature/lesions/processes");
        H5::Group lesion_source_ownership_group =
            file.createGroup("/vasculature/lesions/source_ownership");
        H5::Group nodes_group = file.createGroup("/vasculature/nodes");
        H5::Group tips_group = file.createGroup("/vasculature/tips");

        write_scalar_attribute(meta, "schema_version",
                               kCheckpointDeltaSchemaVersion);
        write_scalar_attribute(meta, "dimension", kCheckpointDimension);
        write_string_attribute(meta, "kind", "field_delta_v2");
        write_string_attribute(meta, "parent_file",
                               parent_path.filename().string());
        write_scalar_attribute(meta, "parent_state_checksum",
                               parent_snapshot.state_checksum);
        write_scalar_attribute(meta, "parent_time_hours",
                               parent_snapshot.clock.time_hours);
        write_scalar_attribute(meta, "chain_length", chain_length);
        write_scalar_attribute(
            meta, "lineage_prefix_count",
            static_cast<std::uint64_t>(parent_snapshot.lineage.size()));
        write_scalar_attribute(meta, "completed_events",
                               snapshot.clock.completed_events);
        write_scalar_attribute(meta, "time_hours", snapshot.clock.time_hours);
        write_scalar_attribute(meta, "next_uid", snapshot.next_uid);
        write_scalar_attribute(meta, "state_checksum", snapshot.state_checksum);
        write_string_attribute(meta, "effective_config_json", config.to_json());
        write_string_attribute(meta, "dynamics_config_json",
                               config.dynamics_json());

        write_stats(stats_group, snapshot.stats);
        write_cell_field_delta(cells_group, delta,
                               snapshot.cell_slot_count,
                               snapshot.cell_free_slots);
        write_lineage(lineage_group,
                      snapshot.lineage.subspan(parent_snapshot.lineage.size()));

        const VasculatureState3D& vasculature = snapshot.vasculature;
        write_scalar_attribute(vasculature_group, "next_vessel_id",
                               vasculature.next_vessel_id);
        write_scalar_attribute(vasculature_group, "next_node_uid",
                               vasculature.next_node_uid);
        write_scalar_attribute(vasculature_group, "next_tip_uid",
                               vasculature.next_tip_uid);
        write_process(process_group, vasculature.process);
        write_scalar_attribute(lesions_group, "next_lesion_id",
                               vasculature.lesions.next_lesion_id);
        write_scalar_attribute(lesions_group, "last_refresh_time_hours",
                               vasculature.lesions.last_refresh_time_hours);
        write_scalar_attribute(lesions_group, "next_refresh_time_hours",
                               vasculature.lesions.next_refresh_time_hours);
        write_scalar_attribute(
            lesions_group, "refresh_schedule_generation",
            vasculature.lesions.refresh_schedule_generation);
        write_lesion_core_identity(lesion_core_group,
                                   vasculature.lesions.core_identity);
        write_lesion_dirty_blocks(lesion_dirty_group,
                                  vasculature.lesions.dirty_blocks);
        write_lesion_processes(lesion_processes_group,
                               vasculature.lesions.processes);
        write_lesion_source_ownership(
            lesion_source_ownership_group,
            vasculature.lesions.source_ownership);
        write_nodes(nodes_group, vasculature.nodes);
        write_tips(tips_group, vasculature.tips);
        write_vector(vasculature_group, "perfused_vessel_ids",
                     vasculature.perfused_vessels);

        file.flush(H5F_SCOPE_GLOBAL);
        file.close();

        {
            H5::H5File verified(temporary.string(), H5F_ACC_RDONLY);
            H5::Group verified_meta = verified.openGroup("/meta");
            if (read_scalar_attribute<std::uint32_t>(
                    verified_meta, "schema_version") !=
                    kCheckpointDeltaSchemaVersion ||
                read_string_attribute(verified_meta, "kind") !=
                    "field_delta_v2") {
                throw std::runtime_error(
                    "incremental checkpoint local schema verification failed");
            }
            H5::Group verified_cells = verified.openGroup("/cells");
            H5::Group verified_births =
                verified_cells.openGroup("births");
            H5::Group verified_updates =
                verified_cells.openGroup("updates");
            if (read_vector<CellUid>(verified_births, "uid").size() !=
                    delta.births.size() ||
                read_vector<CellUid>(verified_updates, "uid").size() !=
                    delta.update_uids.size() ||
                read_vector<CellUid>(verified_cells, "removed_uid").size() !=
                    delta.removed_uids.size()) {
                throw std::runtime_error(
                    "incremental checkpoint local dataset verification failed");
            }
        }
        std::filesystem::rename(temporary, path);
    } catch (...) {
        std::filesystem::remove(temporary);
        throw;
    }
}

void write_hdf5_journal_delta_checkpoint(
    const std::filesystem::path& path,
    const CheckpointJournalSnapshotView3D& snapshot,
    const std::filesystem::path& parent_path,
    std::uint64_t parent_state_checksum,
    double parent_time_hours,
    std::uint64_t chain_length,
    const Model3DConfig& config) {
    if (path.empty() || parent_path.empty() || chain_length == 0 ||
        chain_length > config.checkpoint_max_delta_chain) {
        throw std::invalid_argument(
            "journal checkpoint path or chain is invalid");
    }
    if (!std::filesystem::is_regular_file(parent_path)) {
        throw std::invalid_argument(
            "journal checkpoint parent does not exist");
    }
    const std::filesystem::path output_directory =
        std::filesystem::absolute(path.parent_path().empty()
                                      ? std::filesystem::path{"."}
                                      : path.parent_path())
            .lexically_normal();
    const std::filesystem::path parent_directory =
        std::filesystem::absolute(parent_path.parent_path().empty()
                                      ? std::filesystem::path{"."}
                                      : parent_path.parent_path())
            .lexically_normal();
    if (output_directory != parent_directory ||
        parent_path.filename().empty()) {
        throw std::invalid_argument(
            "journal checkpoint parent must be in the same directory");
    }

    std::vector<CellInit> changed_cells;
    std::vector<Slot> changed_slots;
    std::vector<Slot> removed_slots;
    changed_cells.reserve(snapshot.cell_mutations.size());
    changed_slots.reserve(snapshot.cell_mutations.size());
    removed_slots.reserve(snapshot.cell_mutations.size());
    for (const CheckpointCellMutation3D& mutation :
         snapshot.cell_mutations) {
        if (mutation.slot == kEmptySlot ||
            static_cast<std::size_t>(mutation.slot) >=
                snapshot.cell_slot_count) {
            throw std::invalid_argument(
                "journal checkpoint mutation slot is invalid");
        }
        if (mutation.alive) {
            changed_cells.push_back(mutation.cell);
            changed_slots.push_back(mutation.slot);
        } else {
            removed_slots.push_back(mutation.slot);
        }
    }
    std::vector<std::size_t> order(changed_cells.size());
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::sort(order.begin(), order.end(),
              [&](std::size_t lhs, std::size_t rhs) {
                  return changed_cells[lhs].uid < changed_cells[rhs].uid;
              });
    std::vector<CellInit> sorted_cells;
    std::vector<Slot> sorted_slots;
    sorted_cells.reserve(order.size());
    sorted_slots.reserve(order.size());
    for (const std::size_t index : order) {
        sorted_cells.push_back(changed_cells[index]);
        sorted_slots.push_back(changed_slots[index]);
    }
    std::sort(removed_slots.begin(), removed_slots.end());
    if (std::adjacent_find(removed_slots.begin(), removed_slots.end()) !=
        removed_slots.end()) {
        throw std::invalid_argument(
            "journal checkpoint contains duplicate removed slots");
    }

    std::vector<std::uint8_t> free_kinds;
    std::vector<Slot> free_slots;
    free_kinds.reserve(snapshot.free_list_mutations.size());
    free_slots.reserve(snapshot.free_list_mutations.size());
    for (const FreeListMutation3D mutation :
         snapshot.free_list_mutations) {
        if (mutation.slot == kEmptySlot) {
            throw std::invalid_argument(
                "journal checkpoint free-list slot is invalid");
        }
        free_kinds.push_back(static_cast<std::uint8_t>(mutation.kind));
        free_slots.push_back(mutation.slot);
    }

    if (!path.parent_path().empty()) {
        std::filesystem::create_directories(path.parent_path());
    }
    const std::filesystem::path temporary = path.string() + ".tmp";
    if (std::filesystem::exists(path)) {
        throw std::runtime_error("refusing to overwrite checkpoint: " +
                                 path.string());
    }
    std::filesystem::remove(temporary);
    try {
        const ScopedDatasetWriteOptions dataset_options(config);
        H5::H5File file(temporary.string(), H5F_ACC_TRUNC);
        H5::Group meta = file.createGroup("/meta");
        H5::Group stats_group = file.createGroup("/stats");
        H5::Group cells_group = file.createGroup("/cells");
        H5::Group changed_group = file.createGroup("/cells/changed");
        H5::Group free_group = file.createGroup("/cells/free_list_mutations");
        H5::Group lineage_group = file.createGroup("/lineage");
        H5::Group vasculature_group = file.createGroup("/vasculature");
        H5::Group process_group = file.createGroup("/vasculature/process");
        H5::Group lesions_group = file.createGroup("/vasculature/lesions");
        H5::Group lesion_core_group =
            file.createGroup("/vasculature/lesions/core_identity");
        H5::Group lesion_dirty_group =
            file.createGroup("/vasculature/lesions/dirty_blocks");
        H5::Group lesion_processes_group =
            file.createGroup("/vasculature/lesions/processes");
        H5::Group lesion_source_ownership_group =
            file.createGroup("/vasculature/lesions/source_ownership");
        H5::Group nodes_group = file.createGroup("/vasculature/nodes");
        H5::Group tips_group = file.createGroup("/vasculature/tips");

        write_scalar_attribute(meta, "schema_version",
                               kCheckpointJournalDeltaSchemaVersion);
        write_scalar_attribute(meta, "dimension", kCheckpointDimension);
        write_string_attribute(meta, "kind", "slot_journal_v1");
        write_string_attribute(meta, "parent_file",
                               parent_path.filename().string());
        write_scalar_attribute(meta, "parent_state_checksum",
                               parent_state_checksum);
        write_scalar_attribute(meta, "parent_time_hours",
                               parent_time_hours);
        write_scalar_attribute(meta, "chain_length", chain_length);
        write_scalar_attribute(
            meta, "lineage_prefix_count",
            static_cast<std::uint64_t>(snapshot.lineage_prefix_count));
        write_scalar_attribute(meta, "completed_events",
                               snapshot.clock.completed_events);
        write_scalar_attribute(meta, "time_hours",
                               snapshot.clock.time_hours);
        write_scalar_attribute(meta, "next_uid", snapshot.next_uid);
        write_scalar_attribute(meta, "state_checksum",
                               snapshot.state_checksum);
        write_string_attribute(meta, "effective_config_json",
                               config.to_json());
        write_string_attribute(meta, "dynamics_config_json",
                               config.dynamics_json());

        write_stats(stats_group, snapshot.stats);
        write_scalar_attribute(
            cells_group, "slot_count",
            static_cast<std::uint64_t>(snapshot.cell_slot_count));
        write_vector<Slot>(cells_group, "removed_slot", removed_slots);
        write_cells(changed_group, sorted_cells, sorted_slots,
                    snapshot.cell_slot_count, {});
        write_vector<std::uint8_t>(free_group, "kind", free_kinds);
        write_vector<Slot>(free_group, "slot", free_slots);
        write_lineage(lineage_group, snapshot.lineage_tail);

        const VasculatureState3D& vasculature = snapshot.vasculature;
        write_scalar_attribute(vasculature_group, "next_vessel_id",
                               vasculature.next_vessel_id);
        write_scalar_attribute(vasculature_group, "next_node_uid",
                               vasculature.next_node_uid);
        write_scalar_attribute(vasculature_group, "next_tip_uid",
                               vasculature.next_tip_uid);
        write_process(process_group, vasculature.process);
        write_scalar_attribute(lesions_group, "next_lesion_id",
                               vasculature.lesions.next_lesion_id);
        write_scalar_attribute(lesions_group, "last_refresh_time_hours",
                               vasculature.lesions.last_refresh_time_hours);
        write_scalar_attribute(lesions_group, "next_refresh_time_hours",
                               vasculature.lesions.next_refresh_time_hours);
        write_scalar_attribute(
            lesions_group, "refresh_schedule_generation",
            vasculature.lesions.refresh_schedule_generation);
        write_lesion_core_identity(lesion_core_group,
                                   vasculature.lesions.core_identity);
        write_lesion_dirty_blocks(lesion_dirty_group,
                                  vasculature.lesions.dirty_blocks);
        write_lesion_processes(lesion_processes_group,
                               vasculature.lesions.processes);
        write_lesion_source_ownership(
            lesion_source_ownership_group,
            vasculature.lesions.source_ownership);
        write_nodes(nodes_group, vasculature.nodes);
        write_tips(tips_group, vasculature.tips);
        write_vector(vasculature_group, "perfused_vessel_ids",
                     vasculature.perfused_vessels);

        file.flush(H5F_SCOPE_GLOBAL);
        file.close();
        {
            H5::H5File verified(temporary.string(), H5F_ACC_RDONLY);
            H5::Group verified_meta = verified.openGroup("/meta");
            H5::Group verified_cells = verified.openGroup("/cells");
            if (read_scalar_attribute<std::uint32_t>(
                    verified_meta, "schema_version") !=
                    kCheckpointJournalDeltaSchemaVersion ||
                read_string_attribute(verified_meta, "kind") !=
                    "slot_journal_v1" ||
                read_vector<Slot>(verified_cells, "removed_slot").size() !=
                    removed_slots.size()) {
                throw std::runtime_error(
                    "journal checkpoint local verification failed");
            }
        }
        std::filesystem::rename(temporary, path);
    } catch (...) {
        std::filesystem::remove(temporary);
        throw;
    }
}

CheckpointData3D read_hdf5_checkpoint(const std::filesystem::path& path,
                                      const Model3DConfig& expected_config) {
    const CheckpointReadGuard read_guard(path);
    try {
        H5::H5File file(path.string(), H5F_ACC_RDONLY);
        H5::Group meta = file.openGroup("/meta");
        const std::uint32_t schema =
            read_scalar_attribute<std::uint32_t>(meta, "schema_version");
        if (schema == kCheckpointJournalDeltaSchemaVersion) {
            if (read_scalar_attribute<std::uint32_t>(meta, "dimension") !=
                    kCheckpointDimension ||
                read_string_attribute(meta, "kind") !=
                    "slot_journal_v1") {
                throw std::runtime_error(
                    "unsupported journal checkpoint schema or dimension");
            }
            validate_effective_config_provenance(
                read_string_attribute(meta, "effective_config_json"),
                expected_config);
            if (read_string_attribute(meta, "dynamics_config_json") !=
                expected_config.dynamics_json()) {
                throw std::runtime_error(
                    "checkpoint biological/numerical configuration does not match");
            }
            const std::string parent_file =
                read_string_attribute(meta, "parent_file");
            const std::filesystem::path parent_name(parent_file);
            if (parent_file.empty() ||
                parent_name != parent_name.filename() ||
                parent_name == "." || parent_name == "..") {
                throw std::runtime_error(
                    "journal checkpoint parent path is unsafe");
            }
            const std::uint64_t chain_length =
                read_scalar_attribute<std::uint64_t>(meta, "chain_length");
            if (chain_length == 0 || chain_length > 4096) {
                throw std::runtime_error(
                    "journal checkpoint chain exceeds the reader hard limit");
            }
            const std::filesystem::path parent_path =
                path.parent_path() / parent_name;
            CheckpointData3D parent =
                read_hdf5_checkpoint(parent_path, expected_config);
            {
                H5::H5File parent_file_handle(
                    parent_path.string(), H5F_ACC_RDONLY);
                H5::Group parent_meta =
                    parent_file_handle.openGroup("/meta");
                const std::uint32_t parent_schema =
                    read_scalar_attribute<std::uint32_t>(
                        parent_meta, "schema_version");
                const bool parent_is_base =
                    parent_schema == kLegacyCheckpointSchemaVersion ||
                    parent_schema == kCheckpointSchemaVersion;
                const bool parent_is_delta =
                    parent_schema == kLegacyCheckpointDeltaSchemaVersion ||
                    parent_schema == kCheckpointDeltaSchemaVersion ||
                    parent_schema ==
                        kCheckpointJournalDeltaSchemaVersion;
                const std::uint64_t parent_chain = parent_is_base
                    ? 0
                    : parent_is_delta
                        ? read_scalar_attribute<std::uint64_t>(
                              parent_meta, "chain_length")
                        : std::numeric_limits<std::uint64_t>::max();
                if ((!parent_is_base && !parent_is_delta) ||
                    chain_length != parent_chain + 1) {
                    throw std::runtime_error(
                        "journal checkpoint chain length does not match parent");
                }
            }
            if (read_scalar_attribute<std::uint64_t>(
                    meta, "parent_state_checksum") !=
                    parent.state_checksum ||
                read_scalar_attribute<double>(meta, "parent_time_hours") !=
                    parent.clock.time_hours) {
                throw std::runtime_error(
                    "journal checkpoint parent identity does not match");
            }
            const std::uint64_t lineage_prefix_count =
                read_scalar_attribute<std::uint64_t>(
                    meta, "lineage_prefix_count");
            if (lineage_prefix_count != parent.lineage.size()) {
                throw std::runtime_error(
                    "journal checkpoint lineage prefix does not match parent");
            }

            CheckpointData3D result;
            result.next_uid =
                read_scalar_attribute<std::uint64_t>(meta, "next_uid");
            result.clock.completed_events =
                read_scalar_attribute<std::uint64_t>(
                    meta, "completed_events");
            result.clock.time_hours =
                read_scalar_attribute<double>(meta, "time_hours");
            result.state_checksum =
                read_scalar_attribute<std::uint64_t>(
                    meta, "state_checksum");
            if (result.next_uid < parent.next_uid ||
                !(result.clock.time_hours > parent.clock.time_hours) ||
                result.clock.completed_events <
                    parent.clock.completed_events) {
                throw std::runtime_error(
                    "journal checkpoint does not advance its parent");
            }
            H5::Group stats_group = file.openGroup("/stats");
            result.stats = read_stats(stats_group, true);
            H5::Group cells_group = file.openGroup("/cells");
            const std::uint64_t stored_slot_count =
                read_scalar_attribute<std::uint64_t>(
                    cells_group, "slot_count");
            if (stored_slot_count >
                    static_cast<std::uint64_t>(kEmptySlot) ||
                stored_slot_count < parent.cell_slot_count) {
                throw std::runtime_error(
                    "journal checkpoint slot count is invalid");
            }
            result.cell_slot_count =
                static_cast<std::size_t>(stored_slot_count);
            const std::vector<Slot> removed_slots =
                read_vector<Slot>(cells_group, "removed_slot");
            if (!std::is_sorted(removed_slots.begin(),
                                removed_slots.end()) ||
                std::adjacent_find(removed_slots.begin(),
                                   removed_slots.end()) !=
                    removed_slots.end()) {
                throw std::runtime_error(
                    "journal removed slots are not unique and sorted");
            }
            H5::Group changed_group =
                file.openGroup("/cells/changed");
            std::vector<Slot> changed_slots;
            std::vector<Slot> ignored_free_slots;
            std::size_t changed_slot_count{};
            std::vector<CellInit> changed_cells = read_cells(
                changed_group, result.next_uid, expected_config,
                changed_slots, changed_slot_count, ignored_free_slots, true);
            if (changed_slot_count != result.cell_slot_count ||
                !ignored_free_slots.empty()) {
                throw std::runtime_error(
                    "journal changed-cell metadata is inconsistent");
            }

            std::vector<std::int64_t> slot_to_index(
                result.cell_slot_count, -1);
            for (std::size_t index = 0;
                 index < parent.cell_slots.size(); ++index) {
                const Slot slot = parent.cell_slots[index];
                if (static_cast<std::size_t>(slot) >=
                        result.cell_slot_count ||
                    slot_to_index[slot] != -1) {
                    throw std::runtime_error(
                        "journal parent cell slots are invalid");
                }
                slot_to_index[slot] =
                    static_cast<std::int64_t>(index);
            }
            std::vector<std::uint8_t> remove_parent(
                parent.cells.size(), 0);
            for (const Slot slot : removed_slots) {
                if (static_cast<std::size_t>(slot) >=
                        result.cell_slot_count ||
                    slot_to_index[slot] < 0) {
                    throw std::runtime_error(
                        "journal removes a missing cell slot");
                }
                remove_parent[static_cast<std::size_t>(
                    slot_to_index[slot])] = 1;
                slot_to_index[slot] = -1;
            }
            std::vector<CellInit> appended_cells;
            std::vector<Slot> appended_slots;
            for (std::size_t index = 0;
                 index < changed_cells.size(); ++index) {
                const Slot slot = changed_slots[index];
                if (static_cast<std::size_t>(slot) >=
                    result.cell_slot_count) {
                    throw std::runtime_error(
                        "journal changes an out-of-range cell slot");
                }
                const std::int64_t parent_index = slot_to_index[slot];
                if (parent_index >= 0 &&
                    parent.cells[static_cast<std::size_t>(
                        parent_index)].uid == changed_cells[index].uid) {
                    parent.cells[static_cast<std::size_t>(parent_index)] =
                        changed_cells[index];
                } else {
                    if (parent_index >= 0) {
                        remove_parent[static_cast<std::size_t>(
                            parent_index)] = 1;
                    }
                    appended_cells.push_back(changed_cells[index]);
                    appended_slots.push_back(slot);
                }
            }
            result.cells.reserve(parent.cells.size() +
                                 appended_cells.size());
            result.cell_slots.reserve(parent.cell_slots.size() +
                                      appended_slots.size());
            for (std::size_t index = 0;
                 index < parent.cells.size(); ++index) {
                if (remove_parent[index] != 0) continue;
                result.cells.push_back(parent.cells[index]);
                result.cell_slots.push_back(parent.cell_slots[index]);
            }
            result.cells.insert(result.cells.end(),
                                appended_cells.begin(),
                                appended_cells.end());
            result.cell_slots.insert(result.cell_slots.end(),
                                     appended_slots.begin(),
                                     appended_slots.end());

            result.cell_free_slots = std::move(parent.cell_free_slots);
            H5::Group free_group =
                file.openGroup("/cells/free_list_mutations");
            const std::vector<std::uint8_t> free_kinds =
                read_vector<std::uint8_t>(free_group, "kind");
            const std::vector<Slot> free_slots =
                read_vector<Slot>(free_group, "slot");
            require_equal_sizes(
                "journal free-list mutation", free_kinds.size(),
                {{"slot", free_slots.size()}});
            for (std::size_t index = 0;
                 index < free_kinds.size(); ++index) {
                const Slot slot = free_slots[index];
                if (free_kinds[index] ==
                    static_cast<std::uint8_t>(
                        FreeListMutationKind3D::push)) {
                    result.cell_free_slots.push_back(slot);
                } else if (free_kinds[index] ==
                           static_cast<std::uint8_t>(
                               FreeListMutationKind3D::pop)) {
                    if (result.cell_free_slots.empty() ||
                        result.cell_free_slots.back() != slot) {
                        throw std::runtime_error(
                            "journal free-list pop does not match parent");
                    }
                    result.cell_free_slots.pop_back();
                } else {
                    throw std::runtime_error(
                        "journal free-list mutation kind is invalid");
                }
            }

            H5::Group lineage_group = file.openGroup("/lineage");
            result.lineage = std::move(parent.lineage);
            std::vector<LineageEdge> lineage_tail =
                read_lineage(lineage_group);
            result.lineage.insert(result.lineage.end(),
                                  lineage_tail.begin(),
                                  lineage_tail.end());
            read_vasculature_payload(file, result, expected_config);
            if (read_guard.outermost()) {
                validate_complete_checkpoint(result, expected_config);
            }
            return result;
        }
        if (schema == kCheckpointDeltaSchemaVersion) {
            if (read_scalar_attribute<std::uint32_t>(meta, "dimension") !=
                    kCheckpointDimension ||
                read_string_attribute(meta, "kind") !=
                    "field_delta_v2") {
                throw std::runtime_error(
                    "unsupported field checkpoint schema or dimension");
            }
            validate_effective_config_provenance(
                read_string_attribute(meta, "effective_config_json"),
                expected_config);
            if (read_string_attribute(meta, "dynamics_config_json") !=
                expected_config.dynamics_json()) {
                throw std::runtime_error(
                    "checkpoint biological/numerical configuration does not match");
            }
            const std::string parent_file =
                read_string_attribute(meta, "parent_file");
            const std::filesystem::path parent_name(parent_file);
            if (parent_file.empty() ||
                parent_name != parent_name.filename() ||
                parent_name == "." || parent_name == "..") {
                throw std::runtime_error(
                    "incremental checkpoint parent path is unsafe");
            }
            const std::uint64_t chain_length =
                read_scalar_attribute<std::uint64_t>(meta, "chain_length");
            if (chain_length == 0 || chain_length > 4096) {
                throw std::runtime_error(
                    "incremental checkpoint chain exceeds the reader hard limit");
            }
            const std::filesystem::path parent_path =
                path.parent_path() / parent_name;
            CheckpointData3D parent =
                read_hdf5_checkpoint(parent_path, expected_config);
            {
                H5::H5File parent_file_handle(
                    parent_path.string(), H5F_ACC_RDONLY);
                H5::Group parent_meta =
                    parent_file_handle.openGroup("/meta");
                const std::uint32_t parent_schema =
                    read_scalar_attribute<std::uint32_t>(
                        parent_meta, "schema_version");
                const bool parent_is_base =
                    parent_schema == kLegacyCheckpointSchemaVersion ||
                    parent_schema == kCheckpointSchemaVersion;
                const bool parent_is_delta =
                    parent_schema ==
                        kLegacyCheckpointDeltaSchemaVersion ||
                    parent_schema == kCheckpointDeltaSchemaVersion ||
                    parent_schema ==
                        kCheckpointJournalDeltaSchemaVersion;
                const std::uint64_t parent_chain = parent_is_base
                    ? 0
                    : parent_is_delta
                        ? read_scalar_attribute<std::uint64_t>(
                              parent_meta, "chain_length")
                        : std::numeric_limits<std::uint64_t>::max();
                if ((!parent_is_base && !parent_is_delta) ||
                    chain_length != parent_chain + 1) {
                    throw std::runtime_error(
                        "incremental checkpoint chain length does not match parent");
                }
            }
            if (read_scalar_attribute<std::uint64_t>(
                    meta, "parent_state_checksum") !=
                    parent.state_checksum ||
                read_scalar_attribute<double>(meta, "parent_time_hours") !=
                    parent.clock.time_hours) {
                throw std::runtime_error(
                    "incremental checkpoint parent identity does not match");
            }
            const std::uint64_t lineage_prefix_count =
                read_scalar_attribute<std::uint64_t>(
                    meta, "lineage_prefix_count");
            if (lineage_prefix_count != parent.lineage.size()) {
                throw std::runtime_error(
                    "incremental checkpoint lineage prefix does not match parent");
            }

            CheckpointData3D result;
            result.next_uid =
                read_scalar_attribute<std::uint64_t>(meta, "next_uid");
            result.clock.completed_events =
                read_scalar_attribute<std::uint64_t>(
                    meta, "completed_events");
            result.clock.time_hours =
                read_scalar_attribute<double>(meta, "time_hours");
            result.state_checksum =
                read_scalar_attribute<std::uint64_t>(
                    meta, "state_checksum");
            if (result.next_uid < parent.next_uid ||
                !(result.clock.time_hours > parent.clock.time_hours) ||
                result.clock.completed_events <
                    parent.clock.completed_events) {
                throw std::runtime_error(
                    "incremental checkpoint does not advance its parent");
            }

            H5::Group stats_group = file.openGroup("/stats");
            result.stats = read_stats(stats_group, true);
            H5::Group cells_group = file.openGroup("/cells");
            CellFieldDelta3D delta = read_cell_field_delta(
                cells_group, result.next_uid, expected_config,
                result.cell_slot_count, result.cell_free_slots);
            for (const CellInit& birth : delta.births) {
                if (birth.uid < parent.next_uid) {
                    throw std::runtime_error(
                        "incremental checkpoint birth reuses an old uid");
                }
            }
            {
                std::vector<CellUid> overlap;
                std::set_intersection(
                    delta.update_uids.begin(), delta.update_uids.end(),
                    delta.removed_uids.begin(), delta.removed_uids.end(),
                    std::back_inserter(overlap));
                if (!overlap.empty()) {
                    throw std::runtime_error(
                        "incremental checkpoint both removes and updates a cell");
                }
            }
            apply_cell_field_updates(
                parent.cells, parent.cell_slots, delta);

            result.cells.reserve(parent.cells.size() +
                                 delta.births.size());
            result.cell_slots.reserve(parent.cell_slots.size() +
                                      delta.birth_slots.size());
            std::size_t removed_index = 0;
            for (std::size_t index = 0; index < parent.cells.size();
                 ++index) {
                const CellUid uid = parent.cells[index].uid;
                if (removed_index < delta.removed_uids.size() &&
                    delta.removed_uids[removed_index] < uid) {
                    throw std::runtime_error(
                        "incremental checkpoint removes a missing cell uid");
                }
                if (removed_index < delta.removed_uids.size() &&
                    delta.removed_uids[removed_index] == uid) {
                    ++removed_index;
                    continue;
                }
                result.cells.push_back(parent.cells[index]);
                result.cell_slots.push_back(parent.cell_slots[index]);
            }
            if (removed_index != delta.removed_uids.size()) {
                throw std::runtime_error(
                    "incremental checkpoint removes a missing cell uid");
            }
            result.cells.insert(result.cells.end(), delta.births.begin(),
                                delta.births.end());
            result.cell_slots.insert(result.cell_slots.end(),
                                     delta.birth_slots.begin(),
                                     delta.birth_slots.end());

            H5::Group lineage_group = file.openGroup("/lineage");
            result.lineage = std::move(parent.lineage);
            std::vector<LineageEdge> lineage_tail =
                read_lineage(lineage_group);
            result.lineage.insert(result.lineage.end(),
                                  lineage_tail.begin(),
                                  lineage_tail.end());
            read_vasculature_payload(file, result, expected_config);
            if (read_guard.outermost()) {
                validate_complete_checkpoint(result, expected_config);
            }
            return result;
        }
        if (schema == kLegacyCheckpointDeltaSchemaVersion) {
            if (read_scalar_attribute<std::uint32_t>(meta, "dimension") !=
                    kCheckpointDimension ||
                read_string_attribute(meta, "kind") != "delta_v1") {
                throw std::runtime_error(
                    "unsupported incremental checkpoint schema or dimension");
            }
            validate_effective_config_provenance(
                read_string_attribute(meta, "effective_config_json"),
                expected_config);
            if (read_string_attribute(meta, "dynamics_config_json") !=
                expected_config.dynamics_json()) {
                throw std::runtime_error(
                    "checkpoint biological/numerical configuration does not match");
            }
            const std::string parent_file =
                read_string_attribute(meta, "parent_file");
            const std::filesystem::path parent_name(parent_file);
            if (parent_file.empty() || parent_name != parent_name.filename() ||
                parent_name == "." || parent_name == "..") {
                throw std::runtime_error(
                    "incremental checkpoint parent path is unsafe");
            }
            const std::uint64_t chain_length =
                read_scalar_attribute<std::uint64_t>(meta, "chain_length");
            if (chain_length == 0 || chain_length > 4096) {
                throw std::runtime_error(
                    "incremental checkpoint chain exceeds the reader hard limit");
            }
            const std::filesystem::path parent_path =
                path.parent_path() / parent_name;
            CheckpointData3D parent =
                read_hdf5_checkpoint(parent_path, expected_config);
            {
                H5::H5File parent_file(parent_path.string(), H5F_ACC_RDONLY);
                H5::Group parent_meta = parent_file.openGroup("/meta");
                const std::uint32_t parent_schema =
                    read_scalar_attribute<std::uint32_t>(
                        parent_meta, "schema_version");
                const std::uint64_t parent_chain =
                    parent_schema == kLegacyCheckpointSchemaVersion
                        ? 0
                        : read_scalar_attribute<std::uint64_t>(
                              parent_meta, "chain_length");
                if ((parent_schema != kLegacyCheckpointSchemaVersion &&
                     parent_schema !=
                         kLegacyCheckpointDeltaSchemaVersion) ||
                    chain_length != parent_chain + 1) {
                    throw std::runtime_error(
                        "incremental checkpoint chain length does not match parent");
                }
            }
            if (read_scalar_attribute<std::uint64_t>(
                    meta, "parent_state_checksum") != parent.state_checksum ||
                read_scalar_attribute<double>(meta, "parent_time_hours") !=
                    parent.clock.time_hours) {
                throw std::runtime_error(
                    "incremental checkpoint parent identity does not match");
            }
            const std::uint64_t lineage_prefix_count =
                read_scalar_attribute<std::uint64_t>(
                    meta, "lineage_prefix_count");
            if (lineage_prefix_count != parent.lineage.size()) {
                throw std::runtime_error(
                    "incremental checkpoint lineage prefix does not match parent");
            }

            CheckpointData3D result;
            result.next_uid =
                read_scalar_attribute<std::uint64_t>(meta, "next_uid");
            result.clock.completed_events =
                read_scalar_attribute<std::uint64_t>(meta,
                                                     "completed_events");
            result.clock.time_hours =
                read_scalar_attribute<double>(meta, "time_hours");
            result.state_checksum =
                read_scalar_attribute<std::uint64_t>(meta,
                                                     "state_checksum");
            if (result.next_uid < parent.next_uid ||
                !(result.clock.time_hours > parent.clock.time_hours) ||
                result.clock.completed_events <
                    parent.clock.completed_events) {
                throw std::runtime_error(
                    "incremental checkpoint does not advance its parent");
            }

            H5::Group stats_group = file.openGroup("/stats");
            result.stats = read_stats(stats_group, false);
            H5::Group cells_group = file.openGroup("/cells");
            std::vector<Slot> changed_slots;
            std::size_t stored_slot_count = 0;
            std::vector<Slot> stored_free_slots;
            std::vector<CellInit> changed_cells = read_cells(
                cells_group, result.next_uid, expected_config, changed_slots,
                stored_slot_count, stored_free_slots, false);
            std::vector<CellUid> removed_uids =
                read_vector<CellUid>(cells_group, "removed_uid");
            require_unique_nonzero_ids(removed_uids, "removed cell uid");
            if (!std::is_sorted(removed_uids.begin(), removed_uids.end())) {
                throw std::runtime_error(
                    "incremental checkpoint removed cell uids are not sorted");
            }
            result.cell_slot_count = stored_slot_count;
            result.cell_free_slots = std::move(stored_free_slots);

            result.cells.reserve(parent.cells.size() + changed_cells.size());
            result.cell_slots.reserve(parent.cell_slots.size() +
                                      changed_slots.size());
            std::size_t changed_index = 0;
            std::size_t removed_index = 0;
            for (std::size_t parent_index = 0;
                 parent_index < parent.cells.size(); ++parent_index) {
                const CellUid uid = parent.cells[parent_index].uid;
                if (changed_index < changed_cells.size() &&
                    changed_cells[changed_index].uid < uid) {
                    throw std::runtime_error(
                        "incremental checkpoint introduces a non-monotonic cell uid");
                }
                if (removed_index < removed_uids.size() &&
                    removed_uids[removed_index] < uid) {
                    throw std::runtime_error(
                        "incremental checkpoint removes a missing cell uid");
                }
                const bool removed = removed_index < removed_uids.size() &&
                                     removed_uids[removed_index] == uid;
                const bool changed = changed_index < changed_cells.size() &&
                                     changed_cells[changed_index].uid == uid;
                if (removed && changed) {
                    throw std::runtime_error(
                        "incremental checkpoint both removes and updates a cell");
                }
                if (removed) {
                    ++removed_index;
                    continue;
                }
                if (changed) {
                    if (changed_slots[changed_index] !=
                        parent.cell_slots[parent_index]) {
                        throw std::runtime_error(
                            "incremental checkpoint changes a stable cell slot");
                    }
                    result.cells.push_back(changed_cells[changed_index]);
                    result.cell_slots.push_back(changed_slots[changed_index]);
                    ++changed_index;
                } else {
                    result.cells.push_back(parent.cells[parent_index]);
                    result.cell_slots.push_back(
                        parent.cell_slots[parent_index]);
                }
            }
            if (removed_index != removed_uids.size()) {
                throw std::runtime_error(
                    "incremental checkpoint removes a missing cell uid");
            }
            for (; changed_index < changed_cells.size(); ++changed_index) {
                if (changed_cells[changed_index].uid < parent.next_uid) {
                    throw std::runtime_error(
                        "incremental checkpoint update references a missing old uid");
                }
                result.cells.push_back(changed_cells[changed_index]);
                result.cell_slots.push_back(changed_slots[changed_index]);
            }

            H5::Group lineage_group = file.openGroup("/lineage");
            result.lineage = std::move(parent.lineage);
            std::vector<LineageEdge> lineage_tail =
                read_lineage(lineage_group);
            result.lineage.insert(result.lineage.end(), lineage_tail.begin(),
                                  lineage_tail.end());
            read_vasculature_payload(file, result, expected_config);
            if (read_guard.outermost()) {
                validate_complete_checkpoint(result, expected_config, false);
            }
            return result;
        }
        if (schema == 1 || schema == 2 || schema == 3) {
            throw std::runtime_error(
                "checkpoint schema v" + std::to_string(schema) +
                " is explicitly unsupported by the schema-v4 reader; older "
                "checkpoints do not contain density-modulated angiogenesis hazard state");
        }
        if ((schema != kLegacyCheckpointSchemaVersion &&
             schema != kCheckpointSchemaVersion) ||
            read_scalar_attribute<std::uint32_t>(meta, "dimension") !=
                kCheckpointDimension) {
            throw std::runtime_error("unsupported checkpoint schema or dimension");
        }
        validate_effective_config_provenance(
            read_string_attribute(meta, "effective_config_json"),
            expected_config);
        if (read_string_attribute(meta, "dynamics_config_json") !=
            expected_config.dynamics_json()) {
            throw std::runtime_error(
                "checkpoint biological/numerical configuration does not match");
        }

        CheckpointData3D result;
        result.next_uid = read_scalar_attribute<std::uint64_t>(meta, "next_uid");
        if (result.next_uid == 0) throw std::runtime_error("checkpoint next cell uid is zero");
        result.clock.completed_events =
            read_scalar_attribute<std::uint64_t>(meta, "completed_events");
        result.clock.time_hours = read_scalar_attribute<double>(meta, "time_hours");
        require_finite_nonnegative(result.clock.time_hours, "meta.time_hours");
        result.state_checksum =
            read_scalar_attribute<std::uint64_t>(meta, "state_checksum");

        H5::Group stats_group = file.openGroup("/stats");
        result.stats = read_stats(
            stats_group, schema == kCheckpointSchemaVersion);

        H5::Group cells_group = file.openGroup("/cells");
        result.cells = read_cells(
            cells_group, result.next_uid, expected_config, result.cell_slots,
            result.cell_slot_count, result.cell_free_slots,
            schema == kCheckpointSchemaVersion);
        for (const CellInit& cell : result.cells) {
            const double tolerance = float_storage_time_tolerance(
                result.clock.time_hours, cell.last_update_time);
            if (cell.last_update_time > result.clock.time_hours + tolerance) {
                throw std::runtime_error(
                    "checkpoint cell was updated after the checkpoint clock");
            }
            for (const double event_time : {cell.next_migration_time,
                                            cell.migration_activation_end_time,
                                            cell.next_division_time,
                                            cell.death_deadline}) {
                if (event_time > 0.0 &&
                    event_time < result.clock.time_hours) {
                    throw std::runtime_error(
                        "checkpoint contains a living cell event in the past");
                }
            }
        }

        H5::Group lineage_group = file.openGroup("/lineage");
        result.lineage = read_lineage(lineage_group);
        for (const LineageEdge& edge : result.lineage) {
            if (edge.birth_time > result.clock.time_hours || edge.child_uid >= result.next_uid) {
                throw std::runtime_error("checkpoint lineage time or child uid is invalid");
            }
        }

        H5::Group vasculature_group = file.openGroup("/vasculature");
        result.vasculature.next_vessel_id =
            read_scalar_attribute<std::uint64_t>(vasculature_group, "next_vessel_id");
        result.vasculature.next_node_uid =
            read_scalar_attribute<std::uint64_t>(vasculature_group, "next_node_uid");
        result.vasculature.next_tip_uid =
            read_scalar_attribute<std::uint64_t>(vasculature_group, "next_tip_uid");
        result.vasculature.perfused_vessels =
            read_vector<std::uint64_t>(vasculature_group, "perfused_vessel_ids");

        H5::Group process_group = file.openGroup("/vasculature/process");
        result.vasculature.process = read_process(process_group);
        H5::Group lesions_group = file.openGroup("/vasculature/lesions");
        result.vasculature.lesions.next_lesion_id =
            read_scalar_attribute<std::uint64_t>(lesions_group,
                                                 "next_lesion_id");
        result.vasculature.lesions.last_refresh_time_hours =
            read_scalar_attribute<double>(lesions_group,
                                          "last_refresh_time_hours");
        result.vasculature.lesions.next_refresh_time_hours =
            read_scalar_attribute<double>(lesions_group,
                                          "next_refresh_time_hours");
        result.vasculature.lesions.refresh_schedule_generation =
            read_scalar_attribute<std::uint32_t>(
                lesions_group, "refresh_schedule_generation");
        H5::Group lesion_core_group =
            file.openGroup("/vasculature/lesions/core_identity");
        result.vasculature.lesions.core_identity =
            read_lesion_core_identity(
                lesion_core_group,
                result.vasculature.lesions.next_lesion_id);
        H5::Group lesion_dirty_group =
            file.openGroup("/vasculature/lesions/dirty_blocks");
        result.vasculature.lesions.dirty_blocks =
            read_lesion_dirty_blocks(lesion_dirty_group,
                                     expected_config);
        H5::Group lesion_processes_group =
            file.openGroup("/vasculature/lesions/processes");
        result.vasculature.lesions.processes =
            read_lesion_processes(
                lesion_processes_group,
                result.vasculature.lesions.next_lesion_id,
                result.clock.time_hours);
        H5::Group lesion_source_ownership_group =
            file.openGroup("/vasculature/lesions/source_ownership");
        result.vasculature.lesions.source_ownership =
            read_lesion_source_ownership(
                lesion_source_ownership_group,
                result.vasculature.lesions.next_lesion_id);
        H5::Group nodes_group = file.openGroup("/vasculature/nodes");
        result.vasculature.nodes =
            read_nodes(nodes_group, result.vasculature.next_node_uid,
                       result.vasculature.next_vessel_id, result.clock.time_hours);
        H5::Group tips_group = file.openGroup("/vasculature/tips");
        result.vasculature.tips =
            read_tips(tips_group, result.vasculature.nodes,
                      result.vasculature.next_tip_uid,
                      result.vasculature.next_vessel_id, result.clock.time_hours);
        if (read_guard.outermost()) {
            validate_complete_checkpoint(
                result, expected_config,
                schema == kCheckpointSchemaVersion);
        }
        return result;
    } catch (const H5::Exception& error) {
        throw std::runtime_error("invalid or unreadable HDF5 checkpoint: " +
                                 error.getDetailMsg());
    }
}

}  // namespace atcg3d
