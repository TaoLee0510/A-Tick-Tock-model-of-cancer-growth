"""Materialize ATCG3D checkpoint deltas as server-side vtkPolyData.

The browser never receives these arrays.  They remain in the ParaView Python
process and are rendered through the normal remote-view pipeline.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path
import struct

import h5py
import numpy as np


_BASE_CHECKPOINT_SCHEMAS = frozenset((4, 6))
_LEGACY_ROW_DELTA_SCHEMA = 5
_FIELD_DELTA_SCHEMA = 7
_SLOT_JOURNAL_SCHEMA = 8

# Keep this table in the same order as CheckpointCellField3D.  Every field
# column is compacted independently: a value is present only when its bit is
# set for an update row.
_CELL_FIELD_LAYOUT = (
    ("anchor", 0, ("x", "y", "z")),
    ("parent_uid", 1, ("parent_uid",)),
    ("clone_id", 2, ("clone_id",)),
    ("type", 3, ("type",)),
    ("stage", 4, ("stage",)),
    ("viability", 5, ("viability",)),
    ("flags", 6, ("flags",)),
    ("last_direction", 7, ("last_direction",)),
    ("inherent_growth_rate", 8, ("inherent_growth_rate",)),
    ("density_growth_rate", 9, ("density_growth_rate",)),
    ("migration_rate", 10, ("migration_rate",)),
    ("normal_migration_rate", 11, ("normal_migration_rate",)),
    ("migration_activation_end_time", 12,
     ("migration_activation_end_time",)),
    ("division_work_remaining", 13, ("division_work_remaining",)),
    ("next_migration_time", 14, ("next_migration_time",)),
    ("next_division_time", 15, ("next_division_time",)),
    ("death_deadline", 16, ("death_deadline",)),
    ("last_update_time", 17, ("last_update_time",)),
    ("swap_ready_time", 18, ("swap_ready_time",)),
    ("swap_wait_state", 19, ("swap_wait_state",)),
    ("pending_swap_direction", 20, ("pending_swap_direction",)),
    ("event_sequence", 21, ("event_sequence",)),
    ("migration_schedule_generation", 22,
     ("migration_schedule_generation",)),
    ("division_schedule_generation", 23,
     ("division_schedule_generation",)),
    ("death_schedule_generation", 24, ("death_schedule_generation",)),
)
_KNOWN_CELL_FIELD_MASK = (1 << len(_CELL_FIELD_LAYOUT)) - 1


@dataclass
class CellDisplayState:
    uid: np.ndarray
    slot: np.ndarray
    points: np.ndarray
    lesion_id: np.ndarray
    clone_id: np.ndarray
    cell_type: np.ndarray
    stage: np.ndarray
    viability: np.ndarray
    display_radius: np.ndarray

    def take(self, indices):
        return CellDisplayState(**{
            name: getattr(self, name)[indices]
            for name in self.__dataclass_fields__
        })


def _time_from_checkpoint_name(path: Path) -> float:
    time_hex = path.name.removesuffix(".h5").rsplit("_time_", 1)[-1]
    if len(time_hex) != 16:
        raise ValueError(f"checkpoint has no encoded time: {path}")
    return struct.unpack(">d", int(time_hex, 16).to_bytes(8, "big"))[0]


def _checkpoint_files(run_directory: Path):
    files = []
    for path in (run_directory / "checkpoints").glob(
            "checkpoint_*_time_*.h5"):
        if path.name.endswith(".tmp"):
            continue
        try:
            value = _time_from_checkpoint_name(path)
        except (ValueError, OverflowError, struct.error):
            continue
        if math.isfinite(value) and value >= 0.0:
            files.append((value, path))
    return sorted(files)


def _series_frames(run_directory: Path, name: str):
    path = run_directory / name
    if not path.exists():
        return []
    document = json.loads(path.read_text(encoding="utf-8"))
    return sorted(
        (float(entry["time"]), run_directory / entry["name"])
        for entry in document.get("files", [])
        if not str(entry["name"]).endswith(".tmp")
    )


class CheckpointFrameMaterializer:
    """Reconstruct display fields from a bounded keyframe+delta chain."""

    def __init__(self, run_directory: Path, radii=(1.0, 0.5, 0.25)):
        self.run_directory = Path(run_directory).resolve()
        self.radii = np.asarray(radii, dtype=np.float32)
        if self.radii.shape != (3,) or np.any(~np.isfinite(self.radii)) or \
                np.any(self.radii <= 0.0):
            raise ValueError("display radii must be three finite positive values")
        self.state: CellDisplayState | None = None
        self.time_hours: float | None = None
        self.previous_checkpoint: Path | None = None
        self._polydata_arrays = None

    @staticmethod
    def _sorted(state: CellDisplayState):
        if state.uid.size < 2:
            return state
        if np.any(state.uid[:-1] == state.uid[1:]):
            raise ValueError("display state contains duplicate cell IDs")
        if np.all(state.uid[:-1] < state.uid[1:]):
            return state
        order = np.argsort(state.uid, kind="stable")
        result = state.take(order)
        if np.any(result.uid[:-1] == result.uid[1:]):
            raise ValueError("display state contains duplicate cell IDs")
        return result

    @staticmethod
    def _read_vector(group, name, dtype, context):
        if name not in group:
            raise ValueError(f"{context} is missing dataset {name}")
        values = np.asarray(group[name])
        if values.ndim != 1:
            raise ValueError(f"{context} dataset {name} is not one-dimensional")
        return values.astype(dtype, copy=False)

    @staticmethod
    def _slot_count(group, context, required=False):
        if "slot_count" not in group.attrs:
            if required:
                raise ValueError(f"{context} has no slot_count")
            return None
        value = int(group.attrs["slot_count"])
        if value < 0 or value > np.iinfo(np.uint32).max:
            raise ValueError(f"{context} has invalid slot_count")
        return value

    @staticmethod
    def _validate_strict_ids(values, label, context):
        if np.any(values == 0):
            raise ValueError(f"{context} contains a zero {label}")
        if values.size > 1 and np.any(values[1:] <= values[:-1]):
            raise ValueError(
                f"{context} {label}s are duplicate or not strictly sorted")

    def _validate_state(self, state, context, slot_count=None):
        count = state.uid.size
        for name in state.__dataclass_fields__:
            values = getattr(state, name)
            expected = (count, 3) if name == "points" else (count,)
            if values.shape != expected:
                raise ValueError(
                    f"{context} column {name} has shape {values.shape}, "
                    f"expected {expected}")
        self._validate_strict_ids(state.uid, "UID", context)
        if np.unique(state.slot).size != count:
            raise ValueError(f"{context} contains duplicate cell slots")
        if slot_count is not None and (
                np.any(state.slot >= slot_count) if count else False):
            raise ValueError(f"{context} contains a cell slot outside slot_count")
        if np.any(state.stage >= self.radii.size):
            raise ValueError(f"{context} contains an invalid cell stage")
        if np.any(~np.isfinite(state.points)) or np.any(
                ~np.isfinite(state.display_radius)):
            raise ValueError(f"{context} contains non-finite display values")

    def _state_from_vtkhdf(self, path: Path):
        with h5py.File(path, "r") as handle:
            root = handle["/VTKHDF"]
            point_data = root["PointData"]
            uid = np.asarray(point_data["cell_id"], dtype=np.uint64)
            slot = (np.asarray(point_data["cell_slot"], dtype=np.uint32)
                    if "cell_slot" in point_data else
                    np.arange(uid.size, dtype=np.uint32))
            state = CellDisplayState(
                uid=uid,
                slot=slot,
                points=np.asarray(root["Points"], dtype=np.float32),
                lesion_id=np.asarray(point_data["lesion_id"], dtype=np.uint64),
                clone_id=np.asarray(point_data["clone_id"], dtype=np.uint32),
                cell_type=np.asarray(point_data["cell_type"], dtype=np.uint8),
                stage=np.asarray(point_data["stage"], dtype=np.uint8),
                viability=np.asarray(point_data["viability"], dtype=np.uint8),
                display_radius=np.asarray(
                    point_data["display_radius"], dtype=np.float32),
            )
        state = self._sorted(state)
        self._validate_state(state, f"VTK-HDF frame {path}")
        return state

    def _checkpoint_rows(self, handle: h5py.File, path="/cells",
                         slot_count=None):
        if path not in handle:
            raise ValueError(f"checkpoint is missing group {path}")
        cells = handle[path]
        context = f"checkpoint group {path}"
        uid = self._read_vector(cells, "uid", np.uint64, context)
        slot = self._read_vector(cells, "slot", np.uint32, context)
        x = self._read_vector(cells, "x", np.int32, context)
        y = self._read_vector(cells, "y", np.int32, context)
        z = self._read_vector(cells, "z", np.int32, context)
        clone_id = self._read_vector(
            cells, "clone_id", np.uint32, context)
        cell_type = self._read_vector(cells, "type", np.uint8, context)
        stage = self._read_vector(cells, "stage", np.uint8, context)
        viability = self._read_vector(
            cells, "viability", np.uint8, context)
        count = uid.size
        for name, values in (
                ("slot", slot), ("x", x), ("y", y), ("z", z),
                ("clone_id", clone_id), ("type", cell_type),
                ("stage", stage), ("viability", viability)):
            if values.size != count:
                raise ValueError(
                    f"{context} dataset {name} has {values.size} rows; "
                    f"expected {count}")
        self._validate_strict_ids(uid, "UID", context)
        if np.any(stage >= self.radii.size):
            raise ValueError(f"{context} contains an invalid cell stage")
        anchors = np.column_stack((x, y, z)).astype(np.float32, copy=False)
        offsets = np.where(stage == 0, 1.0, 0.5).astype(np.float32)
        points = anchors + offsets[:, None]
        state = CellDisplayState(
            uid=uid,
            slot=slot,
            points=points,
            # Lesion IDs are derived runtime state and are not required by the
            # current viewer.  Changed rows use 0 until the next VTK keyframe.
            lesion_id=np.zeros(uid.size, dtype=np.uint64),
            clone_id=clone_id,
            cell_type=cell_type,
            stage=stage,
            viability=viability,
            display_radius=self.radii[stage],
        )
        if slot_count is None:
            slot_count = self._slot_count(cells, context)
        self._validate_state(state, context, slot_count)
        return state

    @staticmethod
    def _require_disjoint(left, right, left_name, right_name, path):
        if left.size and right.size and np.intersect1d(
                left, right, assume_unique=True).size:
            raise ValueError(
                f"{path} contains a UID in both {left_name} and {right_name}")

    def _validate_parent(self, meta, path):
        if self.previous_checkpoint is None:
            raise ValueError(f"delta checkpoint has no materialized parent: {path}")
        if "parent_file" not in meta.attrs:
            raise ValueError(f"delta checkpoint has no parent_file: {path}")
        parent = meta.attrs["parent_file"]
        if isinstance(parent, bytes):
            parent = parent.decode("utf-8")
        parent = str(parent)
        parent_path = Path(parent)
        if (not parent or parent in (".", "..") or
                parent_path.name != parent or parent_path.is_absolute()):
            raise ValueError(f"checkpoint delta has unsafe parent_file: {path}")
        if parent != self.previous_checkpoint.name:
            raise ValueError(
                f"checkpoint delta chain is discontinuous at {path}")

    def _combine_after_delta(self, current, removed, births, slot_count,
                             context):
        keep = np.ones(current.uid.size, dtype=np.bool_)
        if removed.size:
            positions = np.searchsorted(current.uid, removed)
            valid = positions < current.uid.size
            if np.any(~valid) or np.any(
                    current.uid[positions[valid]] != removed[valid]):
                raise ValueError(f"{context} removes a missing UID")
            keep[positions] = False
        kept = current.take(keep)
        result = self._sorted(CellDisplayState(**{
            name: np.concatenate((getattr(kept, name),
                                  getattr(births, name)))
            for name in current.__dataclass_fields__
        }))
        self._validate_state(result, context, slot_count)
        return result

    def _apply_legacy_row_delta(self, handle, path):
        if self.state is None:
            raise ValueError("delta checkpoint has no materialized parent")
        cells = handle["/cells"]
        slot_count = self._slot_count(
            cells, f"legacy delta {path}")
        changes = self._checkpoint_rows(
            handle, "/cells", slot_count=slot_count)
        removed = self._read_vector(
            cells, "removed_uid", np.uint64, f"legacy delta {path}")
        self._validate_strict_ids(removed, "removed UID", f"legacy delta {path}")
        self._require_disjoint(
            removed, changes.uid, "removals", "changed rows", path)

        current = self.state.take(np.arange(self.state.uid.size))
        positions = np.searchsorted(current.uid, changes.uid)
        existing = positions < current.uid.size
        existing[existing] &= (
            current.uid[positions[existing]] == changes.uid[existing])
        existing_positions = positions[existing]
        if np.any(current.slot[existing_positions] != changes.slot[existing]):
            raise ValueError(f"legacy delta changes a stable cell slot: {path}")
        for name in current.__dataclass_fields__:
            getattr(current, name)[existing_positions] = \
                getattr(changes, name)[existing]
        births = changes.take(~existing)
        self.state = self._combine_after_delta(
            current, removed, births, slot_count, f"legacy delta {path}")

    def _read_field_updates(self, updates, path):
        context = f"field delta updates {path}"
        update_uid = self._read_vector(
            updates, "uid", np.uint64, context)
        update_slot = self._read_vector(
            updates, "slot", np.uint32, context)
        masks = self._read_vector(
            updates, "field_mask", np.uint64, context)
        if update_slot.size != update_uid.size or masks.size != update_uid.size:
            raise ValueError(f"{context} UID/slot/field_mask sizes differ")
        self._validate_strict_ids(update_uid, "updated UID", context)
        unknown_mask = ((1 << 64) - 1) ^ _KNOWN_CELL_FIELD_MASK
        if np.any(masks == 0) or np.any(
                np.bitwise_and(masks, np.uint64(unknown_mask))):
            raise ValueError(f"{context} contains an invalid field_mask")

        values = {}
        for field_name, bit, dataset_names in _CELL_FIELD_LAYOUT:
            expected = int(np.count_nonzero(
                np.bitwise_and(masks, np.uint64(1 << bit))))
            columns = []
            for dataset_name in dataset_names:
                if dataset_name not in updates:
                    raise ValueError(
                        f"{context} is missing dataset {dataset_name}")
                column = np.asarray(updates[dataset_name])
                if column.ndim != 1 or column.size != expected:
                    raise ValueError(
                        f"{context} dataset {dataset_name} has "
                        f"{column.size if column.ndim == 1 else 'invalid'} "
                        f"rows; expected {expected}")
                columns.append(column)
            values[field_name] = tuple(columns)
        return update_uid, update_slot, masks, values

    def _apply_field_delta(self, handle, path):
        if self.state is None:
            raise ValueError("delta checkpoint has no materialized parent")
        cells = handle["/cells"]
        slot_count = self._slot_count(
            cells, f"field delta {path}", required=True)
        removed = self._read_vector(
            cells, "removed_uid", np.uint64, f"field delta {path}")
        self._validate_strict_ids(removed, "removed UID", f"field delta {path}")
        births = self._checkpoint_rows(
            handle, "/cells/births", slot_count=slot_count)
        birth_slot_count = self._slot_count(
            cells["births"], f"field delta births {path}")
        if birth_slot_count is not None and birth_slot_count != slot_count:
            raise ValueError(
                f"{path} birth slot_count differs from delta slot_count")
        if "updates" not in cells:
            raise ValueError(f"field delta {path} has no updates group")
        update_uid, update_slot, masks, values = self._read_field_updates(
            cells["updates"], path)

        self._require_disjoint(
            removed, update_uid, "removals", "updates", path)
        self._require_disjoint(
            removed, births.uid, "removals", "births", path)
        self._require_disjoint(
            update_uid, births.uid, "updates", "births", path)
        if births.uid.size and np.intersect1d(
                self.state.uid, births.uid, assume_unique=True).size:
            raise ValueError(f"{path} reuses an existing cell UID for a birth")

        current = self.state.take(np.arange(self.state.uid.size))
        positions = np.searchsorted(current.uid, update_uid)
        valid = positions < current.uid.size
        if np.any(~valid) or np.any(
                current.uid[positions[valid]] != update_uid[valid]):
            raise ValueError(f"{path} updates a missing cell UID")
        if np.any(current.slot[positions] != update_slot):
            raise ValueError(f"{path} changes a stable cell slot")

        cursors = {name: 0 for name, _, _ in _CELL_FIELD_LAYOUT}
        for update_index, position in enumerate(positions):
            mask = int(masks[update_index])
            old_stage = int(current.stage[position])
            anchor = current.points[position] - (
                1.0 if old_stage == 0 else 0.5)
            for field_name, bit, _ in _CELL_FIELD_LAYOUT:
                if not (mask & (1 << bit)):
                    continue
                cursor = cursors[field_name]
                columns = values[field_name]
                field_value = tuple(column[cursor] for column in columns)
                cursors[field_name] = cursor + 1
                if field_name == "anchor":
                    anchor = np.asarray(field_value, dtype=np.float32)
                elif field_name == "clone_id":
                    current.clone_id[position] = field_value[0]
                elif field_name == "type":
                    current.cell_type[position] = field_value[0]
                elif field_name == "stage":
                    current.stage[position] = field_value[0]
                elif field_name == "viability":
                    current.viability[position] = field_value[0]
                # All non-display fields deliberately consume their compact
                # column entry so following rows remain aligned.
            stage = int(current.stage[position])
            if stage < 0 or stage >= self.radii.size:
                raise ValueError(f"{path} updates a cell to an invalid stage")
            current.points[position] = anchor + (
                1.0 if stage == 0 else 0.5)
            current.display_radius[position] = self.radii[stage]

        for field_name, _, _ in _CELL_FIELD_LAYOUT:
            expected = values[field_name][0].size
            if cursors[field_name] != expected:
                raise ValueError(
                    f"{path} did not consume all {field_name} update values")
        self.state = self._combine_after_delta(
            current, removed, births, slot_count, f"field delta {path}")

    def _apply_slot_journal(self, handle, path):
        """Apply schema-8 full rows keyed by stable slot.

        A changed slot may update the same UID, introduce a new slot, or reuse
        a slot whose previous cell died during the journal interval. Slots
        remaining empty at the end of the interval are listed separately.
        """
        if self.state is None:
            raise ValueError("slot journal has no materialized parent")
        cells = handle["/cells"]
        context = f"slot journal {path}"
        slot_count = self._slot_count(cells, context, required=True)
        changed = self._checkpoint_rows(
            handle, "/cells/changed", slot_count=slot_count)
        removed_slots = self._read_vector(
            cells, "removed_slot", np.uint32, context)
        if np.unique(removed_slots).size != removed_slots.size:
            raise ValueError(f"{context} contains duplicate removed slots")
        if removed_slots.size and np.any(removed_slots >= slot_count):
            raise ValueError(f"{context} removes a slot outside slot_count")
        if removed_slots.size and changed.slot.size and np.intersect1d(
                removed_slots, changed.slot, assume_unique=True).size:
            raise ValueError(
                f"{context} contains a slot in both changes and removals")

        current = self.state.take(np.arange(self.state.uid.size))
        slot_order = np.argsort(current.slot, kind="stable")
        sorted_slots = current.slot[slot_order]

        changed_positions = np.searchsorted(sorted_slots, changed.slot)
        changed_exists = changed_positions < sorted_slots.size
        changed_exists[changed_exists] &= (
            sorted_slots[changed_positions[changed_exists]] ==
            changed.slot[changed_exists])
        existing_rows = slot_order[changed_positions[changed_exists]]
        for name in current.__dataclass_fields__:
            getattr(current, name)[existing_rows] = getattr(
                changed, name)[changed_exists]
        births = changed.take(~changed_exists)

        keep = np.ones(current.uid.size, dtype=np.bool_)
        if removed_slots.size:
            removed_positions = np.searchsorted(sorted_slots, removed_slots)
            valid = removed_positions < sorted_slots.size
            valid[valid] &= (
                sorted_slots[removed_positions[valid]] ==
                removed_slots[valid])
            if np.any(~valid):
                raise ValueError(f"{context} removes a missing slot")
            keep[slot_order[removed_positions]] = False
        kept = current.take(keep)
        result = self._sorted(CellDisplayState(**{
            name: np.concatenate((
                getattr(kept, name), getattr(births, name)))
            for name in current.__dataclass_fields__
        }))
        self._validate_state(result, context, slot_count)
        self.state = result

    def _load_checkpoint_base(self, path: Path):
        with h5py.File(path, "r") as handle:
            schema = int(handle["/meta"].attrs["schema_version"])
            if schema not in _BASE_CHECKPOINT_SCHEMAS:
                raise ValueError(
                    f"expected checkpoint base schema 4 or 6: {path}")
            self.state = self._sorted(self._checkpoint_rows(handle))
        self.previous_checkpoint = path
        self.time_hours = _time_from_checkpoint_name(path)

    def _apply_checkpoint(self, path: Path):
        with h5py.File(path, "r") as handle:
            meta = handle["/meta"]
            schema = int(meta.attrs["schema_version"])
            if schema in _BASE_CHECKPOINT_SCHEMAS:
                self.state = self._sorted(self._checkpoint_rows(handle))
            elif schema == _LEGACY_ROW_DELTA_SCHEMA:
                self._validate_parent(meta, path)
                self._apply_legacy_row_delta(handle, path)
            elif schema == _FIELD_DELTA_SCHEMA:
                self._validate_parent(meta, path)
                self._apply_field_delta(handle, path)
            elif schema == _SLOT_JOURNAL_SCHEMA:
                self._validate_parent(meta, path)
                self._apply_slot_journal(handle, path)
            else:
                raise ValueError(f"unsupported checkpoint schema {schema}: {path}")
        self.previous_checkpoint = path
        self.time_hours = _time_from_checkpoint_name(path)

    def materialize(self, target_checkpoint: Path):
        target_checkpoint = Path(target_checkpoint).resolve()
        target_time = _time_from_checkpoint_name(target_checkpoint)
        checkpoints = _checkpoint_files(self.run_directory)
        checkpoint_at = {time: path for time, path in checkpoints}

        can_advance = (self.state is not None and self.time_hours is not None and
                       self.time_hours <= target_time)
        if not can_advance:
            keyframes = [entry for entry in _series_frames(
                self.run_directory, "full.vtkhdf.series")
                         if entry[0] <= target_time]
            if keyframes:
                base_time, base_path = keyframes[-1]
                self.state = self._state_from_vtkhdf(base_path)
                self.time_hours = base_time
                self.previous_checkpoint = checkpoint_at.get(base_time)
            else:
                bases = []
                for time_hours, path in checkpoints:
                    if time_hours > target_time:
                        break
                    with h5py.File(path, "r") as handle:
                        if int(handle["/meta"].attrs["schema_version"]) in \
                                _BASE_CHECKPOINT_SCHEMAS:
                            bases.append((time_hours, path))
                if not bases:
                    raise ValueError(
                        f"no full visualization/checkpoint base before {target_time}")
                self._load_checkpoint_base(bases[-1][1])

        for time_hours, path in checkpoints:
            if self.time_hours < time_hours <= target_time:
                self._apply_checkpoint(path)
        if self.time_hours != target_time or target_checkpoint != self.previous_checkpoint:
            raise ValueError(
                f"checkpoint chain cannot reconstruct exact time {target_time}")
        return self._to_polydata()

    def _to_polydata(self):
        from vtkmodules.util.numpy_support import numpy_to_vtk
        from vtkmodules.vtkCommonCore import vtkPoints, vtkUnsignedLongLongArray
        from vtkmodules.vtkCommonDataModel import vtkPolyData

        state = self.state
        if state is None:
            raise RuntimeError("no materialized checkpoint state")
        polydata = vtkPolyData()
        points = vtkPoints()
        points.SetData(numpy_to_vtk(state.points, deep=False))
        polydata.SetPoints(points)
        arrays = {
            "cell_id": state.uid,
            "cell_slot": state.slot,
            "lesion_id": state.lesion_id,
            "clone_id": state.clone_id,
            "cell_type": state.cell_type,
            "stage": state.stage,
            "viability": state.viability,
            "display_radius": state.display_radius,
        }
        vtk_arrays = []
        for name, values in arrays.items():
            vtk_array = numpy_to_vtk(values, deep=False)
            vtk_array.SetName(name)
            polydata.GetPointData().AddArray(vtk_array)
            vtk_arrays.append(vtk_array)
        total = vtkUnsignedLongLongArray()
        total.SetName("total_cell_count")
        total.SetNumberOfTuples(1)
        total.SetValue(0, int(state.uid.size))
        polydata.GetFieldData().AddArray(total)
        self._polydata_arrays = (state, points, vtk_arrays, total)
        return polydata


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Materialize one incremental ATCG3D checkpoint as VTK-HDF")
    parser.add_argument("run_directory", type=Path)
    parser.add_argument("checkpoint", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--compression-level", type=int, default=1)
    args = parser.parse_args(argv)

    materializer = CheckpointFrameMaterializer(args.run_directory)
    data = materializer.materialize(args.checkpoint)
    from vtkmodules.vtkIOHDF import vtkHDFWriter
    writer = vtkHDFWriter()
    writer.SetFileName(str(args.output))
    writer.SetInputData(data)
    writer.SetCompressionLevel(args.compression_level)
    if writer.Write() != 1:
        raise RuntimeError(f"failed to write {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
