import tempfile
from pathlib import Path
import struct
import sys
import unittest

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from visualization.viewer.checkpoint_materializer import (
    _CELL_FIELD_LAYOUT,
    CheckpointFrameMaterializer,
)

_FIELD_DTYPES = {
    "x": np.int32,
    "y": np.int32,
    "z": np.int32,
    "parent_uid": np.uint64,
    "clone_id": np.uint32,
    "type": np.uint8,
    "stage": np.uint8,
    "viability": np.uint8,
    "flags": np.uint8,
    "last_direction": np.uint8,
    "inherent_growth_rate": np.float32,
    "density_growth_rate": np.float32,
    "migration_rate": np.float32,
    "normal_migration_rate": np.float32,
    "migration_activation_end_time": np.float64,
    "division_work_remaining": np.float32,
    "next_migration_time": np.float64,
    "next_division_time": np.float64,
    "death_deadline": np.float32,
    "last_update_time": np.float64,
    "swap_ready_time": np.float64,
    "swap_wait_state": np.uint8,
    "pending_swap_direction": np.uint8,
    "event_sequence": np.uint64,
    "migration_schedule_generation": np.uint32,
    "division_schedule_generation": np.uint32,
    "death_schedule_generation": np.uint32,
}


def checkpoint_name(index, time_hours):
    bits = struct.unpack(">Q", struct.pack(">d", time_hours))[0]
    return f"checkpoint_{index:016d}_time_{bits:016x}.h5"


def write_cell_rows(cells, rows, slot_count=None):
    columns = list(zip(*rows)) if rows else [[] for _ in range(9)]
    names_and_types = (
        ("uid", np.uint64), ("slot", np.uint32),
        ("x", np.int32), ("y", np.int32), ("z", np.int32),
        ("clone_id", np.uint32), ("type", np.uint8),
        ("stage", np.uint8), ("viability", np.uint8),
    )
    for values, (name, dtype) in zip(columns, names_and_types):
        cells.create_dataset(name, data=np.asarray(values, dtype=dtype))
    if slot_count is not None:
        cells.attrs["slot_count"] = np.uint64(slot_count)


def write_cells(handle, rows, removed=(), slot_count=None):
    cells = handle.create_group("cells")
    write_cell_rows(cells, rows, slot_count)
    cells.create_dataset("removed_uid", data=np.asarray(removed, dtype=np.uint64))


def write_field_delta(handle, parent_file, slot_count, births=(), removed=(),
                      updates=()):
    meta = handle.create_group("meta")
    meta.attrs["schema_version"] = np.uint32(7)
    meta.attrs["parent_file"] = parent_file
    cells = handle.create_group("cells")
    cells.attrs["slot_count"] = np.uint64(slot_count)
    cells.create_dataset(
        "free_slots", data=np.asarray([], dtype=np.uint32))
    cells.create_dataset(
        "removed_uid", data=np.asarray(removed, dtype=np.uint64))
    births_group = cells.create_group("births")
    write_cell_rows(births_group, births, slot_count)
    births_group.create_dataset(
        "free_slots", data=np.asarray([], dtype=np.uint32))

    update_group = cells.create_group("updates")
    update_group.create_dataset(
        "uid", data=np.asarray([entry[0] for entry in updates],
                               dtype=np.uint64))
    update_group.create_dataset(
        "slot", data=np.asarray([entry[1] for entry in updates],
                                dtype=np.uint32))
    update_group.create_dataset(
        "field_mask", data=np.asarray([entry[2] for entry in updates],
                                      dtype=np.uint64))
    compact_columns = {name: [] for name in _FIELD_DTYPES}
    for _, _, mask, fields in updates:
        for field_name, bit, dataset_names in _CELL_FIELD_LAYOUT:
            if not (mask & (1 << bit)):
                continue
            value = fields[field_name]
            if len(dataset_names) == 1:
                compact_columns[dataset_names[0]].append(value)
            else:
                for dataset_name, component in zip(dataset_names, value):
                    compact_columns[dataset_name].append(component)
    for name, dtype in _FIELD_DTYPES.items():
        update_group.create_dataset(
            name, data=np.asarray(compact_columns[name], dtype=dtype))


def write_slot_journal(handle, parent_file, slot_count, changed=(),
                       removed_slots=()):
    meta = handle.create_group("meta")
    meta.attrs["schema_version"] = np.uint32(8)
    meta.attrs["kind"] = "slot_journal_v1"
    meta.attrs["parent_file"] = parent_file
    cells = handle.create_group("cells")
    cells.attrs["slot_count"] = np.uint64(slot_count)
    cells.create_dataset(
        "removed_slot", data=np.asarray(removed_slots, dtype=np.uint32))
    changed_group = cells.create_group("changed")
    write_cell_rows(changed_group, changed, slot_count)
    changed_group.create_dataset(
        "free_slots", data=np.asarray([], dtype=np.uint32))
    mutations = cells.create_group("free_list_mutations")
    mutations.create_dataset(
        "kind", data=np.asarray([], dtype=np.uint8))
    mutations.create_dataset(
        "slot", data=np.asarray([], dtype=np.uint32))


class CheckpointMaterializerTest(unittest.TestCase):
    def test_update_birth_death_and_slot_reuse(self):
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            checkpoints = run / "checkpoints"
            checkpoints.mkdir()
            base = checkpoints / checkpoint_name(0, 0.0)
            delta = checkpoints / checkpoint_name(1, 1.0)
            with h5py.File(base, "w") as handle:
                meta = handle.create_group("meta")
                meta.attrs["schema_version"] = np.uint32(4)
                write_cells(handle, [
                    (1, 0, 0, 0, 0, 11, 1, 1, 1),
                    (2, 1, 1, 0, 0, 12, 2, 0, 1),
                ])
            with h5py.File(delta, "w") as handle:
                meta = handle.create_group("meta")
                meta.attrs["schema_version"] = np.uint32(5)
                meta.attrs["parent_file"] = base.name
                write_cells(handle, [
                    (2, 1, 5, 0, 0, 12, 2, 1, 1),
                    (3, 0, 2, 0, 0, 13, 1, 2, 1),
                ], removed=(1,))

            materializer = CheckpointFrameMaterializer(run, (1.0, 0.5, 0.25))
            materializer._load_checkpoint_base(base)
            materializer._apply_checkpoint(delta)
            state = materializer.state
            np.testing.assert_array_equal(state.uid, [2, 3])
            np.testing.assert_array_equal(state.slot, [1, 0])
            np.testing.assert_allclose(state.points, [
                [5.5, 0.5, 0.5], [2.5, 0.5, 0.5],
            ])
            np.testing.assert_allclose(state.display_radius, [0.5, 0.25])
            self.assertEqual(materializer.previous_checkpoint, delta)
            self.assertEqual(materializer.time_hours, 1.0)
            polydata = materializer._to_polydata()
            self.assertEqual(polydata.GetNumberOfPoints(), 2)
            self.assertEqual(polydata.GetNumberOfCells(), 0)
            self.assertIsNotNone(polydata.GetPointData().GetArray("cell_id"))
            self.assertIsNotNone(polydata.GetPointData().GetArray("cell_slot"))
            self.assertEqual(
                polydata.GetFieldData().GetArray("total_cell_count").GetValue(0),
                2,
            )
            from paraview import simple
            producer = simple.TrivialProducer()
            producer.GetClientSideObject().SetOutput(polydata)
            producer.UpdatePipeline()
            self.assertEqual(
                producer.GetDataInformation().GetNumberOfPoints(), 2)
            simple.Delete(producer)

    def test_v6_v7_sparse_fields_birth_death_and_slot_reuse(self):
        """This reconstruction path requires only h5py/numpy, not VTK."""
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            checkpoints = run / "checkpoints"
            checkpoints.mkdir()
            base = checkpoints / checkpoint_name(0, 0.0)
            delta = checkpoints / checkpoint_name(1, 1.0)
            with h5py.File(base, "w") as handle:
                meta = handle.create_group("meta")
                meta.attrs["schema_version"] = np.uint32(6)
                write_cells(handle, [
                    (1, 0, 0, 0, 0, 11, 1, 1, 1),
                    (2, 1, 1, 2, 3, 12, 2, 0, 1),
                    (4, 2, 3, 0, 0, 14, 1, 1, 1),
                ], slot_count=3)

            update_2_mask = (
                (1 << 0) |   # anchor
                (1 << 1) |   # parent_uid, not displayed
                (1 << 2) |   # clone_id
                (1 << 3) |   # type
                (1 << 4) |   # stage
                (1 << 5)     # viability
            )
            update_4_mask = (
                (1 << 4) |   # independently compacted second stage value
                (1 << 6) |   # flags, not displayed
                (1 << 21)    # event_sequence, not displayed
            )
            with h5py.File(delta, "w") as handle:
                write_field_delta(
                    handle,
                    base.name,
                    slot_count=3,
                    births=[
                        (3, 0, 2, 2, 2, 13, 1, 2, 1),
                    ],
                    removed=(1,),
                    updates=[
                        (2, 1, update_2_mask, {
                            "anchor": (5, 6, 7),
                            "parent_uid": 1,
                            "clone_id": 22,
                            "type": 1,
                            "stage": 1,
                            "viability": 0,
                        }),
                        (4, 2, update_4_mask, {
                            "stage": 2,
                            "flags": 7,
                            "event_sequence": 99,
                        }),
                    ],
                )

            materializer = CheckpointFrameMaterializer(
                run, (1.0, 0.5, 0.25))
            materializer._load_checkpoint_base(base)
            materializer._apply_checkpoint(delta)
            state = materializer.state
            np.testing.assert_array_equal(state.uid, [2, 3, 4])
            # UID 1 dies and UID 3 reuses its now-free stable slot 0.
            np.testing.assert_array_equal(state.slot, [1, 0, 2])
            np.testing.assert_allclose(state.points, [
                [5.5, 6.5, 7.5],
                [2.5, 2.5, 2.5],
                [3.5, 0.5, 0.5],
            ])
            np.testing.assert_array_equal(state.clone_id, [22, 13, 14])
            np.testing.assert_array_equal(state.cell_type, [1, 1, 1])
            np.testing.assert_array_equal(state.stage, [1, 2, 2])
            np.testing.assert_array_equal(state.viability, [0, 1, 1])
            np.testing.assert_allclose(
                state.display_radius, [0.5, 0.25, 0.25])
            self.assertEqual(materializer.previous_checkpoint, delta)
            self.assertEqual(materializer.time_hours, 1.0)

    def test_v7_rejects_chain_mask_and_identity_conflicts(self):
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            checkpoints = run / "checkpoints"
            checkpoints.mkdir()
            base = checkpoints / checkpoint_name(0, 0.0)
            with h5py.File(base, "w") as handle:
                meta = handle.create_group("meta")
                meta.attrs["schema_version"] = np.uint32(6)
                write_cells(handle, [
                    (1, 0, 0, 0, 0, 11, 1, 1, 1),
                    (2, 1, 1, 0, 0, 12, 2, 1, 1),
                ], slot_count=2)

            cases = (
                ("wrong_parent", "not_the_parent.h5", (), (), (
                    (2, 1, 1 << 4, {"stage": 2}),
                ), "discontinuous"),
                ("unknown_mask", base.name, (), (), (
                    (2, 1, 1 << 25, {}),
                ), "field_mask"),
                ("remove_update", base.name, (), (2,), (
                    (2, 1, 1 << 4, {"stage": 2}),
                ), "both removals and updates"),
                ("unstable_slot", base.name, (), (), (
                    (2, 0, 1 << 4, {"stage": 2}),
                ), "stable cell slot"),
                ("birth_slot_collision", base.name, (
                    (3, 1, 2, 0, 0, 13, 1, 1, 1),
                ), (), (), "duplicate cell slots"),
            )
            for index, (name, parent, births, removed, updates,
                        message) in enumerate(cases, start=1):
                with self.subTest(name=name):
                    delta = checkpoints / checkpoint_name(index, float(index))
                    with h5py.File(delta, "w") as handle:
                        write_field_delta(
                            handle, parent, 2, births, removed, updates)
                    materializer = CheckpointFrameMaterializer(run)
                    materializer._load_checkpoint_base(base)
                    with self.assertRaisesRegex(ValueError, message):
                        materializer._apply_checkpoint(delta)

    def test_v8_slot_journal_updates_reuses_and_removes_slots(self):
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            checkpoints = run / "checkpoints"
            checkpoints.mkdir()
            base = checkpoints / checkpoint_name(0, 0.0)
            delta = checkpoints / checkpoint_name(1, 1.0)
            with h5py.File(base, "w") as handle:
                meta = handle.create_group("meta")
                meta.attrs["schema_version"] = np.uint32(6)
                write_cells(handle, [
                    (1, 0, 0, 0, 0, 11, 1, 1, 1),
                    (2, 1, 1, 0, 0, 12, 2, 0, 1),
                    (4, 2, 4, 0, 0, 14, 1, 1, 1),
                ], slot_count=3)
            with h5py.File(delta, "w") as handle:
                write_slot_journal(
                    handle,
                    base.name,
                    slot_count=3,
                    changed=[
                        # UID 2 updates in its stable slot.
                        (2, 1, 5, 6, 7, 22, 1, 1, 0),
                        # UID 3 reuses slot 0 after UID 1 dies.
                        (3, 0, 2, 2, 2, 13, 1, 2, 1),
                    ],
                    # UID 4 dies and slot 2 remains empty.
                    removed_slots=(2,),
                )

            materializer = CheckpointFrameMaterializer(run)
            materializer._load_checkpoint_base(base)
            materializer._apply_checkpoint(delta)
            state = materializer.state
            np.testing.assert_array_equal(state.uid, [2, 3])
            np.testing.assert_array_equal(state.slot, [1, 0])
            np.testing.assert_allclose(state.points, [
                [5.5, 6.5, 7.5],
                [2.5, 2.5, 2.5],
            ])
            np.testing.assert_array_equal(state.clone_id, [22, 13])
            np.testing.assert_array_equal(state.cell_type, [1, 1])
            np.testing.assert_array_equal(state.stage, [1, 2])
            np.testing.assert_array_equal(state.viability, [0, 1])
            np.testing.assert_allclose(state.display_radius, [0.5, 0.25])
            self.assertEqual(materializer.previous_checkpoint, delta)


if __name__ == "__main__":
    unittest.main()
