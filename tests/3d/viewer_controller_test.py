import ast
import json
from pathlib import Path
import re
import struct
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from visualization.viewer.controller import FrameCounts, PreviewFullController, SeriesCatalog
from visualization.viewer.app import (
    CELL_TYPE_ANNOTATIONS,
    CELL_TYPE_INDEXED_COLORS,
    VESSEL_SOLID_COLOR,
    VIEW_SLICE_ORIENTATIONS,
    camera_aligned_plane,
    camera_basis,
    camera_oriented_plane,
    clipping_plane_specs,
    localized_viewer_status,
    normalize_viewer_language,
    normalized_plane_normal,
    oriented_clipping_plane_specs,
    projected_bounds,
    viewer_config,
    viewer_text,
)


class FakeBackend:
    def __init__(self):
        self.loads = []
        self.vessel_loads = []
        self.cancelled = []
        self.camera = ("camera", 1)

    def capture_camera(self):
        return self.camera

    def restore_camera(self, camera):
        self.camera = camera

    def load_frame(self, frame, quality, token):
        self.loads.append((frame.time_hours, quality, token))
        return FrameCounts(10_000_000, 100 if quality == "preview" else 10_000_000)

    def load_vessel_frame(self, frame, token):
        self.vessel_loads.append((None if frame is None else frame.time_hours, token))
        return 0 if frame is None else 25

    def cancel_full(self, token):
        self.cancelled.append(token)


def write_series(path: Path, subdirectory: str, times):
    path.write_text(json.dumps({
        "file-series-version": "1.0",
        "files": [
            {"name": f"viz/{subdirectory}/frame_{index:08d}.vtkhdf", "time": value}
            for index, value in enumerate(times)
        ],
    }), encoding="utf-8")


class ViewerControllerTest(unittest.TestCase):
    def test_studio_languages_are_complete_and_default_to_english(self):
        root = Path(__file__).resolve().parents[2]
        html_source = (
            root / "visualization/studio/web/index.html"
        ).read_text(encoding="utf-8")
        javascript = (
            root / "visualization/studio/web/app.js"
        ).read_text(encoding="utf-8")
        self.assertIn('<html lang="en">', html_source)
        self.assertIn('currentLanguage = "en"', javascript)

        english = javascript.split("  en: {", 1)[1].split(
            '  "zh-CN": {', 1)[0]
        simplified_chinese = javascript.split(
            '  "zh-CN": {', 1)[1].split("\n  },\n};", 1)[0]
        key_pattern = re.compile(r'^\s+"([^"]+)":', re.MULTILINE)
        english_keys = set(key_pattern.findall(english))
        chinese_keys = set(key_pattern.findall(simplified_chinese))
        self.assertEqual(english_keys, chinese_keys)

        markup_keys = set(re.findall(
            r'data-i18n(?:-placeholder|-title)?="([^"]+)"', html_source))
        self.assertTrue(markup_keys)
        self.assertEqual(markup_keys - english_keys, set())

    def test_studio_allows_local_viewer_and_viewer_opens_visible_live_tail(self):
        root = Path(__file__).resolve().parents[2]
        plist = (
            root / "visualization/studio/desktop/src-tauri/Info.plist"
        ).read_text(encoding="utf-8")
        self.assertIn("NSAllowsArbitraryLoadsInWebContent", plist)
        self.assertIn("NSAllowsLocalNetworking", plist)
        source = (
            root / "visualization/viewer/app.py"
        ).read_text(encoding="utf-8")
        self.assertIn("load_index(len(state.times) - 1)", source)
        self.assertIn("state.radius_scale = 4.0", source)
        self.assertIn(
            "backend.set_cell_radius_scale(state.radius_scale)", source)

    def test_viewer_language_defaults_and_translations(self):
        self.assertEqual(normalize_viewer_language("en"), "en")
        self.assertEqual(viewer_text("en", "camera"), "Camera")
        self.assertEqual(viewer_text("zh-CN", "camera"), "相机")
        self.assertEqual(
            localized_viewer_status("zh-CN", "quality", "full"), "完整帧")
        self.assertEqual(
            localized_viewer_status("zh-CN", "quality", "future status"),
            "future status")
        with self.assertRaises(ValueError):
            normalize_viewer_language("fr")

    def test_fixed_semantic_colors(self):
        self.assertEqual(CELL_TYPE_ANNOTATIONS, ["1", "r", "2", "K"])
        self.assertEqual(CELL_TYPE_INDEXED_COLORS[:3], [0.0, 0.72, 0.0])
        self.assertEqual(CELL_TYPE_INDEXED_COLORS[3:], [0.90, 0.0, 0.0])
        self.assertEqual(VESSEL_SOLID_COLOR, [0.05, 0.25, 1.0])

    def test_whole_cut_and_slab_plane_specs(self):
        self.assertEqual(clipping_plane_specs("whole", "Z", 3.0, 2.0), [])
        self.assertEqual(
            clipping_plane_specs("cut", "X", 3.0, 2.0, True),
            [([3.0, 0.0, 0.0], [1.0, 0.0, 0.0], 1)],
        )
        self.assertEqual(
            clipping_plane_specs("slab", "Z", 3.0, 2.0),
            [
                ([0.0, 0.0, 2.0], [0.0, 0.0, 1.0], 0),
                ([0.0, 0.0, 4.0], [0.0, 0.0, 1.0], 1),
            ],
        )
        with self.assertRaises(ValueError):
            clipping_plane_specs("slab", "Z", 0.0, 0.0)

    def test_view_aligned_plane_is_normalized_and_stays_frozen(self):
        normal, position = camera_aligned_plane(
            (10.0, 10.0, 10.0), (1.0, 1.0, 1.0))
        self.assertAlmostEqual(sum(value * value for value in normal), 1.0)
        self.assertAlmostEqual(position, sum(normal))
        initial = oriented_clipping_plane_specs(
            "slab", normal, position, 4.0)
        # Moving the slice changes only its offset; an unrelated later camera
        # pose is deliberately absent from the generalized plane API.
        moved = oriented_clipping_plane_specs(
            "slab", normal, position + 2.0, 4.0)
        self.assertEqual(initial[0][1], moved[0][1])
        self.assertEqual(initial[1][1], moved[1][1])
        self.assertNotEqual(initial[0][0], moved[0][0])
        with self.assertRaises(ValueError):
            normalized_plane_normal((0.0, 0.0, 0.0))

    def test_view_basis_supports_axial_sagittal_and_coronal_planes(self):
        forward, right, up = camera_basis(
            (0.0, 0.0, 10.0), (0.0, 0.0, 0.0), (0.0, 1.0, 0.0))
        self.assertEqual(forward, (0.0, 0.0, -1.0))
        self.assertEqual(right, (1.0, -0.0, 0.0))
        self.assertEqual(up, (0.0, 1.0, 0.0))
        expected = {
            "VIEW_AXIAL": forward,
            "VIEW_SAGITTAL": right,
            "VIEW_CORONAL": up,
        }
        for orientation, expected_normal in expected.items():
            normal, offset = camera_oriented_plane(
                (0.0, 0.0, 10.0), (0.0, 0.0, 0.0),
                (0.0, 1.0, 0.0), orientation)
            self.assertEqual(normal, expected_normal)
            self.assertEqual(offset, 0.0)
        self.assertEqual(set(expected), set(VIEW_SLICE_ORIENTATIONS))

    def test_oblique_sagittal_plane_is_orthogonal_and_frozen(self):
        position = (9.0, -4.0, 7.0)
        focal = (1.0, 2.0, -3.0)
        view_up = (0.2, 1.0, 0.3)
        forward, right, up = camera_basis(position, focal, view_up)
        for left, other in ((forward, right), (forward, up), (right, up)):
            self.assertAlmostEqual(sum(a * b for a, b in zip(left, other)), 0.0)
        normal, anchor = camera_oriented_plane(
            position, focal, view_up, "VIEW_SAGITTAL")
        self.assertEqual(normal, right)
        initial = oriented_clipping_plane_specs("cut", normal, anchor, 1.0)
        moved = oriented_clipping_plane_specs("cut", normal, anchor + 5.0, 1.0)
        self.assertEqual(initial[0][1], moved[0][1])
        self.assertNotEqual(initial[0][0], moved[0][0])
        pole_forward, pole_right, pole_up = camera_basis(
            (0, 0, 1), (0, 0, 0), (0, 0, 2))
        self.assertAlmostEqual(
            sum(a * b for a, b in zip(pole_forward, pole_right)), 0.0)
        self.assertAlmostEqual(
            sum(a * b for a, b in zip(pole_forward, pole_up)), 0.0)

    def test_parameter_controls_are_in_the_side_drawer(self):
        app_path = Path(__file__).resolve().parents[2] / "visualization/viewer/app.py"
        source = app_path.read_text(encoding="utf-8")
        drawer = source.index("with layout.drawer:")
        content = source.index("with layout.content:")
        for state_name in (
                "navigation_mode", "radius_scale", "vessel_radius_scale",
                "show_r", "show_k", "show_vessels", "show_influence",
                "slice_mode", "slice_axis", "slice_position",
                "slice_thickness", "slice_invert"):
            control = source.index(f'v_model=("{state_name}"', drawer)
            self.assertLess(control, content)

    def test_projected_bounds_cover_all_corners(self):
        normal = normalized_plane_normal((1.0, 1.0, 1.0))
        minimum, maximum = projected_bounds((-1, 2, -3, 4, -5, 6), normal)
        self.assertAlmostEqual(minimum, (-1 - 3 - 5) / (3 ** 0.5))
        self.assertAlmostEqual(maximum, (2 + 4 + 6) / (3 ** 0.5))

    def test_viewer_reads_influence_cutoff_from_run_metadata(self):
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            (run / "run.json").write_text(json.dumps({
                "effective_config": {
                    "angiogenesis_influence_cutoff_radius_voxels": 12.0,
                },
            }), encoding="utf-8")
            self.assertEqual(
                viewer_config(run)["influence_cutoff_radius_voxels"], 12.0)

    def test_server_ready_callback_accepts_trame_state_keywords(self):
        app_path = Path(__file__).resolve().parents[2] / "visualization/viewer/app.py"
        tree = ast.parse(app_path.read_text(encoding="utf-8"))
        callbacks = [
            node for node in ast.walk(tree)
            if isinstance(node, ast.AsyncFunctionDef) and node.name == "update_loop"
        ]
        self.assertEqual(len(callbacks), 1)
        self.assertIsNotNone(callbacks[0].args.kwarg)

    def test_fast_slider_only_loads_last_full(self):
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            times = [float(value) for value in range(100)]
            write_series(run / "preview.vtkhdf.series", "preview", times)
            write_series(run / "full.vtkhdf.series", "full", times)
            write_series(run / "vessels.vtkhdf.series", "vessels", times)
            backend = FakeBackend()
            controller = PreviewFullController(SeriesCatalog(run), backend, 250)
            for index, value in enumerate(times):
                controller.on_slider(value, index)
            controller.tick(99 + 249)
            self.assertEqual(sum(q == "full" for _, q, _ in backend.loads), 0)
            controller.tick(99 + 250)
            full = [entry for entry in backend.loads if entry[1] == "full"]
            self.assertEqual(len(full), 1)
            self.assertEqual(full[0][0], 99.0)
            self.assertEqual(controller.current_count, 10_000_000)
            self.assertEqual(controller.current_displayed_count, 10_000_000)
            self.assertEqual(controller.current_vessel_count, 25)
            self.assertEqual(backend.vessel_loads[-1][0], 99.0)
            self.assertEqual(backend.camera, ("camera", 1))

    def test_slider_during_full_discards_stale_completion(self):
        class ReentrantBackend(FakeBackend):
            controller = None
            triggered = False

            def load_frame(self, frame, quality, token):
                result = super().load_frame(frame, quality, token)
                if quality == "full" and not self.triggered:
                    self.triggered = True
                    self.controller.on_slider(1.0, 300)
                return result

        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            write_series(run / "preview.vtkhdf.series", "preview", [0.0, 1.0])
            write_series(run / "full.vtkhdf.series", "full", [0.0, 1.0])
            write_series(run / "vessels.vtkhdf.series", "vessels", [0.0, 1.0])
            backend = ReentrantBackend()
            controller = PreviewFullController(SeriesCatalog(run), backend, 250)
            backend.controller = controller
            controller.on_slider(0.0, 0)
            self.assertEqual(controller.tick(250), "preview")
            self.assertEqual(controller.current_time, 1.0)
            self.assertEqual(controller.current_count, 10_000_000)
            self.assertEqual(controller.current_displayed_count, 100)
            self.assertEqual(backend.vessel_loads[-1][0], 1.0)
            self.assertEqual(controller.tick(550), "full")
            self.assertEqual(controller.current_count, 10_000_000)
            self.assertEqual(controller.current_displayed_count, 10_000_000)

    def test_missing_exact_full_keeps_preview(self):
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            write_series(run / "preview.vtkhdf.series", "preview", [0.0, 1.0])
            write_series(run / "full.vtkhdf.series", "full", [0.0])
            write_series(run / "vessels.vtkhdf.series", "vessels", [0.0])
            backend = FakeBackend()
            controller = PreviewFullController(SeriesCatalog(run), backend, 250)
            controller.on_slider(1.0, 0)
            self.assertEqual(controller.tick(250), "preview (no exact full frame)")
            self.assertEqual(sum(q == "full" for _, q, _ in backend.loads), 0)
            self.assertEqual(controller.current_vessel_count, 0)
            self.assertEqual(controller.current_vessel_status, "vessels unavailable")

    def test_atomic_series_refresh_adds_times(self):
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            preview = run / "preview.vtkhdf.series"
            full = run / "full.vtkhdf.series"
            vessels = run / "vessels.vtkhdf.series"
            write_series(preview, "preview", [0.0])
            write_series(full, "full", [])
            write_series(vessels, "vessels", [0.0])
            catalog = SeriesCatalog(run)
            self.assertEqual(catalog.times, [0.0])
            self.assertEqual(catalog.exact_vessel(0.0).time_hours, 0.0)
            write_series(preview, "preview", [0.0, 1.0])
            write_series(vessels, "vessels", [0.0, 1.0])
            catalog.refresh()
            self.assertEqual(catalog.times, [0.0, 1.0])
            self.assertEqual(catalog.exact_vessel(1.0).time_hours, 1.0)

    def test_live_series_replaces_only_the_live_tail(self):
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            write_series(
                run / "preview.vtkhdf.series", "preview", [0.0, 1.0])
            write_series(
                run / "full.vtkhdf.series", "full", [0.0])
            write_series(
                run / "vessels.vtkhdf.series", "vessels", [0.0, 1.0])
            write_series(
                run / "live.vtkhdf.series", "live", [1.5])
            write_series(
                run / "live-vessels.vtkhdf.series", "live", [1.5])
            catalog = SeriesCatalog(run)
            self.assertEqual(catalog.times, [0.0, 1.0, 1.5])
            self.assertEqual(catalog.exact_vessel(1.5).time_hours, 1.5)

            # The simulator atomically overwrites one live series entry.  The
            # viewer must drop the old live instant without disturbing the
            # archived timeline.
            write_series(
                run / "live.vtkhdf.series", "live", [2.0])
            write_series(
                run / "live-vessels.vtkhdf.series", "live", [2.0])
            self.assertTrue(catalog.refresh())
            self.assertEqual(catalog.times, [0.0, 1.0, 2.0])
            self.assertIsNone(catalog.exact_vessel(1.5))
            self.assertEqual(catalog.exact_vessel(2.0).time_hours, 2.0)

    def test_exact_full_falls_back_to_atomic_checkpoint_delta(self):
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary)
            write_series(run / "preview.vtkhdf.series", "preview", [1.0])
            write_series(run / "full.vtkhdf.series", "full", [])
            write_series(run / "vessels.vtkhdf.series", "vessels", [1.0])
            checkpoints = run / "checkpoints"
            checkpoints.mkdir()
            bits = struct.unpack(">Q", struct.pack(">d", 1.0))[0]
            checkpoint = checkpoints / (
                f"checkpoint_{1:016d}_time_{bits:016x}.h5")
            checkpoint.touch()
            catalog = SeriesCatalog(run)
            frame = catalog.exact_full(1.0)
            self.assertIsNotNone(frame)
            self.assertEqual(frame.path, checkpoint.resolve())
            self.assertEqual(frame.storage, "checkpoint")


if __name__ == "__main__":
    unittest.main()
