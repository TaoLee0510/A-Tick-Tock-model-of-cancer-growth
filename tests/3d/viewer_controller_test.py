import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from visualization.viewer.controller import FrameCounts, PreviewFullController, SeriesCatalog


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


if __name__ == "__main__":
    unittest.main()
