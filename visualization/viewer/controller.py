"""Dependency-free timeline state machine used by the trame viewer.

The module deliberately has no trame/VTK imports so debounce and live-series
behavior remain unit-testable on build machines without visualization packages.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Protocol


@dataclass(frozen=True)
class Frame:
    path: Path
    time_hours: float


@dataclass(frozen=True)
class FrameCounts:
    total: int
    displayed: int


class FrameBackend(Protocol):
    def capture_camera(self): ...
    def restore_camera(self, camera) -> None: ...
    def load_frame(self, frame: Frame, quality: str, token: int) -> FrameCounts: ...
    def load_vessel_frame(self, frame: Frame | None, token: int) -> int: ...
    def cancel_full(self, token: int) -> None: ...


class SeriesCatalog:
    def __init__(self, run_directory: Path | str):
        self.run_directory = Path(run_directory).resolve()
        self.preview: list[Frame] = []
        self.full: list[Frame] = []
        self.vessels: list[Frame] = []
        self._signature: tuple[tuple[int, int], tuple[int, int], tuple[int, int]] | None = None
        self.refresh()

    def _file_signature(self, path: Path) -> tuple[int, int]:
        try:
            stat = path.stat()
            return stat.st_mtime_ns, stat.st_size
        except FileNotFoundError:
            return 0, 0

    def _read(self, name: str) -> list[Frame]:
        path = self.run_directory / name
        if not path.exists():
            return []
        document = json.loads(path.read_text(encoding="utf-8"))
        if document.get("file-series-version") != "1.0":
            raise ValueError(f"unsupported ParaView series version in {path}")
        frames: list[Frame] = []
        for entry in document.get("files", []):
            relative = str(entry["name"])
            if relative.endswith(".tmp") or ".tmp/" in relative:
                raise ValueError(f"series references temporary output: {relative}")
            frames.append(Frame(self.run_directory / relative, float(entry["time"])))
        frames.sort(key=lambda frame: frame.time_hours)
        return frames

    def refresh(self) -> bool:
        preview_path = self.run_directory / "preview.vtkhdf.series"
        full_path = self.run_directory / "full.vtkhdf.series"
        vessel_path = self.run_directory / "vessels.vtkhdf.series"
        signature = (self._file_signature(preview_path), self._file_signature(full_path),
                     self._file_signature(vessel_path))
        if signature == self._signature:
            return False
        self.preview = self._read("preview.vtkhdf.series")
        self.full = self._read("full.vtkhdf.series")
        self.vessels = self._read("vessels.vtkhdf.series")
        self._signature = signature
        return True

    def nearest_preview(self, time_hours: float) -> Frame | None:
        if not self.preview:
            return None
        return min(self.preview, key=lambda frame: (abs(frame.time_hours - time_hours),
                                                     frame.time_hours))

    def exact_full(self, time_hours: float, tolerance: float = 1e-9) -> Frame | None:
        for frame in self.full:
            scale = max(1.0, abs(frame.time_hours), abs(time_hours))
            if abs(frame.time_hours - time_hours) <= tolerance * scale:
                return frame
        return None

    def exact_vessel(self, time_hours: float, tolerance: float = 1e-9) -> Frame | None:
        for frame in self.vessels:
            scale = max(1.0, abs(frame.time_hours), abs(time_hours))
            if abs(frame.time_hours - time_hours) <= tolerance * scale:
                return frame
        return None

    @property
    def times(self) -> list[float]:
        return [frame.time_hours for frame in self.preview]


class PreviewFullController:
    def __init__(self, catalog: SeriesCatalog, backend: FrameBackend,
                 debounce_ms: int = 250):
        if debounce_ms < 0:
            raise ValueError("debounce_ms must be non-negative")
        self.catalog = catalog
        self.backend = backend
        self.debounce_ms = debounce_ms
        self.current_time: float | None = None
        self.current_quality = "none"
        self.current_count = 0
        self.current_displayed_count = 0
        self.current_vessel_count = 0
        self.current_vessel_status = "vessels unavailable"
        self._token = 0
        self._last_slider_ms = 0
        self._full_token_done: int | None = None

    @property
    def token(self) -> int:
        return self._token

    def on_slider(self, time_hours: float, now_ms: int) -> str:
        previous = self._token
        self._token += 1
        if previous:
            self.backend.cancel_full(previous)
        self.current_time = float(time_hours)
        self._last_slider_ms = int(now_ms)
        self._full_token_done = None
        self.catalog.refresh()
        frame = self.catalog.nearest_preview(self.current_time)
        if frame is None:
            self.current_quality = "preview unavailable"
            self.current_count = 0
            self.current_displayed_count = 0
            self.current_vessel_count = self.backend.load_vessel_frame(None, self._token)
            self.current_vessel_status = "vessels unavailable"
            return self.current_quality
        camera = self.backend.capture_camera()
        counts = self.backend.load_frame(frame, "preview", self._token)
        self.current_count = counts.total
        self.current_displayed_count = counts.displayed
        vessel = self.catalog.exact_vessel(frame.time_hours)
        self.current_vessel_count = self.backend.load_vessel_frame(vessel, self._token)
        self.current_vessel_status = "vessels" if vessel is not None else "vessels unavailable"
        self.backend.restore_camera(camera)
        self.current_quality = "preview"
        return self.current_quality

    def tick(self, now_ms: int) -> str:
        self.catalog.refresh()
        if self.current_time is None or self._full_token_done == self._token:
            return self.current_quality
        if int(now_ms) - self._last_slider_ms < self.debounce_ms:
            return self.current_quality
        frame = self.catalog.exact_full(self.current_time)
        self._full_token_done = self._token
        if frame is None:
            self.current_quality = "preview (no exact full frame)"
            return self.current_quality
        camera = self.backend.capture_camera()
        token = self._token
        counts = self.backend.load_frame(frame, "full", token)
        if token != self._token:
            return self.current_quality
        vessel = self.catalog.exact_vessel(frame.time_hours)
        vessel_count = self.backend.load_vessel_frame(vessel, token)
        if token != self._token:
            return self.current_quality
        if token == self._token:
            self.current_count = counts.total
            self.current_displayed_count = counts.displayed
            self.current_vessel_count = vessel_count
            self.current_vessel_status = (
                "vessels" if vessel is not None else "vessels unavailable"
            )
            self.backend.restore_camera(camera)
            self.current_quality = "full"
        return self.current_quality
