#!/usr/bin/env pvpython
"""Minimal server-side ParaView + trame viewer for an ATCG3D run directory."""

from __future__ import annotations

import argparse
import asyncio
import math
from pathlib import Path
import sys
import time

try:
    from .controller import Frame, FrameCounts, PreviewFullController, SeriesCatalog
except ImportError:  # Allows direct `pvpython app.py` execution.
    from controller import Frame, FrameCounts, PreviewFullController, SeriesCatalog


CELL_TYPE_ANNOTATIONS = ["1", "r", "2", "K"]
CELL_TYPE_INDEXED_COLORS = [
    0.0, 0.72, 0.0,  # r: green
    0.90, 0.0, 0.0,  # K: red
]
VESSEL_SOLID_COLOR = [0.05, 0.25, 1.0]
SLICE_AXES = {
    "X": (1.0, 0.0, 0.0),
    "Y": (0.0, 1.0, 0.0),
    "Z": (0.0, 0.0, 1.0),
}


def clipping_plane_specs(mode: str, axis: str, position: float,
                         thickness: float, invert: bool = False):
    """Return ParaView Clip plane specifications for whole/cut/slab views."""
    mode = str(mode).lower()
    axis = str(axis).upper()
    position = float(position)
    thickness = float(thickness)
    if mode not in ("whole", "cut", "slab"):
        raise ValueError(f"unsupported slice mode: {mode}")
    if axis not in SLICE_AXES:
        raise ValueError(f"unsupported slice axis: {axis}")
    if not math.isfinite(position) or not math.isfinite(thickness):
        raise ValueError("slice position and thickness must be finite")
    if mode == "whole":
        return []
    normal = SLICE_AXES[axis]
    origin = [normal[index] * position for index in range(3)]
    if mode == "cut":
        return [(origin, list(normal), int(bool(invert)))]
    if thickness <= 0.0:
        raise ValueError("slice thickness must be positive")
    half = thickness * 0.5
    lower = [normal[index] * (position - half) for index in range(3)]
    upper = [normal[index] * (position + half) for index in range(3)]
    return [
        (lower, list(normal), 0),
        (upper, list(normal), 1),
    ]


class ParaViewBackend:
    def __init__(self):
        from paraview import simple

        self.pv = simple
        self.view = simple.GetActiveViewOrCreate("RenderView")
        self.cell_source = None
        self.cell_radius_calculator = None
        self.cell_clip_filters = []
        self.cell_render_source = None
        self.cell_display = None
        self.cell_radius_scale = 1.0
        self.vessel_source = None
        self.vessel_radius_calculator = None
        self.vessel_clip_filters = []
        self.vessel_tube = None
        self.vessel_render_source = None
        self.vessel_display = None
        self.vessel_radius_scale = 1.0
        self.slice_mode = "whole"
        self.slice_axis = "Z"
        self.slice_position = 0.0
        self.slice_thickness = 2.0
        self.slice_invert = False
        self.cancelled_full_tokens: set[int] = set()
        self.cell_count = 0
        self.displayed_cell_count = 0

    @staticmethod
    def _read_total_cell_count(source, displayed_count: int) -> int:
        field = source.GetDataInformation().GetFieldDataInformation()
        info = field.GetArrayInformation("total_cell_count")
        if info is None:
            raise RuntimeError(
                "cell frame has no total_cell_count field metadata; "
                "regenerate it with the current ATCG3D writer"
            )
        value_range = info.GetComponentRange(0)
        if len(value_range) != 2 or value_range[0] != value_range[1]:
            raise RuntimeError("cell frame total_cell_count metadata is invalid")
        total = int(value_range[0])
        if total < displayed_count:
            raise RuntimeError("cell frame total_cell_count is below its point count")
        return total

    def capture_camera(self):
        return (
            list(self.view.CameraPosition),
            list(self.view.CameraFocalPoint),
            list(self.view.CameraViewUp),
            self.view.CameraParallelScale,
        )

    def restore_camera(self, camera) -> None:
        self.view.CameraPosition = camera[0]
        self.view.CameraFocalPoint = camera[1]
        self.view.CameraViewUp = camera[2]
        self.view.CameraParallelScale = camera[3]
        self.pv.Render(self.view)

    def cancel_full(self, token: int) -> None:
        # vtkHDFReader.UpdatePipeline() is a synchronous ParaView call and
        # cannot be interrupted safely.  The token still cancels a read before
        # it starts and prevents a completed stale read from being published.
        self.cancelled_full_tokens.add(int(token))

    def _build_clipped_pipeline(self, input_proxy):
        current = input_proxy
        filters = []
        for origin, normal, invert in clipping_plane_specs(
                self.slice_mode, self.slice_axis, self.slice_position,
                self.slice_thickness, self.slice_invert):
            clip = self.pv.Clip(Input=current)
            clip.ClipType = "Plane"
            clip.ClipType.Origin = origin
            clip.ClipType.Normal = normal
            clip.Invert = invert
            clip.UpdatePipeline()
            filters.append(clip)
            current = clip
        return current, filters

    def _clear_cell_render(self) -> None:
        if self.cell_render_source is not None:
            self.pv.Hide(self.cell_render_source, self.view)
        for proxy in reversed(self.cell_clip_filters):
            self.pv.Delete(proxy)
        self.cell_clip_filters = []
        self.cell_render_source = None
        self.cell_display = None

    def _build_cell_render(self) -> None:
        if self.cell_radius_calculator is None:
            return
        self.cell_render_source, self.cell_clip_filters = \
            self._build_clipped_pipeline(self.cell_radius_calculator)
        self.cell_display = self.pv.Show(self.cell_render_source, self.view)
        self.cell_display.Representation = "Point Gaussian"
        self.cell_display.GaussianRadius = 0.5
        self.pv.ColorBy(self.cell_display, ("POINTS", "cell_type"))
        lookup = self.pv.GetColorTransferFunction("cell_type")
        lookup.InterpretValuesAsCategories = 1
        lookup.Annotations = CELL_TYPE_ANNOTATIONS
        lookup.IndexedColors = CELL_TYPE_INDEXED_COLORS
        self.cell_display.LookupTable = lookup
        self.cell_display.SetScaleArray = ["POINTS", "viewer_display_radius"]
        self.cell_display.ScaleByArray = 1
        self.cell_display.UseScaleFunction = 0
        self.cell_display.SetScalarBarVisibility(self.view, True)

    def data_bounds(self):
        if self.cell_source is None:
            return None
        bounds = tuple(float(value) for value in
                       self.cell_source.GetDataInformation().GetBounds())
        if len(bounds) != 6 or any(not math.isfinite(value) for value in bounds):
            return None
        if any(bounds[2 * axis] > bounds[2 * axis + 1] for axis in range(3)):
            return None
        return bounds

    def set_slice(self, mode: str, axis: str, position: float,
                  thickness: float, invert: bool = False) -> None:
        # Validate before mutating or tearing down the current presentation.
        clipping_plane_specs(mode, axis, position, thickness, invert)
        camera = self.capture_camera()
        self.slice_mode = str(mode).lower()
        self.slice_axis = str(axis).upper()
        self.slice_position = float(position)
        self.slice_thickness = float(thickness)
        self.slice_invert = bool(invert)
        self._clear_cell_render()
        self._build_cell_render()
        self._clear_vessel_render()
        self._build_vessel_render()
        self.restore_camera(camera)

    def load_frame(self, frame: Frame, quality: str, token: int) -> FrameCounts:
        token = int(token)
        if quality == "full" and token in self.cancelled_full_tokens:
            return FrameCounts(self.cell_count, self.displayed_cell_count)

        new_source = self.pv.OpenDataFile(str(frame.path))
        new_source.UpdatePipeline()
        displayed_count = int(new_source.GetDataInformation().GetNumberOfPoints())
        total_count = self._read_total_cell_count(new_source, displayed_count)
        if quality == "full" and token in self.cancelled_full_tokens:
            self.pv.Delete(new_source)
            return FrameCounts(self.cell_count, self.displayed_cell_count)

        new_calculator = self.pv.Calculator(Input=new_source)
        new_calculator.ResultArrayName = "viewer_display_radius"
        new_calculator.Function = (
            f"display_radius*{self.cell_radius_scale:.17g}"
        )
        new_calculator.UpdatePipeline()

        old_source = self.cell_source
        old_calculator = self.cell_radius_calculator
        self._clear_cell_render()
        self.cell_source = new_source
        self.cell_radius_calculator = new_calculator
        self._build_cell_render()
        if old_calculator is not None:
            self.pv.Delete(old_calculator)
        if old_source is not None:
            self.pv.Delete(old_source)
        self.pv.Render(self.view)
        self.cell_count = total_count
        self.displayed_cell_count = displayed_count
        if quality == "full":
            self.cancelled_full_tokens.discard(token)
        return FrameCounts(total_count, displayed_count)

    def set_cell_radius_scale(self, scale: float) -> None:
        scale = float(scale)
        if not scale > 0.0:
            raise ValueError("cell radius scale must be positive")
        self.cell_radius_scale = scale
        if self.cell_radius_calculator is not None:
            self.cell_radius_calculator.Function = (
                f"display_radius*{self.cell_radius_scale:.17g}"
            )
            self.cell_radius_calculator.UpdatePipeline()
            self.pv.Render(self.view)

    def _clear_vessel_render(self) -> None:
        if self.vessel_render_source is not None:
            self.pv.Hide(self.vessel_render_source, self.view)
        for proxy in reversed(self.vessel_clip_filters):
            self.pv.Delete(proxy)
        self.vessel_clip_filters = []
        self.vessel_render_source = None
        self.vessel_display = None
        if self.vessel_tube is not None:
            self.pv.Delete(self.vessel_tube)
            self.vessel_tube = None

    def _clear_vessels(self) -> None:
        self._clear_vessel_render()
        if self.vessel_radius_calculator is not None:
            self.pv.Delete(self.vessel_radius_calculator)
            self.vessel_radius_calculator = None
        if self.vessel_source is not None:
            self.pv.Delete(self.vessel_source)
            self.vessel_source = None

    def _build_vessel_render(self) -> None:
        if self.vessel_source is None:
            return
        tube_input = self.vessel_radius_calculator or self.vessel_source
        radius_array = (
            "viewer_radius_voxels" if self.vessel_radius_calculator is not None
            else "radius_voxels"
        )
        self.vessel_tube = self.pv.Tube(Input=tube_input)
        self.vessel_tube.NumberofSides = 12
        self.vessel_tube.Capping = 1
        self.vessel_tube.Radius = 0.5
        # Radius is a point array on the centerline. ParaView versions differ
        # slightly in Tube property exposure, so fail gracefully to the fixed
        # half-voxel radius while preserving the centerline overlay.
        try:
            self.vessel_tube.VaryRadius = "By Absolute Scalar"
            self.vessel_tube.Scalars = ["POINTS", radius_array]
        except Exception:  # pragma: no cover - ParaView-version dependent
            pass
        self.vessel_tube.UpdatePipeline()
        self.vessel_render_source, self.vessel_clip_filters = \
            self._build_clipped_pipeline(self.vessel_tube)
        self.vessel_display = self.pv.Show(self.vessel_render_source, self.view)
        self.vessel_display.DiffuseColor = VESSEL_SOLID_COLOR
        self.vessel_display.AmbientColor = VESSEL_SOLID_COLOR

    def load_vessel_frame(self, frame: Frame | None, token: int) -> int:
        del token
        self._clear_vessels()
        if frame is None:
            self.pv.Render(self.view)
            return 0
        self.vessel_source = self.pv.OpenDataFile(str(frame.path))
        self.vessel_source.UpdatePipeline()
        try:
            self.vessel_radius_calculator = self.pv.Calculator(Input=self.vessel_source)
            self.vessel_radius_calculator.ResultArrayName = "viewer_radius_voxels"
            self.vessel_radius_calculator.Function = (
                f"radius_voxels*{self.vessel_radius_scale:.17g}"
            )
            self.vessel_radius_calculator.UpdatePipeline()
        except Exception:  # pragma: no cover - ParaView-version dependent
            if self.vessel_radius_calculator is not None:
                self.pv.Delete(self.vessel_radius_calculator)
                self.vessel_radius_calculator = None
        self._build_vessel_render()
        self.pv.Render(self.view)
        return int(self.vessel_source.GetDataInformation().GetNumberOfPoints())

    def set_vessel_radius_scale(self, scale: float) -> None:
        self.vessel_radius_scale = float(scale)
        if self.vessel_radius_calculator is not None:
            self.vessel_radius_calculator.Function = (
                f"radius_voxels*{self.vessel_radius_scale:.17g}"
            )
            self.vessel_radius_calculator.UpdatePipeline()
        if self.vessel_tube is not None:
            if self.vessel_radius_calculator is None:
                self.vessel_tube.Radius = 0.5 * self.vessel_radius_scale
            self.vessel_tube.UpdatePipeline()
            self.pv.Render(self.view)


def build_app(run_directory: Path, debounce_ms: int = 250):
    try:
        from trame.app import get_server
        from trame.ui.vuetify3 import SinglePageLayout
        from trame.widgets import html, vtk as vtk_widgets, vuetify3
    except ImportError as error:
        raise RuntimeError(
            "viewer requires ParaView's Python environment plus trame, "
            "trame-vtk and trame-vuetify"
        ) from error

    catalog = SeriesCatalog(run_directory)
    backend = ParaViewBackend()
    timeline = PreviewFullController(catalog, backend, debounce_ms)
    server = get_server(client_type="vue3")
    state, ctrl = server.state, server.controller
    state.times = catalog.times
    state.time_index = 0
    state.time_hours = state.times[0] if state.times else 0.0
    state.quality = "none"
    state.cell_count = 0
    state.displayed_cell_count = 0
    state.vessel_count = 0
    state.vessel_status = "vessels unavailable"
    state.playing = False
    state.radius_scale = 1.0
    state.vessel_radius_scale = 1.0
    state.slice_mode = "whole"
    state.slice_mode_options = ["whole", "cut", "slab"]
    state.slice_axis = "Z"
    state.slice_axis_options = ["X", "Y", "Z"]
    state.slice_position = 0.0
    state.slice_min = -1.0
    state.slice_max = 1.0
    state.slice_thickness = 2.0
    state.slice_invert = False

    def sync_slice_bounds(recenter=False):
        bounds = backend.data_bounds()
        if bounds is None:
            return
        axis_index = {"X": 0, "Y": 1, "Z": 2}[str(state.slice_axis)]
        minimum = float(bounds[2 * axis_index])
        maximum = float(bounds[2 * axis_index + 1])
        state.slice_min = minimum
        state.slice_max = maximum
        if (recenter or float(state.slice_position) < minimum or
                float(state.slice_position) > maximum):
            state.slice_position = 0.5 * (minimum + maximum)

    def load_index(index):
        catalog.refresh()
        state.times = catalog.times
        if not state.times:
            return
        index = max(0, min(int(index), len(state.times) - 1))
        state.time_index = index
        state.time_hours = state.times[index]
        state.quality = timeline.on_slider(state.time_hours, int(time.monotonic() * 1000))
        state.cell_count = timeline.current_count
        state.displayed_cell_count = timeline.current_displayed_count
        state.vessel_count = timeline.current_vessel_count
        state.vessel_status = timeline.current_vessel_status
        sync_slice_bounds()
        ctrl.view_update()

    @state.change("time_index")
    def _time_changed(time_index, **_):
        load_index(time_index)

    @state.change("radius_scale")
    def _radius_changed(radius_scale, **_):
        backend.set_cell_radius_scale(float(radius_scale))
        ctrl.view_update()

    @state.change("vessel_radius_scale")
    def _vessel_radius_changed(vessel_radius_scale, **_):
        backend.set_vessel_radius_scale(float(vessel_radius_scale))
        ctrl.view_update()

    def apply_slice():
        backend.set_slice(
            str(state.slice_mode), str(state.slice_axis),
            float(state.slice_position), float(state.slice_thickness),
            bool(state.slice_invert))
        ctrl.view_update()

    @state.change("slice_mode", "slice_position", "slice_thickness", "slice_invert")
    def _slice_changed(**_):
        apply_slice()

    @state.change("slice_axis")
    def _slice_axis_changed(**_):
        sync_slice_bounds(recenter=True)
        apply_slice()

    def toggle_play():
        state.playing = not state.playing

    ctrl.toggle_play = toggle_play

    # trame >= 3.13 forwards the current server state as keyword arguments to
    # on_server_ready tasks.  The polling loop does not need that state, but it
    # must accept it to remain compatible with the callback contract.
    async def update_loop(**_):
        while True:
            await asyncio.sleep(0.05)
            changed = catalog.refresh()
            if changed:
                state.times = catalog.times
            previous_quality = state.quality
            state.quality = timeline.tick(int(time.monotonic() * 1000))
            if state.quality != previous_quality:
                sync_slice_bounds()
            state.cell_count = timeline.current_count
            state.displayed_cell_count = timeline.current_displayed_count
            state.vessel_count = timeline.current_vessel_count
            state.vessel_status = timeline.current_vessel_status
            if state.playing and state.times:
                load_index((int(state.time_index) + 1) % len(state.times))
                await asyncio.sleep(0.05)
            ctrl.view_update()

    ctrl.on_server_ready.add_task(update_loop)

    with SinglePageLayout(server) as layout:
        layout.title.set_text("ATCG3D Viewer")
        with layout.toolbar:
            vuetify3.VBtn("{{ playing ? 'Pause' : 'Play' }}", click=ctrl.toggle_play)
            vuetify3.VSlider(
                v_model=("time_index", 0), min=0,
                max=("Math.max(times.length - 1, 0)",), step=1, hide_details=True,
                style="min-width: 360px",
            )
            html.Span("t={{ Number(time_hours).toFixed(3) }} h")
            html.Span("cells={{ cell_count.toLocaleString() }}")
            html.Span("displayed={{ displayed_cell_count.toLocaleString() }}")
            html.Span("vessel nodes={{ vessel_count.toLocaleString() }}")
            html.Span("{{ quality }}")
            html.Span("{{ vessel_status }}")
            vuetify3.VSlider(
                v_model=("radius_scale", 1.0), min=0.1, max=4.0, step=0.1,
                label="Point size", hide_details=True, style="max-width: 220px",
            )
            vuetify3.VSlider(
                v_model=("vessel_radius_scale", 1.0), min=0.1, max=4.0, step=0.1,
                label="Vessel radius", hide_details=True, style="max-width: 220px",
            )
            vuetify3.VSelect(
                v_model=("slice_mode", "whole"), items=("slice_mode_options",),
                label="View", hide_details=True, density="compact",
                style="max-width: 120px",
            )
            vuetify3.VSelect(
                v_model=("slice_axis", "Z"), items=("slice_axis_options",),
                label="Slice axis", hide_details=True, density="compact",
                disabled=("slice_mode === 'whole'",), style="max-width: 110px",
            )
            vuetify3.VSlider(
                v_model=("slice_position", 0.0), min=("slice_min",),
                max=("slice_max",), step=0.5, label="Slice position",
                hide_details=True, disabled=("slice_mode === 'whole'",),
                style="min-width: 220px",
            )
            vuetify3.VSlider(
                v_model=("slice_thickness", 2.0), min=0.5,
                max=("Math.max(slice_max - slice_min, 1)",), step=0.5,
                label="Slab thickness", hide_details=True,
                disabled=("slice_mode !== 'slab'",), style="max-width: 200px",
            )
            vuetify3.VSwitch(
                v_model=("slice_invert", False), label="Invert cut",
                hide_details=True, density="compact",
                disabled=("slice_mode !== 'cut'",),
            )
        with layout.content:
            view = vtk_widgets.VtkRemoteView(backend.view, interactive_ratio=1)
            ctrl.view_update = view.update

    if state.times:
        load_index(0)
        backend.pv.ResetCamera(backend.view)
    return server


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_directory", type=Path)
    parser.add_argument("--debounce-ms", type=int, default=250)
    args, server_args = parser.parse_known_args(argv)
    try:
        server = build_app(args.run_directory, args.debounce_ms)
    except (RuntimeError, ValueError, OSError) as error:
        print(f"atcg3d-viewer: {error}", file=sys.stderr)
        return 1
    server.start(argv=server_args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
