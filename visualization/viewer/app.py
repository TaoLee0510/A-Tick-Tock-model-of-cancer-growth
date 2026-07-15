#!/usr/bin/env pvpython
"""Minimal server-side ParaView + trame viewer for an ATCG3D run directory."""

from __future__ import annotations

import argparse
import asyncio
from pathlib import Path
import sys
import time

try:
    from .controller import Frame, FrameCounts, PreviewFullController, SeriesCatalog
except ImportError:  # Allows direct `pvpython app.py` execution.
    from controller import Frame, FrameCounts, PreviewFullController, SeriesCatalog


class ParaViewBackend:
    def __init__(self):
        from paraview import simple

        self.pv = simple
        self.view = simple.GetActiveViewOrCreate("RenderView")
        self.cell_source = None
        self.cell_radius_calculator = None
        self.cell_display = None
        self.cell_radius_scale = 1.0
        self.vessel_source = None
        self.vessel_radius_calculator = None
        self.vessel_tube = None
        self.vessel_display = None
        self.vessel_color = "perfused"
        self.vessel_radius_scale = 1.0
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
        if old_calculator is not None:
            self.pv.Hide(old_calculator, self.view)
        new_display = self.pv.Show(new_calculator, self.view)
        new_display.Representation = "Point Gaussian"
        new_display.GaussianRadius = 0.5
        self.pv.ColorBy(new_display, ("POINTS", "cell_type"))
        new_display.SetScaleArray = ["POINTS", "viewer_display_radius"]
        new_display.ScaleByArray = 1
        new_display.UseScaleFunction = 0
        self.pv.Render(self.view)
        if old_calculator is not None:
            self.pv.Delete(old_calculator)
        if old_source is not None:
            self.pv.Delete(old_source)
        self.cell_source = new_source
        self.cell_radius_calculator = new_calculator
        self.cell_display = new_display
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

    def _clear_vessels(self) -> None:
        if self.vessel_tube is not None:
            self.pv.Hide(self.vessel_tube, self.view)
            self.pv.Delete(self.vessel_tube)
            self.vessel_tube = None
            self.vessel_display = None
        if self.vessel_radius_calculator is not None:
            self.pv.Delete(self.vessel_radius_calculator)
            self.vessel_radius_calculator = None
        if self.vessel_source is not None:
            self.pv.Delete(self.vessel_source)
            self.vessel_source = None

    def load_vessel_frame(self, frame: Frame | None, token: int) -> int:
        del token
        self._clear_vessels()
        if frame is None:
            self.pv.Render(self.view)
            return 0
        self.vessel_source = self.pv.OpenDataFile(str(frame.path))
        self.vessel_source.UpdatePipeline()
        tube_input = self.vessel_source
        radius_array = "radius_voxels"
        try:
            self.vessel_radius_calculator = self.pv.Calculator(Input=self.vessel_source)
            self.vessel_radius_calculator.ResultArrayName = "viewer_radius_voxels"
            self.vessel_radius_calculator.Function = (
                f"radius_voxels*{self.vessel_radius_scale:.17g}"
            )
            self.vessel_radius_calculator.UpdatePipeline()
            tube_input = self.vessel_radius_calculator
            radius_array = "viewer_radius_voxels"
        except Exception:  # pragma: no cover - ParaView-version dependent
            if self.vessel_radius_calculator is not None:
                self.pv.Delete(self.vessel_radius_calculator)
                self.vessel_radius_calculator = None
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
        self.vessel_display = self.pv.Show(self.vessel_tube, self.view)
        self.set_vessel_color(self.vessel_color)
        self.pv.Render(self.view)
        return int(self.vessel_source.GetDataInformation().GetNumberOfPoints())

    def set_vessel_color(self, array_name: str) -> None:
        if array_name not in ("perfused", "branch_role"):
            raise ValueError(f"unsupported vessel color array: {array_name}")
        self.vessel_color = array_name
        if self.vessel_display is not None:
            self.pv.ColorBy(self.vessel_display, ("POINTS", array_name))
            self.vessel_display.RescaleTransferFunctionToDataRange(True, False)

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
    state.vessel_color = "perfused"
    state.vessel_color_options = ["perfused", "branch_role"]

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

    @state.change("vessel_color")
    def _vessel_color_changed(vessel_color, **_):
        backend.set_vessel_color(str(vessel_color))
        backend.pv.Render(backend.view)
        ctrl.view_update()

    def toggle_play():
        state.playing = not state.playing

    ctrl.toggle_play = toggle_play

    async def update_loop():
        while True:
            await asyncio.sleep(0.05)
            changed = catalog.refresh()
            if changed:
                state.times = catalog.times
            state.quality = timeline.tick(int(time.monotonic() * 1000))
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
                v_model=("vessel_color", "perfused"),
                items=("vessel_color_options",), label="Vessel color",
                hide_details=True, density="compact", style="max-width: 180px",
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
