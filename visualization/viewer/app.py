#!/usr/bin/env pvpython
"""Minimal server-side ParaView + trame viewer for an ATCG3D run directory."""

from __future__ import annotations

import argparse
import asyncio
from pathlib import Path
import sys
import time

try:
    from .controller import Frame, PreviewFullController, SeriesCatalog
except ImportError:  # Allows direct `pvpython app.py` execution.
    from controller import Frame, PreviewFullController, SeriesCatalog


class ParaViewBackend:
    def __init__(self):
        from paraview import simple

        self.pv = simple
        self.view = simple.GetActiveViewOrCreate("RenderView")
        self.source = None
        self.display = None

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
        # Full reads are never queued: the dependency-free controller retains
        # only the newest token until idle. This hook supports future async I/O.
        del token

    def load_frame(self, frame: Frame, quality: str, token: int) -> int:
        del quality, token
        if self.source is not None:
            self.pv.Delete(self.source)
        self.source = self.pv.OpenDataFile(str(frame.path))
        self.source.UpdatePipeline()
        self.display = self.pv.Show(self.source, self.view)
        self.display.Representation = "Point Gaussian"
        self.display.GaussianRadius = 0.5
        self.pv.ColorBy(self.display, ("POINTS", "cell_type"))
        self.display.SetScaleArray = ["POINTS", "display_radius"]
        self.display.ScaleTransferFunction = "PiecewiseFunction"
        self.pv.Render(self.view)
        return int(self.source.GetDataInformation().GetNumberOfPoints())


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
    state.playing = False
    state.radius_scale = 1.0

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
        ctrl.view_update()

    @state.change("time_index")
    def _time_changed(time_index, **_):
        load_index(time_index)

    @state.change("radius_scale")
    def _radius_changed(radius_scale, **_):
        if backend.display is not None:
            backend.display.GaussianRadius = 0.5 * float(radius_scale)
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
            html.Span("{{ quality }}")
            vuetify3.VSlider(
                v_model=("radius_scale", 1.0), min=0.1, max=4.0, step=0.1,
                label="Point size", hide_details=True, style="max-width: 220px",
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
