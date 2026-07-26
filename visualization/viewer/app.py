#!/usr/bin/env pvpython
"""Minimal server-side ParaView + trame viewer for an ATCG3D run directory."""

from __future__ import annotations

import argparse
import asyncio
import json
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
VESSEL_INFLUENCE_COLOR = [0.2, 0.55, 1.0]
SLICE_AXES = {
    "X": (1.0, 0.0, 0.0),
    "Y": (0.0, 1.0, 0.0),
    "Z": (0.0, 0.0, 1.0),
}
VIEW_SLICE_ORIENTATIONS = {
    "VIEW_AXIAL": "forward",
    "VIEW_SAGITTAL": "right",
    "VIEW_CORONAL": "up",
}

VIEWER_TRANSLATIONS = {
    "en": {
        "title": "ATCG3D Viewer",
        "play": "Play",
        "pause": "Pause",
        "time": "time",
        "cells": "cells",
        "displayed": "displayed",
        "vessel_nodes": "vessel nodes",
        "camera": "Camera",
        "primary_drag": "Primary drag",
        "rotate": "Rotate",
        "pan": "Pan",
        "camera_help": (
            "Rotate: left drag. Pan: left drag. Mouse shortcuts: "
            "middle/Shift+left pan; wheel/right drag zoom."
        ),
        "reset_camera": "Reset camera",
        "objects": "Objects",
        "r_cells": "r cells",
        "k_cells": "K cells",
        "vessels": "Vessels",
        "influence": "Influence",
        "voxels": "vox",
        "cell_point_size": "Cell point size",
        "vessel_radius": "Vessel radius",
        "slice": "Slice",
        "mode": "Mode",
        "orientation": "Orientation",
        "whole": "Whole",
        "cut": "One-sided cut",
        "slab": "Slab",
        "world_x": "World X",
        "world_y": "World Y",
        "world_z": "World Z",
        "view_axial": "View axial",
        "view_sagittal": "View sagittal",
        "view_coronal": "View coronal",
        "recapture_view": "Recapture current view",
        "slice_offset": "Offset from captured center",
        "slab_thickness": "Slab thickness",
        "invert_cut": "Invert cut",
        "quality.none": "none",
        "quality.preview unavailable": "preview unavailable",
        "quality.preview": "preview",
        "quality.preview (no exact full frame)": "preview (no exact full frame)",
        "quality.full": "full",
        "vessel.vessels": "vessels",
        "vessel.vessels unavailable": "vessels unavailable",
    },
    "zh-CN": {
        "title": "ATCG3D 三维查看器",
        "play": "播放",
        "pause": "暂停",
        "time": "时间",
        "cells": "细胞",
        "displayed": "已显示",
        "vessel_nodes": "血管节点",
        "camera": "相机",
        "primary_drag": "主拖动操作",
        "rotate": "旋转",
        "pan": "平移",
        "camera_help": (
            "旋转：左键拖动。平移：平移模式下左键拖动。快捷操作："
            "中键或 Shift+左键平移；滚轮或右键拖动缩放。"
        ),
        "reset_camera": "重置相机",
        "objects": "显示对象",
        "r_cells": "r 细胞",
        "k_cells": "K 细胞",
        "vessels": "血管",
        "influence": "血管影响范围",
        "voxels": "体素",
        "cell_point_size": "细胞点大小",
        "vessel_radius": "血管半径",
        "slice": "切片",
        "mode": "模式",
        "orientation": "方向",
        "whole": "完整实体",
        "cut": "单侧切面",
        "slab": "薄层",
        "world_x": "世界坐标 X",
        "world_y": "世界坐标 Y",
        "world_z": "世界坐标 Z",
        "view_axial": "视角轴向面",
        "view_sagittal": "视角矢状面",
        "view_coronal": "视角冠状面",
        "recapture_view": "按当前视角重新捕获",
        "slice_offset": "相对捕获中心的偏移",
        "slab_thickness": "薄层厚度",
        "invert_cut": "反转切面",
        "quality.none": "无数据",
        "quality.preview unavailable": "预览不可用",
        "quality.preview": "预览",
        "quality.preview (no exact full frame)": "预览（无精确完整帧）",
        "quality.full": "完整帧",
        "vessel.vessels": "血管",
        "vessel.vessels unavailable": "血管不可用",
    },
}


def normalize_viewer_language(language: str) -> str:
    """Return a supported UI language; fresh sessions default to English."""
    normalized = str(language or "en")
    if normalized not in VIEWER_TRANSLATIONS:
        raise ValueError(f"unsupported viewer language: {normalized}")
    return normalized


def viewer_text(language: str, key: str) -> str:
    """Look up localized viewer text with an English fallback."""
    language = normalize_viewer_language(language)
    return VIEWER_TRANSLATIONS[language].get(
        key, VIEWER_TRANSLATIONS["en"].get(key, key))


def localized_viewer_status(language: str, category: str, value: str) -> str:
    """Translate known runtime labels while preserving unknown diagnostics."""
    key = f"{category}.{value}"
    translated = viewer_text(language, key)
    return str(value) if translated == key else translated


def normalized_plane_normal(normal):
    """Return a finite unit normal suitable for a persistent view plane."""
    values = tuple(float(value) for value in normal)
    if len(values) != 3 or any(not math.isfinite(value) for value in values):
        raise ValueError("slice normal must contain three finite values")
    length = math.sqrt(sum(value * value for value in values))
    if not length > 0.0:
        raise ValueError("slice normal must be nonzero")
    return tuple(value / length for value in values)


def _cross(left, right):
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def camera_basis(camera_position, focal_point, view_up):
    """Return an orthonormal forward/right/up basis for an arbitrary camera."""
    position = tuple(float(value) for value in camera_position)
    focal = tuple(float(value) for value in focal_point)
    forward = normalized_plane_normal(
        [focal[index] - position[index] for index in range(3)])
    raw_up = normalized_plane_normal(view_up)
    raw_right = _cross(forward, raw_up)
    if sum(value * value for value in raw_right) <= 1.0e-24:
        # A pole-on camera can make its nominal view-up parallel to forward.
        # Choose the least-aligned world axis so every finite camera pose still
        # has a stable screen-relative sagittal plane.
        fallback_up = min(
            ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
            key=lambda axis: abs(sum(a * b for a, b in zip(forward, axis))),
        )
        raw_right = _cross(forward, fallback_up)
    right = normalized_plane_normal(raw_right)
    up = normalized_plane_normal(_cross(right, forward))
    return forward, right, up


def camera_oriented_plane(camera_position, focal_point, view_up,
                          orientation: str):
    """Capture an axial, sagittal, or coronal plane in the current view basis."""
    orientation = str(orientation).upper()
    basis_name = VIEW_SLICE_ORIENTATIONS.get(orientation)
    if basis_name is None:
        raise ValueError(f"unsupported view slice orientation: {orientation}")
    forward, right, up = camera_basis(camera_position, focal_point, view_up)
    normal = {"forward": forward, "right": right, "up": up}[basis_name]
    focal = tuple(float(value) for value in focal_point)
    offset = sum(normal[index] * focal[index] for index in range(3))
    return normal, offset


def oriented_clipping_plane_specs(mode: str, normal, position: float,
                                  thickness: float, invert: bool = False):
    """Build cut/slab planes at signed distance ``position`` along ``normal``."""
    mode = str(mode).lower()
    position = float(position)
    thickness = float(thickness)
    if mode not in ("whole", "cut", "slab"):
        raise ValueError(f"unsupported slice mode: {mode}")
    normal = normalized_plane_normal(normal)
    if not math.isfinite(position) or not math.isfinite(thickness):
        raise ValueError("slice position and thickness must be finite")
    if mode == "whole":
        return []
    origin = [normal[index] * position for index in range(3)]
    if mode == "cut":
        return [(origin, list(normal), int(bool(invert)))]
    if thickness <= 0.0:
        raise ValueError("slice thickness must be positive")
    half = thickness * 0.5
    lower = [normal[index] * (position - half) for index in range(3)]
    upper = [normal[index] * (position + half) for index in range(3)]
    return [(lower, list(normal), 0), (upper, list(normal), 1)]


def clipping_plane_specs(mode: str, axis: str, position: float,
                         thickness: float, invert: bool = False):
    """Return ParaView Clip plane specifications for whole/cut/slab views."""
    axis = str(axis).upper()
    if axis not in SLICE_AXES:
        raise ValueError(f"unsupported slice axis: {axis}")
    return oriented_clipping_plane_specs(
        mode, SLICE_AXES[axis], position, thickness, invert)


def projected_bounds(bounds, normal):
    """Project all eight corners of an axis-aligned box onto a plane normal."""
    if len(bounds) != 6:
        raise ValueError("bounds must contain six values")
    values = tuple(float(value) for value in bounds)
    if any(not math.isfinite(value) for value in values):
        raise ValueError("bounds must be finite")
    normal = normalized_plane_normal(normal)
    projections = [
        x * normal[0] + y * normal[1] + z * normal[2]
        for x in values[0:2]
        for y in values[2:4]
        for z in values[4:6]
    ]
    return min(projections), max(projections)


def camera_aligned_plane(camera_position, focal_point):
    """Capture a plane normal and offset from the current view without tracking it."""
    position = tuple(float(value) for value in camera_position)
    focal = tuple(float(value) for value in focal_point)
    normal = normalized_plane_normal(
        [focal[index] - position[index] for index in range(3)])
    return normal, sum(normal[index] * focal[index] for index in range(3))


def viewer_config(run_directory: Path):
    """Read only visualization-relevant values from atomic run metadata."""
    path = Path(run_directory) / "run.json"
    if not path.exists():
        return {
            "influence_cutoff_radius_voxels": 0.0,
            "display_radii": (1.0, 0.5, 0.25),
        }
    document = json.loads(path.read_text(encoding="utf-8"))
    effective = document.get("effective_config", {})
    cutoff = float(effective.get(
        "angiogenesis_influence_cutoff_radius_voxels", 0.0))
    if not math.isfinite(cutoff) or cutoff < 0.0:
        raise ValueError("run.json contains an invalid vascular influence cutoff")
    radii = tuple(float(effective.get(name, default)) for name, default in (
        ("display_radius_large", 1.0),
        ("display_radius_small", 0.5),
        ("display_radius_ultrasmall", 0.25),
    ))
    if any(not math.isfinite(value) or value <= 0.0 for value in radii):
        raise ValueError("run.json contains invalid cell display radii")
    return {
        "influence_cutoff_radius_voxels": cutoff,
        "display_radii": radii,
    }


class ParaViewBackend:
    def __init__(self, influence_cutoff_radius_voxels: float = 0.0,
                 run_directory: Path | None = None,
                 display_radii=(1.0, 0.5, 0.25)):
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
        self.vessel_halo_calculator = None
        self.vessel_halo_tube = None
        self.vessel_halo_clip_filters = []
        self.vessel_halo_render_source = None
        self.vessel_halo_display = None
        self.vessel_radius_scale = 1.0
        self.influence_cutoff_radius_voxels = float(influence_cutoff_radius_voxels)
        self.checkpoint_materializer = None
        if run_directory is not None:
            try:
                from .checkpoint_materializer import CheckpointFrameMaterializer
            except ImportError:
                from checkpoint_materializer import CheckpointFrameMaterializer
            self.checkpoint_materializer = CheckpointFrameMaterializer(
                run_directory, display_radii)
        self.show_r = True
        self.show_k = True
        self.show_vessels = True
        self.show_influence = False
        self.slice_mode = "whole"
        self.slice_axis = "Z"
        self.slice_normal = SLICE_AXES["Z"]
        self.slice_anchor_position = 0.0
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

    def set_navigation_mode(self, mode: str) -> None:
        """Switch the primary drag between 3D orbit and 2D camera pan."""
        mode = str(mode).lower()
        interaction_mode = {"rotate": "3D", "pan": "2D"}.get(mode)
        if interaction_mode is None:
            raise ValueError(f"unsupported navigation mode: {mode}")
        camera = self.capture_camera()
        self.view.InteractionMode = interaction_mode
        self.restore_camera(camera)

    def reset_camera(self) -> None:
        self.pv.ResetCamera(self.view)
        self.pv.Render(self.view)

    def cancel_full(self, token: int) -> None:
        # vtkHDFReader.UpdatePipeline() is a synchronous ParaView call and
        # cannot be interrupted safely.  The token still cancels a read before
        # it starts and prevents a completed stale read from being published.
        self.cancelled_full_tokens.add(int(token))

    def _build_clipped_pipeline(self, input_proxy):
        current = input_proxy
        filters = []
        for origin, normal, invert in oriented_clipping_plane_specs(
                self.slice_mode, self.slice_normal,
                self.slice_anchor_position + self.slice_position,
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
        lookup.IndexedOpacities = [1.0 if self.show_r else 0.0,
                                   1.0 if self.show_k else 0.0]
        lookup.EnableOpacityMapping = 1
        self.cell_display.LookupTable = lookup
        self.cell_display.SetScaleArray = ["POINTS", "viewer_display_radius"]
        self.cell_display.ScaleByArray = 1
        self.cell_display.UseScaleFunction = 0
        self.cell_display.SetScalarBarVisibility(self.view, True)
        self.cell_display.Visibility = int(self.show_r or self.show_k)

    @staticmethod
    def _source_bounds(source):
        if source is None:
            return None
        bounds = tuple(float(value) for value in
                       source.GetDataInformation().GetBounds())
        if len(bounds) != 6 or any(not math.isfinite(value) for value in bounds):
            return None
        if any(bounds[2 * axis] > bounds[2 * axis + 1] for axis in range(3)):
            return None
        return bounds

    def _vessel_bounds_expansion(self) -> float:
        radius = 0.5 * self.vessel_radius_scale
        if self.vessel_source is not None:
            try:
                point_data = self.vessel_source.GetDataInformation().GetPointDataInformation()
                radius_info = point_data.GetArrayInformation("radius_voxels")
                if radius_info is not None:
                    radius = max(radius, float(radius_info.GetComponentRange(0)[1]) *
                                 self.vessel_radius_scale)
            except Exception:  # pragma: no cover - ParaView-version dependent
                pass
        return max(0.0, radius + self.influence_cutoff_radius_voxels)

    def data_bounds(self):
        candidates = []
        cell_bounds = self._source_bounds(self.cell_source)
        if cell_bounds is not None:
            candidates.append(cell_bounds)
        vessel_bounds = self._source_bounds(self.vessel_source)
        if vessel_bounds is not None:
            expansion = self._vessel_bounds_expansion()
            candidates.append(tuple(
                value + (-expansion if index % 2 == 0 else expansion)
                for index, value in enumerate(vessel_bounds)
            ))
        if not candidates:
            return None
        return tuple(
            min(bounds[index] for bounds in candidates)
            if index % 2 == 0 else
            max(bounds[index] for bounds in candidates)
            for index in range(6)
        )

    def slice_projection_bounds(self):
        bounds = self.data_bounds()
        return None if bounds is None else projected_bounds(bounds, self.slice_normal)

    def set_slice_orientation(self, axis: str):
        axis = str(axis).upper()
        if axis in SLICE_AXES:
            normal = SLICE_AXES[axis]
            projected = self.data_bounds()
            if projected is None:
                position = 0.0
            else:
                minimum, maximum = projected_bounds(projected, normal)
                position = 0.5 * (minimum + maximum)
        elif axis in VIEW_SLICE_ORIENTATIONS:
            normal, position = camera_oriented_plane(
                self.view.CameraPosition, self.view.CameraFocalPoint,
                self.view.CameraViewUp, axis)
        else:
            raise ValueError(f"unsupported slice axis: {axis}")
        self.slice_axis = axis
        self.slice_normal = normal
        self.slice_anchor_position = float(position)
        self.slice_position = 0.0
        return normal, position

    def capture_view_slice(self, orientation: str):
        orientation = str(orientation).upper()
        if orientation not in VIEW_SLICE_ORIENTATIONS:
            raise ValueError("capture requires a view-relative slice orientation")
        return self.set_slice_orientation(orientation)

    def set_slice(self, mode: str, axis: str, position: float,
                  thickness: float, invert: bool = False) -> None:
        # Validate before mutating or tearing down the current presentation.
        axis = str(axis).upper()
        normal = SLICE_AXES.get(
            axis, self.slice_normal if axis in VIEW_SLICE_ORIENTATIONS else None)
        if normal is None:
            raise ValueError(f"unsupported slice axis: {axis}")
        if axis != self.slice_axis:
            normal, _ = self.set_slice_orientation(axis)
        absolute_position = self.slice_anchor_position + float(position)
        oriented_clipping_plane_specs(
            mode, normal, absolute_position, thickness, invert)
        camera = self.capture_camera()
        self.slice_mode = str(mode).lower()
        self.slice_axis = axis
        if self.slice_axis in SLICE_AXES:
            self.slice_normal = SLICE_AXES[self.slice_axis]
        elif self.slice_axis not in VIEW_SLICE_ORIENTATIONS:
            raise ValueError(f"unsupported slice axis: {axis}")
        self.slice_position = float(position)
        self.slice_thickness = float(thickness)
        self.slice_invert = bool(invert)
        self._clear_cell_render()
        self._build_cell_render()
        self._clear_vessel_render()
        self._build_vessel_render()
        self.restore_camera(camera)

    def set_visibility(self, show_r: bool, show_k: bool,
                       show_vessels: bool, show_influence: bool) -> None:
        self.show_r = bool(show_r)
        self.show_k = bool(show_k)
        self.show_vessels = bool(show_vessels)
        self.show_influence = bool(show_influence)
        if self.cell_display is not None:
            lookup = self.cell_display.LookupTable
            lookup.IndexedOpacities = [1.0 if self.show_r else 0.0,
                                       1.0 if self.show_k else 0.0]
            lookup.EnableOpacityMapping = 1
            self.cell_display.Visibility = int(self.show_r or self.show_k)
        if self.vessel_display is not None:
            self.vessel_display.Visibility = int(self.show_vessels)
        if self.vessel_halo_display is not None:
            self.vessel_halo_display.Visibility = int(self.show_influence)
        self.pv.Render(self.view)

    def load_frame(self, frame: Frame, quality: str, token: int) -> FrameCounts:
        token = int(token)
        if quality == "full" and token in self.cancelled_full_tokens:
            return FrameCounts(self.cell_count, self.displayed_cell_count)

        if frame.storage == "checkpoint":
            if self.checkpoint_materializer is None:
                raise RuntimeError("checkpoint materialization is not configured")
            polydata = self.checkpoint_materializer.materialize(frame.path)
            new_source = self.pv.TrivialProducer()
            new_source.GetClientSideObject().SetOutput(polydata)
        else:
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
        if self.vessel_halo_render_source is not None:
            self.pv.Hide(self.vessel_halo_render_source, self.view)
        for proxy in reversed(self.vessel_halo_clip_filters):
            self.pv.Delete(proxy)
        self.vessel_halo_clip_filters = []
        self.vessel_halo_render_source = None
        self.vessel_halo_display = None
        if self.vessel_halo_tube is not None:
            self.pv.Delete(self.vessel_halo_tube)
            self.vessel_halo_tube = None
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
        if self.vessel_halo_calculator is not None:
            self.pv.Delete(self.vessel_halo_calculator)
            self.vessel_halo_calculator = None
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
        self.vessel_display.Visibility = int(self.show_vessels)

        if (self.vessel_halo_calculator is not None and
                self.influence_cutoff_radius_voxels > 0.0):
            self.vessel_halo_tube = self.pv.Tube(Input=self.vessel_halo_calculator)
            self.vessel_halo_tube.NumberofSides = 16
            self.vessel_halo_tube.Capping = 1
            self.vessel_halo_tube.Radius = 0.5
            try:
                self.vessel_halo_tube.VaryRadius = "By Absolute Scalar"
                self.vessel_halo_tube.Scalars = [
                    "POINTS", "viewer_influence_radius_voxels"]
            except Exception:  # pragma: no cover - ParaView-version dependent
                self.vessel_halo_tube.Radius = (
                    0.5 + self.influence_cutoff_radius_voxels)
            self.vessel_halo_tube.UpdatePipeline()
            (self.vessel_halo_render_source,
             self.vessel_halo_clip_filters) = self._build_clipped_pipeline(
                 self.vessel_halo_tube)
            self.vessel_halo_display = self.pv.Show(
                self.vessel_halo_render_source, self.view)
            self.vessel_halo_display.DiffuseColor = VESSEL_INFLUENCE_COLOR
            self.vessel_halo_display.AmbientColor = VESSEL_INFLUENCE_COLOR
            self.vessel_halo_display.Opacity = 0.12
            self.vessel_halo_display.Visibility = int(self.show_influence)

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
        if self.influence_cutoff_radius_voxels > 0.0:
            try:
                self.vessel_halo_calculator = self.pv.Calculator(
                    Input=self.vessel_source)
                self.vessel_halo_calculator.ResultArrayName = (
                    "viewer_influence_radius_voxels")
                self.vessel_halo_calculator.Function = (
                    "radius_voxels+"
                    f"{self.influence_cutoff_radius_voxels:.17g}")
                self.vessel_halo_calculator.UpdatePipeline()
            except Exception:  # pragma: no cover - ParaView-version dependent
                if self.vessel_halo_calculator is not None:
                    self.pv.Delete(self.vessel_halo_calculator)
                    self.vessel_halo_calculator = None
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


def build_app(run_directory: Path, debounce_ms: int = 250,
              language: str = "en"):
    language = normalize_viewer_language(language)
    tr = lambda key: viewer_text(language, key)
    try:
        from trame.app import get_server
        from trame.ui.vuetify3 import SinglePageWithDrawerLayout
        from trame.widgets import html, vtk as vtk_widgets, vuetify3
    except ImportError as error:
        raise RuntimeError(
            "viewer requires ParaView's Python environment plus trame, "
            "trame-vtk and trame-vuetify"
        ) from error

    catalog = SeriesCatalog(run_directory)
    display_config = viewer_config(run_directory)
    backend = ParaViewBackend(
        display_config["influence_cutoff_radius_voxels"], run_directory,
        display_config["display_radii"])
    timeline = PreviewFullController(catalog, backend, debounce_ms)
    server = get_server(client_type="vue3")
    state, ctrl = server.state, server.controller
    state.times = catalog.times
    state.time_index = 0
    state.time_hours = state.times[0] if state.times else 0.0
    state.quality = "none"
    state.quality_label = localized_viewer_status(
        language, "quality", state.quality)
    state.cell_count = 0
    state.displayed_cell_count = 0
    state.vessel_count = 0
    state.vessel_status = "vessels unavailable"
    state.vessel_status_label = localized_viewer_status(
        language, "vessel", state.vessel_status)
    state.playing = False
    state.play_label = tr("play")
    state.pause_label = tr("pause")
    state.navigation_mode = "rotate"
    state.navigation_mode_options = [
        {"title": tr("rotate"), "value": "rotate"},
        {"title": tr("pan"), "value": "pan"},
    ]
    # A physical one-voxel radius is sub-pixel at whole-tumour camera scales.
    # Start with a presentation-only multiplier that keeps the sampled tumour
    # visible; users can still reduce or enlarge it without changing the model.
    state.radius_scale = 4.0
    backend.set_cell_radius_scale(state.radius_scale)
    state.vessel_radius_scale = 1.0
    state.show_r = True
    state.show_k = True
    state.show_vessels = True
    state.show_influence = False
    state.influence_cutoff_radius_voxels = (
        display_config["influence_cutoff_radius_voxels"])
    state.slice_mode = "whole"
    state.slice_mode_options = [
        {"title": tr("whole"), "value": "whole"},
        {"title": tr("cut"), "value": "cut"},
        {"title": tr("slab"), "value": "slab"},
    ]
    state.slice_axis = "Z"
    state.slice_axis_options = [
        {"title": tr("world_x"), "value": "X"},
        {"title": tr("world_y"), "value": "Y"},
        {"title": tr("world_z"), "value": "Z"},
        {"title": tr("view_axial"), "value": "VIEW_AXIAL"},
        {"title": tr("view_sagittal"), "value": "VIEW_SAGITTAL"},
        {"title": tr("view_coronal"), "value": "VIEW_CORONAL"},
    ]
    state.slice_position = 0.0
    state.slice_min = -1.0
    state.slice_max = 1.0
    state.slice_thickness = 2.0
    state.slice_invert = False

    def sync_slice_bounds(recenter=False):
        bounds = backend.data_bounds()
        if bounds is None:
            return
        projected = backend.slice_projection_bounds()
        if projected is None:
            return
        minimum, maximum = projected
        minimum -= backend.slice_anchor_position
        maximum -= backend.slice_anchor_position
        state.slice_min = minimum
        state.slice_max = maximum
        if (recenter or float(state.slice_position) < minimum or
                float(state.slice_position) > maximum):
            state.slice_position = 0.0

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
        state.quality_label = localized_viewer_status(
            language, "quality", state.quality)
        state.vessel_status_label = localized_viewer_status(
            language, "vessel", state.vessel_status)
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

    @state.change("navigation_mode")
    def _navigation_mode_changed(navigation_mode, **_):
        backend.set_navigation_mode(str(navigation_mode))
        ctrl.view_update()

    def reset_camera():
        backend.reset_camera()
        ctrl.view_update()

    ctrl.reset_camera = reset_camera

    @state.change("show_r", "show_k", "show_vessels", "show_influence")
    def _visibility_changed(**_):
        backend.set_visibility(
            bool(state.show_r), bool(state.show_k),
            bool(state.show_vessels), bool(state.show_influence))
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
        backend.set_slice_orientation(str(state.slice_axis))
        sync_slice_bounds(recenter=True)
        apply_slice()

    def capture_view_slice():
        axis = str(state.slice_axis).upper()
        if axis not in VIEW_SLICE_ORIENTATIONS:
            return
        backend.capture_view_slice(axis)
        sync_slice_bounds(recenter=True)
        apply_slice()

    ctrl.capture_view_slice = capture_view_slice

    def toggle_play():
        state.playing = not state.playing

    ctrl.toggle_play = toggle_play

    # trame >= 3.13 forwards the current server state as keyword arguments to
    # on_server_ready tasks.  The polling loop does not need that state, but it
    # must accept it to remain compatible with the callback contract.
    async def update_loop(**_):
        while True:
            await asyncio.sleep(0.05)
            old_times = list(state.times)
            old_index = int(state.time_index)
            was_at_tail = (
                not old_times or old_index >= len(old_times) - 1)
            changed = catalog.refresh()
            if changed:
                state.times = catalog.times
                new_tail = state.times[-1] if state.times else None
                old_tail = old_times[-1] if old_times else None
                if (was_at_tail and new_tail is not None and
                        (old_tail is None or
                         not math.isclose(new_tail, old_tail,
                                          rel_tol=1e-12, abs_tol=1e-12))):
                    load_index(len(state.times) - 1)
            previous_quality = state.quality
            state.quality = timeline.tick(int(time.monotonic() * 1000))
            if state.quality != previous_quality:
                sync_slice_bounds()
            state.cell_count = timeline.current_count
            state.displayed_cell_count = timeline.current_displayed_count
            state.vessel_count = timeline.current_vessel_count
            state.vessel_status = timeline.current_vessel_status
            state.quality_label = localized_viewer_status(
                language, "quality", state.quality)
            state.vessel_status_label = localized_viewer_status(
                language, "vessel", state.vessel_status)
            if state.playing and state.times:
                load_index((int(state.time_index) + 1) % len(state.times))
                await asyncio.sleep(0.05)
            ctrl.view_update()

    ctrl.on_server_ready.add_task(update_loop)

    with SinglePageWithDrawerLayout(server, width=320) as layout:
        layout.title.set_text(tr("title"))
        with layout.toolbar:
            vuetify3.VBtn(
                "{{ playing ? pause_label : play_label }}",
                click=ctrl.toggle_play)
            vuetify3.VSlider(
                v_model=("time_index", 0), min=0,
                max=("Math.max(times.length - 1, 0)",), step=1, hide_details=True,
                style="min-width: 360px",
            )
            html.Span(
                f"{tr('time')}={{{{ Number(time_hours).toFixed(3) }}}} h")
            html.Span(
                f"{tr('cells')}={{{{ cell_count.toLocaleString() }}}}")
            html.Span(
                f"{tr('displayed')}="
                "{{ displayed_cell_count.toLocaleString() }}")
            html.Span(
                f"{tr('vessel_nodes')}="
                "{{ vessel_count.toLocaleString() }}")
            html.Span("{{ quality_label }}")
            html.Span("{{ vessel_status_label }}")

        with layout.drawer:
            with vuetify3.VContainer(fluid=True, classes="pa-3"):
                html.Div(tr("camera"), classes="text-subtitle-2 mb-2")
                vuetify3.VSelect(
                    v_model=("navigation_mode", "rotate"),
                    items=("navigation_mode_options",),
                    item_title="title", item_value="value",
                    label=tr("primary_drag"), hide_details=True,
                    density="compact", classes="mb-2",
                )
                html.Div(
                    tr("camera_help"),
                    classes="text-caption mb-2",
                )
                vuetify3.VBtn(
                    tr("reset_camera"), click=ctrl.reset_camera,
                    size="small", variant="outlined", classes="mb-3",
                )

                vuetify3.VDivider(classes="mb-3")
                html.Div(tr("objects"), classes="text-subtitle-2 mb-1")
                vuetify3.VSwitch(
                    v_model=("show_r", True), label=tr("r_cells"), color="green",
                    hide_details=True, density="compact",
                )
                vuetify3.VSwitch(
                    v_model=("show_k", True), label=tr("k_cells"), color="red",
                    hide_details=True, density="compact",
                )
                vuetify3.VSwitch(
                    v_model=("show_vessels", True), label=tr("vessels"),
                    color="blue",
                    hide_details=True, density="compact",
                )
                vuetify3.VSwitch(
                    v_model=("show_influence", False),
                    label=(
                        f"{tr('influence')} "
                        f"({display_config['influence_cutoff_radius_voxels']:g} "
                        f"{tr('voxels')})"
                    ),
                    color="light-blue", hide_details=True, density="compact",
                    classes="mb-2",
                )
                vuetify3.VSlider(
                    v_model=("radius_scale", 4.0), min=0.1, max=8.0, step=0.1,
                    label=tr("cell_point_size"), hide_details=True,
                    classes="mb-2",
                )
                vuetify3.VSlider(
                    v_model=("vessel_radius_scale", 1.0), min=0.1,
                    max=4.0, step=0.1, label=tr("vessel_radius"),
                    hide_details=True, classes="mb-3",
                )

                vuetify3.VDivider(classes="mb-3")
                html.Div(tr("slice"), classes="text-subtitle-2 mb-2")
                vuetify3.VSelect(
                    v_model=("slice_mode", "whole"),
                    items=("slice_mode_options",),
                    item_title="title", item_value="value",
                    label=tr("mode"), hide_details=True,
                    density="compact", classes="mb-2",
                )
                vuetify3.VSelect(
                    v_model=("slice_axis", "Z"),
                    items=("slice_axis_options",),
                    item_title="title", item_value="value",
                    label=tr("orientation"), hide_details=True,
                    density="compact",
                    disabled=("slice_mode === 'whole'",), classes="mb-2",
                )
                vuetify3.VBtn(
                    tr("recapture_view"), click=ctrl.capture_view_slice,
                    disabled=(
                        "slice_mode === 'whole' || "
                        "!String(slice_axis).startsWith('VIEW_')",
                    ),
                    size="small", variant="outlined", classes="mb-3",
                )
                vuetify3.VSlider(
                    v_model=("slice_position", 0.0), min=("slice_min",),
                    max=("slice_max",), step=0.5,
                    label=tr("slice_offset"),
                    hide_details=True, disabled=("slice_mode === 'whole'",),
                    classes="mb-2",
                )
                vuetify3.VSlider(
                    v_model=("slice_thickness", 2.0), min=0.5,
                    max=("Math.max(slice_max - slice_min, 1)",), step=0.5,
                    label=tr("slab_thickness"), hide_details=True,
                    disabled=("slice_mode !== 'slab'",), classes="mb-2",
                )
                vuetify3.VSwitch(
                    v_model=("slice_invert", False), label=tr("invert_cut"),
                    hide_details=True, density="compact",
                    disabled=("slice_mode !== 'cut'",),
                )
        with layout.content:
            view = vtk_widgets.VtkRemoteView(
                backend.view, interactive_ratio=1,
                style="width: 100%; height: 100%;",
            )
            ctrl.view_update = view.update

    if state.times:
        # A running simulation should open at its latest available preview and
        # remain attached to the live tail. The user can still drag back to any
        # archived frame after the initial render.
        load_index(len(state.times) - 1)
        backend.pv.ResetCamera(backend.view)
    return server


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_directory", type=Path)
    parser.add_argument("--debounce-ms", type=int, default=250)
    parser.add_argument(
        "--language", choices=tuple(VIEWER_TRANSLATIONS), default="en")
    args, server_args = parser.parse_known_args(argv)
    try:
        server = build_app(
            args.run_directory, args.debounce_ms, args.language)
    except (RuntimeError, ValueError, OSError) as error:
        print(f"atcg3d-viewer: {error}", file=sys.stderr)
        return 1
    server.start(argv=server_args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
