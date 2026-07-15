"""Load ATCG3D cells and vascular centerlines in native ParaView.

Run from ParaView's Python shell or with pvpython. This native mode opens one
series directly; automatic preview/full switching is provided by viewer/app.py.
"""

from argparse import ArgumentParser
from pathlib import Path

from paraview.simple import (ColorBy, GetActiveViewOrCreate, GetAnimationScene,
                             Calculator, OpenDataFile, Render, ResetCamera, Show, Tube)


def main():
    parser = ArgumentParser()
    parser.add_argument("run_directory", type=Path)
    parser.add_argument("--quality", choices=("preview", "full"), default="preview")
    parser.add_argument("--vessel-color", choices=("perfused", "branch_role"),
                        default="perfused")
    parser.add_argument("--cell-radius-scale", type=float, default=1.0)
    parser.add_argument("--vessel-radius-scale", type=float, default=1.0)
    args = parser.parse_args()
    if args.cell_radius_scale <= 0.0 or args.vessel_radius_scale <= 0.0:
        parser.error("cell and vessel radius scales must be positive")
    series = args.run_directory / f"{args.quality}.vtkhdf.series"
    source = OpenDataFile(str(series))
    view = GetActiveViewOrCreate("RenderView")
    cell_radius_calculator = Calculator(Input=source)
    cell_radius_calculator.ResultArrayName = "viewer_display_radius"
    cell_radius_calculator.Function = (
        f"display_radius*{args.cell_radius_scale:.17g}"
    )
    display = Show(cell_radius_calculator, view)
    display.Representation = "Point Gaussian"
    display.GaussianRadius = 0.5
    display.SetScaleArray = ["POINTS", "viewer_display_radius"]
    display.ScaleByArray = 1
    display.UseScaleFunction = 0
    ColorBy(display, ("POINTS", "cell_type"))

    vessel_source = None
    vessel_radius_calculator = None
    vessel_tube = None
    vessel_display = None
    vessel_series = args.run_directory / "vessels.vtkhdf.series"
    if vessel_series.exists():
        vessel_source = OpenDataFile(str(vessel_series))
        vessel_radius_calculator = Calculator(Input=vessel_source)
        vessel_radius_calculator.ResultArrayName = "viewer_radius_voxels"
        vessel_radius_calculator.Function = (
            f"radius_voxels*{args.vessel_radius_scale:.17g}"
        )
        vessel_tube = Tube(Input=vessel_radius_calculator)
        vessel_tube.NumberofSides = 12
        vessel_tube.Capping = 1
        vessel_tube.Radius = 0.5 * args.vessel_radius_scale
        try:
            vessel_tube.VaryRadius = "By Absolute Scalar"
            vessel_tube.Scalars = ["POINTS", "viewer_radius_voxels"]
        except Exception:  # ParaView-version dependent property exposure.
            pass
        vessel_display = Show(vessel_tube, view)
        ColorBy(vessel_display, ("POINTS", args.vessel_color))

    GetAnimationScene().UpdateAnimationUsingDataTimeSteps()
    ResetCamera(view)
    Render(view)
    return (source, cell_radius_calculator, display, vessel_source, vessel_radius_calculator,
            vessel_tube, vessel_display, view)


if __name__ == "__main__":
    main()
