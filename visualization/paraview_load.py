"""Load an ATCG3D ParaView series with Point Gaussian rendering.

Run from ParaView's Python shell or with pvpython. This native mode opens one
series directly; automatic preview/full switching is provided by viewer/app.py.
"""

from argparse import ArgumentParser
from pathlib import Path

from paraview.simple import (ColorBy, GetActiveViewOrCreate, OpenDataFile,
                             Render, ResetCamera, Show)


def main():
    parser = ArgumentParser()
    parser.add_argument("run_directory", type=Path)
    parser.add_argument("--quality", choices=("preview", "full"), default="preview")
    args = parser.parse_args()
    series = args.run_directory / f"{args.quality}.vtkhdf.series"
    source = OpenDataFile(str(series))
    view = GetActiveViewOrCreate("RenderView")
    display = Show(source, view)
    display.Representation = "Point Gaussian"
    display.GaussianRadius = 0.5
    display.SetScaleArray = ["POINTS", "display_radius"]
    display.ScaleTransferFunction = "PiecewiseFunction"
    ColorBy(display, ("POINTS", "cell_type"))
    ResetCamera(view)
    Render(view)
    return source, display, view


if __name__ == "__main__":
    main()
