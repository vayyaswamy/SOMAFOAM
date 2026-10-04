#!/usr/bin/env python3
"""Maps and line cuts of a two-dimensional (one cell thick) SOMAFOAM case.

Colour maps of several fields at one time, drawn cell by cell so that
adaptively refined meshes with hanging nodes are shown as they are:

    plot_2d.py CASE --fields electron Arp1 Te --time latest --log electron

With the mesh drawn on top, and the refinement level as a field:

    plot_2d.py CASE --fields electron cellLevel --mesh

A cut along a straight line (coordinates in metres), sampled from the cell
values:

    plot_2d.py CASE --fields electron Te --line 0 4e-4 2e-3 4e-4

Parallel cases are read from the processor directories if the fields have
not been reconstructed. For a vector field the magnitude is plotted;
--component 0, 1 or 2 selects a component.

The figure is shown, or saved with --output file.png.
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import somafoam_io as sio  # noqa: E402


def load(case, field, time, component=None):
    """Polygons, cell centres (in the plane) and values of a field at a
    time, gathered from the case or from its processors."""
    parts = case.parts(field, time)
    axes = None
    polygons, centres, values = [], [], []
    for part in parts:
        mesh = part.mesh(time)
        if axes is None:
            axes = sorted(mesh.axes_by_extent()[:2])
        polygons += mesh.polygons(axes)
        centres.append(mesh.centres[:, axes])
        values.append(part.scalar(field, time, component))
    return axes, polygons, np.vstack(centres), np.concatenate(values)


def draw_map(ax, polygons, values, log, mesh_lines, cmap):
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import LogNorm, Normalize

    scaled = [p*1e3 for p in polygons]
    if log:
        positive = values[values > 0]
        vmin = positive.min() if positive.size else 1e-300
        norm = LogNorm(vmin=vmin, vmax=max(values.max(), vmin*10))
        values = np.maximum(values, vmin)
    else:
        norm = Normalize(vmin=values.min(), vmax=values.max())

    collection = PolyCollection(
        scaled, array=values, cmap=cmap, norm=norm,
        edgecolors="k" if mesh_lines else "face",
        linewidths=0.15 if mesh_lines else 0.2)
    ax.add_collection(collection)
    points = np.vstack(scaled)
    ax.set_xlim(points[:, 0].min(), points[:, 0].max())
    ax.set_ylim(points[:, 1].min(), points[:, 1].max())
    ax.set_aspect("equal")
    return collection


def sample_line(polygons, centres, values, start, end, n):
    """Values along a line, from the cell containing each sample point (the
    nearest cell centre is used, which is exact for the Cartesian and
    refined-Cartesian meshes these cases use, and a good approximation
    otherwise)."""
    s = np.linspace(0.0, 1.0, n)
    points = np.outer(1 - s, start) + np.outer(s, end)
    out = np.empty(n)
    for i, p in enumerate(points):
        out[i] = values[np.argmin(np.sum((centres - p)**2, axis=1))]
    distance = s*np.linalg.norm(np.asarray(end) - np.asarray(start))
    return distance, out


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("case", help="case directory")
    parser.add_argument(
        "--fields", nargs="+", default=["electron", "Arp1", "Te", "Phi"],
        help="field file names, or cellLevel "
             "(default: electron Arp1 Te Phi)")
    parser.add_argument(
        "--time", default="latest",
        help="time in seconds, or 'latest'; the nearest written time is used")
    parser.add_argument("--log", nargs="*", default=[],
                        help="fields to show on a logarithmic scale")
    parser.add_argument("--component", type=int, choices=[0, 1, 2],
                        help="component of vector fields (default: magnitude)")
    parser.add_argument("--mesh", action="store_true",
                        help="draw the cell outlines")
    parser.add_argument(
        "--line", nargs=4, type=float, metavar=("X0", "Y0", "X1", "Y1"),
        help="plot a cut along this line [m] instead of maps")
    parser.add_argument("--samples", type=int, default=400,
                        help="points along the line (default 400)")
    parser.add_argument("--cmap", default="viridis", help="colour map")
    parser.add_argument("--title", help="figure title")
    parser.add_argument("--output", "-o", help="save the figure to this file")
    args = parser.parse_args()

    if args.output:
        import matplotlib
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    case = sio.Case(args.case)
    time = sio.nearest_time(sio.all_times(case), args.time)

    n = len(args.fields)
    names = "xyz"
    title = args.title or os.path.basename(case.path)
    title += f", t = {float(time)*1e6:g} µs"

    if args.line:
        cols = min(n, 3)
        rows = int(np.ceil(n/cols))
        fig, axes = plt.subplots(rows, cols, figsize=(5.5*cols, 3.8*rows),
                                 squeeze=False)
        for ax, field in zip(axes.flat, args.fields):
            _, polygons, centres, values = load(
                case, field, time, args.component)
            d, v = sample_line(polygons, centres, values,
                               args.line[:2], args.line[2:], args.samples)
            ax.plot(d*1e3, v)
            ax.set_xlabel("distance along the line (mm)")
            ax.set_ylabel(sio.label(field, args.component))
            ax.grid(alpha=0.3)
            if field in args.log:
                ax.set_yscale("log")
        for ax in axes.flat[n:]:
            ax.set_visible(False)
        title += (f", line ({args.line[0]:g}, {args.line[1]:g}) to "
                  f"({args.line[2]:g}, {args.line[3]:g}) m")
    else:
        first = load(case, args.fields[0], time, args.component)
        points = np.vstack(first[1])
        size = points.max(axis=0) - points.min(axis=0)
        aspect = size[1]/size[0]
        cols = 1 if aspect < 0.6 else min(n, 3)
        rows = int(np.ceil(n/cols))
        width = 10.0 if cols == 1 else 5.5
        fig, axes = plt.subplots(
            rows, cols,
            figsize=(width*cols, max(width*aspect, 1.6)*rows + 0.6),
            squeeze=False)
        for i, (ax, field) in enumerate(zip(axes.flat, args.fields)):
            axes_used, polygons, _, values = first if i == 0 else load(
                case, field, time, args.component)
            collection = draw_map(ax, polygons, values, field in args.log,
                                  args.mesh, args.cmap)
            fig.colorbar(collection, ax=ax, pad=0.02,
                         label=sio.label(field, args.component))
            ax.set_xlabel(f"{names[axes_used[0]]} (mm)")
            ax.set_ylabel(f"{names[axes_used[1]]} (mm)")
        for ax in axes.flat[n:]:
            ax.set_visible(False)
        title += f", {len(polygons)} cells"

    fig.suptitle(title)
    fig.tight_layout()

    if args.output:
        fig.savefig(args.output, dpi=150)
        print(f"wrote {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
