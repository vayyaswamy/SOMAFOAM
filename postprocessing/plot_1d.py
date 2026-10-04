#!/usr/bin/env python3
"""Profiles and time histories of a one-dimensional SOMAFOAM case.

Profiles of several fields at one or more times:

    plot_1d.py CASE --fields electron Arp1 Arm Te Phi --times 1e-6 latest

Time histories (maximum, mean over the gap and mid-gap value of each field at
every written time):

    plot_1d.py CASE --fields electron Arm Te --history

Fields are the names of the files in the time directories: species densities
(electron, Arp1, Arm, ...), Te, Phi, E, and the period-averaged fields written
by the fieldAverage function object (TeMean, N_electronMean, ...). For a
vector field the component along the gap is plotted.

The figure is shown, or saved with --output file.png. --csv file.csv also
writes the plotted data.
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import somafoam_io as sio  # noqa: E402


def gap_coordinate(mesh):
    """Coordinate along the gap [m] and the order of the cells along it."""
    axis = mesh.axes_by_extent()[0]
    x = mesh.centres[:, axis]
    order = np.argsort(x)
    return axis, x[order], order


def load_profile(case, field, time):
    mesh = case.mesh(time)
    axis, x, order = gap_coordinate(mesh)
    values = case.scalar(field, time, component=axis)
    return x, values[order]


def cell_widths(x):
    """Widths of the cells of a 1D mesh from their centres (for averages on
    a refined mesh)."""
    edges = np.empty(len(x) + 1)
    edges[1:-1] = 0.5*(x[1:] + x[:-1])
    edges[0] = x[0] - (edges[1] - x[0])
    edges[-1] = x[-1] + (x[-1] - edges[-2])
    return np.diff(edges)


def grid(n, width=5.5, height=3.8):
    import matplotlib.pyplot as plt

    cols = min(n, 3)
    rows = int(np.ceil(n/cols))
    fig, axes = plt.subplots(
        rows, cols, figsize=(width*cols, height*rows), squeeze=False)
    for a in axes.flat[n:]:
        a.set_visible(False)
    return fig, axes.flat[:n]


def plot_profiles(case, args):
    times = [case.time_name(t) for t in args.times]
    fig, axes = grid(len(args.fields))
    rows = []

    for ax, field in zip(axes, args.fields):
        for time in times:
            if not case.has_field(field, time):
                print(f"skipping {field} at {time}: not written")
                continue
            x, v = load_profile(case, field, time)
            ax.plot(x*1e3, v, label=f"{float(time)*1e6:g} µs")
            rows += [(field, time, xi, vi) for xi, vi in zip(x, v)]
        ax.set_xlabel("x (mm)")
        ax.set_ylabel(sio.label(field))
        ax.grid(alpha=0.3)
        if field in args.log:
            ax.set_yscale("log")
        if len(times) > 1:
            ax.legend(fontsize=8)

    title = args.title or os.path.basename(case.path)
    if len(times) == 1:
        title += f", t = {float(times[0])*1e6:g} µs"
    fig.suptitle(title)
    fig.tight_layout()

    if args.csv:
        with open(args.csv, "w") as f:
            f.write("field,time_s,x_m,value\n")
            for field, time, xi, vi in rows:
                f.write(f"{field},{time},{xi:.9e},{vi:.9e}\n")
    return fig


def plot_history(case, args):
    fig, axes = grid(len(args.fields))
    rows = []

    for ax, field in zip(axes, args.fields):
        times = case.times_with(field)
        if args.tmin is not None:
            times = [t for t in times if float(t) >= args.tmin]
        if args.tmax is not None:
            times = [t for t in times if float(t) <= args.tmax]
        if not times:
            print(f"skipping {field}: not written")
            continue

        peak, mean, mid = [], [], []
        for time in times:
            x, v = load_profile(case, field, time)
            w = cell_widths(x)
            peak.append(np.abs(v).max())
            mean.append(np.sum(v*w)/np.sum(w))
            mid.append(np.interp(0.5*(x[0] + x[-1]), x, v))

        t = np.array([float(n) for n in times])*1e6
        ax.plot(t, peak, label="maximum")
        ax.plot(t, mean, label="gap average")
        ax.plot(t, mid, label="mid-gap")
        ax.set_xlabel("time (µs)")
        ax.set_ylabel(sio.label(field))
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
        if field in args.log:
            ax.set_yscale("log")
        rows += [(field, n, a, b, c)
                 for n, a, b, c in zip(times, peak, mean, mid)]

    fig.suptitle(args.title or os.path.basename(case.path))
    fig.tight_layout()

    if args.csv:
        with open(args.csv, "w") as f:
            f.write("field,time_s,maximum,gap_average,mid_gap\n")
            for field, n, a, b, c in rows:
                f.write(f"{field},{n},{a:.9e},{b:.9e},{c:.9e}\n")
    return fig


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("case", help="case directory")
    parser.add_argument(
        "--fields", nargs="+",
        default=["electron", "Arp1", "Arm", "Te", "Phi"],
        help="field file names (default: electron Arp1 Arm Te Phi)")
    parser.add_argument(
        "--times", nargs="+", default=["latest"],
        help="times in seconds, or 'latest'; the nearest written time is "
             "used (default: latest)")
    parser.add_argument(
        "--history", action="store_true",
        help="plot time histories instead of profiles")
    parser.add_argument("--tmin", type=float, help="history: first time [s]")
    parser.add_argument("--tmax", type=float, help="history: last time [s]")
    parser.add_argument(
        "--log", nargs="*", default=[],
        help="fields to plot on a logarithmic axis")
    parser.add_argument("--title", help="figure title")
    parser.add_argument("--output", "-o", help="save the figure to this file")
    parser.add_argument("--csv", help="also write the plotted data")
    args = parser.parse_args()

    if args.output:
        import matplotlib
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    case = sio.Case(args.case)
    fig = plot_history(case, args) if args.history else \
        plot_profiles(case, args)

    if args.output:
        fig.savefig(args.output, dpi=150)
        print(f"wrote {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
