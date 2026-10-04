#!/usr/bin/env python3
"""Voltage and current of an electrode from the electrodeVoltageCurrent
function object (1D, 2D or 3D cases).

    plot_electrodes.py CASE --patch electrode --frequency 40e6

Two figures in one: the waveforms over the last periods of the run (voltage,
conduction, displacement and total current), and one value per period of the
applied voltage over the whole run: the RMS of the alternating current, the
period-averaged (DC) current and the period-averaged power V*I. The second
shows whether the discharge has reached a periodic steady state.

Without --frequency only the raw history is plotted.

The function object must be in system/controlDict (see the README of the
repository); its output directory is <case>/<name>/<start time>/<patch>.dat,
and the pieces written by restarts are joined. Currents are positive from the
electrode into the plasma.
"""

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import somafoam_io as sio  # noqa: E402


def per_period(data, frequency):
    """One value per complete period: mid time, RMS of the alternating part
    of the total current, mean total and conduction current, mean V*I."""
    period = 1.0/frequency
    t = data["t"]
    first = int(np.ceil(t[0]/period - 1e-9))
    last = int(np.floor(t[-1]/period + 1e-9))
    edges = np.searchsorted(t, (np.arange(first, last + 1))*period)

    out = {k: [] for k in ("t", "rms", "dc", "dc_conduction", "power")}
    for i0, i1 in zip(edges[:-1], edges[1:]):
        if i1 - i0 < 8:
            continue
        current = data["I_total"][i0:i1]
        mean = current.mean()
        out["t"].append(0.5*(t[i0] + t[i1 - 1]))
        out["rms"].append(np.sqrt(np.mean((current - mean)**2)))
        out["dc"].append(mean)
        out["dc_conduction"].append(data["I_conduction"][i0:i1].mean())
        out["power"].append(np.mean(data["V"][i0:i1]*current))
    return {k: np.array(v) for k, v in out.items()}


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("case", help="case directory")
    parser.add_argument("--patch", default="electrode",
                        help="patch name (default: electrode)")
    parser.add_argument("--name", default="electrodes",
                        help="function object name (default: electrodes)")
    parser.add_argument("--frequency", type=float,
                        help="frequency of the applied voltage [Hz]")
    parser.add_argument("--periods", type=float, default=3,
                        help="periods shown in the waveform panel (default 3)")
    parser.add_argument("--tmin", type=float, help="first time [s]")
    parser.add_argument("--tmax", type=float, help="last time [s]")
    parser.add_argument("--title", help="figure title")
    parser.add_argument("--output", "-o", help="save the figure to this file")
    parser.add_argument("--csv", help="write the per-period values")
    args = parser.parse_args()

    if args.output:
        import matplotlib
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    data = sio.read_electrode(args.case, args.patch, args.name)
    keep = np.ones(len(data["t"]), dtype=bool)
    if args.tmin is not None:
        keep &= data["t"] >= args.tmin
    if args.tmax is not None:
        keep &= data["t"] <= args.tmax
    data = {k: v[keep] for k, v in data.items()}
    if len(data["t"]) < 2:
        sys.exit("no data in the selected time range")

    title = args.title or \
        f"{os.path.basename(os.path.abspath(args.case))}, patch {args.patch}"

    if not args.frequency:
        fig, ax = plt.subplots(2, 1, figsize=(10, 6.5), sharex=True)
        ax[0].plot(data["t"]*1e6, data["V"], color="k")
        ax[0].set_ylabel("voltage (V)")
        ax[1].plot(data["t"]*1e6, data["I_total"]*1e3, label="total")
        ax[1].plot(data["t"]*1e6, data["I_conduction"]*1e3,
                   label="conduction")
        ax[1].set_ylabel("current (mA)")
        ax[1].set_xlabel("time (µs)")
        ax[1].legend()
    else:
        period = 1.0/args.frequency
        pp = per_period(data, args.frequency)
        fig, ax = plt.subplots(2, 3, figsize=(16, 8))

        tail = data["t"] >= data["t"][-1] - args.periods*period
        phase = (data["t"][tail] - data["t"][tail][0])/period
        ax[0, 0].plot(phase, data["V"][tail], color="k")
        ax[0, 0].set_ylabel("voltage (V)")
        ax[0, 1].plot(phase, data["I_total"][tail]*1e3, label="total")
        ax[0, 1].plot(phase, data["I_conduction"][tail]*1e3,
                      label="conduction")
        ax[0, 1].plot(phase, data["I_displacement"][tail]*1e3,
                      label="displacement")
        ax[0, 1].set_ylabel("current (mA)")
        ax[0, 1].legend(fontsize=8)
        ax[0, 2].plot(data["V"][tail], data["I_total"][tail]*1e3)
        ax[0, 2].set_xlabel("voltage (V)")
        ax[0, 2].set_ylabel("total current (mA)")
        for a in ax[0, :2]:
            a.set_xlabel(f"time (periods), last {args.periods:g} periods")

        ax[1, 0].plot(pp["t"]*1e6, pp["rms"]*1e3)
        ax[1, 0].set_ylabel("alternating current, RMS per period (mA)")
        ax[1, 1].plot(pp["t"]*1e6, pp["dc"]*1e3, label="total")
        ax[1, 1].plot(pp["t"]*1e6, pp["dc_conduction"]*1e3,
                      lw=0.9, label="conduction")
        ax[1, 1].set_ylabel("period-averaged current (mA)")
        ax[1, 1].legend(fontsize=8)
        ax[1, 2].plot(pp["t"]*1e6, pp["power"])
        ax[1, 2].set_ylabel("period-averaged V·I (W)")
        for a in ax[1]:
            a.set_xlabel("time (µs)")

        if args.csv:
            np.savetxt(
                args.csv,
                np.column_stack([pp[k] for k in
                                 ("t", "rms", "dc", "dc_conduction", "power")]),
                delimiter=",", comments="",
                header="time_s,I_rms_A,I_dc_A,I_dc_conduction_A,power_W")

    for a in np.ravel(ax):
        a.grid(alpha=0.3)
    fig.suptitle(title)
    fig.tight_layout()

    if args.output:
        fig.savefig(args.output, dpi=150)
        print(f"wrote {args.output}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
