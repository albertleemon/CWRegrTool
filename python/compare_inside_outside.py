"""Simple side-by-side comparison: EEG recorded INSIDE the MRI (any pipeline stage, e.g.
04_cwl_corrected.set) vs EEG recorded OUTSIDE the MRI (raw BrainVision file).

The outside recording gets only what the inside data got anyway (0.33-125 Hz band-pass,
resampling to 500 Hz), so the comparison is like for like. Writes one figure:
  top:    10 s of the same channels, inside and outside, same µV scale
  bottom: mean spectrum over the common EEG channels, inside vs outside
and prints the band power of both.

    python compare_inside_outside.py INSIDE.set OUTSIDE.vhdr OUT.png [--channels Fz Cz Pz Oz O1 O2]
                                     [--start 60] [--seconds 10]
"""
import argparse
import re
import tempfile
from pathlib import Path

import numpy as np
from scipy.signal import resample_poly, welch

from cwregression import bandpass

EXCLUDE = re.compile(r"^(CW.*|ECG.*|EOG.*|65|66)$", re.I)
BANDS = [(1, 4), (4, 8), (8, 13), (13, 30), (30, 60), (60, 125)]


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import mne
    from run_pipeline import resolve_header
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inside"); ap.add_argument("outside"); ap.add_argument("out_png")
    ap.add_argument("--channels", nargs="+", default=["Fz", "Cz", "Pz", "Oz", "O1", "O2"])
    ap.add_argument("--start", type=float, default=60.0, help="start of the example traces (s)")
    ap.add_argument("--seconds", type=float, default=10.0)
    a = ap.parse_args()

    ins = mne.io.read_raw_eeglab(Path(a.inside).expanduser(), preload=True, verbose="error")
    xi, fi = ins.get_data() * 1e6, ins.info["sfreq"]
    out = mne.io.read_raw_brainvision(resolve_header(Path(a.outside).expanduser(),
                                      tempfile.mkdtemp(prefix="cmp_")), preload=True, verbose="error")
    fo = out.info["sfreq"]
    xo = bandpass(out.get_data() * 1e6, fo, 0.33, 125.0)
    if fo != fi:
        g = int(np.gcd(int(fo), int(fi)))
        xo = resample_poly(xo, int(fi) // g, int(fo) // g, axis=1)
    chans = [c for c in ins.ch_names if not EXCLUDE.match(c) and c in out.ch_names]
    show = [c for c in a.channels if c in chans]
    print(f"{len(chans)} common EEG channels; inside {xi.shape[1] / fi:.0f} s, outside {xo.shape[1] / fi:.0f} s")

    pick = lambda names, x, cs: x[[names.index(c) for c in cs]]
    f, pi = welch(pick(ins.ch_names, xi, chans), fi, nperseg=int(4 * fi), axis=-1)
    _, po = welch(pick(out.ch_names, xo, chans), fi, nperseg=int(4 * fi), axis=-1)

    fig, axd = plt.subplot_mosaic([["in", "out"], ["spec", "spec"]], figsize=(14, 9),
                                  layout="constrained", height_ratios=[1, 1])
    s0, n = int(a.start * fi), int(a.seconds * fi)
    t = np.arange(n) / fi
    ti, to = pick(ins.ch_names, xi, show)[:, s0:s0 + n], pick(out.ch_names, xo, show)[:, s0:s0 + n]
    step = 4 * np.median(np.std(np.vstack([ti, to]), 1))
    for k, (tr, title, col) in enumerate([(ti, "inside the MRI", "#2ca02c"), (to, "outside the MRI", "k")]):
        ax = axd["in" if k == 0 else "out"]
        for j, c in enumerate(show):
            ax.plot(t, tr[j] - j * step, color=col, lw=.6)
        ax.set_yticks([-j * step for j in range(len(show))], show)
        ax.set(title=f"{title}, {a.seconds:g} s from {a.start:g} s (offset {step:.0f} µV per trace)",
               xlabel="s", ylim=(-(len(show) - 0.3) * step, 0.9 * step))
    ax = axd["spec"]
    b = (f >= 0.5) & (f <= 125)
    ax.semilogy(f[b], pi.mean(0)[b], color="#2ca02c", lw=1.2, label=f"inside the MRI ({Path(a.inside).name})")
    ax.semilogy(f[b], po.mean(0)[b], "k", lw=1.2, label="outside the MRI (raw, 0.33-125 Hz)")
    ax.set(xlim=(0, 125), xlabel="Hz", ylabel="µV²/Hz", title=f"mean spectrum over {len(chans)} EEG channels")
    ax.grid(alpha=.3, which="both", lw=.4); ax.legend()
    fig.savefig(a.out_png, dpi=110)

    print(f"{'band (Hz)':10s} {'inside':>10s} {'outside':>10s} {'diff dB':>8s}   (µV²/Hz, channel mean)")
    for lo, hi in BANDS:
        w = (f >= lo) & (f < hi)
        vi, vo = pi.mean(0)[w].mean(), po.mean(0)[w].mean()
        print(f"{lo:>4}-{hi:<5} {vi:10.3g} {vo:10.3g} {10 * np.log10(vi / vo):+8.1f}")
    print("wrote", a.out_png)


if __name__ == "__main__":
    main()
