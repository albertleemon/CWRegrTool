"""Mean EEG spectrum at each processing stage of one pipeline run:
  raw (no correction) -> after gradient-artifact correction -> after CWL regression,
plus the null control (regression with random-phase CWL).

The last three come from diagnostics_data.mat (matlab/cwl_diagnose.m). The raw spectrum
is computed from the BrainVision recording on exactly the same samples and EEG
channels, after the same 0.33-125 Hz band-pass and resampling to 500 Hz, so that the
only difference between raw and "after gradient correction" is the correction itself.

    python plot_stage_spectra.py OUT_DIR [--vhdr raw.vhdr] [--title ...]
"""
import argparse
import re
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from scipy.signal import resample_poly, welch

from cwregression import bandpass

EXCLUDE = re.compile(r"^(CW.*|ECG.*|EOG.*|65|66)$", re.I)


def raw_spectrum(vhdr, out_dir, f_target):
    import mne
    from pymatreader import read_mat
    from run_pipeline import resolve_header, volume_onsets
    raw = mne.io.read_raw_brainvision(resolve_header(vhdr, out_dir), preload=True, verbose="error")
    fs = raw.info["sfreq"]
    on = volume_onsets(raw, "R128")
    m = read_mat(str(Path(out_dir) / "04_cwl_corrected.set"))
    m = m.get("EEG", m)
    t0, t1 = (int(v) for v in m["etc"]["cwl_pipeline"]["trimmed_to"])   # 1-based, 500 Hz
    eeg = [i for i, n in enumerate(raw.ch_names) if not EXCLUDE.match(n)]
    # the pipeline cut the data at the first volume marker, then resampled to 500 Hz
    a = on[0]
    b = on[-1] + int(np.median(np.diff(on)))
    up, down = 500, int(fs)
    g = np.gcd(up, down); up, down = up // g, down // g
    P = []
    for i in eeg:
        x = raw.get_data(picks=[i])[0, a:b] * 1e6
        x = resample_poly(bandpass(x[None], fs, 0.33, 125.0)[0], up, down)
        x = x[t0 - 1:t1]
        f, p = welch(x, 500.0, nperseg=2000, noverlap=1000)
        P.append(np.interp(f_target, f, p))
    return np.mean(P, 0), len(eeg)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_dir")
    ap.add_argument("--vhdr", help="raw BrainVision header of the same run (adds the raw curve)")
    ap.add_argument("--title", default="")
    a = ap.parse_args()
    d = Path(a.out_dir).expanduser()
    m = loadmat(d / "diagnostics_data.mat", squeeze_me=True)
    f = m["f"]
    curves = []
    if a.vhdr:
        pr, n = raw_spectrum(Path(a.vhdr).expanduser(), d, f)
        curves.append(("raw (no correction)", pr, "#7f7f7f", "-"))
    curves += [("after gradient correction", m["p0"].mean(1), "#d62728", "-"),
               ("after CWL regression", m["p1"].mean(1), "#2ca02c", "-"),
               ("null: regression with random-phase CWL", m["pn"].mean(1), "#1f77b4", "--")]
    b = (f >= 0.5) & (f <= 125)
    fig, ax = plt.subplots(figsize=(11, 5.5), layout="constrained")
    for lab, p, c, ls in curves:
        ax.semilogy(f[b], p[b], color=c, ls=ls, lw=1.1, label=lab)
    ax.set(xlim=(0, 125), xlabel="Frequency (Hz)", ylabel="Power (µV²/Hz)",
           title=a.title or f"Mean spectrum over {m['p0'].shape[1]} EEG channels, {d.name}")
    ax.legend()
    ax.grid(alpha=.3, which="both", lw=.4)
    fig.savefig(d / "stage_spectra.png", dpi=120)
    rows = [(lo, hi) for lo, hi in [(1, 4), (4, 8), (8, 13), (13, 30), (30, 60), (60, 125)]]
    print("band power (µV²/Hz, mean)       " + "  ".join(f"{lo:>3}-{hi:<3}" for lo, hi in rows))
    for lab, p, *_ in curves:
        print(f"{lab[:32]:32s} " + "  ".join(f"{p[(f >= lo) & (f < hi)].mean():7.3g}" for lo, hi in rows))
    print("wrote", d / "stage_spectra.png")


if __name__ == "__main__":
    main()
