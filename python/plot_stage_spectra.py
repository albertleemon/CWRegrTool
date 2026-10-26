"""Mean EEG spectrum at each processing stage of one pipeline run:
  raw -> 02 gradient template subtraction -> 02b OBS residual removal (if run)
      -> 04 CWL regression (-> 05 BCG AAS, if run), plus the null control
      (regression with random-phase CWL).

"After CWL" and "null" come from diagnostics_data.mat (matlab/cwl_diagnose.m). The raw,
02 and 02b spectra are computed from the 5000 Hz data on exactly the same samples and
EEG channels, after the same 0.33-125 Hz band-pass and resampling to 500 Hz as the
pipeline, so the curves differ only by the correction steps.

    python plot_stage_spectra.py OUT_DIR [--vhdr raw.vhdr] [--title ...]
"""
import argparse
import re
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from scipy.signal import resample_poly, welch

from cwregression import bandpass
from run_manifest import has_stage

EXCLUDE = re.compile(r"^(CW.*|ECG.*|EOG.*|65|66)$", re.I)


def stage_spectrum(raw, start, t0, t1, f_target):
    """Mean EEG spectrum of a 5000 Hz recording, cut at `start` (first volume), band-passed,
    resampled to 500 Hz and restricted to the samples t0..t1 (1-based) kept by the pipeline."""
    fs = raw.info["sfreq"]
    eeg = [i for i, n in enumerate(raw.ch_names) if not EXCLUDE.match(n)]
    up, down = 500, int(fs)
    g = np.gcd(up, down); up, down = up // g, down // g
    P = []
    for i in eeg:
        x = raw.get_data(picks=[i])[0, start:] * 1e6
        x = resample_poly(bandpass(x[None], fs, 0.33, 125.0)[0], up, down)[t0 - 1:t1]
        f, p = welch(x, 500.0, nperseg=2000, noverlap=1000)
        P.append(np.interp(f_target, f, p))
    return np.mean(P, 0)


def stage_spectra(out_dir, vhdr, f_target):
    import mne
    from pymatreader import read_mat
    from run_pipeline import resolve_header, volume_onsets
    m = read_mat(str(Path(out_dir) / "04_cwl_corrected.set"), variable_names=["etc", "EEG"])
    m = m.get("EEG", m)
    t0, t1 = (int(v) for v in m["etc"]["cwl_pipeline"]["trimmed_to"])   # 1-based, 500 Hz
    out = []
    if vhdr:
        raw = mne.io.read_raw_brainvision(resolve_header(vhdr, out_dir), preload=True, verbose="error")
        on = volume_onsets(raw, "R128")
        out.append(("raw (no correction)", stage_spectrum(raw, on[0], t0, t1, f_target), "#7f7f7f", "-"))
        del raw
    for fname, lab, col in (("02_gradient_corrected.set", "after gradient template subtraction", "#d62728"),
                            ("02b_obs_corrected.set", "after OBS residual removal", "#ff7f0e")):
        fp = Path(out_dir) / fname
        if has_stage(out_dir, fp.stem):
            raw = mne.io.read_raw_eeglab(fp, preload=True, verbose="error")
            out.append((lab, stage_spectrum(raw, 0, t0, t1, f_target), col, "-"))
            del raw
    return out


def set_spectrum(fp, f_target):
    """Mean EEG spectrum of a pipeline output that is already 500 Hz and trimmed (05)."""
    import mne
    raw = mne.io.read_raw_eeglab(fp, preload=True, verbose="error")
    eeg = [i for i, n in enumerate(raw.ch_names) if not EXCLUDE.match(n)]
    f, p = welch(raw.get_data(picks=eeg) * 1e6, raw.info["sfreq"], nperseg=2000, noverlap=1000, axis=-1)
    return np.mean([np.interp(f_target, f, q) for q in p], 0)


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
    curves = stage_spectra(d, Path(a.vhdr).expanduser() if a.vhdr else None, f)
    if not any(c[0].startswith("after gradient") for c in curves):
        curves.append(("after gradient correction (03)", m["p0"].mean(1), "#d62728", "-"))
    curves += [("after CWL regression", m["p1"].mean(1), "#2ca02c", "-"),
               ("null: regression with random-phase CWL", m["pn"].mean(1), "#1f77b4", "--")]
    if has_stage(d, "05_bcg_aas"):
        curves.insert(-1, ("after CWL + BCG AAS (05)", set_spectrum(d / "05_bcg_aas.set", f), "#9467bd", "-"))
    b = (f >= 0.5) & (f <= 125)
    fig, ax = plt.subplots(figsize=(11, 5.5), layout="constrained")
    for lab, p, c, ls in curves:
        ax.semilogy(f[b], p[b], color=c, ls=ls, lw=1.1, label=lab)
    ax.set(xlim=(0, 125), xlabel="Frequency (Hz)", ylabel="Power (µV²/Hz)",
           title=a.title or f"Mean spectrum over {m['p0'].shape[1]} EEG channels, {d.name}")
    ax.legend(fontsize=9)
    ax.grid(alpha=.3, which="both", lw=.4)
    fig.savefig(d / "stage_spectra.png", dpi=120)
    rows = [(lo, hi) for lo, hi in [(1, 4), (4, 8), (8, 13), (13, 30), (30, 60), (60, 125)]]
    print("band power (µV²/Hz, mean)             " + "  ".join(f"{lo:>3}-{hi:<3}" for lo, hi in rows))
    for lab, p, *_ in curves:
        print(f"{lab[:38]:38s} " + "  ".join(f"{p[(f >= lo) & (f < hi)].mean():7.3g}" for lo, hi in rows))
    print("wrote", d / "stage_spectra.png")


if __name__ == "__main__":
    main()
