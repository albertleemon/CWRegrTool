"""Compare artifact-corrected EEG recorded INSIDE the scanner with EEG of the same
subject recorded OUTSIDE the scanner (van der Meer et al. 2016, "relative distance of
power spectra", eq. 4, plus band-wise and topographic comparisons).

Stages compared with the outside recording (all 500 Hz, 0.33-125 Hz, same EEG channels):
    "no BCG correction"  = 03_filtered.set   (gradient + OBS corrected, before CWL)
    "CWL"                = 04_cwl_corrected.set
The outside recording (BrainVision) is band-passed 0.33-125 Hz and resampled to 500 Hz
exactly like the pipeline. Optional --crop START END (s) restricts it, e.g. to drop the
part with the participant entering the scanner.

Metrics per stage:
  spectral_distance_pct   eq. 4 on the channel-mean spectrum, 0.5-125 Hz (linear power)
  spectral_distance_1_30Hz_pct  the same on 1-30 Hz only (above 30 Hz an unshielded
                          outside recording usually carries more line/EMG noise)
  band_diff_db            10*log10(inside/outside) of the channel-mean band power
  alpha_topo_distance_pct eq. 6 on the per-channel alpha (8-13 Hz) power, after scaling
                          both topographies to unit sum (shape, not amplitude)
  alpha_topo_r            spatial correlation of the per-channel alpha power
Eq. 5 (Jaccard of the EO/EC alpha surfaces) needs eyes-open AND eyes-closed segments in
both recordings; use metrics.evaluate() for that once those segments are known.

    python evaluate_outside.py OUT_DIR OUTSIDE.vhdr [--crop 10 300]
"""
import argparse
import json
import re
import tempfile
from pathlib import Path

import numpy as np
from scipy.signal import resample_poly, welch

from cwregression import bandpass
from metrics import relative_distance

EXCLUDE = re.compile(r"^(CW.*|ECG.*|EOG.*|65|66)$", re.I)
BANDS = [(1, 4), (4, 8), (8, 13), (13, 30), (30, 60), (60, 90), (90, 125)]


def load_set(path, trim=None):
    import mne
    raw = mne.io.read_raw_eeglab(path, preload=True, verbose="error")
    x = raw.get_data() * 1e6
    if trim is not None and x.shape[1] > trim[1]:
        x = x[:, trim[0] - 1:trim[1]]
    return raw.ch_names, x, raw.info["sfreq"]


def load_outside(path, crop=None):
    import mne
    from run_pipeline import resolve_header
    raw = mne.io.read_raw_brainvision(resolve_header(path, tempfile.mkdtemp(prefix="outside_")),
                                      preload=True, verbose="error")
    fs = raw.info["sfreq"]
    x = raw.get_data() * 1e6
    if crop:
        x = x[:, int(crop[0] * fs):int(crop[1] * fs)]
    x = bandpass(x, fs, 0.33, 125.0)
    g = int(np.gcd(int(fs), 500))
    x = resample_poly(x, 500 // g, int(fs) // g, axis=1)
    return raw.ch_names, x, 500.0


def spectra(names, x, fs, chans):
    f, p = welch(x[[names.index(c) for c in chans]], fs, nperseg=int(4 * fs), axis=-1)
    return f, p


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from pymatreader import read_mat
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_dir"); ap.add_argument("outside")
    ap.add_argument("--crop", type=float, nargs=2, metavar=("START", "END"))
    a = ap.parse_args()
    d = Path(a.out_dir).expanduser()

    m = read_mat(str(d / "04_cwl_corrected.set"), variable_names=["etc", "EEG"])
    m = m.get("EEG", m)
    trim = tuple(int(v) for v in m["etc"]["cwl_pipeline"]["trimmed_to"])
    stages = {"no BCG correction (03)": load_set(d / "03_filtered.set", trim),
              "CWL (04)": load_set(d / "04_cwl_corrected.set")}
    on, ox, ofs = load_outside(Path(a.outside).expanduser(), a.crop)
    chans = [c for c in stages["CWL (04)"][0]
             if not EXCLUDE.match(c) and c in on and all(c in s[0] for s in stages.values())]
    missing = [c for c in stages["CWL (04)"][0] if not EXCLUDE.match(c) and c not in on]
    print(f"{len(chans)} EEG channels in common" + (f"; not in outside recording: {missing}" if missing else ""))
    print(f"outside: {ox.shape[1] / ofs:.0f} s analysed")

    # Is the "outside" recording really outside the magnet? Carbon-wire loops only pick up
    # motion inside B0, so outside the scanner their 1-20 Hz power should be far lower.
    cw_out = [c for c in on if c.upper().startswith("CW")]
    n03, x03, fs03 = stages["no BCG correction (03)"]
    cw_in = [c for c in n03 if c.upper().startswith("CW") and c in cw_out]
    if cw_in:
        fo_, pco = spectra(on, ox, ofs, cw_in)
        fi_, pci = spectra(n03, x03, fs03, cw_in)
        lo = lambda f_, p_: p_[:, (f_ >= 1) & (f_ <= 20)].mean()
        ratio = float(lo(fo_, pco) / lo(fi_, pci))
        print(f"CWL 1-20 Hz power outside/inside: {ratio:.4f}"
              + ("  -> looks like a recording OUTSIDE the magnet" if ratio < 0.05 else
                 "  -> WARNING: loops still see motion; this may be INSIDE the magnet (no scanning)"))
    else:
        ratio = None

    f, po = spectra(on, ox, ofs, chans)
    b = (f >= 0.5) & (f <= 125)
    alpha = (f >= 8) & (f <= 13)
    topo_out = po[:, alpha].mean(1)
    res = {"n_channels": len(chans), "outside_seconds": ox.shape[1] / ofs,
           "cwl_power_1_20Hz_outside_over_inside": ratio, "stages": {}}
    spec = {}
    for name, (n, x, fs) in stages.items():
        _, p = spectra(n, x, fs, chans)
        spec[name] = p
        topo_in = p[:, alpha].mean(1)
        res["stages"][name] = dict(
            spectral_distance_pct=relative_distance(p.mean(0)[b], po.mean(0)[b]),
            spectral_distance_1_30Hz_pct=relative_distance(p.mean(0)[(f >= 1) & (f <= 30)],
                                                           po.mean(0)[(f >= 1) & (f <= 30)]),
            per_channel_spectral_distance_pct_median=float(np.median(
                [relative_distance(p[i, b], po[i, b]) for i in range(len(chans))])),
            band_diff_db={f"{lo}-{hi}": float(10 * np.log10(p.mean(0)[(f >= lo) & (f < hi)].mean() /
                                                            po.mean(0)[(f >= lo) & (f < hi)].mean()))
                          for lo, hi in BANDS},
            alpha_topo_distance_pct=relative_distance(topo_in / topo_in.sum(), topo_out / topo_out.sum()),
            alpha_topo_r=float(np.corrcoef(topo_in, topo_out)[0, 1]))
    (d / "outside_comparison.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))

    col = {"no BCG correction (03)": "#d62728", "CWL (04)": "#2ca02c"}
    fig = plt.figure(figsize=(13, 8), layout="constrained")
    ax = fig.add_subplot(2, 1, 1)
    ax.semilogy(f[b], po.mean(0)[b], "k", lw=1.3, label="outside the scanner")
    for name, p in spec.items():
        r = res["stages"][name]["spectral_distance_pct"]
        ax.semilogy(f[b], p.mean(0)[b], color=col[name], lw=1, label=f"{name}: diff = {r:.1f} %")
    ax.set(xlim=(0, 125), xlabel="Hz", ylabel="µV²/Hz", title=f"mean spectrum over {len(chans)} EEG channels")
    ax.legend()
    ax = fig.add_subplot(2, 2, 3)
    x = np.arange(len(BANDS))
    for k, name in enumerate(spec):
        v = list(res["stages"][name]["band_diff_db"].values())
        ax.bar(x + (k - .5) * .4, v, .4, color=col[name], label=name)
    ax.axhline(0, color="k", lw=.6)
    ax.set_xticks(x, [f"{lo}-{hi}" for lo, hi in BANDS])
    ax.set(ylabel="inside - outside (dB)", xlabel="Hz", title="band power vs outside")
    ax.legend(fontsize=8)
    ax = fig.add_subplot(2, 2, 4)
    for name, p in spec.items():
        ax.loglog(topo_out, p[:, alpha].mean(1), "o", ms=4, color=col[name],
                  label=f"{name}: r = {res['stages'][name]['alpha_topo_r']:.2f}")
    lim = [min(topo_out.min(), min(p[:, alpha].mean(1).min() for p in spec.values())) * .8,
           max(topo_out.max(), max(p[:, alpha].mean(1).max() for p in spec.values())) * 1.2]
    ax.plot(lim, lim, "k--", lw=.6)
    ax.set(xlabel="alpha power outside (µV²/Hz)", ylabel="alpha power inside",
           title="per-channel alpha (8-13 Hz) power")
    ax.legend(fontsize=8)
    fig.savefig(d / "outside_comparison.png", dpi=110)
    print("wrote", d / "outside_comparison.png", "and outside_comparison.json")


if __name__ == "__main__":
    main()
