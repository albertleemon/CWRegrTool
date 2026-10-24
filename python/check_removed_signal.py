"""Did the CWL regression remove brain signal from a channel?

For each EEG channel the removed component is  r = 03_filtered - 04_cwl_corrected.
If r were pure artifact it would be unrelated to the cleaned EEG of the neighbouring
electrodes; if it contains brain activity, it is coherent with them (brain signals are
spatially smooth over a few cm, motion artifacts of different leads are not identical).

Per channel: mean magnitude-squared coherence between its removed component and the
CLEANED signals (04) of its nearest neighbours, in 1-4, 4-8, 8-13 and 13-30 Hz.
Neighbours = the k nearest electrodes on the standard 10-05 montage.

    python check_removed_signal.py OUT_DIR [--channels Cz FC1 FC2 C1 C2 CPz Pz Oz] [--k 4]
                                   [--outside OUTSIDE.vhdr]

With --outside, the coherence of each channel with the same neighbours is also computed
on the outside-the-scanner recording: if a channel is coherent with its neighbours
outside but not inside, it lost brain signal (or gained independent noise) inside.
Note: the removed-part coherence only detects a loss that is a sizeable fraction of the
removed component; a small loss hidden in a large artifact stays below ~0.05.
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np
from scipy.signal import coherence

EXCLUDE = re.compile(r"^(CW.*|ECG.*|EOG.*|65|66)$", re.I)
BANDS = [(1, 4), (4, 8), (8, 13), (13, 30)]


def main():
    import mne
    from pymatreader import read_mat
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_dir")
    ap.add_argument("--channels", nargs="+", default=None, help="default: all EEG channels")
    ap.add_argument("--k", type=int, default=4)
    ap.add_argument("--outside")
    a = ap.parse_args()
    d = Path(a.out_dir).expanduser()
    m = read_mat(str(d / "04_cwl_corrected.set"), variable_names=["etc", "EEG"])
    m = m.get("EEG", m)
    t0, t1 = (int(v) for v in m["etc"]["cwl_pipeline"]["trimmed_to"])
    r3 = mne.io.read_raw_eeglab(d / "03_filtered.set", preload=True, verbose="error")
    r4 = mne.io.read_raw_eeglab(d / "04_cwl_corrected.set", preload=True, verbose="error")
    fs = r4.info["sfreq"]
    x3 = r3.get_data()[:, t0 - 1:t1] * 1e6
    x4 = r4.get_data() * 1e6
    chans = [c for c in r4.ch_names if not EXCLUDE.match(c) and c in r3.ch_names]
    import warnings
    with warnings.catch_warnings():                 # name deprecated in newer MNE, kept for older ones
        warnings.simplefilter("ignore")
        pos = mne.channels.make_standard_montage("standard_1005").get_positions()["ch_pos"]
    chans = [c for c in chans if c in pos]
    P = np.array([pos[c] for c in chans])
    i3 = [r3.ch_names.index(c) for c in chans]; i4 = [r4.ch_names.index(c) for c in chans]
    removed = x3[i3] - x4[i4]
    clean = x4[i4]
    xo = None
    if a.outside:
        import tempfile
        from scipy.signal import resample_poly
        from cwregression import bandpass
        from run_pipeline import resolve_header
        ro = mne.io.read_raw_brainvision(resolve_header(Path(a.outside).expanduser(),
                                         tempfile.mkdtemp(prefix="rm_")), preload=True, verbose="error")
        fo = ro.info["sfreq"]
        g = int(np.gcd(int(fo), int(fs)))
        xo_all = resample_poly(bandpass(ro.get_data() * 1e6, fo, 0.33, 125.0), int(fs) // g, int(fo) // g, axis=1)
        xo = np.array([xo_all[ro.ch_names.index(c)] for c in chans])
    targets = a.channels or chans
    res = {}
    for c in targets:
        if c not in chans:
            print(f"{c}: not available"); continue
        j = chans.index(c)
        dist = np.linalg.norm(P - P[j], axis=1)
        nb = [i for i in np.argsort(dist) if i != j][:a.k]
        coh = []
        for i in nb:
            f, C = coherence(removed[j], clean[i], fs=fs, nperseg=int(4 * fs))
            coh.append(C)
        coh = np.mean(coh, 0)
        # reference level: coherence of the channel's own CLEANED signal with the same
        # neighbours (what genuine EEG looks like)
        ref = np.mean([coherence(clean[j], clean[i], fs=fs, nperseg=int(4 * fs))[1] for i in nb], 0)
        out_coh = {}
        if xo is not None:
            co = np.mean([coherence(xo[j], xo[i], fs=fs, nperseg=int(4 * fs))[1] for i in nb], 0)
            out_coh = {f"coh_outside_{lo}-{hi}": float(co[(f >= lo) & (f < hi)].mean()) for lo, hi in BANDS}
        res[c] = dict(neighbours=[chans[i] for i in nb], **out_coh,
                      removed_rms_uV=float(removed[j].std()), clean_rms_uV=float(clean[j].std()),
                      **{f"coh_removed_{lo}-{hi}": float(coh[(f >= lo) & (f < hi)].mean()) for lo, hi in BANDS},
                      **{f"coh_clean_{lo}-{hi}": float(ref[(f >= lo) & (f < hi)].mean()) for lo, hi in BANDS})
    order = sorted(res, key=lambda c: -res[c]["coh_removed_1-4"])
    print(f"{'chan':5s} {'removed':>8s} {'clean':>7s} | coherence of REMOVED part with cleaned neighbours"
          f" | (cleaned signal itself, for reference)")
    print(f"{'':5s} {'rms µV':>8s} {'rms µV':>7s} | " + " ".join(f"{lo}-{hi}Hz".rjust(7) for lo, hi in BANDS)
          + " | " + " ".join(f"{lo}-{hi}Hz".rjust(7) for lo, hi in BANDS)
          + (" | outside: " + " ".join(f"{lo}-{hi}Hz".rjust(7) for lo, hi in BANDS) if xo is not None else ""))
    for c in order:
        r = res[c]
        print(f"{c:5s} {r['removed_rms_uV']:8.2f} {r['clean_rms_uV']:7.2f} | "
              + " ".join(f"{r[f'coh_removed_{lo}-{hi}']:7.2f}" for lo, hi in BANDS) + " | "
              + " ".join(f"{r[f'coh_clean_{lo}-{hi}']:7.2f}" for lo, hi in BANDS)
              + (" |          " + " ".join(f"{r[f'coh_outside_{lo}-{hi}']:7.2f}" for lo, hi in BANDS)
                 if xo is not None else ""))
    (d / "removed_signal.json").write_text(json.dumps(res, indent=2))
    print("wrote", d / "removed_signal.json")


if __name__ == "__main__":
    main()
