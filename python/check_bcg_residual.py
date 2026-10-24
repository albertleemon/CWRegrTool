"""How much ballistocardiogram (BCG) is left after CWL regression?

R peaks are detected in the ECG channel; the EEG is averaged time-locked to them
(-0.2 .. +0.8 s). In 03_filtered (before CWL) this average is the BCG; in
04_cwl_corrected it is the residual BCG. Brain activity is not locked to the heart beat,
so it averages out (except a small heartbeat-evoked potential, a few µV at most).
The outside-the-scanner recording, if given, provides that floor: no BCG exists there.

Per stage and channel: RMS of the heartbeat-locked average over 0-0.6 s after the R peak.

    python check_bcg_residual.py OUT_DIR [--outside OUTSIDE.vhdr] [--ecg ECG]
"""
import argparse
import json
import re
import tempfile
from pathlib import Path

import numpy as np
from scipy.signal import butter, find_peaks, resample_poly, sosfiltfilt

from cwregression import bandpass

EXCLUDE = re.compile(r"^(CW.*|ECG.*|EOG.*|65|66)$", re.I)


def r_peaks(ecg, fs):
    """R-peak sample indices (band-pass 5-30 Hz, polarity from the largest deflections)."""
    y = sosfiltfilt(butter(3, [5, 30], "band", fs=fs, output="sos"), ecg)
    if np.percentile(-y, 99.9) > np.percentile(y, 99.9):
        y = -y
    pk, _ = find_peaks(y, height=0.5 * np.percentile(y, 99), distance=int(0.33 * fs))
    return pk


def locked_average(x, pk, fs, pre=0.2, post=0.8):
    a, b = int(pre * fs), int(post * fs)
    pk = pk[(pk - a >= 0) & (pk + b <= x.shape[1])]
    E = np.stack([x[:, p - a:p + b] for p in pk], 1)         # chan x beats x time
    return E.mean(1), len(pk), np.arange(-a, b) / fs


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import mne
    from pymatreader import read_mat
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_dir"); ap.add_argument("--outside"); ap.add_argument("--ecg", default="ECG")
    a = ap.parse_args()
    d = Path(a.out_dir).expanduser()
    m = read_mat(str(d / "04_cwl_corrected.set"), variable_names=["etc", "EEG"])
    m = m.get("EEG", m)
    t0, t1 = (int(v) for v in m["etc"]["cwl_pipeline"]["trimmed_to"])

    def load_set(fn, trim):
        raw = mne.io.read_raw_eeglab(d / fn, preload=True, verbose="error")
        x = raw.get_data() * 1e6
        return raw.ch_names, (x[:, t0 - 1:t1] if trim else x), raw.info["sfreq"]

    stages = {"before CWL (03)": load_set("03_filtered.set", True),
              "after CWL (04)": load_set("04_cwl_corrected.set", False)}
    if a.outside:
        from run_pipeline import resolve_header
        raw = mne.io.read_raw_brainvision(resolve_header(Path(a.outside).expanduser(),
                                          tempfile.mkdtemp(prefix="bcg_")), preload=True, verbose="error")
        fs = raw.info["sfreq"]
        x = bandpass(raw.get_data() * 1e6, fs, 0.33, 125.0)
        g = int(np.gcd(int(fs), 500))
        stages["outside the scanner"] = (raw.ch_names, resample_poly(x, 500 // g, int(fs) // g, axis=1), 500.0)

    names04 = stages["after CWL (04)"][0]
    chans = [c for c in names04 if not EXCLUDE.match(c) and all(c in s[0] for s in stages.values())]
    res, avg = {}, {}
    for name, (n, x, fs) in stages.items():
        ecg = x[n.index(a.ecg)]
        pk = r_peaks(ecg, fs)
        A, nb, t = locked_average(x[[n.index(c) for c in chans]], pk, fs)
        w = (t >= 0) & (t <= 0.6)
        rms = np.sqrt((A[:, w] ** 2).mean(1))
        avg[name] = (t, A)
        res[name] = dict(n_beats=nb, heart_rate_bpm=float(60 * fs / np.median(np.diff(pk))),
                         locked_rms_uV_median=float(np.median(rms)), locked_rms_uV_max=float(rms.max()),
                         per_channel=dict(zip(chans, map(float, rms))))
        print(f"{name:22s} beats {nb:4d}, HR {res[name]['heart_rate_bpm']:.0f}/min, "
              f"heartbeat-locked RMS median {np.median(rms):.2f} uV, max {rms.max():.2f} uV "
              f"({chans[int(np.argmax(rms))]})")
    if "before CWL (03)" in res:
        r = np.array(list(res["after CWL (04)"]["per_channel"].values())) / \
            np.array(list(res["before CWL (03)"]["per_channel"].values()))
        res["residual_fraction_median"] = float(np.median(r))
        print(f"residual BCG after CWL: median {100 * np.median(r):.1f} % of the BCG amplitude "
              f"(range {100 * r.min():.1f}-{100 * r.max():.1f} %)")
    (d / "bcg_residual.json").write_text(json.dumps(res, indent=2))

    fig, ax = plt.subplots(1, len(avg) + 1, figsize=(5 * (len(avg) + 1), 4), layout="constrained")
    lim = max(np.abs(A).max() for _, A in avg.values())
    for i, (name, (t, A)) in enumerate(avg.items()):
        ax[i].plot(t, A.T, lw=.5, color="0.4"); ax[i].plot(t, A.mean(0), "r", lw=1.2)
        ax[i].set(title=f"{name}\n(heartbeat-locked average)", xlabel="s from R peak", ylabel="µV")
        ax[i].axvline(0, color="k", lw=.5)
    ax[1].set_ylim(ax[0].get_ylim())
    names = list(avg)
    x = np.arange(len(chans))
    for k, name in enumerate(names):
        ax[-1].bar(x + k * 0.8 / len(names), list(res[name]["per_channel"].values()), 0.8 / len(names), label=name)
    ax[-1].set_xticks(x, chans, rotation=90, fontsize=5)
    ax[-1].set(ylabel="RMS 0-0.6 s (µV)", title="heartbeat-locked RMS per channel"); ax[-1].legend(fontsize=7)
    fig.savefig(d / "bcg_residual.png", dpi=110)
    print("wrote", d / "bcg_residual.png", "and bcg_residual.json")


if __name__ == "__main__":
    main()
