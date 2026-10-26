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
from scipy.signal import resample_poly

from bcg_aas import r_peaks
from cwregression import bandpass
from run_manifest import has_stage

EXCLUDE = re.compile(r"^(CW.*|ECG.*|EOG.*|65|66)$", re.I)


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
    ap.add_argument("--hr", type=float, help="expected heart rate inside (bpm), e.g. from the Siemens .puls log")
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
    if has_stage(d, "05_bcg_aas"):
        stages["CWL + BCG AAS (05)"] = load_set("05_bcg_aas.set", False)
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
        pk = r_peaks(ecg, fs, a.hr if name != "outside the scanner" else None)
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
        if "CWL + BCG AAS (05)" in res:
            r5 = np.array(list(res["CWL + BCG AAS (05)"]["per_channel"].values())) / \
                np.array(list(res["before CWL (03)"]["per_channel"].values()))
            res["residual_fraction_after_aas_median"] = float(np.median(r5))
            print(f"residual BCG after CWL + AAS: median {100 * np.median(r5):.1f} % of the BCG amplitude")
    (d / "bcg_residual.json").write_text(json.dumps(res, indent=2))

    fig, ax = plt.subplots(1, len(avg) + 1, figsize=(5 * (len(avg) + 1), 4), layout="constrained")
    lim = max(np.abs(A).max() for _, A in avg.values())
    for i, (name, (t, A)) in enumerate(avg.items()):
        ax[i].plot(t, A.T, lw=.5, color="0.4"); ax[i].plot(t, A.mean(0), "r", lw=1.2)
        ax[i].set(title=f"{name}\n(heartbeat-locked average)", xlabel="s from R peak", ylabel="µV")
        ax[i].axvline(0, color="k", lw=.5)
    for i in range(1, len(avg)):
        ax[i].set_ylim(ax[0].get_ylim())
    names = list(avg)
    x = np.arange(len(chans))
    for k, name in enumerate(names):
        ax[-1].bar(x + k * 0.8 / len(names), list(res[name]["per_channel"].values()), 0.8 / len(names), label=name)
    ax[-1].set_xticks(x, chans, rotation=90, fontsize=5)
    ax[-1].set(ylabel="RMS 0-0.6 s (µV)", title="heartbeat-locked RMS per channel"); ax[-1].legend(fontsize=7)
    fig.savefig(d / "bcg_residual.png", dpi=110)
    # scalp maps (standard 10-05 positions; channels without a position are skipped)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mont = mne.channels.make_standard_montage("standard_1005")
    tc = [c for c in chans if c in mont.ch_names]
    if len(tc) >= 16:
        info = mne.create_info(tc, 500.0, "eeg"); info.set_montage(mont)
        fig, ax = plt.subplots(1, len(res_names := [n for n in names]) + 1, figsize=(4.3 * (len(names) + 1), 4.3),
                               layout="constrained")
        vmax = max(res[n]["per_channel"][c] for n in names[1:] for c in tc) if len(names) > 1 else None
        for i, n in enumerate(names):
            v = np.array([res[n]["per_channel"][c] for c in tc])
            im, _ = mne.viz.plot_topomap(v, info, axes=ax[i], show=False, cmap="Reds",
                                         vlim=(0, v.max() if i == 0 else vmax))
            ax[i].set_title(f"{n}\nmedian {np.median(v):.1f} µV", fontsize=10)
            plt.colorbar(im, ax=ax[i], shrink=.7, label="µV")
        fr = np.array([res[names[1]]["per_channel"][c] / res[names[0]]["per_channel"][c] for c in tc])
        im, _ = mne.viz.plot_topomap(fr, info, axes=ax[-1], show=False, cmap="viridis", vlim=(0, 0.8))
        ax[-1].set_title("residual fraction (04 / 03)", fontsize=10)
        plt.colorbar(im, ax=ax[-1], shrink=.7)
        fig.savefig(d / "bcg_topomap.png", dpi=110)
    print("wrote", d / "bcg_residual.png", ", bcg_topomap.png and bcg_residual.json")


if __name__ == "__main__":
    main()
