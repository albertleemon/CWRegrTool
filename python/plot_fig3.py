"""Figure in the style of van der Meer et al. (2016) Fig. 3: channel-averaged EEG power
spectra (dB, 0-125 Hz) for each artifact-correction method, each overlaid on the
spectrum recorded OUTSIDE the scanner (black), with the relative distance
diff(%) = ||S_in - S_out|| / ||S_out|| * 100  (eq. 4, linear power, 0.5-125 Hz).

Datasets are given as  "Column/Row=path". Columns are conditions (e.g. "He off",
"He on"), rows are methods. Paths may be EEGLAB .set (pipeline output) or BrainVision
.vhdr. The outside recording (optional) is band-passed 0.33-125 Hz and resampled to
500 Hz like the pipeline does.

Example (MATLAB pipeline output of one run, no outside recording yet):
    python plot_fig3.py fig3.png \
        --data "Inside/No BCG correction=~/cwl_test/matlab_out/03_filtered.set" \
        --data "Inside/CWL=~/cwl_test/matlab_out/04_cwl_corrected.set" \
        --outside ~/BIDS/sub-BOR5A/eeg/sub-BOR5A_task-outside_eeg.vhdr
"""
import argparse
import re
from pathlib import Path

import numpy as np
from scipy.signal import resample_poly, welch

from cwregression import bandpass

COLORS = {"cwl": "#2ca02c", "obs-ica": "#1f3fbf", "obs": "#d62bd6", "aas": "#17becf",
          "no bcg": "#d62728"}
EXCLUDE = re.compile(r"^(CW.*|ECG.*|EOG.*|65|66)$", re.I)


def trimmed_range(set_path):
    """EEG.etc.cwl_pipeline.trimmed_to from a pipeline .set, else None."""
    try:
        from pymatreader import read_mat
        m = read_mat(str(set_path))
        m = m.get("EEG", m)
        t = m["etc"]["cwl_pipeline"]["trimmed_to"]
        return int(t[0]), int(t[1])
    except Exception:
        return None


def load(path, fs_out=500.0, raw_input=False):
    """-> (channel names, data in uV at fs_out)."""
    import mne
    path = Path(path).expanduser()
    if path.suffix == ".vhdr":
        from run_pipeline import resolve_header
        raw = mne.io.read_raw_brainvision(resolve_header(path, path.parent / "_fig3_tmp"),
                                          preload=True, verbose="error")
        raw_input = True
    else:
        raw = mne.io.read_raw_eeglab(path, preload=True, verbose="error")
    x = raw.get_data() * 1e6
    fs = raw.info["sfreq"]
    if raw_input:
        x = bandpass(x, fs, 0.33, 125.0)
    if fs != fs_out:
        g = int(np.gcd(int(fs), int(fs_out)))
        x = resample_poly(x, int(fs_out) // g, int(fs) // g, axis=1)
    return raw.ch_names, x


def spectrum(x, fs):
    f, p = welch(x, fs=fs, nperseg=int(2 * fs), axis=-1)
    return f, p.mean(0)


def color_for(row):
    r = row.lower()
    for k in ("obs-ica", "obs", "aas", "no bcg", "cwl"):
        if r.startswith(k) or k in r:
            return COLORS[k]
    return "#ff7f0e"


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_png")
    ap.add_argument("--data", action="append", default=[], help='"Column/Row=path"')
    ap.add_argument("--mat", action="append", default=[],
                    help='"Column=diagnostics_data.mat" from matlab/cwl_diagnose.m: adds the rows '
                         '"No BCG correction" (03_filtered) and "CWL" (04_cwl_corrected)')
    ap.add_argument("--outside", help="EEG recorded outside the scanner (.vhdr or .set)")
    ap.add_argument("--title", default="Power spectra")
    a = ap.parse_args()

    entries = []
    for d in a.data:
        key, path = d.split("=", 1)
        col, row = key.split("/", 1) if "/" in key else ("Inside", key)
        entries.append((col.strip(), row.strip(), path.strip()))
    mats = {}
    for m in a.mat:
        col, path = m.split("=", 1) if "=" in m else ("Inside", m)
        from scipy.io import loadmat
        mm = loadmat(str(Path(path).expanduser()), squeeze_me=True)
        mats[(col.strip(), "No BCG correction")] = (mm["f"], mm["p0"].mean(1))
        mats[(col.strip(), "CWL")] = (mm["f"], mm["p1"].mean(1))
    keys = [(e[0], e[1]) for e in entries] + list(mats)
    if not keys:
        ap.error("give --data and/or --mat")
    cols = list(dict.fromkeys(k[0] for k in keys))
    rows = list(dict.fromkeys(k[1] for k in keys))

    # common EEG channel set: intersection of all inputs minus CWL/ECG/EOG/65/66
    loaded = {}
    for col, row, path in entries:
        names, x = load(path)
        tr = trimmed_range(Path(path).expanduser()) if path.endswith(".set") else None
        loaded[(col, row)] = (names, x, tr)
    # 03_filtered is longer than 04: cut it to the samples the CWL run kept
    ref_tr = next((v[2] for v in loaded.values() if v[2] is not None), None)
    out = load(a.outside) if a.outside else None
    sets = [set(v[0]) for v in loaded.values()] + ([set(out[0])] if out else [])
    chans = sorted(set.intersection(*sets) - {c for c in set.union(*sets) if EXCLUDE.match(c)}) if sets else []
    print(f"{len(chans)} EEG channels in common")

    def pick(names, x):
        return x[[names.index(c) for c in chans]]

    spec = {}
    for k, (names, x, tr) in loaded.items():
        if tr is None and ref_tr is not None and x.shape[1] > ref_tr[1]:
            x = x[:, ref_tr[0] - 1:ref_tr[1]]
        spec[k] = spectrum(pick(names, x), 500.0)
    spec.update(mats)                              # each entry: (freqs, mean spectrum)
    if out:
        fo, so = spectrum(pick(*out), 500.0)
    db = lambda p: 10 * np.log10(p)
    # y range from 0.5-125 Hz (the high-pass leaves near-zero bins at DC)
    vals = [db(p[(fk >= 0.5) & (fk <= 125)]) for fk, p in spec.values()]
    if out:
        vals.append(db(so[(fo >= 0.5) & (fo <= 125)]))
    vals = np.concatenate(vals)
    ylim = (np.floor(vals.min() / 5) * 5 - 5, np.ceil(vals.max() / 5) * 5 + 5)

    nr = len(rows) + (1 if out else 0)
    fig, ax = plt.subplots(nr, len(cols), figsize=(4.6 * len(cols) + 1.2, 1.9 * nr + 0.6),
                           sharex=True, sharey=True, squeeze=False)
    r0 = 0
    if out:
        for j in range(len(cols)):
            ax[0, j].plot(fo, db(so), "k", lw=1)
            ax[0, j].set_title("Outside" if j == 0 else "")
        r0 = 1
    for j, col in enumerate(cols):
        ax[r0, j].set_title(col if not out else "", fontsize=11)
        if out:
            ax[0, j].set_title(f"Outside  |  {col}")
        for i, row in enumerate(rows):
            a_ = ax[r0 + i, j]
            if (col, row) not in spec:
                a_.axis("off")
                continue
            f, p = spec[(col, row)]
            if out:
                a_.plot(fo, db(so), "k", lw=1)
                band = (f >= 0.5) & (f <= 125)
                s_out = np.interp(f, fo, so)
                d = 100 * np.linalg.norm(p[band] - s_out[band]) / np.linalg.norm(s_out[band])
                a_.text(0.97, 0.88, f"diff(%) = {d:.1f}", transform=a_.transAxes, ha="right", fontsize=9)
            a_.plot(f, db(p), color=color_for(row), lw=1)
            if j == 0:
                a_.set_ylabel(f"{row}\nPower (dB)", color=color_for(row), fontsize=9)
    for a_ in ax.flat:
        a_.set_xlim(0, 125)
        a_.set_ylim(*ylim)
    for a_ in ax[-1]:
        a_.set_xlabel("Frequency (Hz)")
    fig.suptitle(a.title)
    fig.tight_layout()
    fig.savefig(a.out_png, dpi=130)
    print("wrote", a.out_png)


if __name__ == "__main__":
    main()
