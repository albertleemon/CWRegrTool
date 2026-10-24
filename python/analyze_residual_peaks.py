"""Where do the high-frequency peaks that survive CWL regression come from?

Tests on the pipeline output (03_filtered.set = before CWL, 04_cwl_corrected.set = after):

 1. Volume locking. A residual of the MRI gradient artifact repeats with every volume
    (R128 marker), a helium-pump/vibration artifact does not. The data are cut into
    volume epochs; the volume-locked part (mean over all volumes, and a local mean over
    25 neighbouring volumes as in the gradient correction) is subtracted and the
    peaks are measured again. Peaks that vanish -> scanner-locked (gradient residual).
 2. Fine structure. A volume-locked signal only has power on multiples of 1/TR. A
    high-resolution spectrum (18 s segments, 0.056 Hz) around each peak is plotted
    with the 1/TR grid.
 3. Time course and topography of the peak power (per volume, per channel), plus
    correlation with head motion if an fMRIPrep confounds.tsv is given.

    python analyze_residual_peaks.py OUT_DIR [--confounds sub-X_desc-confounds_timeseries.tsv]

OUT_DIR is a cwl_run_pipeline output folder. Writes residual_peaks.png / .json there.
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks, welch

EXCLUDE = re.compile(r"^(CW.*|ECG.*|EOG.*|65|66)$", re.I)


def load_set(path):
    import mne
    raw = mne.io.read_raw_eeglab(path, preload=True, verbose="error")
    fs = raw.info["sfreq"]
    on = np.array(sorted(int(round((a["onset"] - raw.first_time) * fs))
                         for a in raw.annotations if a["description"].strip().endswith("R128")))
    return raw.ch_names, raw.get_data() * 1e6, fs, on


def epochs(x, on, L):
    on = on[(on >= 0) & (on + L <= x.shape[1])]
    return on, np.stack([x[:, o:o + L] for o in on], 1)          # (chan, vol, L)


def subtract_locked(x, on, L, k=None):
    """Remove the volume-locked component: mean over all volumes (k=None) or over the k
    nearest volumes (moving average, like the gradient correction)."""
    on, E = epochs(x, on, L)
    y = x.copy()
    n = len(on)
    for i, o in enumerate(on):
        if k is None:
            tmpl = E.mean(1)
        else:
            lo = int(np.clip(i - k // 2, 0, n - k)); tmpl = E[:, lo:lo + k].mean(1)
        y[:, o:o + L] = x[:, o:o + L] - tmpl
    return y, on


def peak_heights(f, p, freqs):
    """dB above the local median (±1..3 Hz) at each frequency."""
    L = 10 * np.log10(p)
    out = []
    for fr in freqs:
        j = np.argmin(abs(f - fr))
        nb = (abs(f - fr) > 1) & (abs(f - fr) < 3)
        out.append(float(L[max(j - 1, 0):j + 2].max() - np.median(L[nb])))
    return out


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_dir")
    ap.add_argument("--confounds")
    ap.add_argument("--fmin", type=float, default=40, help="look for peaks above this frequency")
    a = ap.parse_args()
    d = Path(a.out_dir).expanduser()

    names, x, fs, on = load_set(d / "03_filtered.set")
    eeg = [i for i, n in enumerate(names) if not EXCLUDE.match(n)]
    cwl = [i for i, n in enumerate(names) if n.upper().startswith("CW")]
    tr = int(round(np.median(np.diff(on))))
    print(f"{len(on)} volume markers, TR = {tr} samples = {tr / fs:.3f} s, {len(eeg)} EEG channels")

    # --- 1. peaks in the mean EEG spectrum before CWL --------------------------------
    nper = int(4 * fs)
    f, p = welch(x[eeg], fs, nperseg=nper, axis=-1)
    pm = p.mean(0)
    Ld = 10 * np.log10(pm)
    base = np.convolve(Ld, np.ones(41) / 41, "same")
    idx, _ = find_peaks(Ld - base, height=3)
    idx = idx[(f[idx] >= a.fmin) & (f[idx] < fs / 2 - 2)]
    pk = f[idx]
    f0 = np.median(pk / np.round(pk / np.median(np.diff(pk)))) if len(pk) > 1 else np.nan
    print("peaks (Hz):", np.round(pk, 2), f"-> spacing ~ {f0:.3f} Hz")

    # --- 2. remove the volume-locked part -----------------------------------------
    y_all, _ = subtract_locked(x, on, tr)
    y_loc, _ = subtract_locked(x, on, tr, k=25)
    seg = slice(on[0], on[-1] + tr)          # compare on the volume-covered samples only
    spec = lambda z: welch(z[eeg, seg], fs, nperseg=nper, axis=-1)[1].mean(0)
    p_pre, p_all, p_loc = spec(x), spec(y_all), spec(y_loc)
    res = {"peak_hz": [round(float(v), 2) for v in pk],
           "harmonic_spacing_hz": round(float(f0), 3), "tr_s": tr / fs,
           "height_db_before": peak_heights(f, p_pre, pk),
           "height_db_minus_volume_locked_mean": peak_heights(f, p_all, pk),
           "height_db_minus_volume_locked_25": peak_heights(f, p_loc, pk)}
    after = d / "04_cwl_corrected.set"
    if after.exists():
        n2, x2, _, _ = load_set(after)
        p_cwl = welch(x2[[n2.index(names[i]) for i in eeg]], fs, nperseg=nper, axis=-1)[1].mean(0)
        res["height_db_after_cwl"] = peak_heights(f, p_cwl, pk)

    # --- 3. fine structure: high-resolution spectrum with the 1/TR grid --------------
    fh, ph = welch(x[eeg, seg], fs, nperseg=int(10 * tr), axis=-1)
    ph = ph.mean(0)
    # where do the fine-structure maxima sit relative to the 1/TR grid?
    grid = fs / tr
    offs = []
    for fr in pk:
        w = (fh > fr - 0.6) & (fh < fr + 0.6)
        fm = fh[w][np.argmax(ph[w])]
        offs.append(float((fm / grid - np.round(fm / grid)) * grid))
    res["fine_peak_offset_from_1/TR_grid_hz"] = [round(o, 3) for o in offs]
    res["frequency_resolution_hz"] = round(float(fh[1] - fh[0]), 4)

    # --- 4. per-volume time course and per-channel topography of the peak power -------
    _, E = epochs(x[eeg], on, tr)
    fe = np.fft.rfftfreq(tr, 1 / fs)
    P = np.abs(np.fft.rfft(E, axis=-1)) ** 2                      # chan x vol x freq
    sel = np.zeros_like(fe, bool)
    for fr in pk:
        sel |= abs(fe - fr) < 0.6
    vol_power = P[:, :, sel].sum(-1).mean(0)                     # per volume
    chan_power = P[:, :, sel].sum(-1).mean(1)                    # per channel
    order = np.argsort(chan_power)[::-1]
    res["strongest_channels"] = [names[eeg[i]] for i in order[:8]]
    res["weakest_channels"] = [names[eeg[i]] for i in order[-5:]]
    res["volume_power_cv"] = float(vol_power.std() / vol_power.mean())
    if cwl:
        _, Ec = epochs(x[cwl], on, tr)
        Pc = (np.abs(np.fft.rfft(Ec, axis=-1)) ** 2)[:, :, sel].sum(-1).mean(0)
        res["corr_eeg_vs_cwl_peak_power_over_volumes"] = float(np.corrcoef(vol_power, Pc)[0, 1])
    fd = None
    if a.confounds:
        import csv
        with open(a.confounds) as fh_:
            rows = list(csv.DictReader(fh_, delimiter="\t"))
        fd = np.array([float(r["framewise_displacement"]) if r["framewise_displacement"] not in ("n/a", "") else np.nan
                       for r in rows])
        # epochs start at the first marker kept in 03_filtered; align to the end of the run
        fd = fd[-len(vol_power):] if len(fd) >= len(vol_power) else None
        if fd is not None:
            ok = ~np.isnan(fd)
            res["corr_peak_power_vs_framewise_displacement"] = float(np.corrcoef(vol_power[ok], fd[ok])[0, 1])

    # --- figure --------------------------------------------------------------------
    fig = plt.figure(figsize=(13, 9), layout="constrained")
    ax = fig.add_subplot(3, 1, 1)
    b = (f > 20) & (f < fs / 2)
    ax.semilogy(f[b], p_pre[b], label="before CWL (03_filtered)", lw=1)
    ax.semilogy(f[b], p_loc[b], label="minus volume-locked mean (25 volumes)", lw=1)
    ax.semilogy(f[b], p_all[b], label="minus volume-locked mean (all volumes)", lw=1)
    if after.exists():
        ax.semilogy(f[b], p_cwl[b], label="after CWL (04)", lw=1)
    for fr in pk:
        ax.axvline(fr, color="0.85", lw=.6, zorder=0)
    ax.set(xlabel="Hz", ylabel="uV^2/Hz", title=f"mean EEG spectrum; peaks spaced ~{f0:.2f} Hz")
    ax.legend(fontsize=8)
    show = pk[np.argsort(peak_heights(f, p_pre, pk))[::-1][:3]] if len(pk) else []
    for i, fr in enumerate(sorted(show)):
        ax = fig.add_subplot(3, 3, 4 + i)
        w = (fh > fr - 1.5) & (fh < fr + 1.5)
        ax.semilogy(fh[w], ph[w], lw=.8)
        for g in np.arange(np.ceil((fr - 1.5) / grid), np.floor((fr + 1.5) / grid) + 1) * grid:
            ax.axvline(g, color="r", lw=.4, alpha=.5)
        ax.set(title=f"{fr:.2f} Hz, 0.056 Hz resolution\nred: multiples of 1/TR", xlabel="Hz")
    ax = fig.add_subplot(3, 2, 5)
    ax.plot(vol_power / vol_power.mean(), lw=.8, label="EEG peak power")
    if fd is not None:
        ax2 = ax.twinx(); ax2.plot(fd, color="0.6", lw=.6); ax2.set_ylabel("FD (mm)", color="0.5")
    ax.set(xlabel="volume", ylabel="relative power", title="peak power per volume")
    ax = fig.add_subplot(3, 2, 6)
    ax.bar(range(len(eeg)), 10 * np.log10(chan_power / np.median(chan_power)))
    ax.set_xticks(range(len(eeg)), [names[i] for i in eeg], rotation=90, fontsize=5)
    ax.set(ylabel="dB vs median channel", title="peak power per channel")
    fig.savefig(d / "residual_peaks.png", dpi=110)
    (d / "residual_peaks.json").write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
