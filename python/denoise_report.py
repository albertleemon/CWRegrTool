"""Denoising report for one pipeline run: indicators, mean spectrum of EVERY step, scalp maps.

Stages (all on the same samples and EEG channels, 0.33-125 Hz, 500 Hz, so the curves differ
only by the correction steps):
    raw (needs --vhdr) -> 02 gradient -> 02b OBS (if run) -> 03 filtered -> 04 CWL -> 05 BCG AAS (if run)
    and the outside-the-scanner recording (--outside), the reference without any artifact.

Indicators per stage (median over EEG channels unless noted)
    rms_1_30_uV            RMS 1-30 Hz from the spectrum
    band power             1-4, 4-8, 8-13, 13-30, 30-60, 60-90, 90-125 Hz, µV²/Hz; dB vs outside
    spectral_distance_pct  eq. 4 of van der Meer et al. to the outside spectrum (0.5-125 Hz and 1-30 Hz)
    log_spectrum_r         correlation of log10 spectra with the outside, 1-45 Hz
    slice_peak_db          excess of the slice-frequency harmonics over the local background (dB)
    bcg_locked_rms_uV      RMS 0-0.6 s of the heartbeat-locked average (03, 04, 05, outside; R peaks of
                           the inside ECG are detected once on 03 and used for all inside stages)
    bcg_residual_pct       the same, as % of 03
    step_change_db_1_30    change of 1-30 Hz power relative to the previous stage
    cwl_variance_removed   1 - var(04)/var(03), median over channels
    alpha_topo_r           spatial correlation of 8-13 Hz power with the outside
    alpha_peak_hz / _db_over_1f / alpha_power_vs_outside_db / alpha_peak_shift_hz
                           posterior-channel alpha peak above its 1/f background (fit on 4-6 and 15-30 Hz):
                           frequency (nan if the peak is < 3 dB), height, power relative to the outside (dB), frequency shift
Writes into OUT_DIR: denoise_spectra.png, denoise_alpha.png, denoise_indicators.png, denoise_topomaps.png,
denoise_indicators.json / .csv.

    python denoise_report.py OUT_DIR [--outside OUTSIDE.vhdr] [--vhdr RAW.vhdr] [--hr BPM]
                                     [--slice-freq 16.84] [--crop START END]
"""
import argparse
import csv
import json
import re
import warnings
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from scipy.signal import resample_poly, welch

from bcg_aas import r_peaks
from check_bcg_residual import locked_average
from cwregression import bandpass
from evaluate_outside import load_outside
from metrics import relative_distance
from run_manifest import has_stage

EXCLUDE = re.compile(r"^(CW.*|ECG.*|EOG.*|65|66)$", re.I)
BANDS = [(1, 4), (4, 8), (8, 13), (13, 30), (30, 60), (60, 90), (90, 125)]
POSTERIOR = ["O1", "O2", "Oz", "POz", "PO3", "PO4", "PO7", "PO8", "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "Pz"]
COLORS = {"raw (no correction)": "#7f7f7f", "02 gradient": "#d62728", "02b OBS": "#ff7f0e",
          "03 filtered": "#8c564b", "04 CWL": "#2ca02c", "05 BCG AAS": "#9467bd",
          "outside the scanner": "#000000"}


def load_any(path, t0, t1, start=0):
    """All channels of a pipeline .set (or an mne Raw) -> band-passed, 500 Hz, samples t0..t1 (1-based)."""
    import mne
    raw = path if hasattr(path, "get_data") else mne.io.read_raw_eeglab(path, preload=True, verbose="error")
    fs = raw.info["sfreq"]
    x = raw.get_data()[:, start:] * 1e6
    names = raw.ch_names
    if fs != 500:                                  # 5000 Hz stages: same filter + resampling as the pipeline
        x = bandpass(x, fs, 0.33, 125.0)
        g = int(np.gcd(int(fs), 500))
        x = resample_poly(x, 500 // g, int(fs) // g, axis=1)
    return names, x[:, t0 - 1:t1], 500.0


def slice_peak_db(f, p, f_slice):
    """Median excess (dB) of P at the harmonics of f_slice over the median of the ±3 Hz neighbourhood."""
    out = []
    for h in np.arange(f_slice, 125.0, f_slice):
        near = (np.abs(f - h) <= 0.5)
        bg = (np.abs(f - h) > 1.0) & (np.abs(f - h) <= 3.0)
        if near.any() and bg.any():
            out.append(10 * np.log10(p[near].max() / np.median(p[bg])))
    return float(np.median(out)) if out else float("nan")


def analyse(names, x, fs, chans, pk=None, ecg="ECG", f_slice=16.84):
    ci = [names.index(c) for c in chans]
    f, p = welch(x[ci], fs, nperseg=int(4 * fs), noverlap=int(2 * fs), axis=-1)       # chan x freq
    r = dict(f=f, p=p, var=x[ci].var(1))
    if pk is not None:
        A, _, t = locked_average(x[ci], pk, fs)
        r["locked_rms"] = np.sqrt((A[:, (t >= 0) & (t <= 0.6)] ** 2).mean(1))
    return r


def alpha_peak(f, p):
    """Alpha peak of one spectrum above its 1/f background. The background is a line in log-log
    coordinates fitted on the flanks 4-6 Hz and 15-30 Hz (so a broadband noise difference is not
    mistaken for an alpha difference). Returns peak frequency (7-14 Hz), peak height (dB above the
    background), the alpha power above the background (8-13 Hz, µV²) and the background curve."""
    flank = ((f >= 4) & (f < 6)) | ((f >= 15) & (f <= 30))
    k = np.polyfit(np.log10(f[flank]), np.log10(p[flank]), 1)
    bg = 10 ** np.polyval(k, np.log10(np.maximum(f, f[1])))
    ex = 10 * np.log10(p / bg)
    sel = (f >= 7) & (f <= 14)
    i = int(np.argmax(ex[sel]))
    a = (f >= 8) & (f <= 13)
    # a "peak" lower than 3 dB above the background is noise, not an alpha peak: no frequency is reported
    return dict(peak_hz=float(f[sel][i]) if ex[sel][i] >= 3 else float("nan"), peak_db=float(ex[sel][i]),
                power_above_bg=float(np.clip(p[a] - bg[a], 0, None).sum() * (f[1] - f[0])), bg=bg)


def band(f, p, lo, hi):
    return p[..., (f >= lo) & (f < hi)].mean(-1)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import mne
    from pymatreader import read_mat
    from run_pipeline import resolve_header, volume_onsets
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_dir")
    ap.add_argument("--outside", help="outside-the-scanner BrainVision header (reference)")
    ap.add_argument("--vhdr", help="raw BrainVision header of this run (adds the 'raw' stage)")
    ap.add_argument("--hr", type=float, help="expected heart rate inside (bpm)")
    ap.add_argument("--ecg", default="ECG")
    ap.add_argument("--marker", default="R128", help="volume marker (for the raw stage)")
    ap.add_argument("--slice-freq", type=float, default=16.84, help="slice frequency (Hz) for the harmonic peak indicator")
    ap.add_argument("--crop", type=float, nargs=2, metavar=("START", "END"), help="seconds of the outside recording to use")
    a = ap.parse_args()
    d = Path(a.out_dir).expanduser()
    m = read_mat(str(d / "04_cwl_corrected.set"), variable_names=["etc", "EEG"])
    m = m.get("EEG", m)
    t0, t1 = (int(v) for v in m["etc"]["cwl_pipeline"]["trimmed_to"])

    # --- load the stages one at a time (the 5000 Hz ones are large) --------------------------------
    srcs = []                                       # (name, loader)
    if a.vhdr:
        def _raw():
            from run_pipeline import resolve_header as rh
            raw = mne.io.read_raw_brainvision(rh(Path(a.vhdr).expanduser(), d), preload=True, verbose="error")
            return load_any(raw, t0, t1, start=int(volume_onsets(raw, a.marker)[0]))
        srcs.append(("raw (no correction)", _raw, False))
    for stage, label, lock in (("02_gradient_corrected", "02 gradient", False), ("02b_obs_corrected", "02b OBS", False),
                               ("03_filtered", "03 filtered", True), ("04_cwl_corrected", "04 CWL", True),
                               ("05_bcg_aas", "05 BCG AAS", True)):
        if has_stage(d, stage):
            full = stage in ("02_gradient_corrected", "02b_obs_corrected", "03_filtered")   # untrimmed stages
            srcs.append((label, (lambda s=stage, fl=full: load_any(d / f"{s}.set", t0, t1) if fl
                                 else load_any(d / f"{s}.set", 1, 10**12)), lock))

    loaded = {}
    chans = None
    for label, loader, lock in srcs:
        n, x, fs = loader()
        loaded[label] = (n, x, fs)
        e = [c for c in n if not EXCLUDE.match(c)]
        chans = e if chans is None else [c for c in chans if c in e]
    # inside R peaks, detected once on the filtered stage
    pk_in = None
    n3, x3, fs3 = loaded["03 filtered"]
    if a.ecg in n3:
        pk_in = r_peaks(x3[n3.index(a.ecg)], fs3, a.hr)
    outside = None
    if a.outside:
        outside = load_outside(Path(a.outside).expanduser(), a.crop)
        chans = [c for c in chans if c in outside[0]]
    print(f"{len(chans)} EEG channels in common; stages: {', '.join(loaded)}" + (", outside" if outside else ""))

    R = {}
    for label, (n, x, fs) in loaded.items():
        lock = label in ("03 filtered", "04 CWL", "05 BCG AAS")
        R[label] = analyse(n, x, fs, chans, pk_in if lock else None, a.ecg)
    loaded.clear()
    if outside:
        n, x, fs = outside
        pk_o = r_peaks(x[n.index(a.ecg)], fs) if a.ecg in n else None
        R["outside the scanner"] = analyse(n, x, fs, chans, pk_o, a.ecg)
    f = R["04 CWL"]["f"]
    ref = R.get("outside the scanner")
    b_all, b_130 = (f >= 0.5) & (f <= 125), (f >= 1) & (f <= 30)
    df = f[1] - f[0]

    # --- indicators ----------------------------------------------------------------------------------
    ind, prev = {}, None
    for label, r in R.items():
        p = r["p"]; pm = p.mean(0)
        row = dict(rms_1_30_uV=float(np.median(np.sqrt(p[:, b_130].sum(1) * df))),
                   band_power_uV2_per_Hz={f"{lo}-{hi}": float(band(f, pm, lo, hi)) for lo, hi in BANDS},
                   slice_peak_db=slice_peak_db(f, pm, a.slice_freq))
        if ref is not None:
            po = ref["p"].mean(0)
            row["band_diff_vs_outside_db"] = {f"{lo}-{hi}": float(10 * np.log10(band(f, pm, lo, hi) / band(f, po, lo, hi)))
                                              for lo, hi in BANDS}
            row["spectral_distance_pct"] = relative_distance(pm[b_all], po[b_all])
            row["spectral_distance_1_30Hz_pct"] = relative_distance(pm[b_130], po[b_130])
            bc = (f >= 1) & (f <= 45)
            row["log_spectrum_r"] = float(np.corrcoef(np.log10(pm[bc]), np.log10(po[bc]))[0, 1])
            al = (f >= 8) & (f <= 13)
            row["alpha_topo_r"] = float(np.corrcoef(p[:, al].mean(1), ref["p"][:, al].mean(1))[0, 1])
        pi = [i for i, c in enumerate(chans) if c in POSTERIOR]
        if len(pi) >= 3:
            ap_ = alpha_peak(f, p[pi].mean(0))
            r["alpha"] = ap_
            row.update(alpha_peak_hz=ap_["peak_hz"], alpha_peak_db_over_1f=ap_["peak_db"],
                       alpha_power_over_1f_uV2=ap_["power_above_bg"])
        if prev is not None and "outside" not in label:
            row["step_change_db_1_30"] = float(10 * np.log10(pm[b_130].sum() / prev[b_130].sum()))
        prev = pm
        if "locked_rms" in r:
            row["bcg_locked_rms_uV"] = float(np.median(r["locked_rms"]))
        ind[label] = row
    if "locked_rms" in R.get("03 filtered", {}):
        for label, r in R.items():
            if "locked_rms" in r and label != "outside the scanner":
                ind[label]["bcg_residual_pct"] = float(100 * np.median(r["locked_rms"] / R["03 filtered"]["locked_rms"]))
    if ref is not None and "alpha" in ref:
        for label, r in R.items():
            if "alpha" in r and "outside" not in label and ref["alpha"]["power_above_bg"] > 0 and r["alpha"]["power_above_bg"] > 0:
                ind[label]["alpha_power_vs_outside_db"] = float(10 * np.log10(r["alpha"]["power_above_bg"] / ref["alpha"]["power_above_bg"]))
                ind[label]["alpha_peak_shift_hz"] = r["alpha"]["peak_hz"] - ref["alpha"]["peak_hz"]   # nan if either has no peak
    ind["04 CWL"]["cwl_variance_removed"] = float(np.median(1 - R["04 CWL"]["var"] / R["03 filtered"]["var"]))
    (d / "denoise_indicators.json").write_text(json.dumps(ind, indent=2))
    cols = ["rms_1_30_uV", "spectral_distance_pct", "spectral_distance_1_30Hz_pct", "log_spectrum_r", "alpha_topo_r",
            "alpha_peak_hz", "alpha_peak_db_over_1f", "alpha_power_vs_outside_db", "alpha_peak_shift_hz", "slice_peak_db", "bcg_locked_rms_uV", "bcg_residual_pct", "step_change_db_1_30", "cwl_variance_removed"]
    with open(d / "denoise_indicators.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["stage"] + cols + [f"dB_vs_outside_{lo}-{hi}" for lo, hi in BANDS])
        for s, row in ind.items():
            w.writerow([s] + [row.get(c, "") for c in cols] +
                       [row.get("band_diff_vs_outside_db", {}).get(f"{lo}-{hi}", "") for lo, hi in BANDS])
    print(f"{'stage':22s}" + "".join(f"{c[:14]:>15s}" for c in cols))
    for s, row in ind.items():
        print(f"{s:22s}" + "".join(f"{row[c]:15.3g}" if c in row else f"{'-':>15s}" for c in cols))

    # --- figure 1: mean spectrum of every step + difference to the outside -----------------------------
    nul = None
    if (d / "diagnostics_data.mat").exists():
        mm = loadmat(d / "diagnostics_data.mat", squeeze_me=True)
        if "pn" in mm and len(mm["f"]) == len(f):
            nul = mm["pn"].mean(1)
    fig, ax = plt.subplots(1, 2 if ref is not None else 1, figsize=(15 if ref is not None else 9, 5.5),
                           layout="constrained", squeeze=False)
    ax = ax[0]
    for label, r in R.items():
        ax[0].semilogy(f[b_all], r["p"].mean(0)[b_all], color=COLORS.get(label), lw=1.6 if "outside" in label else 1.1,
                       ls="--" if "outside" in label else "-", label=label)
    if nul is not None:
        ax[0].semilogy(f[b_all], nul[b_all], color="#1f77b4", ls=":", lw=1, label="null: random-phase CWL")
    ax[0].set(xlim=(0, 125), xlabel="Frequency (Hz)", ylabel="Power (µV²/Hz)",
              title=f"Mean spectrum of every step ({len(chans)} EEG channels), {d.name}")
    ax[0].legend(fontsize=8); ax[0].grid(alpha=.3, which="both", lw=.4)
    if ref is not None:
        po = ref["p"].mean(0)
        for label, r in R.items():
            if "outside" in label:
                continue
            ax[1].plot(f[b_all], 10 * np.log10(r["p"].mean(0)[b_all] / po[b_all]), color=COLORS.get(label), lw=1, label=label)
        ax[1].axhline(0, color="k", lw=.8)
        ax[1].set(xlim=(0, 125), xlabel="Frequency (Hz)", ylabel="dB relative to outside", title="Distance to the outside spectrum (0 dB = identical)")
        ax[1].legend(fontsize=8); ax[1].grid(alpha=.3, lw=.4)
    fig.savefig(d / "denoise_spectra.png", dpi=120)

    # --- figure: posterior alpha, with and without the 1/f background -------------------------------------
    if all("alpha" in r for r in R.values()):
        pi = [i for i, c in enumerate(chans) if c in POSTERIOR]
        z = (f >= 3) & (f <= 20)
        fig, axa = plt.subplots(1, 2, figsize=(14, 5), layout="constrained")
        for label, r in R.items():
            pp = r["p"][pi].mean(0)
            kw = dict(color=COLORS.get(label), lw=1.6 if "outside" in label else 1.1, ls="--" if "outside" in label else "-")
            axa[0].semilogy(f[z], pp[z], label=label, **kw)
            axa[1].plot(f[z], 10 * np.log10(pp[z] / r["alpha"]["bg"][z]), **kw)
            axa[1].plot(r["alpha"]["peak_hz"], r["alpha"]["peak_db"], "o", color=COLORS.get(label), ms=5)
        axa[0].set(xlabel="Frequency (Hz)", ylabel="Power (µV²/Hz)", title=f"Posterior channels ({len(pi)}), 3-20 Hz")
        axa[1].set(xlabel="Frequency (Hz)", ylabel="dB above the 1/f background", title="Alpha peak above 1/f (background fitted on 4-6 and 15-30 Hz)")
        for axx in axa:
            axx.axvspan(8, 13, color="0.9", zorder=0); axx.grid(alpha=.3, lw=.4)
        axa[0].legend(fontsize=8)
        fig.savefig(d / "denoise_alpha.png", dpi=120)

    # --- figure 2: indicators per stage ---------------------------------------------------------------
    stages = list(R)
    panels = [("spectral_distance_pct", "Spectral distance to outside (%)\nlower = better"),
              ("rms_1_30_uV", "RMS 1-30 Hz (µV, median)\noutside = reference"),
              ("bcg_locked_rms_uV", "Heartbeat-locked RMS (µV)\nlower = better"),
              ("slice_peak_db", "Slice-harmonic peaks (dB over background)\n0 = none"),
              ("log_spectrum_r", "Log-spectrum correlation with outside (1-45 Hz)\nhigher = better"),
              ("alpha_topo_r", "Alpha (8-13 Hz) topography correlation\nhigher = better"),
              ("alpha_power_vs_outside_db", "Posterior alpha power above 1/f vs outside (dB)\n0 = same as outside")]
    panels = [(k, t) for k, t in panels if any(k in ind[s] for s in stages)]
    fig, axs = plt.subplots(3, 3, figsize=(15, 11), layout="constrained")
    for axp, (k, t) in zip(axs.ravel(), panels):
        vals = [ind[s].get(k, np.nan) for s in stages]
        if k in ("spectral_distance_pct", "rms_1_30_uV"):
            vals = [v if v > 0 else np.nan for v in vals]       # log axis: 0 (the reference itself) cannot be drawn
        axp.bar(range(len(stages)), vals, color=[COLORS.get(s) for s in stages])
        axp.set_xticks(range(len(stages)), [s.replace(" (no correction)", "").replace("outside the scanner", "outside") for s in stages],
                       rotation=40, ha="right", fontsize=8)
        axp.set_title(t, fontsize=9)
        if k in ("spectral_distance_pct", "rms_1_30_uV"):
            axp.set_yscale("log")                    # the uncorrected stages are orders of magnitude larger
        for i, v in enumerate(vals):
            if np.isfinite(v):
                axp.text(i, v, f"{v:.3g}", ha="center", va="bottom", fontsize=7)
    for axp in axs.ravel()[len(panels):]:
        axp.axis("off")
    fig.savefig(d / "denoise_indicators.png", dpi=120)

    # --- figure 3: scalp maps -------------------------------------------------------------------------
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mont = mne.channels.make_standard_montage("standard_1005")
    tc = [c for c in chans if c in mont.ch_names]
    if len(tc) >= 16:
        ti = [chans.index(c) for c in tc]
        info = mne.create_info(tc, 500.0, "eeg"); info.set_montage(mont)
        rows = [("1-30 Hz RMS (µV)", lambda r: np.sqrt(r["p"][ti][:, b_130].sum(1) * df), "viridis", False)]
        if any("locked_rms" in r for r in R.values()):
            rows.append(("heartbeat-locked RMS (µV)", lambda r: r["locked_rms"][ti] if "locked_rms" in r else None, "Reds", False))
        if ref is not None:
            rows.append(("1-30 Hz power vs outside (dB)", lambda r: 10 * np.log10(r["p"][ti][:, b_130].sum(1) / ref["p"][ti][:, b_130].sum(1)),
                         "RdBu_r", True))
        fig, axm = plt.subplots(len(rows), len(R), figsize=(3.2 * len(R), 3.2 * len(rows)), layout="constrained", squeeze=False)
        for i, (title, fn, cmap, div) in enumerate(rows):
            vs = {s: fn(r) for s, r in R.items()}
            have = [v for v in vs.values() if v is not None]
            if div:
                lim = max(np.abs(v).max() for v in have); vlim = (-lim, lim)
            else:
                # shared scale; the uncorrected stages would flatten it, so cap at the 95th percentile of the cleaned ones
                ref_vals = [v for s, v in vs.items() if v is not None and s in ("04 CWL", "05 BCG AAS", "outside the scanner")] or have
                vlim = (0, float(np.percentile(np.concatenate(ref_vals), 99)) * 1.5)
            for j, s in enumerate(R):
                v = vs[s]
                if v is None:
                    axm[i, j].axis("off"); continue
                im, _ = mne.viz.plot_topomap(v, info, axes=axm[i, j], show=False, cmap=cmap, vlim=vlim, contours=0)
                axm[i, j].set_title(f"{s}\n{title if j == 0 else ''}\nmedian {np.median(v):.2f}", fontsize=8)
            plt.colorbar(im, ax=axm[i, :].tolist(), shrink=.7)
        fig.savefig(d / "denoise_topomaps.png", dpi=110)
        print("wrote denoise_topomaps.png")
    else:
        print(f"only {len(tc)} channels with a standard position: scalp maps skipped")
    print("wrote", d / "denoise_spectra.png", ",", "denoise_alpha.png, denoise_indicators.png/.json/.csv")


if __name__ == "__main__":
    main()
