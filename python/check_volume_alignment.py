"""Why does the gradient-artifact residual survive in some volumes?

Works on the RAW BrainVision recording (5000 Hz, before any correction). For a few
channels it cuts one epoch per volume (R128 marker, one TR long), builds the same
moving-average template as the pipeline (25 nearest volumes), and per volume reports
  - residual high-frequency power (60-250 Hz) after subtracting the template,
  - the time shift of the volume against its template (sub-sample, from the
    cross-correlation peak): timing jitter makes the residual grow with frequency,
  - clipping: samples at the amplifier limit (INT_16 x resolution, +-16383.5 uV here).
Then it tells which of the two explains the bursts of residual power.

    python check_volume_alignment.py sub-X_eeg.vhdr out_dir [--channels T8 C6 FT10 Cz CW1]
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np
from scipy.signal import butter, sosfiltfilt


def frac_lag(a, b, maxlag=5):
    """Lag of a relative to b in samples (parabolic peak of the cross-correlation)."""
    n = len(a)
    c = np.fft.irfft(np.fft.rfft(a, 2 * n) * np.conj(np.fft.rfft(b, 2 * n)))
    c = np.r_[c[-maxlag:], c[:maxlag + 1]]
    k = int(np.argmax(c))
    if 0 < k < len(c) - 1:
        y0, y1, y2 = c[k - 1], c[k], c[k + 1]
        d = 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
    else:
        d = 0.0
    return k - maxlag + d


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import mne
    from run_pipeline import resolve_header, volume_onsets
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("vhdr"); ap.add_argument("out_dir")
    ap.add_argument("--channels", nargs="+",
                    default=["T8", "C6", "FT10", "TP10", "C2", "FC1", "Cz", "CPz", "CW1", "CW3"])
    ap.add_argument("--marker", default="R128")
    ap.add_argument("--k", type=int, default=25)
    a = ap.parse_args()
    out = Path(a.out_dir).expanduser(); out.mkdir(parents=True, exist_ok=True)

    vhdr = resolve_header(Path(a.vhdr).expanduser(), out)
    h = Path(vhdr).read_text(encoding="utf-8", errors="replace")
    raw = mne.io.read_raw_brainvision(vhdr, preload=False, verbose="error")
    chans = [c for c in a.channels if c in raw.ch_names]
    raw.pick(chans).load_data(verbose="error")
    fs = raw.info["sfreq"]
    x = raw.get_data() * 1e6
    on = volume_onsets(raw, a.marker)
    tr = int(np.median(np.diff(on)))
    on = on[on + tr <= x.shape[1]]
    nv = len(on)
    # amplifier limit from the header resolution (INT_16)
    res = {}
    for c in chans:
        m = re.search(rf"^Ch\d+={re.escape(c)},[^,]*,([\d.]+)", h, re.M)
        res[c] = float(m.group(1)) * 32767 if m else np.inf
    print(f"{nv} volumes, TR {tr} samples; channels {chans}")

    sos = butter(4, 60, "high", fs=fs, output="sos")
    idx = on[:, None] + np.arange(tr)[None, :]
    report, P, LAG, CLIP = {}, [], [], []
    for ci, c in enumerate(chans):
        E = x[ci][idx]                                            # vol x tr
        clip = (np.abs(E) >= 0.999 * res[c]).sum(1)
        T = np.empty_like(E)
        for i in range(nv):
            lo = int(np.clip(i - a.k // 2, 0, nv - a.k))
            T[i] = E[lo:lo + a.k].mean(0)
        R = E - T
        Rh = sosfiltfilt(sos, R, axis=1)
        p = (Rh ** 2).mean(1)
        Eh, Th = sosfiltfilt(sos, E, axis=1), sosfiltfilt(sos, T, axis=1)
        lag = np.array([frac_lag(Eh[i], Th[i]) for i in range(nv)])
        P.append(p); LAG.append(lag); CLIP.append(clip)
        top = np.argsort(p)[::-1][:10]
        report[c] = dict(
            amplifier_limit_uV=res[c], max_abs_uV=float(np.abs(E).max()),
            volumes_with_clipping=int((clip > 0).sum()),
            residual_power_cv=float(p.std() / p.mean()),
            lag_sd_samples=float(lag.std()), max_abs_lag_samples=float(np.abs(lag).max()),
            corr_residual_vs_abs_lag=float(np.corrcoef(p, np.abs(lag))[0, 1]),
            corr_residual_vs_clipping=(float(np.corrcoef(p, clip)[0, 1]) if clip.any() else None),
            top_volumes=[dict(volume=int(i + 1), residual_rel=float(p[i] / np.median(p)),
                              lag_samples=round(float(lag[i]), 3), clipped_samples=int(clip[i]))
                         for i in top])
        r = report[c]
        print(f"{c:5s} max|x| {r['max_abs_uV']:8.0f} uV (limit {r['amplifier_limit_uV']:.0f}), "
              f"clipped volumes {r['volumes_with_clipping']:3d}, lag sd {r['lag_sd_samples']:.3f} smp, "
              f"corr(resid,|lag|) {r['corr_residual_vs_abs_lag']:+.2f}, "
              f"corr(resid,clip) {r['corr_residual_vs_clipping'] if r['corr_residual_vs_clipping'] is not None else 'n/a'}")
    P, LAG, CLIP = map(np.array, (P, LAG, CLIP))
    (out / "volume_alignment.json").write_text(json.dumps(report, indent=2))

    fig, ax = plt.subplots(3, 1, figsize=(13, 8), sharex=True, layout="constrained")
    v = np.arange(1, nv + 1)
    for ci, c in enumerate(chans):
        ax[0].semilogy(v, P[ci] / np.median(P[ci]), lw=.7, label=c)
        ax[1].plot(v, LAG[ci], lw=.7)
        ax[2].plot(v, CLIP[ci], lw=.7)
    ax[0].set(ylabel="residual 60-250 Hz\n(x median)", title="per-volume residual after 25-volume template")
    ax[0].legend(fontsize=7, ncol=5)
    ax[1].set(ylabel="shift vs template\n(samples at %g Hz)" % fs)
    ax[2].set(ylabel="clipped samples", xlabel="volume")
    fig.savefig(out / "volume_alignment.png", dpi=110)
    print("wrote", out / "volume_alignment.png", "and volume_alignment.json")


if __name__ == "__main__":
    main()
