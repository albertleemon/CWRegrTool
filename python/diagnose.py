"""Diagnostics for a run_pipeline.py output folder: before/after CWL spectra, an
example trace, variance removed per channel, clipping report. Writes
diagnostics.png and diagnostics.json -- send these two files back for review.

    python diagnose.py out/
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.signal import welch


def main(out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out_dir = Path(out_dir)
    z = np.load(out_dir / "cwl_corrected.npz", allow_pickle=True)
    log = json.loads((out_dir / "pipeline_log.json").read_text())
    fs = float(z["fs"]); names = list(z["ch_names"])
    eeg_idx = z["eeg_idx"]; cwl_idx = z["cwl_idx"]
    pre = z["pre_cwl"].astype(float)            # EEG after gradient correction + filter
    post = z["data"][eeg_idx]
    cwl = z["data"][cwl_idx]
    enames = [names[i] for i in eeg_idx]

    f, p_pre = welch(pre, fs, nperseg=int(4 * fs), axis=-1)
    _, p_post = welch(post, fs, nperseg=int(4 * fs), axis=-1)
    _, p_cwl = welch(cwl, fs, nperseg=int(4 * fs), axis=-1)
    var_removed = 100 * (1 - post.var(1) / pre.var(1))
    band = lambda p, lo, hi: float(p[:, (f >= lo) & (f <= hi)].mean())

    fig, ax = plt.subplots(2, 2, figsize=(13, 8))
    ax[0, 0].semilogy(f, p_pre.mean(0), label="before CWL")
    ax[0, 0].semilogy(f, p_post.mean(0), label="after CWL")
    ax[0, 0].set(xlim=(0, 125), xlabel="Hz", ylabel="power (uV^2/Hz)", title="EEG mean spectrum")
    ax[0, 0].legend()
    ax[0, 1].semilogy(f, p_cwl.T, lw=.7)
    ax[0, 1].set(xlim=(0, 125), xlabel="Hz", title="CWL spectra")
    ch = enames.index("Oz") if "Oz" in enames else 0
    t0 = int(60 * fs); n = int(10 * fs); t = np.arange(n) / fs
    ax[1, 0].plot(t, pre[ch, t0:t0 + n], lw=.7, label="before")
    ax[1, 0].plot(t, post[ch, t0:t0 + n], lw=.7, label="after")
    ax[1, 0].set(xlabel="s", ylabel="uV", title=f"{enames[ch]}, 10 s example"); ax[1, 0].legend()
    ax[1, 1].bar(range(len(enames)), var_removed)
    ax[1, 1].set(xlabel="EEG channel #", ylabel="% variance removed", title="variance removed by CWL")
    fig.tight_layout(); fig.savefig(out_dir / "diagnostics.png", dpi=110)

    summary = dict(
        mean_var_removed_pct=float(var_removed.mean()),
        power_ratio_after_before={k: band(p_post, lo, hi) / band(p_pre, lo, hi)
                                  for k, (lo, hi) in dict(delta_theta=(1, 7), alpha=(8, 13),
                                                         beta=(14, 30), gamma=(30, 80)).items()},
        rms_before_uV=float(np.sqrt((pre ** 2).mean())), rms_after_uV=float(np.sqrt((post ** 2).mean())),
        clipped_channels_fraction=log.get("clipped_channels_fraction"),
        n_volumes=log["n_volumes"], method=log["method"], samples_out=log["samples_out"])
    (out_dir / "diagnostics.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main(sys.argv[1])
