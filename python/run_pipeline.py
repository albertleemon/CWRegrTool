"""End-to-end offline pipeline for simultaneous EEG-fMRI with carbon-wire loops.

    BrainVision (.vhdr) -> gradient-artifact correction (all channels, incl. CWL)
    -> band-pass 0.33-125 Hz (all channels) -> resample to 500 Hz
    -> CWL windowed regression (EEG channels only) -> corrected EEG (.fif / .npz)

Example (the NB recording of this project):
    python run_pipeline.py NB.vhdr out/ --rp rp_run01.txt \
        --exclude 65 66 --ecg ECG --cwl-regex '^CW'
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np
from scipy.signal import resample_poly

from cwregression import bandpass, cw_regression
from gradient_correction import (gradient_correct, moving_average_matrix,
                                 rp_weighting_matrix)


def volume_onsets(raw, marker):
    """Sample indices of all annotations whose description ends with `marker`
    (BrainVision 'Response,R128' is read by MNE as 'Response/R128')."""
    ann = raw.annotations
    on = [int(round((a["onset"] - raw.first_time) * raw.info["sfreq"]))
          for a in ann if a["description"].split("/")[-1] == marker]
    return np.array(sorted(on), int)


def run(vhdr, out_dir, rp=None, marker="R128", cwl_regex=r"^CW", ecg="ECG",
        exclude=(), n_templates=25, threshold=0.3, span=None,
        band=(0.33, 125.0), target_fs=500, window_s=4.0, delay_s=0.021):
    import mne
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    raw = mne.io.read_raw_brainvision(vhdr, preload=True, verbose="error")
    fs = raw.info["sfreq"]
    names = raw.ch_names
    data = raw.get_data() * 1e6                          # volts -> microvolts

    cwl_idx = [i for i, n in enumerate(names) if re.match(cwl_regex, n)]
    ecg_idx = [i for i, n in enumerate(names) if n == ecg]
    eeg_idx = [i for i, n in enumerate(names)
               if i not in cwl_idx + ecg_idx and n not in exclude]
    if not cwl_idx:
        raise SystemExit(f"no CWL channels match {cwl_regex!r}: {names}")

    onsets = volume_onsets(raw, marker)
    if len(onsets) < 3:
        raise SystemExit(f"found {len(onsets)} '{marker}' markers; need the volume markers")
    tr = int(np.median(np.diff(onsets)))
    pre, post = span if span else (0, tr - 1)            # artifact span relative to marker

    # 1. gradient artifact: every channel, CWL included
    if rp is not None:
        W = rp_weighting_matrix(np.loadtxt(rp), len(onsets), threshold, n_templates)
        method = "rp_info"
    else:
        W = moving_average_matrix(len(onsets), n_templates)
        method = "moving_average"
    if onsets[-1] + post >= data.shape[1]:
        raise SystemExit("last artifact exceeds the recording: drop the last marker or trim")
    data = gradient_correct(data, fs, onsets, pre, post, W)

    # data outside the corrected span still contains gradient artifacts -> cut it
    a, b = onsets[0] + pre, onsets[-1] + post + 1
    data = data[:, a:b]

    # 2. identical zero-phase band-pass on every channel, then resample
    data = bandpass(data, fs, *band)
    if target_fs != fs:
        g = int(np.gcd(int(fs), int(target_fs)))
        data = resample_poly(data, int(target_fs) // g, int(fs) // g, axis=1)
        fs = float(target_fs)
    if np.isnan(data).any():
        raise SystemExit("NaN in data before CWL regression")

    # 3. CWL regression on the EEG channels only
    eeg, art = cw_regression(data[eeg_idx], data[cwl_idx], fs, window_s=window_s,
                             delay_s=delay_s, filter_band=None, return_artifact=True)

    out = data.copy()
    out[eeg_idx] = eeg
    info = mne.create_info(names, fs, ch_types="eeg")
    mne.io.RawArray(out * 1e-6, info, verbose="error").save(
        out_dir / "cwl_corrected_raw.fif", overwrite=True)
    np.savez_compressed(out_dir / "cwl_corrected.npz", data=out, ch_names=names, fs=fs,
                        eeg_idx=eeg_idx, cwl_idx=cwl_idx, artifact=art)
    (out_dir / "pipeline_log.json").write_text(json.dumps(dict(
        vhdr=str(vhdr), method=method, n_volumes=len(onsets), tr_samples=tr,
        span=[pre, post], n_templates=n_templates, threshold_mm=threshold, band=band,
        fs_out=fs, window_s=window_s, delay_s=delay_s,
        eeg_channels=[names[i] for i in eeg_idx], cwl_channels=[names[i] for i in cwl_idx],
        ecg=[names[i] for i in ecg_idx], samples_out=int(out.shape[1])), indent=2))
    return out, names


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("vhdr"); p.add_argument("out_dir")
    p.add_argument("--rp", help="SPM rp_*.txt (omit -> moving average)")
    p.add_argument("--marker", default="R128")
    p.add_argument("--cwl-regex", default=r"^CW")
    p.add_argument("--ecg", default="ECG")
    p.add_argument("--exclude", nargs="*", default=[], help="channel names to leave out of the EEG targets")
    p.add_argument("--templates", type=int, default=25)
    p.add_argument("--threshold", type=float, default=0.3, help="mm")
    p.add_argument("--span", type=int, nargs=2, metavar=("PRE", "POST"),
                   help="artifact span in samples relative to the marker (default 0 TR-1)")
    p.add_argument("--target-fs", type=float, default=500)
    a = p.parse_args()
    run(a.vhdr, a.out_dir, a.rp, a.marker, a.cwl_regex, a.ecg, a.exclude, a.templates,
        a.threshold, tuple(a.span) if a.span else None, target_fs=a.target_fs)
