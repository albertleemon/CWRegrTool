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
from gradient_correction import (gradient_correct, load_motion, motion_steps,
                                 moving_average_matrix, rp_weighting_matrix)


def volume_onsets(raw, marker):
    """Sample indices of all annotations whose description ends with `marker`
    (BrainVision 'Response,R128' is read by MNE as 'Response/R128')."""
    ann = raw.annotations
    on = [int(round((a["onset"] - raw.first_time) * raw.info["sfreq"]))
          for a in ann if a["description"].split("/")[-1] == marker]
    return np.array(sorted(on), int)


def check_volume_markers(onsets, tr):
    """A missing marker leaves a stretch without gradient correction that the trimming
    would keep; duplicated or extra markers corrupt the templates."""
    tol = max(5, int(round(0.01 * tr)))
    bad = np.flatnonzero(np.abs(np.diff(onsets) - tr) > tol)
    if len(bad):
        raise SystemExit(
            f"{len(bad)} volume interval(s) deviate from the TR of {tr} samples by more than {tol} "
            f"(first at markers {(bad[:5] + 1).tolist()}, intervals {np.diff(onsets)[bad[:5]].tolist()} "
            "samples): missing, duplicated or extra markers. Fix them or pass --no-marker-check.")


def resolve_header(vhdr, out_dir):
    """If the header points to files that no longer exist (e.g. DataFile=RS1.eeg after
    renaming the recording), write a patched header/marker copy that points to the
    sibling .eeg/.vmrk with the same stem."""
    vhdr = Path(vhdr)
    h = vhdr.read_text(encoding="utf-8", errors="replace")
    data = re.search(r"DataFile=(.+)", h).group(1).strip()
    mrk = re.search(r"MarkerFile=(.+)", h).group(1).strip()
    ok = (vhdr.parent / data).exists() and (vhdr.parent / mrk).exists()
    if ok:
        return vhdr
    eeg, vmrk = vhdr.with_suffix(".eeg"), vhdr.with_suffix(".vmrk")
    if not (eeg.exists() and vmrk.exists()):
        raise SystemExit(f"header points to {data}/{mrk}, which are missing, and no "
                         f"{eeg.name}/{vmrk.name} next to the header")
    tmp = Path(out_dir) / "_patched"
    tmp.mkdir(parents=True, exist_ok=True)
    m = vmrk.read_text(encoding="utf-8", errors="replace")
    m = re.sub(r"DataFile=.+", f"DataFile={eeg.resolve()}", m)
    (tmp / vmrk.name).write_text(m, encoding="utf-8")
    h = re.sub(r"DataFile=.+", f"DataFile={eeg.resolve()}", h)
    h = re.sub(r"MarkerFile=.+", f"MarkerFile={vmrk.name}", h)
    (tmp / vhdr.name).write_text(h, encoding="utf-8")
    print(f"header pointed to {data}; using {eeg.name} instead")
    return tmp / vhdr.name


def run(vhdr, out_dir, rp=None, marker="R128", cwl_regex=r"^CW", ecg="ECG",
        exclude=(), n_templates=25, threshold=0.3, span=None,
        band=(0.33, 125.0), target_fs=500, window_s=4.0, delay_s=0.021, obs_npc=0,
        obs_hp=70.0, artifact_lowpass=None, bcg_aas_on=False, bcg_hr=None,
        check_markers=True):
    import mne
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    vhdr = resolve_header(vhdr, out_dir)
    raw = mne.io.read_raw_brainvision(vhdr, preload=True, verbose="error")
    fs = raw.info["sfreq"]
    fs_in = fs
    names = raw.ch_names
    data = raw.get_data() * 1e6                          # volts -> microvolts

    cwl_idx = [i for i, n in enumerate(names) if re.match(cwl_regex, n)]
    ecg_idx = [i for i, n in enumerate(names) if n == ecg]
    eeg_idx = [i for i, n in enumerate(names)
               if i not in cwl_idx + ecg_idx and n not in exclude]
    if not cwl_idx:
        raise SystemExit(f"no CWL channels match {cwl_regex!r}: {names}")

    clip = (np.abs(data) >= 0.999 * np.abs(data).max(axis=1, keepdims=True)).mean(axis=1)
    clipped = {names[i]: float(c) for i, c in enumerate(clip) if c > 1e-4}

    onsets = volume_onsets(raw, marker)
    if len(onsets) < 3:
        raise SystemExit(f"found {len(onsets)} '{marker}' markers; need the volume markers")
    tr = int(np.median(np.diff(onsets)))
    if check_markers:
        check_volume_markers(onsets, tr)
    pre, post = span if span else (0, tr - 1)            # artifact span relative to marker

    # 1. gradient artifact: every channel, CWL included
    n_events = 0
    if rp is not None:
        motion = load_motion(rp)
        if len(motion) != len(onsets):
            print(f"WARNING: {len(motion)} motion rows but {len(onsets)} volume markers")
        steps = motion_steps(motion[:len(onsets)])
        n_events = int((steps > threshold).sum())
        print(f"head motion: {n_events} volume steps > {threshold} mm (max step {steps.max():.3f} mm)")
        if n_events == 0:
            print("  -> no motion events: templates fall back to moving average; "
                  "lower --threshold to use the motion information")
        W = rp_weighting_matrix(motion, len(onsets), threshold, n_templates)
        method = "rp_info"
    else:
        W = moving_average_matrix(len(onsets), n_templates)
        method = "moving_average"
    if onsets[-1] + post >= data.shape[1]:
        raise SystemExit("last artifact exceeds the recording: drop the last marker or trim")
    data = gradient_correct(data, fs, onsets, pre, post, W, obs_npc=obs_npc, obs_hp=obs_hp)

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
    pre_cwl = data[eeg_idx].astype(np.float32)
    eeg, art = cw_regression(data[eeg_idx], data[cwl_idx], fs, window_s=window_s,
                             delay_s=delay_s, filter_band=None, return_artifact=True,
                             artifact_lowpass=artifact_lowpass)

    aas_info = None
    if bcg_aas_on:                      # optional: heartbeat-locked AAS after CWL (Allen 1998)
        from bcg_aas import bcg_aas
        if not ecg_idx:
            raise SystemExit(f"--bcg-aas needs an ECG channel ({ecg!r} not found)")
        eeg, pk, hr = bcg_aas(eeg, data[ecg_idx[0]], fs, hr_hint=bcg_hr)
        aas_info = dict(n_beats=int(len(pk)), heart_rate_bpm=float(hr))
        print(f"BCG AAS after CWL: {len(pk)} beats, heart rate {hr:.1f} bpm")
    out = data.copy()
    out[eeg_idx] = eeg
    # keep what the input had: channel roles, positions, and the events inside the kept span
    types = ["ecg" if i in ecg_idx else "eeg" if i in eeg_idx else "misc" for i in range(len(names))]
    info = mne.create_info(names, fs, ch_types=types)
    out_raw = mne.io.RawArray(out * 1e-6, info, verbose="error")
    montage = raw.get_montage()
    if montage is not None:
        out_raw.set_montage(montage, on_missing="ignore", verbose="error")
    ann = raw.annotations
    t = np.array([an["onset"] - raw.first_time for an in ann]) - a / fs_in    # s from the first output sample
    dur = np.array([an["duration"] for an in ann])
    keep = (t >= 0) & (t < out.shape[1] / fs)
    out_raw.set_annotations(mne.Annotations(t[keep], dur[keep],
                                            [an["description"] for an, k in zip(ann, keep) if k]))
    out_raw.save(out_dir / "cwl_corrected_raw.fif", overwrite=True)
    np.savez_compressed(out_dir / "cwl_corrected.npz", data=out, ch_names=names, fs=fs,
                        eeg_idx=eeg_idx, cwl_idx=cwl_idx, artifact=art,
                        pre_cwl=pre_cwl)
    (out_dir / "pipeline_log.json").write_text(json.dumps(dict(
        vhdr=str(vhdr), method=method, obs_npc=obs_npc, obs_hp=obs_hp, artifact_lowpass=artifact_lowpass, bcg_aas=aas_info, n_volumes=len(onsets), tr_samples=tr,
        span=[pre, post], n_templates=n_templates, threshold_mm=threshold, band=band,
        fs_out=fs, window_s=window_s, delay_s=delay_s,
        motion_events=(n_events if rp is not None else None),
        eeg_channels=[names[i] for i in eeg_idx], cwl_channels=[names[i] for i in cwl_idx],
        ecg=[names[i] for i in ecg_idx], clipped_channels_fraction=clipped, samples_out=int(out.shape[1])), indent=2))
    return out, names


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("vhdr"); p.add_argument("out_dir")
    p.add_argument("--rp", help="SPM rp_*.txt or fMRIPrep *_desc-confounds_timeseries.tsv (omit -> moving average)")
    p.add_argument("--marker", default="R128")
    p.add_argument("--no-marker-check", action="store_true",
                   help="do not stop on missing/duplicated/uneven volume markers")
    p.add_argument("--cwl-regex", default=r"^CW")
    p.add_argument("--ecg", default="ECG")
    p.add_argument("--exclude", nargs="*", default=[], help="channel names to leave out of the EEG targets")
    p.add_argument("--templates", type=int, default=25)
    p.add_argument("--threshold", type=float, default=0.3, help="mm")
    p.add_argument("--span", type=int, nargs=2, metavar=("PRE", "POST"),
                   help="artifact span in samples relative to the marker (default 0 TR-1)")
    p.add_argument("--target-fs", type=float, default=500)
    p.add_argument("--obs-npc", type=int, default=0,
                   help="FASTR residual PCA after the gradient template (Niazy 2005); 0 = off, e.g. 4")
    p.add_argument("--obs-hp", type=float, default=70.0, help="high-pass (Hz) of the residual for OBS")
    p.add_argument("--bcg-aas", action="store_true",
                   help="optional heartbeat-locked AAS after the CWL regression (default off)")
    p.add_argument("--bcg-hr", type=float, help="expected heart rate (bpm) for the R-peak detection")
    p.add_argument("--artifact-lowpass", type=float, default=None,
                   help="low-pass (Hz) of the CWL artifact estimate before subtraction (e.g. 90)")
    a = p.parse_args()
    run(a.vhdr, a.out_dir, a.rp, a.marker, a.cwl_regex, a.ecg, a.exclude, a.templates,
        a.threshold, tuple(a.span) if a.span else None, target_fs=a.target_fs,
        obs_npc=a.obs_npc, obs_hp=a.obs_hp, artifact_lowpass=a.artifact_lowpass,
        bcg_aas_on=a.bcg_aas, bcg_hr=a.bcg_hr, check_markers=not a.no_marker_check)
