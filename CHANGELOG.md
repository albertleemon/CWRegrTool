# Changelog

## v0.3.0
- `python/denoise_report.py` (run by `evaluate_subject.sh`): indicators per step (spectral distance to outside,
  RMS, band powers and dB vs outside, slice-harmonic peaks, heartbeat-locked residual, step-wise change, CWL variance
  removed, log-spectrum and alpha-topography correlation), mean spectrum of EVERY step (raw/02/02b/03/04/05/outside),
  scalp maps per step. Outputs `denoise_*.png/.json/.csv`.
- Run-safety fixes from the review: run manifest (stale 02b/05 no longer picked up), volume-marker gap check,
  Python output keeps events/channel types/positions, BrainVision channel units, CWL valid span from the tool's
  window weights, stronger `cwl_crosscheck` (8 checks, errors on FAIL), fuller `EEG.etc.cwl_pipeline` provenance.

## v0.2.0
- Optional heartbeat-locked BCG AAS (05), stand-alone MATLAB package (no EEGLAB), MATLAB-vs-Python cross-check.
