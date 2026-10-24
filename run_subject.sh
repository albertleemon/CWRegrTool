#!/usr/bin/env bash
# Complete EEG-fMRI artifact correction for one run, with all diagnostics.
#
#   ./run_subject.sh RAW.vhdr OUT_DIR [excluded labels, default "65,66"]
#
# Stages (each saved in OUT_DIR):
#   02_gradient_corrected.set  gradient template subtraction, 5000 Hz (Allen 2000 / Moosmann 2009)
#   02b_obs_corrected.set      OBS removal of the gradient residual, 5000 Hz (Niazy 2005)
#   03_filtered.set            0.33-125 Hz band-pass + 500 Hz
#   04_cwl_corrected.set       CWL regression (van der Meer 2016)  <- final result
# Diagnostics:
#   stage_spectra.png          mean spectrum raw / 02 / 02b / 04 / random-phase-CWL null
#   diagnostics_py.png, diagnostics.json, residual_peaks.png/.json, matlab_run.log
#
# Environment (optional):
#   MATLAB_BIN     MATLAB executable (default: matlab)
#   PYTHON         Python with numpy/scipy/mne/pymatreader/matplotlib (default: python)
#   CWL_OBS        OBS components (default 4; 0 = skip the OBS stage)
#   CWL_OBS_HP     high-pass (Hz) of the residual for OBS (default 70; 60 also treats the 67 Hz harmonic)
#   CWL_LOWPASS    low-pass (Hz) of the estimated CWL artifact before subtraction (default off;
#                  e.g. 90: above 90 Hz the EEG stays as after OBS - no loop noise added, nothing removed)
#   CWL_RP         SPM rp_*.txt or fMRIPrep confounds .tsv (default: moving-average templates)
#   CWL_THRESHOLD  head-motion threshold in mm (default 0.3)
set -euo pipefail
VHDR="$(realpath "$1")"; OUT="$(realpath -m "$2")"; EXCL="${3:-65,66}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PYTHON:-python}"
mkdir -p "$OUT"
export CWL_OBS="${CWL_OBS:-4}" CWL_OBS_HP="${CWL_OBS_HP:-}" CWL_LOWPASS="${CWL_LOWPASS:-}"
echo "settings: CWL_OBS=$CWL_OBS CWL_OBS_HP=${CWL_OBS_HP:-70} CWL_LOWPASS=${CWL_LOWPASS:-off} CWL_RP=${CWL_RP:-none} CWL_THRESHOLD=${CWL_THRESHOLD:-0.3}" | tee "$OUT/settings.txt"

echo "== 1/4 correction (MATLAB): gradient -> OBS ($CWL_OBS components) -> filter -> CWL"
"$HERE/matlab/run_headless.sh" "$VHDR" "$OUT" "$EXCL"

echo "== 2/4 CWL diagnostics incl. random-phase null (MATLAB)"
"${MATLAB_BIN:-matlab}" -nodisplay -nosplash -batch \
  "cd('$HERE/matlab'); setup_cwl_pipeline; cwl_diagnose('$OUT')" 2>&1 | tee -a "$OUT/matlab_run.log"

echo "== 3/4 stage spectra (raw / gradient / OBS / CWL / null)"
"$PY" "$HERE/python/plot_stage_spectra.py" "$OUT" --vhdr "$VHDR"

echo "== 4/4 figures and residual-peak check"
"$PY" "$HERE/python/plot_diagnostics_mat.py" "$OUT"
"$PY" "$HERE/python/analyze_residual_peaks.py" "$OUT" > "$OUT/residual_peaks.txt"
echo "done: $OUT"
