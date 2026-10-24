#!/usr/bin/env bash
# Run the pipeline without the MATLAB GUI and keep a full log.
#   ./run_headless.sh /path/to/file.vhdr /path/to/out_dir [excluded labels, e.g. "65,66"]
# Optional environment variables:
#   MATLAB_BIN     MATLAB executable (default: matlab), e.g. /usr/local/MATLAB/R2024b/bin/matlab
#   CWL_RP         rp_*.txt or fMRIPrep *_desc-confounds_timeseries.tsv
#   CWL_THRESHOLD  head-motion threshold in mm (default 0.3)
#   CWL_OBS        FASTR residual PCA components after the gradient template (0 = off, e.g. 4)
#   CWL_OBS_HP     high-pass (Hz) of the residual for OBS (default 70, e.g. 60)
#   CWL_LOWPASS    low-pass (Hz) of the CWL regressors before the regression (default off, e.g. 90)
# Log: <out_dir>/matlab_run.log ; non-zero exit code on error.
set -euo pipefail
VHDR="$1"; OUT="$2"; EXCL="${3:-65,66}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$OUT"
MATLAB_CMD="${MATLAB_BIN:-matlab}"
echo "MATLAB used: $(command -v "$MATLAB_CMD" || echo "$MATLAB_CMD")" | tee "$OUT/matlab_run.log"
export CWL_VHDR="$VHDR" CWL_OUT="$OUT" CWL_EXCLUDE="$EXCL" CWL_RP="${CWL_RP:-}" CWL_THRESHOLD="${CWL_THRESHOLD:-}" CWL_OBS="${CWL_OBS:-}" CWL_OBS_HP="${CWL_OBS_HP:-}" CWL_LOWPASS="${CWL_LOWPASS:-}"
"$MATLAB_CMD" -nodisplay -nosplash -nodesktop -batch \
  "cd('$HERE'); setup_cwl_pipeline; debug_cwl_env; run_headless_job" 2>&1 | tee -a "$OUT/matlab_run.log"
exit "${PIPESTATUS[0]}"
