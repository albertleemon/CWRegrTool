#!/usr/bin/env bash
# Run the pipeline without the MATLAB GUI and keep a full log.
#   ./run_headless.sh /path/to/file.vhdr /path/to/out_dir [excluded labels, e.g. "65,66"]
# Optional environment variables:
#   MATLAB_BIN     MATLAB executable (default: matlab), e.g. /usr/local/MATLAB/R2024b/bin/matlab
#   CWL_RP         rp_*.txt or fMRIPrep *_desc-confounds_timeseries.tsv
#   CWL_THRESHOLD  head-motion threshold in mm (default 0.3)
# Log: <out_dir>/matlab_run.log ; non-zero exit code on error.
set -euo pipefail
VHDR="$1"; OUT="$2"; EXCL="${3:-65,66}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$OUT"
export CWL_VHDR="$VHDR" CWL_OUT="$OUT" CWL_EXCLUDE="$EXCL" CWL_RP="${CWL_RP:-}" CWL_THRESHOLD="${CWL_THRESHOLD:-}"
"${MATLAB_BIN:-matlab}" -nodisplay -nosplash -nodesktop -batch \
  "cd('$HERE'); setup_cwl_pipeline; debug_cwl_env; run_headless_job" 2>&1 | tee "$OUT/matlab_run.log"
exit "${PIPESTATUS[0]}"
