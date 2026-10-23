#!/usr/bin/env bash
# Run the pipeline without the MATLAB GUI and keep a full log.
#   ./run_headless.sh /path/to/file.vhdr /path/to/out_dir [excluded labels, e.g. "65,66"]
# Needs `matlab` on PATH. Log: <out_dir>/matlab_run.log ; non-zero exit code on error.
set -euo pipefail
VHDR="$1"; OUT="$2"; EXCL="${3:-65,66}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$OUT"
export CWL_VHDR="$VHDR" CWL_OUT="$OUT" CWL_EXCLUDE="$EXCL"
matlab -nodisplay -nosplash -nodesktop -batch \
  "cd('$HERE'); setup_cwl_pipeline; debug_cwl_env; run_headless_job" 2>&1 | tee "$OUT/matlab_run.log"
