#!/usr/bin/env bash
# All inside-vs-outside evaluations for one processed run.
#
#   ./evaluate_subject.sh OUT_DIR OUTSIDE.vhdr
#
# OUT_DIR is a run_subject.sh output folder (needs 03_filtered.set and 04_cwl_corrected.set).
# Writes into OUT_DIR:
#   inside_vs_outside.png / .channels.csv   traces, mean spectrum, per-channel RMS 1-30 Hz
#   outside_comparison.png / .json          spectral distance (eq. 4), band differences, alpha topography
#   bcg_residual.png / .json                heartbeat-locked residual (03, 04, outside)
#   removed_signal.json                     coherence of the removed part / of the channels with neighbours
#   evaluation.txt                          all terminal output
set -euo pipefail
OUT="$(realpath "$1")"; OUTSIDE="$(realpath "$2")"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="${PYTHON:-python}"
export MPLBACKEND=Agg
{
echo "== inside vs outside (traces, spectrum, per-channel RMS)"
"$PY" "$HERE/python/compare_inside_outside.py" "$OUT/04_cwl_corrected.set" "$OUTSIDE" \
      "$OUT/inside_vs_outside.png" --also "$OUT/03_filtered.set" --channels Fz FC1 FC2 Cz C1 C2 CPz Pz Oz
echo; echo "== paper metrics vs outside"
"$PY" "$HERE/python/evaluate_outside.py" "$OUT" "$OUTSIDE"
echo; echo "== heartbeat-locked residual"
"$PY" "$HERE/python/check_bcg_residual.py" "$OUT" --outside "$OUTSIDE"
echo; echo "== removed component / neighbour coherence"
"$PY" "$HERE/python/check_removed_signal.py" "$OUT" --outside "$OUTSIDE" \
      --channels Fz FC1 FC2 Cz C1 C2 CPz Pz Oz
} 2>&1 | grep -v "matplotlibrc\|Bad value in file" | tee "$OUT/evaluation.txt"
