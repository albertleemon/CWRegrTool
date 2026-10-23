"""Convert an fMRIPrep confounds .tsv into an SPM-style rp_*.txt (6 columns, mm / rad),
for tools that only read rp files (the MATLAB Bergen functions).

    python fmriprep_to_rp.py sub-X_task-rest_run-01_desc-confounds_timeseries.tsv rp_run01.txt
"""
import sys
import numpy as np
from gradient_correction import load_motion, motion_steps

rp = load_motion(sys.argv[1])
np.savetxt(sys.argv[2], rp, fmt="%.8f")
st = motion_steps(rp)
print(f"{len(rp)} volumes; frame-to-frame translation: max {st.max():.3f} mm, "
      f"steps > 0.1 mm: {(st > 0.1).sum()}, > 0.2 mm: {(st > 0.2).sum()}, > 0.3 mm: {(st > 0.3).sum()}")
