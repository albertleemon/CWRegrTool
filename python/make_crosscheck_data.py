"""Synthetic test data + Python results for matlab/cwl_crosscheck.m (MATLAB vs Python).

    python python/make_crosscheck_data.py [OUT.mat]     (default: crosscheck_data.mat here)
"""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from scipy.io import savemat
from gradient_correction import gradient_correct, obs_residual, moving_average_matrix
from bcg_aas import bcg_aas
from cwregression import cw_regression
from scipy.signal import resample_poly
rng = np.random.default_rng(1)
# gradient: 3 channels, fs 5000, 40 volumes TR 0.5 s with slowly changing shape
fs = 5000; TR = 2500; nv = 40
n = TR * (nv + 2)
t = np.arange(TR) / fs
tmpl = 300 * np.sin(2 * np.pi * 16.84 * 5 * t) * np.exp(-t * 3)
x = rng.standard_normal((3, n)) * 5
on = TR + np.arange(nv) * TR
for i, o in enumerate(on):
    x[:, o:o + TR] += tmpl * (1 + 0.01 * i) + 20 * np.sin(2*np.pi*101*t + 0.02*i)
W = moving_average_matrix(nv, 25)
g = gradient_correct(x, fs, on, 0, TR - 1, W)
idx = on[:, None] + np.arange(TR)
o, fr = obs_residual(g[0], fs, idx, 4, 60.0)
# BCG: 500 Hz, ECG with big T wave
fs2 = 500; m = 500 * 300; tt = np.arange(m) / fs2
beats = np.cumsum(60 / 63 + 0.03 * rng.standard_normal(330)); beats = beats[beats < 299]
ecg = np.zeros(m); bcg = np.zeros(m)
for b in beats:
    s = (tt - b)
    ecg += 1000 * np.exp(-(s / 0.01) ** 2) + 1500 * np.exp(-((s - 0.3) / 0.06) ** 2)
    bcg += 15 * np.exp(-((s - 0.25) / 0.08) ** 2) * np.sin(2 * np.pi * 4 * s)
eeg = np.vstack([bcg * a + 5 * rng.standard_normal(m) for a in (1, -0.5, 2)])
ecg += 20 * rng.standard_normal(m)
y, pk, hr = bcg_aas(eeg, ecg, fs2)
# CWL regression: 500 Hz, 3 EEG + 2 loops, artifact = lagged mix of the loops
fs3 = 500; n3 = fs3 * 30; t3 = np.arange(n3) / fs3
cw3 = np.vstack([np.sin(2*np.pi*1.2*t3 + p) + 0.3*np.sin(2*np.pi*7.3*t3 + 2*p) for p in (0, 1)]) \
    + 0.05 * rng.standard_normal((2, n3))
eeg3 = rng.standard_normal((3, n3)) * 3 + np.outer([4, -2, 3], cw3[0]) + np.outer([1, 2, -1], np.roll(cw3[1], 3))
cw_out = cw_regression(eeg3, cw3, fs3, filter_band=None)
# resampling 5000 -> 500 Hz of a smooth signal
t4 = np.arange(5000 * 8) / 5000
x4 = np.sin(2*np.pi*10*t4) + 0.5*np.sin(2*np.pi*40*t4 + 1) + 0.2*np.sin(2*np.pi*90*t4)
rs_y = resample_poly(x4, 1, 10)
savemat(sys.argv[1] if len(sys.argv) > 1 else "crosscheck_data.mat", dict(x=x, on=on + 1, W=W, g=g, o=o, fr=fr, eeg=eeg, ecg=ecg,
        y=y, pk=pk + 1, hr=hr, nbeats=len(beats),
        cw_eeg=eeg3, cw_cwl=cw3, cw_out=cw_out, rs_x=x4, rs_y=rs_y))
print("python: beats true", len(beats), "found", len(pk), "hr", round(hr, 1))
