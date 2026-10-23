"""Draw diagnostics.png from the diagnostics_data.mat written by matlab/cwl_diagnose.m.

    python plot_diagnostics_mat.py ~/cwl_test/matlab_crop
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

d = Path(sys.argv[1]).expanduser()
m = loadmat(d / "diagnostics_data.mat", squeeze_me=True)
f = m["f"]
fig, ax = plt.subplots(2, 2, figsize=(13, 8))
ax[0, 0].semilogy(f, m["p0"].mean(1), label="before CWL")
ax[0, 0].semilogy(f, m["p1"].mean(1), label="after CWL")
ax[0, 0].semilogy(f, m["pn"].mean(1), "--", label="null (random-phase CWL)")
ax[0, 0].set(xlim=(0, 125), xlabel="Hz", ylabel="uV^2/Hz", title="EEG mean spectrum"); ax[0, 0].legend()
ax[0, 1].semilogy(f, m["pc"], lw=.8)
ax[0, 1].set(xlim=(0, 125), xlabel="Hz", title="CWL spectra"); ax[0, 1].legend(list(m["cwl_names"]), fontsize=8)
ax[1, 0].plot(m["t"], m["tr"][0], lw=.7, label="before")
ax[1, 0].plot(m["t"], m["tr"][1], lw=.7, label="after")
ax[1, 0].set(xlabel="s", ylabel="uV", title=f"{m['oz']}, first 10 s"); ax[1, 0].legend()
k = np.arange(len(m["vr"]))
ax[1, 1].bar(k - .2, m["vr"], .4, label="CWL"); ax[1, 1].bar(k + .2, m["vrn"], .4, label="null")
ax[1, 1].set(xlabel="EEG channel #", ylabel="% variance removed", title="variance removed"); ax[1, 1].legend()
# ignore the near-zero bins at DC (removed by the high-pass) when scaling the y axis
band = (f > 0.5) & (f <= 125)
for a, p in ((ax[0, 0], m["p0"]), (ax[0, 1], m["pc"])):
    a.set_ylim(bottom=0.5 * np.min(p[band]), top=2 * np.max(p[band]))
fig.tight_layout()
fig.savefig(d / "diagnostics_py.png", dpi=110)
print("wrote", d / "diagnostics_py.png")
