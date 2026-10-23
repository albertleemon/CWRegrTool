import numpy as np
from pybv import write_brainvision
from run_pipeline import run


def test_end_to_end(tmp_path):
    fs, tr, nvol = 5000, 9000, 40
    n = tr * (nvol + 2)
    rng = np.random.default_rng(0)
    t = np.arange(n) / fs
    eeg = rng.normal(scale=3, size=(4, n)) + 15 * np.sin(2 * np.pi * 10.37 * t)
    bcg = np.sin(2 * np.pi * 1.1 * t)
    cw = np.vstack([bcg * (1 + i) + rng.normal(scale=.05, size=n) for i in range(2)])
    art_bcg = 40 * bcg
    grad = np.zeros(n)
    on = tr + tr * np.arange(nvol)
    tt = np.arange(tr) / fs
    for o in on:
        grad[o:o + tr] = 2000 * np.sin(2 * np.pi * 15.4 * tt) * np.exp(-tt)
    data = eeg + art_bcg + grad
    cwd = cw + grad * 0.5
    all_ = np.vstack([data, cwd, np.zeros((1, n))]) * 1e-6          # volts
    names = ["Fz", "Cz", "Pz", "Oz", "CW1", "CW2", "ECG"]
    write_brainvision(data=all_, sfreq=fs, ch_names=names, fname_base="t", folder_out=tmp_path,
                      events=np.array([[o, 1, 128] for o in on]), unit="V")
    out, nm = run(tmp_path / "t.vhdr", tmp_path / "out", marker="S  1")
    assert out.shape[0] == 7
    # artifact-free reference at 500 Hz in the trimmed span
    from cwregression import bandpass
    from scipy.signal import resample_poly
    ref = resample_poly(bandpass((eeg * 1.0)[:, on[0]:on[-1] + tr], fs), 1, 10, axis=1)
    got = out[:4]
    r = lambda x: np.sqrt(np.mean(x ** 2))
    print("residual", r(got - ref[:, :got.shape[1]]), "artifact", r(40 * bcg))
    assert r(got - ref[:, :got.shape[1]]) < 0.15 * r(40 * bcg)
