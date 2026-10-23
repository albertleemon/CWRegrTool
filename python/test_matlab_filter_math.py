"""Python mirror of matlab/cwl_filter.m (same steps) checked against scipy, to validate
the math of the toolbox-free MATLAB filter."""
import numpy as np
from scipy.signal import butter, lfilter, sosfiltfilt


def design(kind, fc, fs, n):
    wc = 2 * fs * np.tan(np.pi * fc / fs)
    p = np.exp(1j * np.pi * (2 * np.arange(1, n + 1) + n - 1) / (2 * n))
    s, z0 = (wc * p, -1.0) if kind == "low" else (wc / p, 1.0)
    zp = (2 * fs + s) / (2 * fs - s)
    secs = []
    zr = 1.0 if kind == "low" else -1.0

    def norm(b, a):
        g = np.sum(a * zr ** np.arange(len(a))) / np.sum(b * zr ** np.arange(len(b)))
        return b * g, a
    for q in zp[zp.imag > 1e-12]:
        secs.append(norm(np.poly([z0, z0]).real, np.poly([q, np.conj(q)]).real))
    for q in zp[np.abs(zp.imag) <= 1e-12]:
        secs.append(norm(np.array([1.0, -z0]), np.array([1.0, -q.real])))
    return secs


def steady_state(b, a):
    m = len(a) - 1
    C = np.array([[-a[1]]]) if m == 1 else np.vstack([-a[1:], np.eye(m - 1, m)])
    return np.linalg.solve(np.eye(m) - C.T, b[1:] - a[1:] * b[0])


def run(x, secs):
    for b, a in secs:
        x, _ = lfilter(b, a, x, zi=steady_state(b, a) * x[0])
    return x


def zerophase(x, secs, padlen):
    n = len(x)
    padlen = min(padlen, n - 1)
    xe = np.concatenate([2 * x[0] - x[padlen:0:-1], x, 2 * x[-1] - x[-2:-padlen - 2:-1]])
    y = run(xe, secs)
    y = run(y[::-1], secs)[::-1]
    return y[padlen:padlen + n]


def test_matches_scipy_butter():
    fs = 5000
    x = np.random.default_rng(1).normal(size=200000)
    for kind, fc, n in [("high", 1.0, 2), ("low", 125.0, 3), ("high", 0.33, 3)]:
        mine = zerophase(x, design(kind, fc, fs, n), round(3 * fs / fc))
        ref = sosfiltfilt(butter(n, fc, btype=kind, fs=fs, output="sos"), x)
        err = np.abs(mine - ref)[60000:140000].max()      # edges differ only by padding convention
        assert err < 1e-5 * np.abs(ref).max(), (kind, fc, n, err)
    # edge behaviour: no big transient from a constant offset
    y = zerophase(np.full(40000, 50.0), design("high", 0.33, fs, 3), round(3 * fs / 0.33))
    assert np.abs(y).max() < 1e-3
