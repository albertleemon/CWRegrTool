import numpy as np
from cwregression import cw_regression, delay_embed


def synth(fs=500, dur=60, seed=0):
    rng = np.random.default_rng(seed)
    n = fs * dur
    t = np.arange(n) / fs
    # "brain": 10 Hz alpha with slow amplitude modulation + pink-ish noise
    alpha = np.sin(2 * np.pi * 10 * t) * (1 + np.sin(2 * np.pi * 0.1 * t))
    brain = 10 * alpha[None] + np.cumsum(rng.normal(size=(8, n)), 1) * 0.05
    # motion sources: BCG (~1.1 Hz pulse-like) and helium pump (40 Hz-ish bursts)
    bcg = np.sin(2 * np.pi * 1.1 * t) + 0.5 * np.sin(2 * np.pi * 2.3 * t + 1)
    he = np.sin(2 * np.pi * 41 * t) * (1 + np.sin(2 * np.pi * 0.3 * t))
    src = np.vstack([bcg, he])
    # loops see the sources with different gains/delays and some own noise
    cwl = np.vstack([np.roll(src[i % 2], 3 * i) * (1 + .3 * i) for i in range(6)])
    cwl += rng.normal(scale=0.05, size=cwl.shape)
    # electrodes see sources through delays (+ slow gain drift)
    art = np.zeros((8, n))
    for c in range(8):
        art[c] = (60 * (1 + 0.2 * np.sin(2 * np.pi * 0.02 * t + c)) * np.roll(bcg, 4 + c)
                  + 15 * np.roll(he, -3 + c))
    return fs, brain, art, cwl


def rmse(a, b):
    return np.sqrt(np.mean((a - b) ** 2))


def test_artifact_removed():
    fs, brain, art, cwl = synth()
    out = cw_regression(brain + art, cwl, fs)
    from cwregression import bandpass
    clean = bandpass(brain, fs)
    before = rmse(bandpass(brain + art, fs), clean)
    after = rmse(out, clean)
    print(f"rmse vs. clean EEG: before={before:.2f}  after={after:.2f}")
    assert after < 0.1 * before


def test_alpha_preserved():
    fs, brain, art, cwl = synth()
    out = cw_regression(brain + art, cwl, fs)
    f = np.fft.rfftfreq(out.shape[1], 1 / fs)
    p = lambda x: np.abs(np.fft.rfft(x[0])) ** 2
    i = np.argmin(abs(f - 10))
    ratio = p(out)[i] / p(brain)[i]
    print("alpha power ratio:", ratio)
    assert 0.8 < ratio < 1.2


def test_unrelated_regressors_do_limited_damage():
    fs, brain, _, cwl = synth()
    out = cw_regression(brain, cwl * 0 + np.random.default_rng(1).normal(size=cwl.shape), fs)
    from cwregression import bandpass
    # 126 regressors (6 CWL x 43 lags) per 4 s window overfit the slow EEG
    # somewhat even when unrelated; the damage should stay well below 50 %.
    assert rmse(out, bandpass(brain, fs)) < 0.5 * np.std(bandpass(brain, fs))
