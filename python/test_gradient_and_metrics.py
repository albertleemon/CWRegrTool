import numpy as np
from gradient_correction import (rp_weighting_matrix, moving_average_matrix,
                                 gradient_correct)
from metrics import (relative_distance, jaccard_overlap, evaluate,
                     eoec_band, power_spectrum)


def make_scan(fs=1000, tr=0.5, n_vol=120, jump_at=60, seed=0):
    rng = np.random.default_rng(seed)
    L = int(fs * tr)
    n = L * n_vol + 2 * L
    t = np.arange(L) / fs
    base = 3000 * np.sin(2 * np.pi * 40 * t) * np.exp(-t * 3)
    shifted = np.roll(base, 7) * 1.15                    # template after head motion
    brain = rng.normal(scale=5, size=(3, n))
    brain += 20 * np.sin(2 * np.pi * 10 * np.arange(n) / fs)
    art = np.zeros((3, n))
    onsets = L + L * np.arange(n_vol)
    for i, o in enumerate(onsets):
        art[:, o:o + L] = (base if i < jump_at else shifted) * np.array([[1], [.7], [.4]])
    rp = np.zeros((n_vol, 6))
    rp[jump_at:, 0] = 1.0                                # 1 mm jump
    return fs, L, brain, art, onsets, rp


def test_rp_matrix_avoids_templates_across_motion():
    *_, rp = make_scan()
    W = rp_weighting_matrix(rp, 120, threshold=0.3, k=25)
    assert (W.sum(1) == 25).all()
    # volumes before the jump never borrow templates from after it
    assert W[10, 61:].sum() == 0 and W[100, :59].sum() == 0
    assert np.allclose(moving_average_matrix(120, 25).sum(1), 25)


def test_gradient_correction_rp_beats_moving_average():
    fs, L, brain, art, onsets, rp = make_scan()
    rms = lambda x: np.sqrt(np.mean(x ** 2))
    sl = slice(onsets[0], onsets[-1] + L)
    res = {}
    for name, W in [("rp", rp_weighting_matrix(rp, 120)), ("ma", moving_average_matrix(120, 25))]:
        out = gradient_correct(brain + art, fs, onsets, 0, L - 1, W)
        res[name] = rms((out - brain)[:, sl])
    print("residual rms", res, "artifact rms", rms(art[:, sl]))
    assert res["rp"] < 0.05 * rms(art[:, sl])
    assert res["rp"] < res["ma"]


def test_metrics_identity_and_scale():
    f = np.linspace(0, 25, 251)
    y = np.exp(-(f - 10) ** 2 / 4)
    assert relative_distance(y, y) == 0
    assert abs(relative_distance(2 * y, y) - 100) < 1e-9
    assert jaccard_overlap(f, y * 2, y, y * 2, y, (8, 12)) == 1.0
    assert jaccard_overlap(f, y * 2 + 5, y + 5, y * 2, y, (8, 12)) < 0.5


def test_eoec_evaluate_reference_vs_noisy():
    fs, n, nch = 250, 250 * 120, 6
    rng = np.random.default_rng(1)
    t = np.arange(n) / fs
    topo = np.array([0, 0, 0, 1, 1, 0.5])[:, None]
    def eeg(alpha, noise):
        return alpha * topo * np.sin(2 * np.pi * 10 * t + rng.uniform(0, 6)) + rng.normal(size=(nch, n)) * noise
    ref_o, ref_c = eeg(1, 1), eeg(4, 1)
    posterior = [3, 4, 5]
    good = evaluate(fs, ref_o, ref_c, eeg(1, 1), eeg(4, 1), posterior)
    bad = evaluate(fs, ref_o, ref_c, eeg(1, 6), eeg(4, 6), posterior)
    print(good, bad)
    assert 8 <= good["band"][0] <= 10 and good["band"][1] >= 12.5
    assert good["jaccard"] > bad["jaccard"]
    assert good["spectral_distance_closed"] < bad["spectral_distance_closed"]
