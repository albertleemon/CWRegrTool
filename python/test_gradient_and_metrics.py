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


def test_load_fmriprep_tsv(tmp_path):
    from gradient_correction import load_motion, motion_steps
    cols = ["trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z", "framewise_displacement"]
    rows = ["\t".join(cols), "\t".join(["0.0"] * 6 + ["n/a"]), "\t".join(["0.5", "0", "0", "0", "0", "0", "0.5"])]
    p = tmp_path / "c_desc-confounds_timeseries.tsv"
    p.write_text("\n".join(rows) + "\n")
    rp = load_motion(p)
    assert rp.shape == (2, 6) and abs(motion_steps(rp)[0] - 0.5) < 1e-12


def test_obs_removes_shape_changes_but_keeps_alpha():
    """Artifact whose high-frequency shape changes in some volumes: the 25-volume
    template leaves a residual, the OBS step removes most of it, alpha stays."""
    from gradient_correction import gradient_correct, moving_average_matrix
    from scipy.signal import welch
    fs, L, nv = 5000, 9000, 80
    rng = np.random.default_rng(3)
    n = L * (nv + 2)
    t = np.arange(n) / fs
    brain = 20 * np.sin(2 * np.pi * 10.3 * t + 1) + rng.normal(scale=2, size=n)
    tt = np.arange(L) / fs
    base = sum(np.sin(2 * np.pi * 16.84 * k * tt) for k in range(4, 9)) * 300
    alt = sum(np.sin(2 * np.pi * 16.84 * k * tt + 0.4) for k in range(4, 9)) * 300
    x = brain.copy()
    on = L + L * np.arange(nv)
    for i, o in enumerate(on):
        w = 0.6 + 0.4 * np.cos(2 * np.pi * i / 17)           # slowly changing shape mix
        x[o:o + L] += w * base + (1 - w) * alt
    W = moving_average_matrix(nv, 25)
    sl = slice(on[0], on[-1] + L)
    f, _ = welch(x[sl], fs, nperseg=fs)
    band = lambda y, lo, hi: welch(y[sl], fs, nperseg=fs)[1][(f >= lo) & (f <= hi)].mean()
    a0 = gradient_correct(x[None], fs, on, 0, L - 1, W)[0]
    a4 = gradient_correct(x[None], fs, on, 0, L - 1, W, obs_npc=4)[0]
    hf0, hf4, hfb = band(a0, 65, 160), band(a4, 65, 160), band(brain, 65, 160)
    al0, al4 = band(a0, 9, 12), band(a4, 9, 12)
    print("65-160 Hz power: template only %.3g, +OBS %.3g, clean %.3g; alpha ratio %.3f"
          % (hf0, hf4, hfb, al4 / al0))
    assert hf4 < 0.1 * hf0
    assert abs(al4 / al0 - 1) < 0.01


def test_bcg_aas_removes_locked_residual_keeps_eeg():
    from bcg_aas import bcg_aas
    fs, n = 500, 500 * 180
    rng = np.random.default_rng(4)
    t = np.arange(n) / fs
    beats = np.cumsum(60 / 63 + rng.normal(scale=0.04, size=200)); beats = beats[beats < 178]
    ecg = rng.normal(size=n) * 20
    resid = np.zeros((3, n))
    k = np.arange(int(0.6 * fs)) / fs
    shape = -15 * np.exp(-0.5 * ((k - 0.25) / 0.05) ** 2)
    for b in beats:
        i = int(b * fs)
        ecg += 800 * np.exp(-0.5 * ((t - b) / 0.008) ** 2) + 3000 * np.exp(-0.5 * ((t - b - 0.33) / 0.03) ** 2)
        resid[:, i:i + len(k)] += shape[: max(0, min(len(k), n - i))] * np.array([[1], [.5], [.2]])
    brain = rng.normal(size=(3, n)) * 5 + 8 * np.sin(2 * np.pi * 10.3 * t)
    out, pk, hr = bcg_aas(brain + resid, ecg, fs)
    from check_bcg_residual import locked_average
    rms = lambda x: np.sqrt(np.mean(x ** 2))
    lock = lambda x: rms(locked_average(x, (beats * fs).astype(int), fs, 0.0, 0.6)[0])
    # heartbeat-locked residual (what AAS targets) and broadband change of the EEG
    print("HR %.1f, beats %d/%d; locked residual before %.2f, after %.2f uV; EEG power x%.3f"
          % (hr, len(pk), len(beats), lock(resid), lock(out - brain),
             (out - resid * 0).var() / (brain + resid).var()))
    assert abs(hr - 63) < 3 and abs(len(pk) - len(beats)) <= 2
    assert lock(out - brain) < 0.15 * lock(resid)
