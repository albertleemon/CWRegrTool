"""MRI gradient-artifact correction with realignment-parameter-informed
template selection (Moosmann et al. 2009), as used before CWL regression in
van der Meer et al. (2016). Python port of the logic in BergenToolboxModified
(m_rp_info.m, m_moving_average.m, my_m_iscloser2.m, CorrectionMatrix.m).

Not ported: the GUI, clipping detection / NaN interpolation, and the optional
baseline correction step.
"""
import numpy as np
from scipy.signal import butter, sosfiltfilt


def moving_average_matrix(n_vol, k):
    """Row j selects the k volumes closest in time to volume j (Allen 2000)."""
    idx = np.arange(n_vol)
    W = np.zeros((n_vol, n_vol))
    for j in range(n_vol):
        dist = np.abs(idx - j).astype(float)
        dist[j] = 0.0                      # same ordering as the MATLAB triangle
        W[j, np.argsort(dist, kind="stable")[:k]] = 1
    return W


def _runs(mask_idx):
    """Split sorted indices into runs of consecutive integers."""
    runs, cur = [], [mask_idx[0]]
    for i in mask_idx[1:]:
        if i == cur[-1] + 1:
            cur.append(i)
        else:
            runs.append(cur)
            cur = [i]
    runs.append(cur)
    return runs


def _drop_minor_neighbours(motion, k):
    """my_m_iscloser2: of suprathreshold runs closer than k volumes, keep the
    longer (then larger-area) one. Returns indices of removed samples."""
    nz = np.flatnonzero(motion)
    if nz.size == 0:
        return np.array([], int)
    runs = _runs(nz)
    removed = []
    while True:
        mark = []
        for i in range(len(runs) - 1):
            if runs[i + 1][0] - runs[i][-1] < k + 1:
                a, b = runs[i], runs[i + 1]
                if len(a) != len(b):
                    mark.append(i if len(a) < len(b) else i + 1)
                else:
                    sa, sb = motion[a].sum(), motion[b].sum()
                    mark.append(i + 1 if sb < sa else i)   # ties: left one
        if not mark:
            break
        mark = sorted(set(mark))
        smallest = min(mark, key=lambda m: (len(runs[m]), motion[runs[m]].sum()))
        removed += runs[smallest]
        del runs[smallest]
    return np.array(sorted(removed), int)


def load_motion(path):
    """Realignment parameters as (n_vol, 6) [trans x,y,z mm, rot x,y,z rad].

    Accepts an SPM rp_*.txt (whitespace separated, 6 columns) or an fMRIPrep
    *_desc-confounds_timeseries.tsv (columns trans_x..rot_z). Both use mm / radians;
    only the translations matter for the weighting matrix."""
    path = str(path)
    if path.endswith(".tsv"):
        import csv
        with open(path, newline="") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        cols = ["trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z"]
        return np.array([[float(r[c]) for c in cols] for r in rows])
    return np.loadtxt(path)


def motion_steps(rp):
    """Frame-to-frame translation (mm) between consecutive volumes, length n_vol-1."""
    return np.linalg.norm(np.diff(np.asarray(rp, float)[:, :3], axis=0), axis=1)


def rp_weighting_matrix(rp, n_vol, threshold=0.3, k=25):
    """Weighting matrix (n_vol x n_vol) from SPM realignment parameters.

    rp : (n_scans, 6) rp_*.txt (3 translations in mm, 3 rotations in rad)
    n_vol : number of volume markers in the EEG. If rp has fewer rows the
        missing ones are assumed to be dummy scans at the start (zeros are
        prepended); if more, rp is truncated.
    threshold : frame-to-frame translation (mm) counted as head motion
        (the paper used 0.3 mm). Rotations are ignored, as in Bergen.
    """
    rp = np.asarray(rp, float)[:n_vol]
    trans = np.linalg.norm(np.diff(rp[:, :3], axis=0), axis=1)
    motion = np.where(trans > threshold, trans, 0.0)
    motion = np.concatenate([np.zeros(1), motion])       # motion[i]: vol i-1 -> i
    motion = np.concatenate([np.zeros(n_vol - motion.size), motion])  # dummies
    if not motion.any():
        return moving_average_matrix(n_vol, k)

    minor = _drop_minor_neighbours(motion, k)
    motion[minor] = 0.0
    if not motion.any():
        return moving_average_matrix(n_vol, k)

    scaling = k / motion[motion > 0].min()
    W = np.zeros((n_vol, n_vol))
    for j in range(n_vol):                               # 0-based centre volume
        dist = np.empty(n_vol)
        dist[: j + 1] = np.arange(j + 1, 0, -1)
        dist[j + 1:] = np.arange(2, n_vol - j + 1)
        # crossing a motion event between volume j and candidate costs `scaling*motion`
        dist = dist + scaling * np.cumsum(np.concatenate([-motion[: j + 1], motion[j + 1:]]))
        dist[motion > 0] = np.nan                        # never use moving volumes
        order = np.argsort(dist, kind="stable")          # NaN sorts last
        W[j, order[:k]] = 1
    W[:, minor] = 0
    return W


def obs_residual(x, fs, idx, npc=4, hp=70.0, exclude_factor=100.0):
    """Optimal-basis-set removal of the gradient residual (Niazy et al. 2005, FASTR).

    x : one channel after template subtraction; idx : (n_vol, L) sample indices.
    The residual is high-passed at `hp` Hz, the volume epochs are decomposed with PCA
    and the first `npc` components are fitted to and subtracted from every volume.
    Only the high-passed residual is fitted, so content below ~hp Hz (EEG rhythms)
    is not touched. Volumes whose residual power exceeds `exclude_factor` x the median
    (e.g. an incomplete last volume) do not enter the PCA but are still corrected.
    Returns the corrected channel and the fraction of high-pass residual variance removed.
    """
    if npc <= 0:
        return x, 0.0
    sos = butter(4, hp, btype="high", fs=fs, output="sos")
    R = sosfiltfilt(sos, x)[idx]                                  # (n_vol, L)
    pw = (R ** 2).mean(1)
    use = pw <= exclude_factor * np.median(pw)
    _, _, Vt = np.linalg.svd(R[use] - R[use].mean(0), full_matrices=False)
    B = Vt[:npc]                                                  # (npc, L), orthonormal rows
    fit = (R @ B.T) @ B
    y = x.copy()
    y[idx] -= fit
    return y, float(1 - ((R - fit) ** 2).sum() / (R ** 2).sum())


def gradient_correct(eeg, fs, onsets, offset_pre, offset_post, W,
                     highpass=1.0, corr_thresh=2.5, obs_npc=0, obs_hp=70.0):
    """Average-artifact subtraction with a weighting matrix.

    eeg : (n_chan, n_samples); all channels, including CWL, should be passed.
    onsets : sample index of every volume marker (e.g. the R128 events)
    offset_pre/post : artifact spans onsets[i]+offset_pre .. +offset_post
        (inclusive) in samples; spans must lie inside the recording
    W : (n_vol, n_vol) weighting matrix; row i selects the templates for volume i
    Templates are averaged from a copy high-passed at `highpass` Hz (keeps drifts
    out of the templates); templates whose dissimilarity to the others exceeds
    `corr_thresh` x the median are dropped. Samples outside the spans are left
    untouched (and therefore still contain gradient artifacts).
    obs_npc > 0 adds the FASTR optimal-basis-set step (see obs_residual) after the
    template subtraction; it is not part of van der Meer et al. (2016).
    """
    eeg = np.asarray(eeg, float)
    onsets = np.asarray(onsets, int)
    L = offset_post - offset_pre + 1
    starts = onsets + offset_pre
    if starts.min() < 0 or starts.max() + L > eeg.shape[1]:
        raise ValueError("artifact span exceeds the recording; trim the data")
    if W.shape != (len(onsets),) * 2:
        raise ValueError("weighting matrix must be n_vol x n_vol")

    out = eeg.copy()
    sos = butter(2, highpass, btype="high", fs=fs, output="sos") if highpass else None
    idx = starts[:, None] + np.arange(L)[None, :]               # (n_vol, L)
    for c in range(eeg.shape[0]):
        A = eeg[c][idx]
        A2 = sosfiltfilt(sos, eeg[c])[idx] if sos is not None else A
        for i in range(len(starts)):
            sel = np.flatnonzero(W[i])
            if sel.size > 2:
                with np.errstate(invalid="ignore", divide="ignore"):
                    dis = np.nan_to_num(1 - np.corrcoef(A2[sel]), nan=0.0)
                bad = dis.mean(0) > corr_thresh * np.median(dis.mean(0))
                if bad.any() and not bad.all():
                    sel = sel[~bad]
            out[c][idx[i]] = A[i] - A2[sel].mean(0)
        if obs_npc:
            out[c], _ = obs_residual(out[c], fs, idx, obs_npc, obs_hp)
    return out
