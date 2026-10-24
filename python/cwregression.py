"""Carbon-wire-loop (CWL) windowed regression for EEG/fMRI artifact removal.

Python re-implementation of the method in van der Meer et al. (2016),
NeuroImage 125:880-894, and of the MATLAB EEGlab plugin in this repository
(m_do_taperedhann.m / tools.fit_regmat_to_signalmat.m).

Method
------
For every overlapping Hann-weighted window the EEG of each channel is
regressed on all CWL signals, each delay-embedded over +/-K samples (K = 21 ms
by default) so that a temporal lag and a linear "transfer function" between
loop and electrode can be fitted. The fitted part is the movement artifact
(BCG + helium pump) in that window. Windows overlap so that their weights sum
to 1; the weighted fits are overlap-added into one artifact trace which is
subtracted from the EEG.
"""
import numpy as np
from scipy.signal import butter, sosfiltfilt


def bandpass(x, fs, lo=0.33, hi=125.0, order=3):
    """Zero-phase band-pass along the last axis (applied to EEG and CWL)."""
    hi = min(hi, 0.99 * fs / 2)
    sos = butter(order, [lo, hi], btype="band", fs=fs, output="sos")
    return sosfiltfilt(sos, x, axis=-1)


def hann(n):
    """Symmetric Hann window of n samples (same as hann.m)."""
    return 0.5 - 0.5 * np.cos(2 * np.pi * np.arange(n) / (n - 1))


def delay_embed(regs, k):
    """Stack regs shifted by -k..k samples -> (n_regs*(2k+1), n_samples).

    Shifts use the true neighbouring samples; only the record edges are
    zero-padded. Row order is delay-major (all regressors at one lag, then next).
    """
    n = regs.shape[1]
    pad = np.pad(regs, ((0, 0), (k, k)))
    return np.vstack([pad[:, k + s : k + s + n] for s in range(-k, k + 1)])


def cw_regression(eeg, cwl, fs, window_s=4.0, delay_s=0.021, taper=1,
                  filter_band=(0.33, 125.0), return_artifact=False, cwl_lowpass=None):
    """Remove motion-related artifacts from `eeg` using carbon-wire signals.

    eeg : (n_chan, n_samples) MR-gradient-corrected EEG
    cwl : (n_cwl,  n_samples) carbon-wire-loop signals, same sampling clock
    window_s : window length in seconds (paper: 4 s)
    delay_s  : half-width of delay embedding in seconds (paper: 21 ms)
    taper    : 2**taper hops per window; taper=1 -> 50 % overlap Hann (sum=1)
    filter_band : band-pass applied to both inputs before regression, or None
    cwl_lowpass : optional low-pass (Hz) applied to the CWL regressors only

    Returns the corrected EEG (and the subtracted artifact if requested).
    """
    eeg = np.asarray(eeg, float)
    cwl = np.asarray(cwl, float)
    if eeg.shape[1] != cwl.shape[1]:
        raise ValueError("eeg and cwl must have the same number of samples")
    if filter_band:
        eeg = bandpass(eeg, fs, *filter_band)
        cwl = bandpass(cwl, fs, *filter_band)
    if cwl_lowpass:
        sos = butter(8, cwl_lowpass, btype="low", fs=fs, output="sos")   # steep: 8th order, zero phase
        cwl = sosfiltfilt(sos, cwl, axis=-1)

    k = int(np.floor(fs * delay_s))
    n_hops = 2 ** taper
    n_win = int(np.floor(fs * window_s + 1))
    while (n_win - 1) % n_hops:          # window edges must land on samples
        n_win += 1
    step = (n_win - 1) // n_hops
    n = eeg.shape[1]
    if n < n_win:
        raise ValueError("recording shorter than one window")

    w = hann(n_win)
    expregs = delay_embed(cwl, k)
    fit = np.zeros_like(eeg)
    wsum = np.zeros(n)

    starts = list(range(0, n - n_win + 1, step))
    if starts[-1] + n_win < n:           # cover the tail with a final window
        starts.append(n - n_win)
    for s in starts:
        sl = slice(s, s + n_win)
        ws = w
        if s == starts[0] or s == starts[-1]:   # keep full weight at the
            ws = w.copy()                       # record edges, where no
            if s == starts[0]:                  # neighbour window overlaps
                ws[: n_win // 2] = 1.0
            if s == starts[-1]:
                ws[n_win // 2 :] = 1.0
        R = expregs[:, sl] * ws          # windowed regressors
        Y = eeg[:, sl] * ws              # windowed data
        beta = np.linalg.lstsq(R.T, Y.T, rcond=None)[0]   # (n_regs, n_chan)
        fit[:, sl] += (beta.T @ R)       # ~ w * artifact in this window
        wsum[sl] += ws

    # Overlapping windows sum to 1 in the interior; dividing by the summed
    # weights (as m_do_taperedhann.m does) also keeps the record edges unbiased.
    ok = wsum > 1e-12
    artifact = np.zeros_like(eeg)
    artifact[:, ok] = fit[:, ok] / wsum[ok]
    corrected = eeg - artifact
    return (corrected, artifact) if return_artifact else corrected
