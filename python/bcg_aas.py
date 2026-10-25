"""Optional step after CWL regression: average artifact subtraction of the remaining
heartbeat-locked BCG (Allen et al. 1998, NeuroImage 8:229-239). Not part of
van der Meer et al. (2016), where AAS was a competing method, not a complement.

R peaks come from the ECG channel (robust to the large T waves of in-scanner ECG). Each
beat owns the samples from `pre` s before its R peak to `pre` s before the next R peak;
from those, the template (mean of the k nearest beats, same window) is subtracted over a
length of one median RR interval. Applied to the EEG channels only.
Removes everything locked to the heart beat, including the heartbeat-evoked potential
and the cardiac field artifact: do not use for heartbeat-evoked-potential studies.
"""
import numpy as np
from scipy.signal import butter, find_peaks, sosfiltfilt


def r_peaks(ecg, fs, hr_hint=None):
    """R-peak sample indices. Pan-Tompkins style: 10-35 Hz band-pass (QRS is sharp, the
    in-scanner T wave is large but slow), squared derivative, 80 ms integration; beat
    period from the autocorrelation of that energy (or hr_hint, bpm); peaks closer than
    0.7 of a period are not allowed; each peak is moved to the largest |ECG| deflection
    within +-50 ms."""
    y = sosfiltfilt(butter(3, [10, 35], "band", fs=fs, output="sos"), ecg)
    e = np.convolve(np.diff(y, prepend=y[0]) ** 2, np.ones(int(0.08 * fs)) / int(0.08 * fs), "same")
    if hr_hint:
        period = 60.0 / hr_hint
    else:
        seg = e[: int(min(len(e), 120 * fs))]
        seg = seg - seg.mean()
        ac = np.correlate(seg, seg, "full")[len(seg) - 1:]
        lo, hi = int(0.33 * fs), int(2.0 * fs)              # 30-180 bpm
        period = (lo + int(np.argmax(ac[lo:hi]))) / fs
    pk, _ = find_peaks(e, height=0.3 * np.percentile(e, 99), distance=int(0.7 * period * fs))
    w = int(0.05 * fs)
    return np.array([p - w + int(np.argmax(np.abs(y[p - w:p + w]))) for p in pk
                     if p - w >= 0 and p + w <= len(y)], int)


def bcg_aas(eeg, ecg, fs, k=21, pre=0.25, hr_hint=None):
    """Return (corrected eeg, R-peak indices, heart rate in bpm)."""
    eeg = np.asarray(eeg, float)
    pk = r_peaks(ecg, fs, hr_hint)
    rr = np.diff(pk)
    L = int(np.median(rr))
    a = int(pre * fs)
    ok = (pk - a >= 0) & (pk - a + L <= eeg.shape[1])
    pk_ok = pk[ok]
    E = np.stack([eeg[:, p - a:p - a + L] for p in pk_ok], 1)   # chan x beats x L
    out = eeg.copy()
    n = len(pk_ok)
    for i, p in enumerate(pk_ok):
        lo = int(np.clip(i - k // 2, 0, max(n - k, 0)))
        tmpl = E[:, lo:lo + k].mean(1)
        # this beat owns the samples up to `pre` before the next R peak
        nxt = pk_ok[i + 1] - a if i + 1 < n else p - a + L
        m = min(L, nxt - (p - a))
        out[:, p - a:p - a + m] -= tmpl[:, :m]
    return out, pk, 60.0 * fs / np.median(rr)
