"""Comparison metrics of van der Meer et al. (2016) for judging an artifact
correction against EEG recorded outside the scanner:
  - relative spectral distance (eq. 4)
  - eyes-open/eyes-closed (EO/EC) alpha band, channel selection and weights
  - Jaccard overlap of the alpha-power surfaces (eq. 5)
  - relative distance of topographic distributions (eq. 6)
"""
import numpy as np
from scipy.signal import welch


def power_spectrum(x, fs, nperseg=None):
    """Welch power spectrum along the last axis -> (freqs, (n_chan, n_freq))."""
    return welch(x, fs=fs, nperseg=nperseg or int(4 * fs), axis=-1)


def relative_distance(a, b):
    """sqrt(sum (a-b)^2) / sqrt(sum b^2) * 100  (b = reference / 'outside')."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    return 100.0 * np.sqrt(np.sum((a - b) ** 2)) / np.sqrt(np.sum(b ** 2))


def spectral_distance(sp_in, sp_out):
    """Eq. 4 on (channel-averaged) spectra."""
    return relative_distance(np.mean(sp_in, axis=0) if np.ndim(sp_in) > 1 else sp_in,
                             np.mean(sp_out, axis=0) if np.ndim(sp_out) > 1 else sp_out)


def eoec_band(freqs, sp_open, sp_closed, posterior, fmax=25.0,
              rel_thresh=0.10, min_lower=7.5, min_upper=12.5):
    """Subject-specific EO/EC frequency band from the reference ('outside') EEG.

    sp_open, sp_closed : (n_chan, n_freq) power spectra
    posterior : indices of posterior channels (paper: Pz, PO3, POz, PO4, Oz, ...)
    The band is where the closed-minus-open increase, averaged over the posterior
    channels, is at least `rel_thresh` of its maximum; contiguous around the peak.
    Limits are constrained to lower >= min_lower and upper >= min_upper Hz, as
    stated in the paper (the sentence is ambiguous; this is the literal reading).
    """
    m = freqs <= fmax
    inc = (sp_closed - sp_open)[posterior].mean(0)
    inc = np.where(m, inc, -np.inf)
    peak = int(np.argmax(inc))
    ok = inc >= rel_thresh * inc[peak]
    lo = hi = peak
    while lo > 0 and ok[lo - 1]:
        lo -= 1
    while hi < len(freqs) - 1 and ok[hi + 1]:
        hi += 1
    f_lo, f_hi = max(freqs[lo], min_lower), max(freqs[hi], min_upper)
    return f_lo, f_hi


def eoec_channel_weights(freqs, sp_open, sp_closed, band):
    """Per-channel EO/EC power increase in `band`, the selected channels (those
    above half the sum of min and max increase), and their weights."""
    m = (freqs >= band[0]) & (freqs <= band[1])
    inc = (sp_closed - sp_open)[:, m].sum(1)
    thr = 0.5 * (inc.min() + inc.max())
    sel = np.flatnonzero(inc > thr)
    w = np.zeros(len(inc))
    w[sel] = inc[sel] / inc[sel].sum()
    return inc, sel, w


def weighted_spectrum(sp, w):
    return w @ sp


def jaccard_overlap(freqs, y_closed, y_open, r_closed, r_open, band):
    """Eq. 5. At every frequency the EO/EC 'surface' is the interval between the
    open and closed spectra; overlap = intersection / union of those intervals
    (summed over the band) between the scanner data Y and the reference R."""
    m = (freqs >= band[0]) & (freqs <= band[1])
    ylo, yhi = np.minimum(y_closed, y_open)[m], np.maximum(y_closed, y_open)[m]
    rlo, rhi = np.minimum(r_closed, r_open)[m], np.maximum(r_closed, r_open)[m]
    inter = np.clip(np.minimum(yhi, rhi) - np.maximum(ylo, rlo), 0, None).sum()
    union = (yhi - ylo).sum() + (rhi - rlo).sum() - inter
    return inter / union if union > 0 else 0.0


def alpha_increase_per_channel(freqs, sp_open, sp_closed, band):
    m = (freqs >= band[0]) & (freqs <= band[1])
    return (sp_closed - sp_open)[:, m].sum(1)


def topographic_distance(d_in, d_out):
    """Eq. 6: relative distance between alpha-power-increase topographies."""
    return relative_distance(d_in, d_out)


def evaluate(fs, ref_open, ref_closed, in_open, in_closed, posterior):
    """Convenience: all three metrics for one corrected data set vs reference.
    Inputs are (n_chan, n_samples) arrays, same channel order."""
    f, ro = power_spectrum(ref_open, fs)
    _, rc = power_spectrum(ref_closed, fs)
    _, io = power_spectrum(in_open, fs)
    _, ic = power_spectrum(in_closed, fs)
    band = eoec_band(f, ro, rc, posterior)
    _, sel, w = eoec_channel_weights(f, ro, rc, band)
    jac = jaccard_overlap(f, weighted_spectrum(ic, w), weighted_spectrum(io, w),
                          weighted_spectrum(rc, w), weighted_spectrum(ro, w), band)
    topo = topographic_distance(alpha_increase_per_channel(f, io, ic, band),
                                alpha_increase_per_channel(f, ro, rc, band))
    return dict(band=band, channels=sel, spectral_distance_closed=spectral_distance(ic, rc),
                jaccard=jac, topographic_distance=topo)
