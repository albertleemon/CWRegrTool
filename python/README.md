# Python port of the van der Meer et al. (2016) pipeline

| file | content |
|---|---|
| `gradient_correction.py` | rp-informed average-artifact subtraction (Moosmann 2009; port of BergenToolboxModified). Run on all channels, CWL included, before filtering |
| `cwregression.py` | CWL windowed regression (`cw_regression`) |
| `metrics.py` | relative spectral distance, EO/EC band + channel weights, Jaccard overlap, topographic distance |

Order: `rp_weighting_matrix` -> `gradient_correct` -> band-pass 0.33-125 Hz (EEG and CWL) -> resample to 500 Hz -> `cw_regression` -> `metrics.evaluate` against the outside-scanner recording.

Not ported: Bergen's clipping/NaN handling and baseline correction, and the other
comparison methods (AAS/OBS/OBS-ICA). Tests: `python -m pytest python`.
