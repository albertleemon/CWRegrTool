# Python port of the van der Meer et al. (2016) pipeline

| file | content |
|---|---|
| `gradient_correction.py` | rp-informed average-artifact subtraction (Moosmann 2009; port of BergenToolboxModified). Run on all channels, CWL included, before filtering |
| `cwregression.py` | CWL windowed regression (`cw_regression`) |
| `metrics.py` | relative spectral distance, EO/EC band + channel weights, Jaccard overlap, topographic distance |

Order: `rp_weighting_matrix` -> `gradient_correct` -> band-pass 0.33-125 Hz (EEG and CWL) -> resample to 500 Hz -> `cw_regression` -> `metrics.evaluate` against the outside-scanner recording.

Not ported: Bergen's clipping/NaN handling and baseline correction, and the other
comparison methods (AAS/OBS/OBS-ICA). Tests: `python -m pytest python`.

## Real data
```
pip install -r python/requirements.txt
python python/run_pipeline.py NB.vhdr out/ --rp rp_run01.txt --exclude 65 66
```
Outputs `cwl_corrected_raw.fif`, `cwl_corrected.npz`, `pipeline_log.json`.
See `docs/PIPELINE.md` for the logic and a pre-flight checklist; the MATLAB version is in `matlab/`.


Full manual (Chinese): `docs/MANUAL_PYTHON.md` — pipeline, all options, evaluation/report scripts (`denoise_report.py` etc., which read the MATLAB `.set` outputs), tests, troubleshooting.
