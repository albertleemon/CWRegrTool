# Handoff: state of the work (for a new session)

Repository: albertleemon/CWRegrTool, branch `claude/awesome-clarke-obslcf`.
User's fork of the Bergen toolbox: albertleemon/BergenToolboxModified (vendored unmodified files in `matlab/vendor/`).
Language with the user: Chinese.

## Goal
Reproduce the carbon-wire-loop (CWL) EEG-fMRI artifact correction of van der Meer et al.
(NeuroImage 2016; Data in Brief 2016). The work has a MATLAB and a Python version, it runs on the
user's real data, and it is evaluated against a recording made outside the magnet.

## Data (on the user's machine; not in the repo)
- `/media/bianca/DataPool/Meng/Neurofeedback/BIDS/sub-XXX/eeg/`
  - `run-01` = inside the scanner: R128 volume markers, TR 1.8 s, 30 slices (slice frequency 16.84 Hz), 5000 Hz, 72 channels.
  - `run-00` = outside the scanner.
- Channels:
  - ECG = Ch32.
  - Ch65/66 are unused: always exclude them.
  - CW6..CW1 = Ch67-72.
  - Cap: standard Brain Products 64-channel cap. The lower-edge electrodes lie outside the CWL coverage.
- Head motion: fMRIPrep `*_desc-confounds_timeseries.tsv`, used as `CWL_RP` / `opts.rp_file`.
- Subjects tested: sub-BOR5A and sub-kvhDR. The user runs MATLAB R2024b locally. The cloud session has no MATLAB, but Octave can be installed.

## Pipeline (`matlab/cwl_run_pipeline.m`; Python `python/run_pipeline.py`)
02 gradient AAS, rp-informed (Bergen `m_rp_info` / `m_moving_average`)
→ 02b optional OBS residual PCA (Niazy 2005; `obs_npc`, `obs_hp`)
→ 03 band-pass 0.33-125 Hz + resample to 500 Hz
→ 04 CWL regression (original code: 4 s Hann window, ±21 ms delays), optional low-pass of the estimated artifact (`cwl_lowpass`), trimmed to the corrected samples
→ 05 optional heartbeat-locked BCG AAS (`bcg_aas`, default OFF).

Recommended settings: `CWL_OBS=4 CWL_OBS_HP=60 CWL_LOWPASS=90 ./run_subject.sh RAW.vhdr OUT 65,66`, then `./evaluate_subject.sh OUT OUTSIDE_run-00.vhdr`.

## The MATLAB version is self-contained (latest change)
- `matlab/` needs only MATLAB + the Signal Processing Toolbox. It does not need EEGLAB.
- Own code: `cwl_load_brainvision`, `cwl_saveset` / `cwl_loadset`, and `cwl_resample` (the `pop_resample` algorithm).
- `cwl_cwregression` wraps the vendored `run_cwregression`.
- Verified in Octave against the old EEGLAB version on synthetic BrainVision data: all five stages are bitwise identical, and so are the events.
- `cwl_crosscheck.m` + `python/make_crosscheck_data.py` compare MATLAB with Python. All PASS in Octave.
- Octave caveats, for testing only: Octave lacks `filtfilt(sos,g,x)`, `import`, and MATLAB `pwelch` semantics. Use shims (sos2tf + filtfilt; a no-op `import`; flattened `+misc`) and do not `cd` into `matlab/`. The shims are not in the repo.

## Main findings so far
- **Gradient residual peaks:** the remaining peaks at 84.5 / 101 / 117.75 Hz are slice harmonics. They come from volume-to-volume changes in the artifact shape, not from timing jitter, clipping or motion. OBS 4 components at 60 Hz removes them.
- **CWL low-pass:** CWL raised the power above 90 Hz. Fixed by low-passing the estimated artifact; low-passing the regressors instead was tried, was ill-conditioned, and is abandoned.
- **Cz:** reduced only in sub-BOR5A. It is normal in sub-kvhDR, so the effect is subject-specific.
- **sub-kvhDR:**
  - 1-8 Hz power is about +7 dB above the outside recording.
  - Heartbeat-locked residual: 7.9 µV before CWL → 3.1 µV after CWL; the outside floor is 0.8 µV.
  - The largest residual is at the left-posterior temporal edge electrodes.
  - The locked residual explains only about 10% of the low-frequency excess.
- **Heart rate / pulse data:** R-peak detection in `python/bcg_aas.py` (`r_peaks`) is robust to the in-scanner T wave; it gives 62 bpm, which matches the PULS log. Pulse data (`.puls`) is not needed for processing.

## Open / next steps
1. The user should run `cwl_crosscheck` in real MATLAB (expect 8× PASS: template matrix, gradient, OBS fraction + waveform, R-peaks, BCG AAS, CWL regression interior, resampling; it errors on any FAIL; events/I-O/band-pass/full pipeline are not covered), then run sub-kvhDR with the standalone version and confirm that `04` matches the earlier EEGLAB-based result.
2. Run sub-kvhDR with `CWL_BCG_AAS=1` (optionally `CWL_BCG_HR=64`) and compare 05 vs 04 in `evaluation.txt` / `bcg_residual.json`. BCG AAS has not yet been run on real data.
3. Open questions for the user:
   - Is the resting state eyes open or closed? The alpha question was deferred: the user said to ignore it for now.
   - Which frequency band will the analysis use?
4. Not done: the paper's OBS / OBS-ICA BCG comparison methods, EEG–BOLD analysis, and the Python pipeline on real data.

## Docs
- `docs/PIPELINE.md`: logic.
- `docs/RUN_LOCALLY.md`: how to run, parameter table.
- `docs/IMPROVEMENTS.md`: per-step comparison with the paper.
- `matlab/README.md`: the standalone MATLAB package.
- `docs/EEG_fMRI_CWL_v2.2.md`: the user's revised GUI SOP.
