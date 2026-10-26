function EEG = cwl_run_pipeline(vhdr, opts)
% CWL_RUN_PIPELINE  Offline EEG-fMRI preprocessing with carbon-wire loops.
%
%   BrainVision -> 02 gradient-artifact template subtraction (all channels incl. CWL)
%   -> 02b optional OBS removal of the gradient residual (opts.obs_npc > 0)
%   -> 03 band-pass 0.33-125 Hz (all channels) + resample to 500 Hz
%   -> 04 CWL windowed regression (EEG channels only); every stage is saved as .set
%
%   EEG = cwl_run_pipeline('D:\data\NB.vhdr', opts)
%
%   opts (all optional)
%     .out_dir      output folder (default: <vhdr folder>/preprocessed)
%     .marker       volume marker type          (default 'R128')
%     .rp_file      SPM rp_*.txt or fMRIPrep *_desc-confounds_timeseries.tsv;
%                   [] -> moving average templates
%     .cwl_regex    CWL channel labels          (default '^CW')
%     .ecg_label    ECG channel label           (default 'ECG')
%     .exclude      cellstr of labels that are neither EEG targets nor CWL
%                   (e.g. {'65','66'})
%     .n_templates  templates per volume        (default 25)
%     .threshold    head-motion threshold in mm (default 0.3)
%     .span         [pre post] samples relative to the marker (default [0 TR-1])
%     .band         band-pass in Hz             (default [0.33 125])
%     .target_fs    sampling rate for regression (default 500)
%     .window_s     regression window, s        (default 4)
%     .delay_s      delay embedding half-width  (default 0.021)
%     .obs_npc      FASTR residual PCA components after the gradient template
%                   (Niazy 2005), 0 = off (default), e.g. 4; not part of the paper
%     .obs_hp       high-pass (Hz) of the residual for the OBS fit (default 70; e.g. 60 to
%                   also treat the 67 Hz harmonic)
%     .bcg_aas      1 = after the CWL regression, subtract the remaining heartbeat-locked
%                   BCG by AAS (Allen 1998), saved as 05_bcg_aas.set (default 0 = off;
%                   not part of the paper; removes heartbeat-evoked potentials too)
%     .bcg_hr       expected heart rate (bpm) for the R-peak detection (optional)
%     .cwl_lowpass  low-pass (Hz) of the ESTIMATED CWL artifact before it is subtracted
%                   (default 0 = off, e.g. 90): above that frequency the EEG stays exactly
%                   as before the regression (no loop noise added, nothing removed)
%
%   Requires only MATLAB + the Signal Processing Toolbox; everything else is in this
%   folder (run setup_cwl_pipeline first to set the path). No EEGLAB needed: the
%   outputs are EEGLAB datasets (.set + .fdt) that EEGLAB and MNE can open.

if nargin < 2, opts = struct(); end
d = @(f, v) setdefault(opts, f, v);
out_dir     = d('out_dir', fullfile(fileparts(vhdr), 'preprocessed'));
marker      = d('marker', 'R128');
rp_file     = d('rp_file', []);
cwl_regex   = d('cwl_regex', '^CW');
ecg_label   = d('ecg_label', 'ECG');
exclude     = d('exclude', {});
n_templates = d('n_templates', 25);
threshold   = d('threshold', 0.3);
band        = d('band', [0.33 125]);
target_fs   = d('target_fs', 500);
window_s    = d('window_s', 4);
delay_s     = d('delay_s', 0.021);
obs_npc     = d('obs_npc', 0);
obs_hp      = d('obs_hp', 70);
cwl_lowpass = d('cwl_lowpass', 0);
bcg_aas     = d('bcg_aas', 0);
bcg_hr      = d('bcg_hr', []);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% 0. import
EEG = cwl_load_brainvision(vhdr);
fs  = EEG.srate;
labels = {EEG.chanlocs.labels};
cwl_idx = find(~cellfun('isempty', regexp(labels, cwl_regex, 'once')));
ecg_idx = find(strcmp(labels, ecg_label));
eeg_idx = setdiff(1:EEG.nbchan, [cwl_idx ecg_idx find(ismember(labels, exclude))]);
assert(~isempty(cwl_idx), 'no CWL channels match %s', cwl_regex);

types  = {EEG.event.type};
onsets = round([EEG.event(strcmp(types, marker)).latency]);
assert(numel(onsets) > 2, 'found %d "%s" markers', numel(onsets), marker);
TR = round(median(diff(onsets)));
span = d('span', [0 TR-1]);

%% 1. gradient artifact on ALL channels (CWL included)
nvol = numel(onsets);
if ~isempty(rp_file)
    rp_file = motion_to_rp(rp_file, out_dir, nvol, threshold);
    [~, W] = m_rp_info(rp_file, nvol, threshold, n_templates);
    method = 'rp_info';
else
    W = m_moving_average(nvol, n_templates);
    method = 'moving_average';
end
data = cwl_gradient_correct(EEG.data, fs, onsets, span(1), span(2), W);
a = onsets(1) + span(1);  b = onsets(end) + span(2);   % outside this: not corrected
EEG.data = data(:, a:b);  clear data
EEG.pnts  = size(EEG.data, 2);
EEG.times = (0:EEG.pnts-1) / fs * 1000;  EEG.xmax = EEG.times(end) / 1000;
EEG.event = EEG.event([EEG.event.latency] >= a & [EEG.event.latency] <= b);
for i = 1:numel(EEG.event), EEG.event(i).latency = EEG.event(i).latency - a + 1; end
EEG = cwl_saveset(EEG, out_dir, '02_gradient_corrected.set');

%% 1b. optional: OBS removal of the gradient residual (Niazy et al. 2005), all channels
obs_frac = [];
if obs_npc > 0
    idx = (onsets(:) - a + 1 + span(1)) + (0:span(2)-span(1));   % volume epochs in the cut data
    obs_frac = zeros(EEG.nbchan, 1);
    for c = 1:EEG.nbchan
        [y, obs_frac(c)] = cwl_obs_residual(double(EEG.data(c, :)), fs, idx, obs_npc, obs_hp);
        EEG.data(c, :) = y;
    end
    fprintf('OBS residual removal (%d components, >%g Hz): median %.0f%% of the residual removed\n', ...
            obs_npc, obs_hp, 100 * median(obs_frac));
    EEG = cwl_saveset(EEG, out_dir, '02b_obs_corrected.set');
end

%% 2. identical zero-phase band-pass on every channel, then resample
[z, p, k] = butter(3, band / (fs / 2), 'bandpass');
[sos, g]  = zp2sos(z, p, k);
EEG.data  = cwl_sosfiltfilt(sos, g, double(EEG.data).').';
EEG = cwl_resample(EEG, target_fs);              % anti-aliased, same as EEGLAB pop_resample
assert(~any(isnan(EEG.data(:))), 'NaN in data before CWL regression');
EEG = cwl_saveset(EEG, out_dir, '03_filtered.set');

%% 3. CWL regression (EEG channels only), original CWRegrTool code
before = EEG.data(eeg_idx, :);
EEG = cwl_cwregression(EEG, window_s, delay_s, cwl_idx, eeg_idx);
% the tool leaves the first ~half window and last ~1-3 windows uncorrected:
% keep only the samples it actually changed (determined before any artifact low-pass,
% which would smear tiny values into the uncorrected edges)
changed = find(any(EEG.data(eeg_idx, :) ~= before, 1));
if cwl_lowpass > 0
    % low-pass the artifact estimate, not the regressors: low-passing the regressors
    % makes the lagged fit ill-conditioned and amplifies the transition band
    [z, p, k] = butter(8, cwl_lowpass / (EEG.srate / 2), 'low');   % steep, zero phase
    [sos, g]  = zp2sos(z, p, k);
    art = double(before) - double(EEG.data(eeg_idx, :));
    art = cwl_sosfiltfilt(sos, g, art.').';
    EEG.data(eeg_idx, :) = double(before) - art;
    fprintf('CWL artifact estimate low-passed at %g Hz before subtraction\n', cwl_lowpass);
end
fprintf('CWL regression changed samples %d..%d of %d\n', changed(1), changed(end), EEG.pnts);
EEG = crop_points(EEG, changed(1), changed(end));

EEG.etc.cwl_pipeline = struct('method', method, 'obs_npc', obs_npc, 'obs_hp', obs_hp, 'obs_frac', obs_frac, ...
    'cwl_lowpass', cwl_lowpass, 'n_volumes', nvol, 'tr_samples', TR, ...
    'span', span, 'band', band, 'eeg_channels', {labels(eeg_idx)}, ...
    'cwl_channels', {labels(cwl_idx)}, 'trimmed_to', [changed(1) changed(end)]);
EEG = cwl_saveset(EEG, out_dir, '04_cwl_corrected.set');

%% 4. optional: heartbeat-locked AAS of the remaining BCG (Allen et al. 1998)
if bcg_aas
    assert(~isempty(ecg_idx), 'opts.bcg_aas needs an ECG channel (%s not found)', ecg_label);
    [y, pk, hr] = cwl_bcg_aas(EEG.data(eeg_idx, :), EEG.data(ecg_idx(1), :), EEG.srate, [], [], bcg_hr);
    EEG.data(eeg_idx, :) = y;
    fprintf('BCG AAS after CWL: %d beats, heart rate %.1f bpm\n', numel(pk), hr);
    EEG.etc.cwl_pipeline.bcg_aas = struct('n_beats', numel(pk), 'heart_rate_bpm', hr);
    EEG = cwl_saveset(EEG, out_dir, '05_bcg_aas.set');
end
end

function v = setdefault(s, f, v)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); end
end

function EEG = crop_points(EEG, a, b)
% keep samples a..b (like EEGLAB pop_select 'point'); boundary events mark the cuts
size_before = EEG.pnts;
EEG.data  = EEG.data(:, a:b);
EEG.pnts  = b - a + 1;
EEG.xmax  = EEG.xmin + (EEG.pnts - 1) / EEG.srate;
EEG.times = linspace(EEG.xmin * 1000, EEG.xmax * 1000, EEG.pnts);
n  = numel(EEG.event);
ev = EEG.event([EEG.event.latency] >= a & [EEG.event.latency] <= b);
for i = 1:numel(ev), ev(i).latency = ev(i).latency - a + 1; end
if n
    bnd = EEG.event(1);  bnd = structfun(@(v) [], bnd, 'UniformOutput', false);
    bnd.type = 'boundary';
    b0 = bnd;  b0.latency = 0.5;               b0.duration = a - 1;
    b1 = bnd;  b1.latency = EEG.pnts + 0.5;    b1.duration = size_before - b;
    ev = [b0 ev b1];
end
EEG.event = ev;
end

function rp_txt = motion_to_rp(f, out_dir, nvol, threshold)
% Accept an fMRIPrep confounds .tsv (trans_x..rot_z, mm/rad) as well as an SPM rp file;
% report how many frame-to-frame translations exceed the threshold.
if endsWith(f, '.tsv')
    T  = readtable(f, 'FileType', 'text', 'Delimiter', '\t', 'TreatAsMissing', 'n/a');
    rp = T{:, {'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'}};
    rp_txt = fullfile(out_dir, 'rp_from_confounds.txt');
    writematrix(rp, rp_txt, 'Delimiter', ' ', 'FileType', 'text');
else
    rp = load(f);  rp_txt = f;
end
if size(rp, 1) ~= nvol
    warning('cwl:rp', '%d motion rows but %d volume markers', size(rp, 1), nvol);
end
steps = sqrt(sum(diff(rp(:, 1:3)).^2, 2));
fprintf('head motion: %d volume steps > %g mm (max step %.3f mm)\n', ...
        sum(steps > threshold), threshold, max(steps));
if ~any(steps > threshold)
    fprintf('  -> no motion events: templates equal moving average; lower opts.threshold to use motion\n');
end
end
