function EEG = cwl_run_pipeline(vhdr, opts)
% CWL_RUN_PIPELINE  Offline EEG-fMRI preprocessing with carbon-wire loops.
%
%   BrainVision -> gradient-artifact correction (all channels incl. CWL)
%   -> band-pass 0.33-125 Hz (all channels) -> resample to 500 Hz
%   -> CWL windowed regression (EEG channels only) -> saved .set
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
%
%   Requires EEGLAB + bva-io (pop_loadbv), the Signal Processing Toolbox, and
%   this repository on the path (run setup_cwl_pipeline first).

if nargin < 2, opts = struct(); end
d = @(f, v) setdefault(opts, f, v);
[vdir, vname, vext] = fileparts(vhdr);
out_dir     = d('out_dir', fullfile(vdir, 'preprocessed'));
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
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% 0. import
[vdir, vfile] = fix_header(vdir, [vname vext], out_dir);
EEG = pop_loadbv(vdir, vfile);
EEG = eeg_checkset(EEG);
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
EEG.data = data(:, a:b);
EEG.pnts  = size(EEG.data, 2);
EEG.times = (0:EEG.pnts-1) / fs * 1000;  EEG.xmax = EEG.times(end) / 1000;
EEG.event = EEG.event([EEG.event.latency] >= a & [EEG.event.latency] <= b);
for i = 1:numel(EEG.event), EEG.event(i).latency = EEG.event(i).latency - a + 1; end
EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = pop_saveset(EEG, 'filename', '02_gradient_corrected.set', 'filepath', out_dir);

%% 2. identical zero-phase band-pass on every channel, then resample
[z, p, k] = butter(3, band / (fs / 2), 'bandpass');
[sos, g]  = zp2sos(z, p, k);
EEG.data  = filtfilt(sos, g, double(EEG.data).').';
EEG = pop_resample(EEG, target_fs);              % anti-aliased resample (Signal Processing Toolbox)
assert(~any(isnan(EEG.data(:))), 'NaN in data before CWL regression');
EEG = pop_saveset(EEG, 'filename', '03_filtered.set', 'filepath', out_dir);

%% 3. CWL regression (EEG channels only), original CWRegrTool code
before = EEG.data(eeg_idx, :);
EEG = pop_cwregression(EEG, EEG.srate, window_s, delay_s, 1, 'hann', ...
                       cwl_idx, eeg_idx, 'taperedhann', 0);
% the tool leaves the first ~half window and last ~1-3 windows uncorrected:
% keep only the samples it actually changed
changed = find(any(EEG.data(eeg_idx, :) ~= before, 1));
fprintf('CWL regression changed samples %d..%d of %d\n', changed(1), changed(end), EEG.pnts);
EEG = pop_select(EEG, 'point', [changed(1) changed(end)]);

EEG.etc.cwl_pipeline = struct('method', method, 'n_volumes', nvol, 'tr_samples', TR, ...
    'span', span, 'band', band, 'eeg_channels', {labels(eeg_idx)}, ...
    'cwl_channels', {labels(cwl_idx)}, 'trimmed_to', [changed(1) changed(end)]);
EEG = pop_saveset(EEG, 'filename', '04_cwl_corrected.set', 'filepath', out_dir);
end

function v = setdefault(s, f, v)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); end
end

function [vdir, vfile] = fix_header(vdir, vfile, out_dir)
% If the header names files that do not exist (e.g. DataFile=RS1.eeg after the
% recording was renamed), build a patched copy in out_dir/_patched that points to the
% sibling .eeg/.vmrk with the same stem (the .eeg is symlinked on Linux/macOS, copied on Windows).
h = fileread(fullfile(vdir, vfile));
data = strtrim(regexp(h, 'DataFile=([^\r\n]+)', 'tokens', 'once'));
mrk  = strtrim(regexp(h, 'MarkerFile=([^\r\n]+)', 'tokens', 'once'));
if exist(fullfile(vdir, data{1}), 'file') && exist(fullfile(vdir, mrk{1}), 'file'), return; end
[~, stem] = fileparts(vfile);
eeg = fullfile(vdir, [stem '.eeg']);  vmrk = fullfile(vdir, [stem '.vmrk']);
assert(exist(eeg, 'file') && exist(vmrk, 'file'), ...
    'header points to %s / %s which are missing, and no %s.eeg/.vmrk next to it', data{1}, mrk{1}, stem);
tmp = fullfile(out_dir, '_patched');  if ~exist(tmp, 'dir'), mkdir(tmp); end
hh = regexprep(h, 'DataFile=[^\r\n]+', ['DataFile=' stem '.eeg']);
hh = regexprep(hh, 'MarkerFile=[^\r\n]+', ['MarkerFile=' stem '.vmrk']);
mm = regexprep(fileread(vmrk), 'DataFile=[^\r\n]+', ['DataFile=' stem '.eeg']);
write_text(fullfile(tmp, [stem '.vhdr']), hh);
write_text(fullfile(tmp, [stem '.vmrk']), mm);
target = fullfile(tmp, [stem '.eeg']);
if ~exist(target, 'file')
    if isunix, system(['ln -s "' eeg '" "' target '"']); else, copyfile(eeg, target); end
end
fprintf('header pointed to %s; using %s.eeg instead\n', data{1}, stem);
vdir = tmp;  vfile = [stem '.vhdr'];
end

function write_text(f, t)
fid = fopen(f, 'w', 'n', 'UTF-8'); fwrite(fid, t, 'char'); fclose(fid);
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
