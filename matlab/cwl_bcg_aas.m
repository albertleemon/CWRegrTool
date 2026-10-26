function [out, pk, hr] = cwl_bcg_aas(eeg, ecg, fs, k, pre, hr_hint)
% CWL_BCG_AAS  Optional step after CWL regression: average artifact subtraction of the
% remaining heartbeat-locked BCG (Allen et al. 1998, NeuroImage 8:229-239). Not part of
% van der Meer et al. (2016), where AAS was a competing method, not a complement.
%
%   [out, pk, hr] = cwl_bcg_aas(eeg, ecg, fs)            % k = 21 beats, pre = 0.25 s
%   [out, pk, hr] = cwl_bcg_aas(eeg, ecg, fs, k, pre, hr_hint)
%
%   eeg  nchan x nsamples (EEG channels only), ecg 1 x nsamples, fs in Hz
%   Each beat owns the samples from `pre` s before its R peak to `pre` s before the next
%   one; the template (mean of the k nearest beats) is subtracted over one median RR.
%   Removes everything locked to the heart beat, including the heartbeat-evoked
%   potential: do not use for heartbeat-evoked-potential studies.
%   Same algorithm as python/bcg_aas.py.

if nargin < 4 || isempty(k),   k = 21;    end
if nargin < 5 || isempty(pre), pre = 0.25; end
if nargin < 6, hr_hint = []; end
pk = cwl_r_peaks(double(ecg(:)).', fs, hr_hint);
rr = diff(pk);
L  = round(median(rr));
a  = round(pre * fs);
pk = pk(pk - a >= 1 & pk - a + L - 1 <= size(eeg, 2));
n  = numel(pk);
eeg = double(eeg);
E  = zeros(size(eeg, 1), n, L);
for i = 1:n
    E(:, i, :) = reshape(eeg(:, pk(i) - a : pk(i) - a + L - 1), [], 1, L);
end
out = eeg;
for i = 1:n
    lo   = max(1, min(i - floor(k / 2), n - k + 1));
    tmpl = squeeze(mean(E(:, lo:min(lo + k - 1, n), :), 2));
    if size(eeg, 1) == 1, tmpl = tmpl(:).'; end
    if i < n, m = min(L, pk(i + 1) - pk(i)); else, m = L; end   % up to `pre` before the next R
    s0 = pk(i) - a;
    out(:, s0 : s0 + m - 1) = out(:, s0 : s0 + m - 1) - tmpl(:, 1:m);
end
hr = 60 * fs / median(rr);
end

function pk = cwl_r_peaks(ecg, fs, hr_hint)
% R peaks robust to the large in-scanner T wave (see python/bcg_aas.py r_peaks)
[z, p, kk] = butter(3, [10 35] / (fs / 2), 'bandpass');
[sos, g] = zp2sos(z, p, kk);
y = cwl_sosfiltfilt(sos, g, ecg(:)).';
d = [0 diff(y)].^2;
w = round(0.08 * fs);
e = conv(d, ones(1, w) / w, 'same');
if ~isempty(hr_hint) && hr_hint > 0
    period = 60 / hr_hint;
else
    seg = e(1:min(numel(e), round(120 * fs)));  seg = seg - mean(seg);
    ac  = xcorr(seg);  ac = ac(numel(seg):end);
    lo = round(0.33 * fs); hi = round(2.0 * fs);
    [~, im] = max(ac(lo + 1:hi + 1));
    period = (lo + im - 1) / fs;
end
es = sort(e);  q = 0.99 * (numel(es) - 1);          % 99th percentile as numpy (no Statistics Toolbox)
p99 = es(floor(q) + 1) + (q - floor(q)) * (es(min(floor(q) + 2, numel(es))) - es(floor(q) + 1));
[~, pk] = findpeaks(e, 'MinPeakHeight', 0.3 * p99, ...
                    'MinPeakDistance', round(0.7 * period * fs));
h = round(0.05 * fs);
pk = pk(pk - h >= 1 & pk + h <= numel(y));
for i = 1:numel(pk)
    [~, j] = max(abs(y(pk(i) - h : pk(i) + h)));
    pk(i) = pk(i) - h + j - 1;
end
end
