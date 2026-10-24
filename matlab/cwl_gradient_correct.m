function [out, obs_frac] = cwl_gradient_correct(data, fs, onsets, pre, post, W, highpass, corr_thresh, obs_npc, obs_hp)
% CWL_GRADIENT_CORRECT  Average-artifact subtraction with a weighting matrix.
%
%   out = cwl_gradient_correct(data, fs, onsets, pre, post, W)
%
%   Headless equivalent of BergenToolboxModified/CorrectionMatrix.m (no GUI,
%   no clipping/NaN handling, no baseline step).
%
%   data    nchan x nsamples (ALL channels, CWL and ECG included)
%   onsets  sample index (1-based) of every volume marker
%   pre,post artifact span relative to the marker, inclusive, in samples
%   W       nvol x nvol weighting matrix; row i selects the templates for volume i
%   highpass (default 1 Hz) templates are built from a high-passed copy
%   corr_thresh (default 2.5) templates with dissimilarity > thresh*median are dropped
%   obs_npc  (default 0 = off) FASTR residual PCA after the template subtraction
%            (Niazy et al. 2005, see cwl_obs_residual); not part of van der Meer 2016
%   obs_hp   (default 70 Hz) high-pass of the residual for that PCA
%
%   Samples outside the corrected spans are returned unchanged (they still
%   contain gradient artifacts) -- trim them afterwards.

if nargin < 7 || isempty(highpass),    highpass = 1;     end
if nargin < 8 || isempty(corr_thresh), corr_thresh = 2.5; end
if nargin < 9 || isempty(obs_npc),     obs_npc = 0;      end
if nargin < 10 || isempty(obs_hp),     obs_hp = 70;      end

data   = double(data);
onsets = round(onsets(:));
L      = post - pre + 1;
starts = onsets + pre;
if min(starts) < 1 || max(starts) + L - 1 > size(data, 2)
    error('cwl:span', 'artifact span exceeds the recording; trim the data first');
end
if ~isequal(size(W), [numel(onsets) numel(onsets)])
    error('cwl:W', 'weighting matrix must be nvol x nvol');
end

idx = starts + (0:L-1);                       % nvol x L sample indices
out = data;
obs_frac = zeros(size(data, 1), 1);
if highpass > 0
    [z, p, k] = butter(2, highpass / (fs / 2), 'high');
    [sos, g]  = zp2sos(z, p, k);
end

for c = 1:size(data, 1)
    x  = data(c, :);
    A  = x(idx);                              % nvol x L
    if highpass > 0
        xf = cwl_sosfiltfilt(sos, g, x(:)).';
        A2 = xf(idx);
    else
        A2 = A;
    end
    for i = 1:numel(onsets)
        sel = find(W(i, :));
        if numel(sel) > 2
            cc  = corrcoef(A2(sel, :).');
            dis = 1 - cc;  dis(isnan(dis)) = 0;
            m   = mean(dis, 1);
            bad = m > corr_thresh * median(m);
            if any(bad) && ~all(bad)
                sel = sel(~bad);
            end
        end
        out(c, idx(i, :)) = A(i, :) - mean(A2(sel, :), 1);
    end
    if obs_npc > 0
        [out(c, :), obs_frac(c)] = cwl_obs_residual(out(c, :), fs, idx, obs_npc, obs_hp);
    end
end
end
