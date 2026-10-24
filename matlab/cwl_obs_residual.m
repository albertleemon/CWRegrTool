function [y, frac] = cwl_obs_residual(x, fs, idx, npc, hp, exclude_factor)
% CWL_OBS_RESIDUAL  Optimal-basis-set removal of the gradient residual
% (Niazy et al. 2005, NeuroImage 28:720-737; the OBS step of FASTR).
%
%   [y, frac] = cwl_obs_residual(x, fs, idx, npc, hp)
%
%   x    one channel (row vector) AFTER the gradient template subtraction
%   idx  nvol x L sample indices of the volume epochs
%   npc  number of principal components to remove (FASTR uses about 4)
%   hp   high-pass (Hz) applied to the residual before the PCA (default 70): only the
%        high-passed residual is fitted, so content below ~hp Hz (EEG rhythms) is kept
%   exclude_factor  volumes whose residual power exceeds this x the median (e.g. an
%        incomplete last volume) do not enter the PCA but are still corrected (default 100)
%
%   frac = fraction of the high-passed residual variance that was removed.
%   Not part of van der Meer et al. (2016).

if nargin < 5 || isempty(hp), hp = 70; end
if nargin < 6 || isempty(exclude_factor), exclude_factor = 100; end
y = x;  frac = 0;
if npc <= 0, return; end
[z, p, k] = butter(4, hp / (fs / 2), 'high');
[sos, g]  = zp2sos(z, p, k);
xh = cwl_sosfiltfilt(sos, g, x(:)).';
R  = xh(idx);                                         % nvol x L
pw = mean(R.^2, 2);
use = pw <= exclude_factor * median(pw);
Ru = R(use, :);
[~, ~, V] = svd(Ru - mean(Ru, 1), 'econ');
B   = V(:, 1:min(npc, size(V, 2))).';                % npc x L, orthonormal rows
fit = (R * B.') * B;
y(idx) = x(idx) - fit;
frac = 1 - sum((R(:) - fit(:)).^2) / sum(R(:).^2);
end
