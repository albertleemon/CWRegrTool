function r = cwl_eoec_metrics(fs, ref_open, ref_closed, in_open, in_closed, posterior)
% CWL_EOEC_METRICS  Comparison metrics of van der Meer et al. (2016).
%
%   r = cwl_eoec_metrics(fs, ref_open, ref_closed, in_open, in_closed, posterior)
%
%   ref_*  nchan x nsamples EEG recorded OUTSIDE the scanner (eyes open/closed)
%   in_*   artifact-corrected EEG recorded INSIDE the scanner, same channel order
%   posterior  indices of posterior channels (Pz, PO3, POz, PO4, Oz, ...)
%
%   Returns r.band, r.channels, r.spectral_distance_closed (eq. 4),
%   r.jaccard (eq. 5) and r.topographic_distance (eq. 6).

nper = round(4 * fs);
[f, ro] = spec(ref_open,  fs, nper);  [~, rc] = spec(ref_closed, fs, nper);
[~, io] = spec(in_open,   fs, nper);  [~, ic] = spec(in_closed,  fs, nper);

% --- subject-specific EO/EC band from the reference data -------------------
inc = mean(rc(posterior, :) - ro(posterior, :), 1);
inc(f > 25) = -Inf;
[mx, pk] = max(inc);
ok = inc >= 0.10 * mx;
lo = pk; hi = pk;
while lo > 1 && ok(lo-1),           lo = lo - 1; end
while hi < numel(f) && ok(hi+1),    hi = hi + 1; end
band = [max(f(lo), 7.5), max(f(hi), 12.5)];   % literal reading of the paper's limits

% --- channel selection and weights ------------------------------------------
m    = f >= band(1) & f <= band(2);
incc = sum(rc(:, m) - ro(:, m), 2);
sel  = find(incc > 0.5 * (min(incc) + max(incc)));
w    = zeros(1, size(rc, 1));  w(sel) = incc(sel) / sum(incc(sel));

yc = w * ic; yo = w * io; rcw = w * rc; row = w * ro;

% --- Jaccard overlap of the EO/EC surfaces (eq. 5) --------------------------
ylo = min(yc, yo); yhi = max(yc, yo); rlo = min(rcw, row); rhi = max(rcw, row);
inter = sum(max(0, min(yhi(m), rhi(m)) - max(ylo(m), rlo(m))));
uni   = sum(yhi(m) - ylo(m)) + sum(rhi(m) - rlo(m)) - inter;

r.band     = band;
r.channels = sel(:).';
r.spectral_distance_closed = reldist(mean(ic, 1), mean(rc, 1));
r.jaccard  = inter / uni;
r.topographic_distance = reldist(sum(ic(:, m) - io(:, m), 2).', sum(rc(:, m) - ro(:, m), 2).');
end

function d = reldist(a, b)                    % eqs. 4 and 6
d = 100 * sqrt(sum((a(:) - b(:)).^2)) / sqrt(sum(b(:).^2));
end

function [f, P] = spec(x, fs, nper)           % Welch, periodic Hann, 50 % overlap
win = 0.5 - 0.5 * cos(2 * pi * (0:nper-1)' / nper);
[P, f] = pwelch(double(x).', win, round(nper / 2), nper, fs);
P = P.'; f = f(:).';
end
