function y = cwl_filter(x, fs, type, fc, order)
% CWL_FILTER  Zero-phase Butterworth filter in plain MATLAB (no toolboxes).
%
%   y = cwl_filter(x, fs, 'high', 1, 2)          % 2nd-order high-pass at 1 Hz
%   y = cwl_filter(x, fs, 'low',  125, 3)
%   y = cwl_filter(x, fs, 'band', [0.33 125], 3) % high-pass(order) then low-pass(order)
%
%   x is channels x samples. The filter is run forward and backward (zero phase,
%   so the effective order is doubled) with odd-extension padding and
%   steady-state initial conditions, like filtfilt. Replaces butter/zp2sos/filtfilt,
%   which need the Signal Processing Toolbox.

if strcmp(type, 'band')
    y = cwl_filter(cwl_filter(x, fs, 'high', fc(1), order), fs, 'low', fc(2), order);
    return
end
secs = design(type, fc, fs, order);
y = zeros(size(x), 'double');
for c = 1:size(x, 1)
    y(c, :) = zerophase(double(x(c, :)), secs, round(3 * fs / fc));
end
end

function secs = design(type, fc, fs, n)
wc = 2 * fs * tan(pi * fc / fs);                     % pre-warped analog cutoff
k  = 1:n;
p  = exp(1i * pi * (2 * k + n - 1) / (2 * n));       % unit Butterworth poles
if strcmp(type, 'low'), s = wc * p; z0 = -1; else, s = wc ./ p; z0 = 1; end
zp = (2 * fs + s) ./ (2 * fs - s);                   % bilinear transform
secs = {};
for q = zp(imag(zp) > 1e-12)                         % one pole of each conjugate pair
    secs{end+1} = norm_section(poly([z0 z0]), real(poly([q conj(q)])), type); %#ok<AGROW>
end
for q = zp(abs(imag(zp)) <= 1e-12)                   % real pole (odd order)
    secs{end+1} = norm_section([1 -z0], [1 -real(q)], type); %#ok<AGROW>
end
end

function sec = norm_section(b, a, type)
if strcmp(type, 'low'), zr = 1; else, zr = -1; end   % unit gain at DC / Nyquist
e = zr .^ (0:numel(b)-1);
g = sum(a .* zr .^ (0:numel(a)-1)) / sum(b .* e);
sec = {b * g, a};
end

function y = zerophase(x, secs, padlen)
n = numel(x);
padlen = min(padlen, n - 1);
xe = [2*x(1) - x(padlen+1:-1:2), x, 2*x(end) - x(end-1:-1:end-padlen)];
y = run(xe, secs);
y = fliplr(run(fliplr(y), secs));
y = y(padlen+1 : padlen+n);
end

function y = run(x, secs)
y = x;
for i = 1:numel(secs)
    b = secs{i}{1};  a = secs{i}{2};
    y = filter(b, a, y, steady_state(b, a) * y(1));
end
end

function zi = steady_state(b, a)                      % filter state for a constant input of 1
m = numel(a) - 1;
if m == 1
    C = -a(2);
else
    C = [-a(2:end); eye(m-1, m)];
end
zi = (eye(m) - C.') \ (b(2:end).' - a(2:end).' * b(1));
end
