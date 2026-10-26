function ok = cwl_crosscheck(mat_file)
% CWL_CROSSCHECK  Run the MATLAB gradient correction, OBS and BCG AAS on the synthetic
% data written by python/make_crosscheck_data.py and compare with the Python results.
%
%   (shell)  python python/make_crosscheck_data.py crosscheck_data.mat
%   (MATLAB) setup_cwl_pipeline; cwl_crosscheck('crosscheck_data.mat')
%
% Covered: template matrix, gradient correction, OBS (fraction AND output waveform), R-peaks,
% BCG AAS, CWL regression (original tool vs the Python port, where the tool applies full window weight),
% 5000 -> 500 Hz resampling. NOT covered: BrainVision/.set I/O, events, band-pass, the
% edge samples of the CWL regression (the tool leaves them uncorrected, Python corrects
% them) and the full pipeline end to end. Errors out if any check FAILs.
%
% Small differences (< 1% of the artifact) come from filtfilt edge padding, which
% differs between MATLAB and scipy.
S  = load(mat_file);
on = double(S.on(:));
W  = m_moving_average(size(S.W, 1), 25);
t  = {};
t(end+1, :) = {'Bergen moving-average matrix', isequal(W, S.W), ''};
g  = cwl_gradient_correct(S.x, 5000, on, 0, 2499, W);
dg = max(abs(g(:) - S.g(:)));
t(end+1, :) = {'gradient template subtraction', dg < 3, sprintf('max diff %.3g (artifact ~300)', dg)};
idx = on + (0:2499);
[o, fr] = cwl_obs_residual(g(1, :), 5000, idx, 4, 60);
t(end+1, :) = {'OBS residual removal', abs(fr - S.fr) < 0.01, sprintf('removed %.4f vs %.4f', fr, S.fr)};
[y, pk, hr] = cwl_bcg_aas(S.eeg, S.ecg, 500);
same = numel(pk) == numel(S.pk) && max(abs(pk(:) - double(S.pk(:)))) <= 2;
dobs = max(abs(o(:) - S.o(:)));
t(end+1, :) = {'OBS output waveform', dobs < 0.05, sprintf('max diff %.3g (residual ~%.3g)', dobs, max(abs(S.o(:))))};
t(end+1, :) = {'R-peak detection', same && abs(hr - S.hr) < 0.5, ...
               sprintf('%d beats (python %d, true %d), HR %.1f', numel(pk), numel(S.pk), S.nbeats, hr)};
dy = max(abs(y(:) - S.y(:)));
t(end+1, :) = {'BCG AAS', same && dy < 0.5, sprintf('max diff %.3g uV', dy)};
% CWL regression: original tool (via cwl_cwregression) vs Python port
E = struct('data', [S.cw_eeg; S.cw_cwl], 'srate', 500);
[E, wts] = cwl_cwregression(E, 4, 0.021, 4:5, 1:3);
v = wts > 0.99;  v(1:2001) = false;   % Python weights the first half window fully, the tool blends two windows
dc = E.data(1:3, v) - S.cw_out(:, v);
ref = S.cw_out(:, v);
rc = sqrt(mean(dc(:).^2)) / sqrt(mean(ref(:).^2));
t(end+1, :) = {'CWL regression (interior)', rc < 1e-4 && nnz(v) > 0.5 * numel(v), ...
               sprintf('relative RMS diff %.2g on %d of %d samples', rc, nnz(v), numel(v))};
% resampling 5000 -> 500 Hz (EEGLAB-style kernel vs scipy resample_poly), away from the edges
R = cwl_resample(struct('data', S.rs_x(:).', 'srate', 5000, 'nbchan', 1, 'pnts', numel(S.rs_x), ...
                        'xmin', 0, 'event', []), 500);
k = 51:numel(S.rs_y) - 50;                       % skip 0.1 s at both ends
rsy = S.rs_y(:).';
dr = max(abs(R.data(k) - rsy(k)));
t(end+1, :) = {'resample 5000->500 Hz', dr < 0.01, sprintf('max diff %.3g (signal amplitude ~1.7)', dr)};
for i = 1:size(t, 1)
    r = 'FAIL'; if t{i, 2}, r = 'PASS'; end
    fprintf('%-32s %s  %s\n', t{i, 1}, r, t{i, 3});
end
ok = all([t{:, 2}]);
if ~ok, error('cwl:crosscheck', 'MATLAB and Python disagree: see the FAIL lines above'); end
end
