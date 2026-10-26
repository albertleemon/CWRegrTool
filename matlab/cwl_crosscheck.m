function ok = cwl_crosscheck(mat_file)
% CWL_CROSSCHECK  Run the MATLAB gradient correction, OBS and BCG AAS on the synthetic
% data written by python/make_crosscheck_data.py and compare with the Python results.
%
%   (shell)  python python/make_crosscheck_data.py crosscheck_data.mat
%   (MATLAB) setup_cwl_pipeline; cwl_crosscheck('crosscheck_data.mat')
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
t(end+1, :) = {'R-peak detection', same && abs(hr - S.hr) < 0.5, ...
               sprintf('%d beats (python %d, true %d), HR %.1f', numel(pk), numel(S.pk), S.nbeats, hr)};
dy = max(abs(y(:) - S.y(:)));
t(end+1, :) = {'BCG AAS', same && dy < 0.5, sprintf('max diff %.3g uV', dy)};
for i = 1:size(t, 1)
    r = 'FAIL'; if t{i, 2}, r = 'PASS'; end
    fprintf('%-32s %s  %s\n', t{i, 1}, r, t{i, 3});
end
ok = all([t{:, 2}]);
end
