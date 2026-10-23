function s = cwl_diagnose(out_dir)
% CWL_DIAGNOSE  Before/after CWL-regression check for a cwl_run_pipeline output folder.
%
%   s = cwl_diagnose('~/cwl_test/matlab_crop')
%
% Compares 03_filtered.set (before CWL regression) with 04_cwl_corrected.set on the same
% samples and writes diagnostics.png and diagnostics.json into out_dir (same content as
% python/diagnose.py): mean EEG spectrum before/after, CWL spectra, a 10 s example
% trace (Oz), % variance removed per channel and band-power ratios, plus a null
% control: the same regression with phase-randomised CWL (over-fitting floor).

pre  = pop_loadset('filename', '03_filtered.set', 'filepath', out_dir);
post = pop_loadset('filename', '04_cwl_corrected.set', 'filepath', out_dir);
info = post.etc.cwl_pipeline;
fs   = post.srate;
lab  = {post.chanlocs.labels};
eeg  = find(ismember(lab, info.eeg_channels));
cwl  = find(ismember(lab, info.cwl_channels));
seg  = info.trimmed_to(1):info.trimmed_to(2);          % same samples in 03_filtered
x0   = double(pre.data(eeg, seg));
x1   = double(post.data(eeg, :));
assert(isequal(size(x0), size(x1)), 'before/after lengths differ');

nper = round(4 * fs);
win  = 0.5 - 0.5 * cos(2 * pi * (0:nper-1)' / nper);
[p0, f] = pwelch(x0.', win, nper/2, nper, fs);
p1 = pwelch(x1.', win, nper/2, nper, fs);
pc = pwelch(double(post.data(cwl, :)).', win, nper/2, nper, fs);
vr = 100 * (1 - var(x1, 0, 2) ./ var(x0, 0, 2));

% --- null control: same regression with phase-randomised CWL ------------------------
% The surrogate CWL keep their spectra and their mutual correlations but lose any
% temporal relation to the EEG. Whatever the regression removes with them is the
% over-fitting floor (126 regressors per window) and must be well below the real run.
rng(1);
c   = double(pre.data(cwl, :)).';
nfr = size(c, 1);  half = floor((nfr - 1) / 2);
ph  = exp(1i * 2 * pi * rand(half, 1));
F   = fft(c);
F(2:half+1, :)       = F(2:half+1, :) .* ph;
F(end-half+1:end, :) = F(end-half+1:end, :) .* conj(flipud(ph));
nul = pre;  nul.data(cwl, :) = real(ifft(F)).';
nul = pop_cwregression(nul, fs, 4, 0.021, 1, 'hann', cwl, eeg, 'taperedhann', 0);
xn  = double(nul.data(eeg, seg));
pn  = pwelch(xn.', win, nper/2, nper, fs);
vrn = 100 * (1 - var(xn, 0, 2) ./ var(x0, 0, 2));


fig = figure('Visible', 'off', 'Position', [100 100 1300 800], 'Renderer', 'painters');
subplot(2, 2, 1); semilogy(f, mean(p0, 2), f, mean(p1, 2), f, mean(pn, 2), '--'); xlim([0 125]);
xlabel('Hz'); ylabel('\muV^2/Hz'); title('EEG mean spectrum');
legend('before CWL', 'after CWL', 'null (random-phase CWL)');
subplot(2, 2, 2); semilogy(f, pc); xlim([0 125]); xlabel('Hz'); title('CWL spectra');
legend(info.cwl_channels, 'Location', 'northeast');
ch = find(strcmp(lab(eeg), 'Oz'), 1); if isempty(ch), ch = 1; end
n  = min(round(10 * fs), size(x0, 2));  t = (0:n-1) / fs;
subplot(2, 2, 3); plot(t, x0(ch, 1:n), t, x1(ch, 1:n));
xlabel('s'); ylabel('\muV'); title([lab{eeg(ch)} ', first 10 s']); legend('before', 'after');
subplot(2, 2, 4); bar([vr vrn]); xlabel('EEG channel #'); ylabel('% variance removed');
legend('CWL', 'null');
title('variance removed by CWL regression');
% painters renderer + print: saveas gives a blank image under -nodisplay on some systems
set(fig, 'Renderer', 'painters');
print(fig, fullfile(out_dir, 'diagnostics.png'), '-dpng', '-r110');  close(fig);
% the plotted data, for python/plot_diagnostics_mat.py (MATLAB figures can come out
% blank under -nodisplay)
tr = [x0(ch, 1:n); x1(ch, 1:n)];  oz = lab{eeg(ch)};  cwl_names = info.cwl_channels;
save(fullfile(out_dir, 'diagnostics_data.mat'), 'f', 'p0', 'p1', 'pn', 'pc', 'vr', 'vrn', ...
     't', 'tr', 'oz', 'cwl_names', '-v7');

bp = @(p, lo, hi) mean(mean(p(f >= lo & f <= hi, :)));

bands = struct('delta_theta', [1 7], 'alpha', [8 13], 'beta', [14 30], 'gamma', [30 80]);
r = struct();
for b = fieldnames(bands)'
    r.(b{1}) = bp(p1, bands.(b{1})(1), bands.(b{1})(2)) / bp(p0, bands.(b{1})(1), bands.(b{1})(2));
end
rn = struct();
for b = fieldnames(bands)'
    rn.(b{1}) = bp(pn, bands.(b{1})(1), bands.(b{1})(2)) / bp(p0, bands.(b{1})(1), bands.(b{1})(2));
end
s = struct('mean_var_removed_pct', mean(vr), 'power_ratio_after_before', r, ...
    'null_mean_var_removed_pct', mean(vrn), 'null_power_ratio_after_before', rn, ...
    'rms_before_uV', sqrt(mean(x0(:).^2)), 'rms_after_uV', sqrt(mean(x1(:).^2)), ...
    'method', info.method, 'n_volumes', info.n_volumes, 'seconds_analysed', size(x1, 2) / fs);
fid = fopen(fullfile(out_dir, 'diagnostics.json'), 'w');
fwrite(fid, jsonencode(s, 'PrettyPrint', true), 'char');  fclose(fid);
disp(s); disp(r); disp('null control (phase-randomised CWL):'); disp(rn);
fprintf('wrote %s and diagnostics.json\n', fullfile(out_dir, 'diagnostics.png'));
end
