function debug_cwl_env()
% DEBUG_CWL_ENV  Print what MATLAB resolves for every function the pipeline needs
% and run a tiny test of each. Run after setup_cwl_pipeline:
%   matlab -nodisplay -nosplash -batch "setup_cwl_pipeline; debug_cwl_env" 2>&1 | tee env.log
fprintf('MATLAB %s | %s | matlabroot %s\n', version, computer, matlabroot);
fprintf('Signal Processing Toolbox licensed: %d\n', license('test', 'Signal_Toolbox'));
names = {'signalwavelet.internal.isodd','butter','filtfilt','zp2sos','pwelch', ...
         'pop_loadbv','pop_resample','pop_saveset','pop_select','eeg_checkset', ...
         'pop_cwregression','run_cwregression','m_taperedhann','hann','m_rp_info', ...
         'm_moving_average','cwl_gradient_correct','cwl_run_pipeline'};
for i = 1:numel(names)
    w = which('-all', names{i});
    if isempty(w), fprintf('  MISSING  %s\n', names{i});
    else, fprintf('  %-34s %s%s\n', names{i}, w{1}, ternary(numel(w) > 1, sprintf('  (+%d more: shadowing?)', numel(w)-1), ''));
    end
end
t = @(name, f) try_run(name, f);
t('butter high',  @() butter(2, 0.01, 'high'));
t('butter band',  @() butter(3, [0.0013 0.5], 'bandpass'));
t('zp2sos',       @() zp2sos(1, 0.5, 1));
t('filtfilt',     @() filtfilt(ones(1,2)/2, 1, randn(1000,1)));
t('pwelch',       @() pwelch(randn(4000,2), hann_periodic(512), 256, 512, 500));
t('hann(repo)',   @() hann(9));
t('resample',     @() resample(randn(1000,1), 1, 10));
t('corrcoef',     @() corrcoef(randn(100, 5)));
end

function try_run(name, f)
try, f(); fprintf('  ok       %s\n', name);
catch e, fprintf('  FAILED   %s -> %s\n', name, strtok(e.message, newline)); end
end
function w = hann_periodic(n), w = 0.5 - 0.5 * cos(2*pi*(0:n-1)'/n); end
function s = ternary(c, a, b), if c, s = a; else, s = b; end, end
