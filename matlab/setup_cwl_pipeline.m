function setup_cwl_pipeline()
% SETUP_CWL_PIPELINE  Put this self-contained pipeline on the MATLAB path.
%
%   cd CWRegrTool/matlab ; setup_cwl_pipeline
%
% Needs only MATLAB and the Signal Processing Toolbox. Everything else is in this folder:
%   .                     the pipeline (cwl_*.m), BrainVision reader, .set reader/writer
%   vendor/CWRegrTool     CWL regression core (van der Meer et al. 2016), MIT
%   vendor/BergenToolboxModified  head-motion template selection (Moosmann 2009), GPL v2
%   vendor/firfilt        anti-aliasing filter design (EEGLAB firfilt), GPL
% No EEGLAB, no download, no internet. The folder can be copied anywhere.

here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(here, 'vendor', 'CWRegrTool'));             % no genpath: +misc/+tools are packages
addpath(fullfile(here, 'vendor', 'BergenToolboxModified'));
addpath(fullfile(here, 'vendor', 'firfilt'));
check_signal_toolbox();
others = setdiff(which('-all', 'run_cwregression'), {fullfile(here, 'vendor', 'CWRegrTool', 'run_cwregression.m')});
if ~isempty(others)
    fprintf('note: other copies on the path (unused, this folder comes first): %s\n', strjoin(others, ', '));
end
fprintf('CWL pipeline ready (%s). Next: EEG = cwl_run_pipeline(''<path>/run.vhdr'', opts);\n', here);
end

function check_signal_toolbox()
% Fail early, with a clear message, if the Signal Processing Toolbox does not work.
% A failure here means the MATLAB path itself is broken (see docs/RUN_LOCALLY.md).
if exist('OCTAVE_VERSION', 'builtin'), return; end      % Octave: only used for testing
assert(license('test', 'Signal_Toolbox') ~= 0, 'Signal Processing Toolbox required');
checks = {'butter',   @() butter(2, 0.1, 'high'); ...
          'filtfilt', @() filtfilt(1, 1, randn(100, 1)); ...
          'resample', @() resample(randn(100, 1), 1, 10); ...
          'pwelch',   @() pwelch(randn(1000, 1), 256, 128, 256, 500)};
bad = {};
for i = 1:size(checks, 1)
    try, checks{i, 2}(); catch e, bad{end+1} = sprintf('  %s: %s', checks{i, 1}, e.message); end %#ok<AGROW>
end
if ~isempty(bad)
    error('cwl:signal', 'Signal Processing Toolbox functions fail in %s (%s) -- MATLAB path broken?\n%s', ...
          version, matlabroot, strjoin(bad, newline));
end
end
