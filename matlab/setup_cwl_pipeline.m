function setup_cwl_pipeline(toolbox_dir)
% SETUP_CWL_PIPELINE  Download EEGLAB + bva-io (once) and set up the MATLAB path.
%
%   setup_cwl_pipeline                 % toolboxes go to <repo>/matlab/toolboxes
%   setup_cwl_pipeline('D:\toolboxes') % or any folder
%
% Needs git on the system path for the first run (or clone the two repositories
% yourself into <toolbox_dir>/eeglab and <toolbox_dir>/eeglab/plugins/bva-io).
% No MATLAB toolbox is required for the pipeline itself.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);                              % repo root = CWRegrTool
if nargin < 1, toolbox_dir = fullfile(here, 'toolboxes'); end
if ~exist(toolbox_dir, 'dir'), mkdir(toolbox_dir); end

eeglab_dir = fullfile(toolbox_dir, 'eeglab');
if ~exist(eeglab_dir, 'dir')
    run_git(['clone --depth 1 https://github.com/sccn/eeglab.git "' eeglab_dir '"']);
end
bva_dir = fullfile(eeglab_dir, 'plugins', 'bva-io');
if ~exist(bva_dir, 'dir')
    run_git(['clone --depth 1 https://github.com/sccn/bva-io.git "' bva_dir '"']);
end

addpath(eeglab_dir);
disable_eeglab_online_check(eeglab_dir);   % no web calls at startup (can crash some MATLAB installs)
eeglab nogui;
addpath(bva_dir);                                    % no genpath: it would add the .git folders
addpath(root);                                       % CWRegrTool (+tools, +misc, pop_cwregression)
addpath(here);                                       % this pipeline
addpath(fullfile(here, 'vendor', 'BergenToolboxModified'));
heal_matlab_path();
assert(license('test', 'Signal_Toolbox') ~= 0, 'Signal Processing Toolbox required');
fprintf('CWL pipeline ready. Next: EEG = cwl_run_pipeline(''<path>/NB.vhdr'', opts);\n');
end

function run_git(args)
[st, msg] = system(['git ' args]);
if st ~= 0, error('git failed: %s\nInstall git or clone manually.', msg); end
end

function disable_eeglab_online_check(eeglab_dir)
% EEGLAB queries sccn.ucsd.edu for plugin versions at every start. On some systems that
% web call corrupts the MATLAB heap ("malloc(): invalid size") and kills the session.
% Set option_checkversion = 0 in ~/eeg_options.m (a .bak copy is kept).
opt = fullfile(getenv('HOME'), 'eeg_options.m');
if ~exist(opt, 'file')
    copyfile(fullfile(eeglab_dir, 'functions', 'adminfunc', 'eeg_optionsbackup.m'), opt);
end
t  = fileread(opt);
t2 = regexprep(t, 'option_checkversion\s*=\s*1', 'option_checkversion = 0');
if ~strcmp(t, t2)
    copyfile(opt, [opt '.bak']);
    fid = fopen(opt, 'w'); fwrite(fid, t2, 'char'); fclose(fid);
    fprintf('EEGLAB online version check disabled in %s (backup: %s.bak)\n', opt, opt);
end
end

function heal_matlab_path()
% A saved pathdef.m can lack MATLAB's own toolbox folders. On this kind of install
% toolbox/shared/signalwavelet was missing, which breaks butter/filtfilt with
% "Unable to resolve signalwavelet.internal.isodd". Add it back for this session only.
if isempty(which('signalwavelet.internal.isodd'))
    sw = fullfile(matlabroot, 'toolbox', 'shared', 'signalwavelet');
    if exist(sw, 'dir')
        addpath(genpath(sw), '-end');           % genpath skips +package/private folders
        fprintf('added missing MATLAB folders under %s\n', sw);
    end
end
try
    butter(2, 0.1, 'high');
catch e
    error('cwl:signal', ['butter() still fails after path repair: %s\n' ...
          'Run debug_cwl_env to see which function is missing.'], e.message);
end
end
