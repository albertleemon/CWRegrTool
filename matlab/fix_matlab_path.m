function fix_matlab_path(mode)
% FIX_MATLAB_PATH  Rebuild a broken toolbox/local/pathdef.m from the installation.
%
%   fix_matlab_path            % dry run: report what would change, write nothing
%   fix_matlab_path('apply')   % back up pathdef.m, set the rebuilt path, savepath
%
% restoredefaultpath computes MATLAB's default path by running
% toolbox/local/getphlpaths.pl (Perl) over the product files in toolbox/local/path;
% if that fails it silently falls back to pathdef.m -- which is the broken file.
% This function runs the Perl script itself (MATLAB's bundled Perl first, then the
% system one), keeps the non-MathWorks folders you added (SPM, FSL, ...), drops this
% repository's folders, and saves the result.

if nargin < 1, mode = 'dryrun'; end
script = fullfile(matlabroot, 'toolbox', 'local', 'getphlpaths.pl');
perls  = {fullfile(matlabroot, 'sys', 'perl', 'glnxa64', 'bin', 'perl'), 'perl'};
def = {};
for i = 1:numel(perls)
    if i == 1 && ~exist(perls{i}, 'file'), continue; end
    [st, out] = system(sprintf('"%s" "%s" "%s"', perls{i}, script, matlabroot));
    fprintf('perl: %s -> status %d\n', perls{i}, st);
    if st == 0
        parts = regexp(strtrim(out), '[;:\n]', 'split');
        parts = parts(~cellfun('isempty', parts));
        parts = parts(startsWith(parts, matlabroot));
        if numel(parts) > 100, def = parts; break; end
    else
        fprintf('%s\n', out);
    end
end
assert(~isempty(def), 'getphlpaths.pl produced no usable path; paste the output above');

must = {'signalwavelet', fullfile('signal', 'sigtools'), fullfile('toolbox', 'eml')};
for k = 1:numel(must)
    fprintf('rebuilt path contains %-18s: %d folders\n', must{k}, sum(contains(def, must{k})));
end

repo = fileparts(fileparts(mfilename('fullpath')));
cur  = strsplit(path, pathsep);
user = cur(~startsWith(cur, matlabroot) & ~startsWith(cur, repo));   % SPM, FSL, ...
oldroot = cur(startsWith(cur, matlabroot));
fprintf('MathWorks folders: now %d, rebuilt %d (%d added, %d removed)\n', numel(oldroot), ...
        numel(def), numel(setdiff(def, oldroot)), numel(setdiff(oldroot, def)));
fprintf('kept %d of your own folders, e.g. %s\n', numel(user), strjoin(user(1:min(3, end)), ', '));

if ~strcmpi(mode, 'apply')
    fprintf('\nDry run only. Run fix_matlab_path(''apply'') to write pathdef.m.\n');
    return
end
pd = fullfile(matlabroot, 'toolbox', 'local', 'pathdef.m');
bak = fullfile(getenv('HOME'), sprintf('pathdef_backup_%s.m', datestr(now, 'yyyymmdd_HHMMSS')));
copyfile(pd, bak);
fprintf('backup: %s\n', bak);
path(strjoin([def, user], pathsep));
rehash toolboxcache
if savepath(pd) ~= 0, error('savepath failed: is %s writable?', pd); end
fprintf('saved %s. Restart MATLAB, then test:\n', pd);
fprintf('  butter(2,0.1,''high''); resample(randn(100,1),1,10); pwelch(randn(1000,1));\n');
end
