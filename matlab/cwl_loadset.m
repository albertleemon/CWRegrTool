function EEG = cwl_loadset(set_file)
% CWL_LOADSET  Load an EEGLAB dataset (.set, data inside or in a .fdt next to it),
% without EEGLAB.
%
%   EEG = cwl_loadset(fullfile(out_dir, '03_filtered.set'))

S = load(set_file, '-mat');
if isfield(S, 'EEG'), EEG = S.EEG; else, EEG = S; end
if ischar(EEG.data)
    fid = fopen(fullfile(fileparts(set_file), EEG.data), 'r', 'ieee-le');
    assert(fid > 0, 'cannot open %s', EEG.data);
    EEG.data = fread(fid, [EEG.nbchan, EEG.pnts * EEG.trials], 'float32=>single');
    fclose(fid);
end
end
