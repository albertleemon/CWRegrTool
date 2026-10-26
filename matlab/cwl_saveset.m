function EEG = cwl_saveset(EEG, out_dir, fname)
% CWL_SAVESET  Save an EEG structure as an EEGLAB dataset (.set header + .fdt float32
% data, little endian), readable by EEGLAB's pop_loadset, MNE (read_raw_eeglab) and
% cwl_loadset, without EEGLAB. Like pop_saveset, the returned data are single precision.
%
%   EEG = cwl_saveset(EEG, out_dir, '03_filtered.set')

[~, stem] = fileparts(fname);
EEG.filename = fname;  EEG.filepath = out_dir;  EEG.datfile = [stem '.fdt'];
EEG.data   = single(reshape(EEG.data, EEG.nbchan, []));
EEG.saved  = 'yes';
fid = fopen(fullfile(out_dir, EEG.datfile), 'w', 'ieee-le');
assert(fid > 0, 'cannot write %s', fullfile(out_dir, EEG.datfile));
fwrite(fid, EEG.data, 'float32');
fclose(fid);
hdr = EEG;  hdr.data = EEG.datfile;
save(fullfile(out_dir, fname), '-v7', '-struct', 'hdr');
fprintf('saved %s\n', fullfile(out_dir, fname));
end
