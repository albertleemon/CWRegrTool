% Called by run_headless.sh; reads CWL_VHDR, CWL_OUT, CWL_EXCLUDE from the environment.
opts = struct();
opts.out_dir = getenv('CWL_OUT');
ex = strsplit(getenv('CWL_EXCLUDE'), ',');
opts.exclude = ex(~cellfun('isempty', ex));
try
    EEG = cwl_run_pipeline(getenv('CWL_VHDR'), opts);
    disp(EEG.etc.cwl_pipeline);
catch err
    fprintf(2, '\nPIPELINE FAILED: %s\n', err.message);
    for k = 1:numel(err.stack)
        fprintf(2, '  at %s (line %d) in %s\n', err.stack(k).name, err.stack(k).line, err.stack(k).file);
    end
    rethrow(err);
end
