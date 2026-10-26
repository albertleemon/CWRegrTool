function EEG = cwl_resample(EEG, freq)
% CWL_RESAMPLE  Anti-aliased resampling of continuous data, without EEGLAB.
% Same algorithm as EEGLAB's pop_resample (with the Signal Processing Toolbox; default
% fc = 0.9, df = 0.2): Kaiser-windowed sinc anti-aliasing filter designed with firfilt's
% firwsord/firws/windows (vendor/firfilt), DC padding at both ends, MATLAB resample()
% with that kernel; event latencies rescaled. Assumes no data breaks (boundary events)
% inside the data, which holds in cwl_run_pipeline (resampling happens on one segment).
%
%   EEG = cwl_resample(EEG, 500)
%
% pop_resample: Copyright (C) 2001 Arnaud Delorme, Salk Institute; BSD license (EEGLAB).

fc = 0.9;  df = 0.2;
[p, q] = rat(freq / EEG.srate, 1e-12);
nyq = 1 / max([p q]);
m = firwsord('kaiser', 2, df * nyq, 0.002);                 % anti-aliasing filter order
b = p * firws(m, fc * nyq, windows('kaiser', m + 1, 5));
nPad = ceil((m / 2) / q) * q;                               % multiple of q, for unpadding
fprintf('resampling %g -> %g Hz\n', EEG.srate, EEG.srate * p / q);
out = [];
for c = 1:EEG.nbchan
    x = double(EEG.data(c, :)).';
    y = resample([repmat(x(1), nPad, 1); x; repmat(x(end), nPad, 1)], p, q, b(:));
    y = y(nPad * p / q + 1 : end - nPad * p / q);
    if isempty(out), out = zeros(EEG.nbchan, numel(y)); end
    out(c, :) = y;
end
EEG.data  = out;
EEG.srate = EEG.srate * p / q;
EEG.pnts  = size(out, 2);
EEG.xmax  = EEG.xmin + (EEG.pnts - 1) / EEG.srate;
EEG.times = linspace(EEG.xmin * 1000, EEG.xmax * 1000, EEG.pnts);
for i = 1:numel(EEG.event)
    EEG.event(i).latency = (EEG.event(i).latency - 1) * p / q + 1;
    if isfield(EEG.event, 'duration') && ~isempty(EEG.event(i).duration)
        EEG.event(i).duration = EEG.event(i).duration * p / q;
    end
end
EEG.urevent = [];
end
