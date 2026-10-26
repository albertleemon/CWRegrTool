function [EEG, weights] = cwl_cwregression(EEG, window_s, delay_s, cwl_idx, eeg_idx)
% CWL_CWREGRESSION  Carbon-wire-loop regression of the original CWRegrTool, without
% EEGLAB: the computation of pop_cwregression(EEG, srate, window_s, delay_s, 1, 'hann',
% cwl_idx, eeg_idx, 'taperedhann', 0) (tapered Hann windows, overlap factor 1).
%
%   [EEG, weights] = cwl_cwregression(EEG, 4, 0.021, cwl_idx, eeg_idx)
%
% weights: the summed window weights of the tool per sample; 0 = not corrected.

[data, cfg] = run_cwregression(EEG.data, EEG.srate, window_s, delay_s, 1, @hann, ...
                               cwl_idx, eeg_idx, 'taperedhann', 0);
EEG.data(eeg_idx, :) = EEG.data(eeg_idx, :) - data.subtracted_data;
cfg.cwregression.taperingfunction = func2str(cfg.cwregression.taperingfunction);
EEG.cwregression.cfg = cfg;
weights = data.subtracted_data_weights(:).';
end
