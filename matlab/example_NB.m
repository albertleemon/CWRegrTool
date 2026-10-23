% Example: the NB recording (72 channels, 5000 Hz, TR 1.8 s, R128 volume markers,
% CW6..CW1 = channels 67-72, ECG = channel 32, channels 65/66 unused).
setup_cwl_pipeline;

opts = struct();
opts.marker   = 'R128';
opts.rp_file  = 'D:\EEG_data\NB\fmri_motion\rp_run01.txt';   % [] -> moving average
opts.exclude  = {'65', '66'};
opts.out_dir  = 'D:\EEG_data\NB\preprocessed';
EEG = cwl_run_pipeline('D:\EEG_data\NB\NB.vhdr', opts);
