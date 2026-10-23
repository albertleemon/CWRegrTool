% Example for sub-BOR5A (72 channels, 5000 Hz, TR 1.8 s, R128 markers,
% CW6..CW1 = channels 67-72, ECG = channel 32, channels 65/66 unused).
% Run from the repository folder:  cd CWRegrTool/matlab ; example_BOR5A
setup_cwl_pipeline;

vhdr = '/media/bianca/DataPool/Meng/Neurofeedback/BIDS/sub-BOR5A/eeg/sub-BOR5A_task-rest_run-01_eeg.vhdr';
opts = struct();
opts.marker  = 'R128';
opts.exclude = {'65', '66'};
opts.out_dir = fullfile(getenv('HOME'), 'cwl_test', 'matlab_out');   % outside the BIDS folder
% opts.rp_file = '/path/to/rp_run01.txt';        % optional; default: moving average
EEG = cwl_run_pipeline(vhdr, opts);
