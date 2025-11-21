function cfg = example_config()
% EXAMPLE_CONFIG  Parameters for run_analysis

cfg.out_dir = './outputs';

% Data
cfg.traces_csv   = '/Users/peyshyuanchin/Documents/Project_AON/Htr2C_bfafcondition/103025_110625_110725_H1_celltraces.csv';  % file with Time_s__CellStatus, accepted, accepted_1...
cfg.dio_bf_path  = '/Users/peyshyuanchin/Documents/Project_AON/Htr2C_bfafcondition/H1_103025_GPIO.csv';
cfg.dio_af_path  = '/Users/peyshyuanchin/Documents/Project_AON/Htr2C_bfafcondition/H1_110625_GPIO.csv';
cfg.odors_bf_csv = '/Users/peyshyuanchin/Documents/Project_AON/Htr2C_bfafcondition/5stim_10 trials.csv';  % vector of labels (can be non-1..K)
cfg.odors_af_csv = '/Users/peyshyuanchin/Documents/Project_AON/Htr2C_bfafcondition/5stim_10 trials.csv';
cfg.dio_day3_path   = '/Users/peyshyuanchin/Documents/Project_AON/Htr2C_bfafcondition/110725_H1_GPIO.csv';     % DIO for day 3
cfg.odors_day3_csv  = '/Users/peyshyuanchin/Documents/Project_AON/Htr2C_bfafcondition/5stim_10 trials.csv';    % labels for day 3


cfg.drop_first_label_bf = true;   % check the number of trials of the odor stim file to see if it matches the times of odor presentations
cfg.drop_first_label_af = true;   % replicate your original "odorsaf(1)=[];" when needed
cfg.drop_first_label_day3 = true;                   
cfg.dio_threshold       = 1000;   % threshold to binarize DIO

% Session parsing
cfg.session_gap_threshold = 1.0;  % seconds; a negative jump < -1 means a new session

% Sampling / alignment windows
cfg.FR = 15;  % Hz
cfg.BS = 3;   % seconds before odor
cfg.RS = 10;  % seconds after odor

% Mineral oil & features
cfg.mo_index = 1;                 % index of MO in odor set, or [] if not used
cfg.decoder_window = [3 5];       % seconds after onset to export features (e.g., 2 s window)

% Reliability
cfg.num_splits = 10;
cfg.random_seed = 42;
% In example_config.m (recommended)
cfg.do_mo_subtraction   = true;   % perform MO subtraction step
cfg.reclassify_after_mo = true;   % recompute RC/NC/IC on MO-subtracted data
cfg.plot_after_mo_sets  = true;   % make plots for those recomputed sets
end
