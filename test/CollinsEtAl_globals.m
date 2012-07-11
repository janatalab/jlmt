function defs = CollinsEtAl_globals(defs)

% Return global parameters for Collins et al 201X tonality model comparisons analyses
% 
%   defs = CollinsEtAl_globals(defs)
% 
% PJ 2011.10.27 
% TC 2011-201X
% FB 2012.07.10 adapted for JLMT release

%% set paths
defs.paths.install_root = fileparts(fileparts(which('jlmt_proc_series')));
defs.paths.data_root = fullfile(defs.paths.install_root,'data');
defs.paths.project_root = fullfile(defs.paths.install_root,'test','CollinsEtAl2012');
defs.paths.analysis_path = fullfile(defs.paths.project_root, 'analyses');
defs.paths.log_path = fullfile(defs.paths.project_root, 'logs');
defs.paths.fig_path = fullfile(defs.paths.project_root, 'figs');
defs.paths.matpath = fullfile(defs.paths.project_root,'matfiles');
defs.paths.tablepath = fullfile(defs.paths.project_root,'tables');

%% Closure parameters
defs.paths.closure_matrices.theoretical = fullfile(defs.paths.data_root,...
    'maps','closure_struct-20120703T111311.mat');
defs.paths.closure_matrices.empirical = fullfile(defs.paths.data_root,...
    'maps','closure_struct-20120710T153125.mat');
% Time window over which penultimate and final event correlations are
% calculated.
defs.closure_matrix.win = [100 300];

%% JLMT parameters
pparams = jlmt_preproc_params('modulation');

% Integration times that we want to use for this project
pparams.inHalfDecayTimes = [0.1 2 4];
names = calc_li('calc_names',pparams.inHalfDecayTimes);

% torus projection weight matrix path
defs.jlmt.paths.map_fname = fullfile(defs.paths.data_root,'maps',...
    'map_10-Dec-2006_16:18.mat');

% pitch class projection weight matrix path
defs.jlmt.paths.pcmap_fname = fullfile(defs.paths.data_root,'maps',...
    'pp2pitchclass_nnet_full_20120611T113342.mat');

% pitch class to torus weight matrix
defs.jlmt.paths.pc_to_torus_fname = fullfile(defs.paths.data_root,'maps',...
    'pc_ci2toract_map_12-Jun-2012_15:47.mat');

% Steps to execute and save
defs.jlmt.glob.process = {'ani','pp','li','toract','pc','li','toract','ani','rp'};
defs.jlmt.glob.save_calc = defs.jlmt.glob.process;

% Parameters for individual steps
defs.jlmt.ani = params_ani('PlotFlag',0, ...
    'DownSampleFactor', pparams.downsample_factor, ...
    'NumOfChannels', pparams.nchan_ani,'prev_steps',[]);
defs.jlmt.ani(2) = params_ani('PlotFlag',0, ...
    'DownSampleFactor', 5, ...
    'NumOfChannels', pparams.nchan_ani,...
    'inDataType','sig_st',...
    'prev_steps',{'ani','pp','li','toract','pc','li','toract'});

defs.jlmt.pp = params_pp('PlotFlag', 0, ...
    'LowFrequency', [], ...
    'FrameWidth', pparams.frame_width, ...
    'FrameStepSize', pparams.frame_stepsize, ...
    'Atten', 'ipem_squash_hf');

defs.jlmt.li = params_li('PlotFlag', 0, ...
    'HalfDecayTimes', pparams.inHalfDecayTimes);

defs.jlmt.pc = params_pc(...
    'nnet',struct('fname',defs.jlmt.paths.pcmap_fname),...
    'inDataType','pp');

defs.jlmt.toract = params_toract('li_siglist',names,...
    'HalfDecayTimes',pparams.inHalfDecayTimes,'calc_spher_harm',1,...
    'som',struct('fname',defs.jlmt.paths.map_fname),...
    'spher_harm',struct('nharm_theta',3,'nharm_phi',4,'min_rsqr',0.95),...
    'prev_steps',{'ani','pp','li'});
defs.jlmt.toract(2) = params_toract('li_siglist',names,...
    'HalfDecayTimes',pparams.inHalfDecayTimes,'calc_spher_harm',1,...
    'som',struct('fname',defs.jlmt.paths.pc_to_torus_fname),...
    'spher_harm',struct('nharm_theta',3,'nharm_phi',4,'min_rsqr',0.95),...
    'prev_steps',{'ani','pp','li','toract','pc','li'});

% set rhythm profiler parameters
defs.jlmt.rp = rp_paramGroups_v2('input_type','ani',...
    'gain_type','beta_distribution','Fs',100,...
    'param_group','reson_filterQSpacing_periodBasedDecay');
defs.jlmt.rp.perform.calcOnsetInfo = 1;
defs.jlmt.rp.onsetInfo.Fs = defs.jlmt.rp.Fs;
% When using resonator output to estimate the ontimes of audio events, this
% parameter sets a sliding window in which only the highest peak is
% retained. For each peak found, if there are other peaks within
% +/- peak_height_window seconds of that peak, all but the highest peak
% will be discarded. When peaks get discarded, this is an indication that
% events are happening too quickly to be processed individually.
defs.jlmt.rp.onsetInfo.peak_height_window = 0.05;
% As a proportion of the range of the input signal, the observed peak must
% be greater than this parameter in order to be output as an onset.
defs.jlmt.rp.onsetInfo.thresh = 0.2;
% Some ANI signals have exhibited large spikes at the beginning of the
% signal where there probably shouldn't have been one.
defs.jlmt.rp.onsetInfo.throwOutFirstOnset = 0;

return

