function defs = CollinsEtAl_globals(defs)

% Return global parameters for Collins et al 2012 analyses
% 
%   defs = CollinsEtAl_globals(defs)
% 
% This script returns global parameters 
% 
% PJ 2011.10.27 
% FB 2012.07.10 adapted for JLMT release

%% set paths
defs.paths.install_root = fileparts(which('jlmt_proc_series'));
defs.paths.data_root = fullfile(defs.paths.install_root,'data');
defs.paths.project_root = fullfile(defs.paths.install_root,'test','CollinsEtAl2012');
defs.paths.analysis_path = fullfile(defs.paths.project_root, 'analyses');
defs.paths.log_path = fullfile(defs.paths.project_root, 'logs');
defs.paths.fig_path = fullfile(defs.paths.project_root, 'figs');
defs.paths.matpath = fullfile(defs.paths.project_root,'matfiles');
defs.paths.tablepath = fullfile(defs.paths.project_root,'tables');

%% JLMT parameters
pparams = jlmt_preproc_params('modulation');

% Integration times that we want to use for this project
pparams.inHalfDecayTimes = [0.1 0.2 0.5 1.5 2 4];
names = calc_context('calc_names',pparams.inHalfDecayTimes);

% torus projection weight matrix path
defs.jlmt.paths.map_fname = fullfile(defs.paths.data_root,...
    'map_10-Dec-2006_16:18.mat');

% pitch class projection weight matrix path
defs.jlmt.paths.pcmap_fname = fullfile(defs.paths.data_root,...
    'pp2pitchclass_nnet_full_20120418T220430.mat');

% Steps to execute and save
defs.jlmt.glob.process = {'ani','pp','li','toract','pc','li','toract'};
defs.jlmt.glob.save_calc = defs.jlmt.glob.process;

% Parameters for individual steps
defs.jlmt.ani = params_ani('PlotFlag',0, ...
    'DownSampleFactor', pparams.downsample_factor, ...
    'NumOfChannels', pparams.nchan_ani); 

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
    'som',struct('fname',defs.jlmt.paths.pcmap_fname),...
    'spher_harm',struct('nharm_theta',3,'nharm_phi',4,'min_rsqr',0.95),...
    'prev_steps',{'ani','pp','li','toract','pc','li'});

% Attach dataset info
defs.datasets = tonmodcomp_datasets;


return

