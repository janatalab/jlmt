function defs = CollinsEtAl_globals(defs)

% Return global parameters for Collins et al 2012 analyses
% 
%   defs = CollinsEtAl_globals(defs)
% 
% This script returns global parameters 
% 
% PJ 2011.10.27 

ensemble_globals;

%% set paths
defs.paths.install_root = fileparts(which('jlmt_proc_series'));
defs.paths.data_root = fileparts(defs.paths.install_root,'data');
defs.paths.project_root = fullfile(defs.paths.install_root,'test','CollinsEtAl2012');
defs.paths.analysis_path = fullfile(defs.paths.project_root, 'analyses');
defs.paths.log_path = fullfile(defs.paths.project_root, 'logs');
defs.paths.fig_path = fullfile(defs.paths.project_root, 'figs');
defs.paths.matpath = fullfile(defs.paths.project_root,'matfiles');
defs.paths.tablepath = fullfile(defs.paths.project_root,'tables');

%% IPEM parameters
pparams = ipem_preproc_params('modulation');

% Integration times that we want to use for this project
pparams.inHalfDecayTimes = [0.1 0.2 0.5 1.5 2 4];
names = calc_context('calc_names',pparams.inHalfDecayTimes);

% torus projection weight matrix path
defs.ipem.paths.map_fname = fullfile(defs.paths.data_root,...
    'map_10-Dec-2006_16:18.mat');

% pitch class projection weight matrix path
defs.ipem.paths.pcmap_fname = fullfile(defs.paths.data_root,...
    'pp2pitchclass_nnet_full_20120418T220430.mat');

% Steps to execute and save
defs.ipem.glob.process = {'ani','pp','ci','pc','pc_ci','toract'};
defs.ipem.glob.save_calc = defs.ipem.glob.process;

% Parameters for individual steps
defs.ipem.ani = params_ani('PlotFlag',0, ...
    'DownSampleFactor', pparams.downsample_factor, ...
    'NumOfChannels', pparams.nchan_ani); 

defs.ipem.pp = params_pp('PlotFlag', 0, ...
    'LowFrequency', [], ...
    'FrameWidth', pparams.frame_width, ...
    'FrameStepSize', pparams.frame_stepsize, ...
    'Atten', 'ipem_squash_hf');

defs.ipem.ci = params_ci('PlotFlag', 0, ...
    'HalfDecayTimes', pparams.inHalfDecayTimes);

defs.ipem.pc = params_pc(...
    'nnet',struct('fname',defs.ipem.paths.pcmap_fname));

defs.ipem.pc_ci = params_ci('PlotFlag', 0, ...
    'HalfDecayTimes', pparams.inHalfDecayTimes);

defs.ipem.toract.som.fname = defs.ipem.paths.map_fname;
defs.ipem.toract.ci_siglist = names;
defs.ipem.toract.HalfDecayTimes = pparams.inHalfDecayTimes;
defs.ipem.toract.norm = 'x./repmat(sum(x),size(x,1),1)';
defs.ipem.toract.calc_spher_harm = 1;
defs.ipem.toract.spher_harm.nharm_theta = 3;
defs.ipem.toract.spher_harm.nharm_phi = 4;
defs.ipem.toract.spher_harm.min_rsqr = 0.95;
hdc = pparams.inHalfDecayTimes;
nhd = length(hdc);
tmetrics = {'toract_corr','toract_kldist'};
ntmetric = length(tmetrics);
for k=1:ntmetric
  ntc = 0;
  defs.ipem.toract.metrics.(tmetrics{k}).tc_pair = zeros(0,2);
  for l=1:nhd
    for m=l+1:nhd
      ntc = ntc+1;
      defs.ipem.toract.metrics.(tmetrics{k}).tc_pair(ntc,:) = [hdc(l) hdc(m)];
    end % for m=l+1:nhd
  end % for l=1:nhd
end % for k=1:nmetric
defs.ipem.toract.metrics.toract_entropy=1;

% Attach dataset info
defs.datasets = tonmodcomp_datasets;


return

