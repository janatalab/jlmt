% BASIC JLMT UNIT TEST SCRIPT
% 
% This script tests jlmt_proc_series and JLMT analyses on a .wav file to
% verify that everything is working. 
%
% Copyright (c) 2007-2012 The Regents of the University of California
% All Rights Reserved.
%
% 13 June 2007 - Stefan Tomic
% 22 May 2012 - Fred Barrett

% Collect test file paths
jlmt_root = fileparts(fileparts(which('jlmt_proc_series')));
test_data_dir = fullfile(jlmt_root,'data','test_jlmt');
flist = listFilesOfType(test_data_dir, {'wav','mp3'});
flist = check_stim_dirs(flist,'srcroot',test_data_dir,'destroot',test_data_dir);

% Specify minimum parameters for torus projections. Defaults will be used
% for each step.
params.glob.process = {{'ani','pp','li','toract','tonreg_fmri'},...
    {'ani','pp','pc','li','toract','tonreg_fmri'}};
params.glob.force_recalc = {};
params.glob.save_calc = params.glob.process;
params.toract = params_toract('som',struct('fname',...
    fullfile(jlmt_root,'data','maps','map_10-Dec-2006_16:18.mat')),...
    'prev_steps',{'ani','pp','li'});
params.toract(2) = params_toract('som',struct('fname',...
    fullfile(jlmt_root,'data','maps','pc_ci2toract_map_12-Jun-2012_15:47.mat')),...
    'prev_steps',{'ani','pp','pc','li'});
outData = jlmt_proc_series(flist,params);
