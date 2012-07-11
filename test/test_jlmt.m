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
test_data_dir = fullfile(fileparts(fileparts(which('jlmt_proc_series'))),...
    'data','test_jlmt');
input_file_pattern = {'*.wav','*.mp3'};

fnames = {};
for k=1:length(input_file_pattern)
  dirStruct = dir(fullfile(test_data_dir,input_file_pattern{k}));
  fnames = [fnames dirStruct.name];
end % for k=1:length(input_file_pattern
flist = check_stim_dirs(fnames,'srcroot',test_data_dir,'destroot',test_data_dir);

% Specify minimum parameters for torus projections. Defaults will be used
% for each step.
params.glob.process = {{'ani','pp','li','toract','tonreg_fmri'},...
    {'ani','pp','pc','li'}};
params.glob.force_recalc = {};
params.glob.save_calc = params.glob.process;
outData = jlmt_proc_series(flist,params);
