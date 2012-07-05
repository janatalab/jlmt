% BASIC TTt UNIT TEST SCRIPT
% 
% This script tests ipem_proc_series and TTt analyses on a .wav file to
% verify that everything is working. 
%
% Copyright (c) 2007-2012 The Regents of the University of California
% All Rights Reserved.
%
% 13 June 2007 - Stefan Tomic
% 22 May 2012 - Fred Barrett

% Collect test file paths
test_data_dir = fullfile(fileparts(which('ipem_proc_series')),'data','test_ipem');
dirStruct = dir(fullfile(test_data_dir,'*.wav'));
fnames = {dirStruct.name};
flist = check_stim_dirs(fnames,'srcroot',test_data_dir,'destroot',test_data_dir);

% Specify minimum parameters. Defaults will be used for each step.
% params.glob.process = ...
%       {{'ani','pp','ci','toract'},{'ani','pp','pc','ci','toract'}};
params.glob.process = {'ani','pp','ci','toract'};
params.glob.force_recalc = {};
params.glob.save_calc = params.glob.process;
ipem_proc_series(flist,params);
