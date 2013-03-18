function pc = calc_pitchclass(indata,params)
% Projects periodicity pitch or context (integrated periodicity pitch) images into pitch class space.
% 
%   pc = calc_pitchclass(indata,params);
%
% This script takes auditory data represented in periodicity pitch or
% integrated periodicity pitch space, and projects it to pitch class space
% using a trained neural network model, and the MATLAB Neural Network
% Toolbox.
% 
% This script will accept both context images and periodicity pitch
% images. It also calculates pitch class metrics, if desired.
%
% IF you want to project context images to pitch class space, you MUST
% specify 'params.ci_siglist', since the output from calc_context.m can
% contain multiple context images (integrated at different timeconstants).
% 
% calc_pitchclass is intended to be called by ipem_proc_series.m. If you
% specify params.ci_siglist, ipem_proc_series.m will assume that you want
% to project context image data into pitch class space, and will therefore
% send CI data to this function. If you do not specify params.ci_siglist,
% it will assume that you want to project periodicity pitch images to pitch
% class space, and will therefore send pp data to this function.
% 
% REQUIRES
%   indata - output from calc_context.m or calc_pp.m
%   params.nnet.fname - path to the neural network model that will be used
%       to project indata into pitch class space. It is assumed that
%       params.nnet.fname holds a variable named "net" that is the output
%       of the MATLAB Neural Network Toolbox, containing a trained neural
%       network model that will project pp or ci-space data into pitch
%       class space.
%   params.HalfDecayTimes (used if indata.type = 'ci') - Half Decay Times
%       used by IPEMLeakyIntegration
%   params.ci_siglist (used if indata.type = 'ci') - names of images in
%       indata to use. If this is not specified, ipem_proc_series will
%       assume that you want to project periodicity pitch images to pitch
%       class space, and it will send PP data rather than CI data to this
%       function.
%   params.norm - equation (using @inline) used to normalize the input
%   params.metrics - this is a struct whose fieldnames will be taken as
%       functions to be executed in order to calculate metrics on
%       pitchclass data. The first argument to the metric function will be
%       'pc' (the output of this function), and the second arguent to the
%       metric function will be params.metrics.(fieldname). Output from
%       each metric function will be stored in a column of the 'pc' data
%       structure with the column name set to the metric function name
% 
% RETURNS
%   pc - an ensemble data struct with the following variables:
%
% Copyright (c) 2011-2013 The Regents of the University of California
% All Rights Reserved.
%
% 2011.05.06 FB - adapted from calc_toract.m

% storage for parameter struct
pc.params = {};

error(nargchk(2,2,nargin))

indata_cols = set_var_col_const(indata.vars);

% copy the sampling rate
params.Fs = indata.data{indata_cols.Fs};

% Load the nnet model
if ~exist(params.nnet.fname,'file')
  fprintf('NNet Model file (%s) does not exist\n',params.nnet.fname);
  return
end
load(params.nnet.fname,'net');

% are we dealing with 'ci' data or 'pp' data?
if strmatch(indata.type,'ci','exact')
  pc.time_constants = params.HalfDecayTimes;
  pc.labels = params.ci_siglist;
  nsig = length(params.ci_siglist);
  for isig = 1:nsig
    signame = params.ci_siglist{isig};
    cidx = find(ismember(indata.data{indata_cols.names},signame));
    if isempty(cidx)
      warning('no context image found for %s, SKIPPING',signame);
      continue
    end
  
    cimg = indata.data{indata_cols.signals}{cidx};
    if isfield(params,'norm') && ~isempty(params.norm);
      normFunc = inline(params.norm);
      cimg = normFunc(cimg);
    end
    pc.sig{isig} = sim(net,cimg);
  end % for isig=
elseif strmatch(indata.type,'pp','exact')
  pc.sig{1} = sim(net,indata.data{indata_cols.sig});
else
  error('unknown indata type: %s',indata.type);
end

pc.Fs = params.Fs;
pc = ensemble_tree2datastruct(pc);
pc_dataCols = set_var_col_const(pc.vars);
pc.type = 'pc';

%store parameters from input data as well as the current calculation
storeParams = indata.data{indata_cols.params};
storeParams.pc = params;
pc.data{pc_dataCols.params} = storeParams;

%% calculate metrics?
if isfield(params,'metrics')
  metfuncs = fieldnames(params.metrics);
  nmetfunc = length(metfuncs);
  for k=1:nmetfunc
    fh = parse_fh(metfuncs{k});
    pc.vars = [pc.vars metfuncs{k}];
    l = length(pc.vars);
    pc.data{l} = fh(pc,params.metrics.(metfuncs{k}));
  end % for k=1:nmetfunc
end % if isfield(params,'metrics

