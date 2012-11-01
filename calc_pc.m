function pc = calc_pc(indata,params)

% Projects periodicity pitch or context (integrated periodicity pitch) images into pitch class space.
% 
%   pc = calc_pc(indata,params);
%
% This script takes auditory data represented in periodicity pitch or
% integrated periodicity pitch space, and projects it to pitch class space
% using a trained neural network model, and the MATLAB Neural Network
% Toolbox.
% 
% This script will accept both context images and periodicity pitch
% images.
%
% IF you want to project context images to pitch class space, you MUST
% specify 'params.li_siglist', since the output from calc_context.m can
% contain multiple context images (integrated at different timeconstants).
% 
% calc_pitchclass is intended to be called by jlmt_proc_series.m. If you
% specify params.li_siglist, jlmt_proc_series.m will assume that you want
% to project context image data into pitch class space, and will therefore
% send LI data to this function. If you do not specify params.li_siglist,
% it will assume that you want to project periodicity pitch images to pitch
% class space, and will therefore send pp data to this function.
% 
% REQUIRES
%   indata - output from calc_context.m or calc_pp.m
%   params.nnet.fname - path to the neural network model that will be used
%       to project indata into pitch class space. It is assumed that
%       params.nnet.fname holds a variable named "net" that is the output
%       of the MATLAB Neural Network Toolbox, containing a trained neural
%       network model that will project pp or li-space data into pitch
%       class space.
%   params.HalfDecayTimes (used if indata.type = 'li') - Half Decay Times
%       used by IPEMLeakyIntegration
%   params.li_siglist (used if indata.type = 'li') - names of images in
%       indata to use. If this is not specified, jlmt_proc_series will
%       assume that you want to project periodicity pitch images to pitch
%       class space, and it will send PP data rather than LI data to this
%       function.
%   params.norm - equation (using @inline) used to normalize the input
% 
% RETURNS
%   pc - an ensemble data struct with the following variables:
%
% Copyright (c) 2011-2012 The Regents of the University of California
% All Rights Reserved.
%

% 2011.05.06 FB - adapted from calc_toract.m
% 2012.10.26 PJ - eliminated unnecessary assignment of empty params field
% 2012.11.01 TC - inserted "pc.params = [];" at line 100, to avoid error at
%                  line 111, where there is a reference to params.

pc = [];
if ischar(indata) && ~isempty(strmatch(indata,'getDefaultParams'))
    if ~exist('params','var'), params = ''; end
    pc = getDefaultParams(params);
    return
end

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

% are we dealing with 'li' data or 'pp' data?
if strmatch(indata.type,'li','exact')
  pc.time_constants = params.HalfDecayTimes;
  pc.labels = params.li_siglist;
  nsig = length(params.li_siglist);
  for isig = 1:nsig
    signame = params.li_siglist{isig};
    lidx = find(ismember(indata.data{indata_cols.names},signame));
    if isempty(lidx)
      warning('no context image found for %s, SKIPPING',signame);
      continue
    end
  
    limg = indata.data{indata_cols.signals}{lidx};
    if isfield(params,'norm') && ~isempty(params.norm);
      normFunc = inline(params.norm);
      limg = normFunc(limg);
    end
    pc.sig{isig} = sim(net,limg);
  end % for isig=
elseif strmatch(indata.type,'pp','exact')
  pc.sig{1} = sim(net,indata.data{indata_cols.sig});
else
  error('unknown indata type: %s',indata.type);
end

pc.Fs = params.Fs;
pc.params = [];
pc = ensemble_tree2datastruct(pc);
pc_dataCols = set_var_col_const(pc.vars);
pc.type = 'pc';

%store parameters from input data as well as the current calculation
storeParams = indata.data{indata_cols.params};
storeParams.pc = params;
pc.data{pc_dataCols.params} = storeParams;

%makes sure that all necessary param fields are populated and
%non-empty. Sets empty or non-existent values to default values
function params = getDefaultParams(params)

def = params_pc;

fnames = fieldnames(def);

if isempty(params) || ~isstruct(params)
  params = def;
  return
else

  for ifld = 1:length(fnames)
    if (~isfield(params,fnames{ifld}) || isempty(params.(fnames{ifld})))
      params.(fnames{ifld}) = def.(fnames{ifld});
    end
  end

end

return
