function li = calc_li(inData,par)
% JLMT wrapper for IPEMLeakyIntegration()
% 
%   li = calc_li(pp,li_params);
%
% Takes an input matrix and related parameters, and creates leaky
% integrated images for the requested timescales. The old method used
% IPEMContextualityIndex, with parameters for HalfDecayChords and
% HalfDecayToneCenters. If these parameters are provided, then the return
% structure will include the 5 li signal matrices returned by
% IPEMContextualityIndex. Otherwise, the return structure will contain one
% li signal matrix for each requested timescale.
% 
% This was first conceived of in terms of contextuality images, which are
% 0.2s or 2s leaky integrated auditory nerve images or periodicity pitch
% images. This has been expanded to allow for multiple other times and
% spaces, such as 10s (or any other length) integration windows and inputs
% such as pitch chroma vectors.
% 
% NOTE: if you pass 'calc_names' in as the first argument, and a vector of
% HalfDecayTimes in as the second argument, this function will return the
% names that it would generate for these times, given no specific
% HalfDecayNames.
% 
% REQUIRES
%   indata
%       .vars = {'sig','Fs','periods'}
%   par
%       .HalfDecayTimes
%       .HalfDecayNames (optional)
%       (.HalfDecayChords - deprecated)
%       (.HalfDecayToneCenters - deprecated)
% 
% RETURNS
%   li data struct
%       .vars = {'Fs','periods','params',...
%           'names','signals'}
%
% Copyright (c) 2007-2012 The Regents of the University of California
% All Rights Reserved.
%
% 12/04/06 Petr Janata
% 4/14/07 Stefan Tomic - changed data to conform to ensemble data struct
% 2010/02/18 FB - adapted to support calculation of multiple (more than 2)
% contextuality images, using IPEMLeakyIntegration rather than
% IPEMContextualityIndex. Added a bit of header documentation.
% 2012/07/03 FB - changed to 'calc_li' to reflect expansion beyound
% contextuality images calculated from ani or pp space

li = [];

if ischar(inData)
  if ~isempty(strmatch(inData,'calc_names','exact'))
    li = HalfDecayNameGen(par);
    return
  elseif ~isempty(strmatch(inData,'getDefaultParams'))
    if ~exist('par','var'), par = ''; end
    li = getDefaultParams(par);
    return
  end % if ~isempty(strmatch(inData,'calc_names...
end

%see if inData was a 'get_default_params' tag
if(ischar(inData) && strcmp(inData,'getDefaultParams'))
else
  aud = inData;
end

li = init_li_struct;
liCols = set_var_col_const(li.vars);
dataCols = set_var_col_const(inData.vars);
inSig = inData.data{dataCols.sig};

if(iscell(inSig)), inSig = inSig{1}; end

switch inData.type
    case 'pc'
      Fs = inData.data{dataCols.params}.pc.Fs;
    otherwise
      Fs = inData.data{dataCols.Fs};
end

% The following two data entries were previously saved as params.
% They are saved into the data struct now to facilitate param
% matching with saved calculations
li.data{liCols.Fs} = Fs;

if isfield(dataCols,'periods')
  periods = inData.data{dataCols.periods}; 
  if (iscell(periods)), periods = periods{1}; end
  li.data{liCols.periods}    = periods;
end

%store parameters from input data as well as the current calculation
storeParams = inData.data{dataCols.params};
storeParams.li = par;
li.data{liCols.params} = storeParams;

if isfield(par,'HalfDecayChords') && ...
        isfield(par,'HalfDecayToneCenters') && ...
        exist('periods','var')

  % Perform the calculation, using IPEMContextualityIndex
  [sig1, sig2, ref_sig1_match, ref_sig2_match, sig1_sig2_match] = ...
      IPEMContextualityIndex( ...
      inSig, ...
      Fs, ...
      periods, ...
      par.SnapShot, ...
      par.HalfDecayChords, ...
      par.HalfDecayToneCenters, ...
      par.Enlargement, ...
      par.PlotFlag);
  
  % save data to output struct
  li.data{liCols.names} = {'sig1','sig2','ref_sig1_match',...
      'ref_sig2_match','sig1_sig2_match'};
  li.data{liCols.signals} = {sig1,sig2,ref_sig1_match,...
      ref_sig2_match,sig1_sig2_match};
  
elseif isfield(par,'HalfDecayTimes')

  % get timescales
  times = par.HalfDecayTimes;
  ntime = length(times);

  % get or create output signal names
  if isfield(par,'HalfDecayNames')
    names = par.HalfDecayNames;
  else
    names = HalfDecayNameGen(times);
  end

  % save signal names
  li.data{liCols.names} = names;

  % Perform the calculation using IPEMLeakyIntegration
  for it=1:ntime
    li.data{liCols.signals}{it} = IPEMLeakyIntegration(inSig,Fs,...
        times(it),par.Enlargement,0);
  end
end

if isfield(par,'metrics')
  metfuncs = fieldnames(par.metrics);
  nmetfunc = length(metfuncs);
  for k=1:nmetfunc
    fh = parse_fh(metfuncs{k});
    li.vars = [li.vars metfuncs{k}];
    l = length(li.vars);
    li.data{l} = fh(li,par.metrics.(metfuncs{k}));
  end % for k=1:nmetfunc
end % if isfield(params,'metrics

function names = HalfDecayNameGen(times)
  ntime = length(times);
  names = cell(1,ntime);
  for j=1:ntime
    t = times(j);
    if mod(t,fix(t))
      names{j} = regexprep(sprintf('tc_%1.4f',t),'\.','p');
    else
      names{j} = sprintf('tc_%d',t);
    end
  end
  
  return
  
%makes sure that all necessary param fields are populated and
%non-empty. Sets empty or non-existent values to default values
function params = getDefaultParams(params)

def = params_li;

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
