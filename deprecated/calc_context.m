function ci = calc_context(inData,par)
% Calculates context images, using IPEMLeakyIntegration()
% 
%   ci = calc_context2(pp,ci_params);
%
% Takes an auditory nerve image, and related parameters, and creates
% contextuality images for the requested timescales. The old method used
% IPEMContextualityIndex, with parameters for HalfDecayChords and
% HalfDecayToneCenters. If these parameters are provided, then the return
% structure will include the 5 ci signal matrices returned by
% IPEMContextualityIndex. Otherwise, the return structure will contain one
% ci signal matrix for each requested timescale.
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
%   ci data struct
%       .vars = {'Fs','periods','params',...
%           'names','signals'}
%
% Copyright (c) 2007-2013 The Regents of the University of California
% All Rights Reserved.
%
% 12/04/06 Petr Janata
% 4/14/07 Stefan Tomic - changed data to conform to ensemble data struct
% 2010/02/18 FB - adapted to support calculation of multiple (more than 2)
% contextuality images, using IPEMLeakyIntegration rather than
% IPEMContextualityIndex. Added a bit of header documentation.

if(ischar(inData) && strmatch(inData,'calc_names','exact'))
    ci = HalfDecayNameGen(par);
    return
end

ci = init_ci_struct;
ciCols = set_var_col_const(ci.vars);
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
ci.data{ciCols.Fs} = Fs;

if isfield(dataCols,'periods')
  periods = inData.data{dataCols.periods}; 
  if (iscell(periods)), periods = periods{1}; end
  ci.data{ciCols.periods}    = periods;
end

%store parameters from input data as well as the current calculation
storeParams = inData.data{dataCols.params};
storeParams.ci = par;
ci.data{ciCols.params} = storeParams;

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
  ci.data{ciCols.names} = {'sig1','sig2','ref_sig1_match',...
      'ref_sig2_match','sig1_sig2_match'};
  ci.data{ciCols.signals} = {sig1,sig2,ref_sig1_match,...
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
  ci.data{ciCols.names} = names;

  % Perform the calculation using IPEMLeakyIntegration
  for it=1:ntime
    ci.data{ciCols.signals}{it} = IPEMLeakyIntegration(inSig,Fs,...
        times(it),par.Enlargement,0);
  end
end

if isfield(par,'metrics')
  metfuncs = fieldnames(par.metrics);
  nmetfunc = length(metfuncs);
  for k=1:nmetfunc
    fh = parse_fh(metfuncs{k});
    ci.vars = [ci.vars metfuncs{k}];
    l = length(ci.vars);
    ci.data{l} = fh(ci,par.metrics.(metfuncs{k}));
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