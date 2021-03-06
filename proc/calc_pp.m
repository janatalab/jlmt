function pp = calc_pp(inData,varargin)
% pp = calc_pp(ani,pp_params);
%
% Calculates a periodicity pitch image from an ANI image using IPEM Toolbox
% routine IPEMPeriodicityPitch()
%
% Copyright (c) 2007-2013 The Regents of the University of California
% All Rights Reserved.
%
% 12/04/06 Petr Janata
% 12/10/06 PJ - Added attenuation function capability
% 4/14/07 Stefan Tomic - data structures now conform to ensemble datastruct
% 30Oct2012 PJ - added support for varargin and handling of default
%                parameters
% 29Jul2019 PJ - Fixed indexing of varargin in getDefaultParams

if nargin > 1 && isstruct(varargin{1})
  par = varargin{1};
end

%see if inData was a 'get_default_params' tag
if ischar(inData) && strcmp(inData,'getDefaultParams')
  if exist('par','var')
    pp = getDefaultParams('params',par);
  else
    pp = getDefaultParams(varargin{:});
  end
  return
else
  aud = inData;
end

pp = init_pp_struct;
ppCols = set_var_col_const(pp.vars);

inDataCols = set_var_col_const(inData.vars);
aniSig = inData.data{inDataCols.sig};
sampleFreq = inData.data{inDataCols.Fs};

% Perform the actual calculation
[outSig, outFs, outPeriods, outFiltANI] = ...
    IPEMPeriodicityPitch( ...
    aniSig, ...
    sampleFreq, ...
    par.LowFrequency, ...    
    par.FrameWidth, ...    
    par.FrameStepSize, ...    
    par.PlotFlag);

% Apply attenuation function if one is specified
if ~isempty(par.Atten)
  attenFunc = pp_atten_func(par.Atten);
  outSig = outSig .* repmat(attenFunc(size(outSig,1))',1,size(outSig,2));
end

pp.data{ppCols.sig} = outSig;
pp.data{ppCols.Fs} = outFs;
pp.data{ppCols.periods} = outPeriods;
pp.data{ppCols.filtANI} = outFiltANI;


%store parameters from input data as well as the current calculation
storeParams = inData.data{inDataCols.params};
storeParams.pp = par;
pp.data{ppCols.params} = storeParams;


return

% makes sure that all necessary param fields are populated and
% non-empty. Sets empty or non-existent values to default values
function params = getDefaultParams(varargin)

for iarg = 1:2:nargin
  switch varargin{iarg}
    case 'params'
      params = varargin{iarg+1};
  end
end

def = params_pp(varargin{:});

fnames = fieldnames(def);

if ~exist('params','var') || ~isstruct(params)
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