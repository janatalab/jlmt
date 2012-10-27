function toract = calc_toract(inData,params)
% Projects leaky integrated images into toroidal space using a projection map specified in params.
% 
%   toract = calc_toract(inData,params);
%
% Also calculates spherical harmonics, if desired.
%
% REQUIRES
%   inData - output from calc_li
%   params.som.fname - path to the self organizing map that transforms
%       context image data into torus activation data
%   params.HalfDecayTimes - Half Decay Times used by IPEMLeakyIntegration
%   params.li_siglist - names of images in inData to use
%   params.norm - equation (using @inline) used to normalize inData
%   params.calc_spher_harm - if specified, these parameters will be used to
%       guide construction of spherical harmonics
%       .nharm_theta
%       .nharm_phi
%       .min_rsqr
% 
% RETURNS
%   toract - an ensemble data struct with the following variables:
%       
%
% Copyright (c) 2007-2012 The Regents of the University of California
% All Rights Reserved.
%
% 12/10/06 Petr Janata
% 01/08/07 PJ - adapted from project2som in the tontraj project directory.  The
%               main change is that the input is now an inData strucure
%               rather than a list of files.  This places the function at the
%               same level as calc_ani, etc.
% 4/14/07 Stefan Tomic - changed data to conform to ensemble data struct
% 2010.02.19 FB - changed params.ci.siglist to params.ci_siglist
% 2012.07.03 FB - changed cntxt_img to inData, ci_siglist to li_siglist,
% added getDefaultParams

toract = [];
if ischar(inData) && ~isempty(strmatch(inData,'getDefaultParams'))
    if ~exist('params','var'), params = ''; end
    toract = getDefaultParams(params);
    return
end

%storage for parameter struct
toract.params = {};

error(nargchk(2,2,nargin))

inData_cols = set_var_col_const(inData.vars);

% copy the sampling rate
params.Fs = inData.data{inData_cols.Fs};

% Load the weight matrix
if ~isfield(params.som,'fname')
  som_fname = '';
else
  som_fname = params.som.fname;
end
if ~exist(som_fname)
  fprintf('SOM file (%s) does not exist\n', som_fname);
  return
end

load(som_fname,'sM')
weights = sM.codebook;  % the weight matrix
toract.time_constants = params.HalfDecayTimes;
toract.labels = params.li_siglist;

nsig = length(params.li_siglist);
for isig = 1:nsig
  signame = params.li_siglist{isig};
  lidx = find(ismember(inData.data{inData_cols.names},signame));
  if isempty(lidx)
    warning('no context image found for %s, SKIPPING',signame);
    continue
  end
  
  limg = inData.data{inData_cols.signals}{lidx};
  npts = size(limg,2);
  
  if ~isempty(params.norm);
    normFunc = inline(params.norm);
    limg = normFunc(limg);
  end
  act = weights*limg;  % calculate activation on SOM map
  toract.activations{isig} = reshape(act,[sM.topol.msize npts]);
  
  % Calculate spherical harmonics if desired
  if params.calc_spher_harm
    nharm_theta = params.spher_harm.nharm_theta;
    nharm_phi = params.spher_harm.nharm_phi;
    
    tmp = [];
    [tmp.data, tmp.theta, tmp.phi, tmp.spher_names, tmp.Rsqr, X] = ...
	toroidal_spect(act', sM.topol.msize, nharm_theta, nharm_phi);

    if isfield('min_rsqr',params.spher_harm) && ...
	  ~isempty(params.spher_harm.min_rsqr)
      if any(tmp.Rsqr < params.spher_harm.min_rsqr)
	warning(['calc_toract: Fit of spherical harmonics less than criterion' ...
	      ' value (%1.3f)'], params.spher_harm.min_rsqr)
      end
    end
    
    % Remove columns that have all zeros, i.e. those that incorporate a sin(0)
    % term
    good_cols = find(sum(X.^2));
    X = X(:,good_cols);
    tmp.spher_names = tmp.spher_names(good_cols);
    tmp.data = tmp.data(:,good_cols);
    
    % Fit the data with splines so that users of this data can interpolate to arbitrary
    % timepoints.
    
    % Generate a time vector
    timevect = (1:size(tmp.data,1))/params.Fs;
    
    % Do the spline interpolation. Note that spherical harmonics columns become
    % rows.
    fprintf('Performing spline interpolation of spherical harmonic timeseries ...\n');
    tmp.spline = spline(timevect,tmp.data');
    
    toract.spherharm{isig} = tmp;
  end % if params.calc_spher_harm
end % for isig=

toract = ensemble_tree2datastruct(toract);
toract_dataCols = set_var_col_const(toract.vars);

toract.type = 'toract';
%store parameters from input data as well as the current calculation
storeParams = inData.data{inData_cols.params};
storeParams.toract = params;
toract.data{toract_dataCols.params} = storeParams;

%makes sure that all necessary param fields are populated and
%non-empty. Sets empty or non-existent values to default values
function params = getDefaultParams(params)

def = params_toract;

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
