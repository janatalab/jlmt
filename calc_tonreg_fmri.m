function [regval,regnames] = calc_tonreg_fmri(toract,params)
% generates fmri tonality tracking regressors for a given stimulus
% 
%   [regval,regnames] = calc_tonreg_fmri(indata,params);
%
% This script uses the splines returned during torus projection to resample
% the spherical harmonics onto the appropriate time-scale.
% 
% REQUIRES
%   indata - the output of calc_toract
%   params.scanner - fMRI study design specs for regressor construction
%       .TR - length of time taken to collect each brain volume
%       .dt - micro-timing slices within each TR
%   params.use_sig - string name of integration time constant [created by
%       calc_context('calc_names',<integration time constant vector>)] to
%       use when extracting torus activation timecourses
% 
% RETURNS
%   regval - a matrix of spherical harmonic timecourses that represent
%       change in the given torus projection over time
%   regnames - a cell array of strings identifying each spherical harmonic
%       timecourse in the 'regval' matrix
% 
% 2012.05.22 FB - adapted from fmri_generate_tonreg.m

% init vars
regval = [];
regnames = {};

%see if toract was a 'get_default_params' tag
if(ischar(toract) && strcmp(toract,'getDefaultParams'))
  if ~exist('params','var'), params = ''; end
  regval = getDefaultParams(params);
  return
end

ensemble_globals;
sdefs = params.scanner;
toract_cols = set_var_col_const(toract.vars);

% get the spline data
try use_sig = params.use_sig; catch use_sig = 'tc_2'; end
sidx = strmatch(use_sig,toract.data{toract_cols.labels});
sphm_st = toract.data{toract_cols.spherharm}{sidx};
regnames = sphm_st.spher_names;
spdata = sphm_st.spline;

if iscell(spdata.form)
  spdata.form = spdata.form{1};
  spdata.breaks = spdata.breaks{1};
  spdata.coefs = spdata.coefs{1};
  regnames = regnames{1};
end

% Create the time-scale for the MRI data
fmri_timevect = min(spdata.breaks):sdefs.TR/sdefs.dt:max(spdata.breaks);
regval = ppval(spdata, fmri_timevect)';

% Remove the DC offset term if it is present and small, i.e. previously removed
dc_col = strmatch('cc00',regnames,'exact');
if ~isempty(dc_col) && abs(mean(regval(:,dc_col))) < 1E-09
  regval(:,dc_col) = [];
  regnames(dc_col) = [];
end

return

%makes sure that all necessary param fields are populated and
%non-empty. Sets empty or non-existent values to default values
function params = getDefaultParams(params)

def = params_tonreg_fmri;

fnames = fieldnames(def);

if(~exist('params','var'))
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