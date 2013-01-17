function outdata = calc_mode_estimate(inData,varargin)

% estimates average and timecourse of mode given a torus activation set
% 
%   outdata = calc_mode_estimate(inData,params)
% 
% REQUIRES
%   inData - the output of jlmt calc_toract()
%   params
%       .toract_mode_map.fname - path to file containing the mode map for
%           the given toract model. This must contain the variables
%           'bound_params' and 'grp_mtx', which can be generated using the
%           utility file 'make_toract_mode_map.m'
%       .li_siglist - names of images in inData to use
% 
% RETURNS
%   outdata.vars/data
%       params - the input parameter structure
%       time_constants - time constants associated with each toract source
%       labels - time constant labels associated with each toract source
%       avg - the average mode estimate of the entire given toract
%           timeseries. > 0 indicates that the source signal is, on
%           average, in major mode, < 0 indicates that the source signal
%           is, on average, in minor mode. 
%       pct_minor - % of the timecourse estimated to lie in minor mode
%       pct_major - % of the timecourse estimated to lie in major mode
%       timeseries - the timeseries of the mode estimate for the entire
%           given toract timeseries
% 
% FB 2013.01.16

if nargin > 1 && isstruct(varargin{1})
  params = varargin{1};
end
  
if ischar(inData) && strcmp(inData,'getDefaultParams')
  if exist('params','var')
    outdata = getDefaultParams('params',params);
  else
    outdata = getDefaultParams(varargin{:});
  end
  return
end

% initialize variables
outdata = ensemble_init_data_struct();
outdata.vars = {'params','time_constants','labels',...
    'avg','pct_minor','pct_major','timeseries'};
oc = set_var_col_const(outdata.vars);
for k=1:length(outdata.vars)
  outdata.data{k} = [];
end
outdata.data{oc.labels} = {};

% load toract mode map information
load(params.toract_mode_map.fname);

tc = set_var_col_const(inData.vars);
nsig = length(params.li_siglist);
for isig=1:nsig
    signame = params.li_siglist{isig};
    lidx = find(ismember(inData.data{tc.labels},signame));
    if isempty(lidx)
      warning('no context image found for %s, SKIPPING',signame);
      continue
    end
    
    planes = inData.data{tc.activations}{lidx};
    new_plane = planes(:,:,1);

    % get the average value for all cells within a given region
    areas = cell(size(bound_params,1)+1,1);
    for k=1:size(planes,1)
      for l=1:size(planes,2)
        areas{grp_mtx(k,l)}(end+1,:) = planes(k,l,:);
      end % for l=1:size(new_plane,2
    end % for k=1:size(new_plane,1

    tcs = cell2mat(cellfun(@nanmean,areas,'UniformOutput',false));
    mode_tc = line_class*tcs;

    outdata = ensemble_add_data_struct_row(outdata,'avg',mean(mode_tc),...
        'pct_minor',sum(mode_tc < 0)/length(mode_tc),...
        'pct_major',sum(mode_tc > 0)/length(mode_tc),...
        'timeseries',mode_tc,...
        'time_constants',inData.data{tc.time_constants}(lidx),...
        'labels',inData.data{tc.labels}{lidx});
end % for isig=1:nsig

% append the parameters structure to the output data
outdata.data{oc.params} = params;

%%%%%%%%%%%%%%%
%makes sure that all necessary param fields are populated and
%non-empty. Sets empty or non-existent values to default values
function params = getDefaultParams(varargin)

for iarg = 1:2:nargin
  switch varargin{iarg}
    case 'params'
      params = varargin{iarg}+1;
  end
end

def = params_mode_estimate(varargin{:});

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
