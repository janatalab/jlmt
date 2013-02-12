function outdata = calc_mode_estimate(inData,varargin)

% estimates average and timecourse of mode given a torus activation set
% 
%   outdata = calc_mode_estimate(inData,params)
% 
% Five different mode estimate techniques are currently implemented:
%   1. ratio of the avg pixel value in major vs minor mode torus regions
%   2. ratio of the avg half-rectified (only accepting positive values)
%      pixel value in major vs minor torus regions
%   3. region (major or minor) containing the max value in the torus
%      projection
%   4. ratio of the max pixel value in major vs minor mode torus regions
% 
% REQUIRES
%   inData - the output of jlmt calc_toract()
%   params
%       .toract_mode_map.fname - path to file containing the mode map for
%           the given toract model. This must contain the variables
%           'bound_params' and 'grp_mtx', which can be generated using the
%           utility file 'make_toract_mode_map.m'
%       .li_siglist - names of images in inData to use
%       .mode_tc_fun - a string or cell array of strings indicating the
%           function or functions containing the algorithm to estimate the
%           mode timecourse (mode_est_max,mode_est_ratio,mode_est_hr_ratio,
%           mode_est_max_ratio,mode_est_hr_ratio_2)
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
    'avg','pct_minor','pct_major','timeseries','method'};
oc = set_var_col_const(outdata.vars);
for k=1:length(outdata.vars)
  outdata.data{k} = [];
end
outdata.data{oc.labels} = {};
outdata.data{oc.method} = {};
outdata.data{oc.timeseries} = {};

% load toract mode map information
load(params.toract_mode_map.fname);
try labelSize = params.labelSize; catch labelSize=14; end
map_coords = som_vis_coords(sM.topol.lattice,sM.topol.msize);
grp_mtx(ismember(grp_mtx(:),find(line_class==-1))) = 0; % mark minor regions
grp_mtx(ismember(grp_mtx(:),find(line_class==1))) = 1; % mark major regions

% iterate over algorithms and signals, calculate mode timecourse
tc = set_var_col_const(inData.vars);
nsig = length(params.li_siglist);
if ~iscell(params.mode_tc_fun),params.mode_tc_fun={params.mode_tc_fun};end
for ifun=1:length(params.mode_tc_fun)
  fh = parse_fh(params.mode_tc_fun{ifun});
  
  for isig=1:nsig
    % get the signal
    signame = params.li_siglist{isig};
    lidx = find(ismember(inData.data{tc.labels},signame));
    if isempty(lidx)
      warning('no context image found for %s, SKIPPING',signame);
      continue
    end
    
    planes = inData.data{tc.activations}{lidx};

    % get the average value for all cells within a given region
    [mode_tc,avg,pct_major,pct_minor] = fh(planes,grp_mtx);
    
    outdata = ensemble_add_data_struct_row(outdata,...
        'avg',avg,'pct_minor',pct_major,'pct_major',pct_minor,...
        'timeseries',mode_tc,...
        'method',func2str(fh),...
        'time_constants',inData.data{tc.time_constants}(lidx),...
        'labels',inData.data{tc.labels}{lidx},...
        'params',params);
  end % for isig=1:nsig
end % for ifun=1:length(params.mode_tc_fun

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



%%%%% old plotting code

%     figure();
%       clf;
%       subplot(1,2,1)
%       imagesc(new_plane);
% 
%       coord_translation = size(new_plane)./sM.topol.msize;
%       label_coords = (map_coords(label_idxs,:)-0.5).*repmat(coord_translation,length(label_idxs),1);
%       t = text(label_coords(:,1),label_coords(:,2),new_label_array, ...
%           'fontweight','bold','fontsize', labelSize, ...
%           'horizontalalignment', 'center', ...
%           'verticalalign', 'middle');
% 
%       for l=1:size(bound_params,1)
%         avg_slope = mean(bound_params(l,1));
%         avg_int = mean(bound_params(l,2));
%         x = [0 -avg_int/avg_slope];
%         y = [avg_int 0];
%         line(x,y,'color','w');
%       end % for l=1:size(bound_params,1
%       
%       subplot(1,2,1);
%       text(mxidx,myidx,'X','color','w','fontsize',15,'fontweight','bold');
%       subplot(1,2,2)
%       imagesc(grp_mtx);
%       text(mxidx,myidx,'X','color','w','fontsize',15,'fontweight','bold');
%       pause(0.25);
