% compile metric readout data for tonmodcomp
% 
% This script calculates metrics for different spaces, windows, and
% calculation types.
% 
% There are parameters below that you may want to experiment with.
% 

% 08Nov2011 PJ - fixed search for incomplete analyses
% 16Nov2011 TC - As I added new readout windows (for 0-600 ms and 201-600
% ms), I changed the indices of params.readoutWindows (line 32) so that
% these are still correct.
% 26Feb2012 TC - Added pc to params.useSpaces, to pursue Petr's suggestion
% of a pool of variables based on chromagrams. I use the term '3 pools' to
% label the output file, as we are now conisdering three pools of
% variables.

%% initialize analysis variables
if ~exist('params','var')
  params = CollinsEtAl_globals;
end

%%
%%  PARAMETERS THAT CAN BE CHANGED
%% 
params.useSpaces = {'li','toract'};
params.timeConstantList = {[0.1 4],[4 4]};

% specify general metrics to calculate
params.calcMetrics = {'jlPost_target_mean','map_max'};

% Specify metric-specific parameters
params.jlCorr.type = 'average';
params.jlCorr.preproc = 'average_time';
params.jlKLdist = params.jlCorr;

analysis_fname = fullfile(params.paths.analysis_path,'tonmodcomp_metrics.mat');
check_dir(fileparts(analysis_fname),0,1);

%% define metric readout windows
rinfo = struct();
nr = 0;

% 0-200 ms
nr = nr+1;
rinfo(nr).name = 'post_0_200ms';
rinfo(nr).type = 'post';
rinfo(nr).pre_ms = [];
rinfo(nr).post_ms = [0 200];
rinfo(nr).event_onset_ms = [];

% 201-600 ms
nr = nr+1;
rinfo(nr).name = 'post_201_600ms';
rinfo(nr).type = 'post';
rinfo(nr).pre_ms = [];
rinfo(nr).post_ms = [201 600];
rinfo(nr).event_onset_ms = [];

%% Pre/post comparisons
nr = nr+1;
rinfo(nr).name = 'prepost_100base_0_200ms';
rinfo(nr).type = 'prepost';
rinfo(nr).pre_ms = [-100 0];
rinfo(nr).post_ms = [0 200];
rinfo(nr).event_onset_ms = [];
rinfo(nr).preIsLongerTimeConstant = 0;

nr = nr+1;
rinfo(nr).name = 'prepost_100base_201_600ms';
rinfo(nr).type = 'prepost';
rinfo(nr).pre_ms = [-100 0];
rinfo(nr).post_ms = [201 600];
rinfo(nr).event_onset_ms = [];
rinfo(nr).preIsLongerTimeConstant = 0;

params.readoutWindows = rinfo;

%%
%%  BEYOND HERE: CHANGE AT YOUR OWN RISK
%%

% check to see if an analysis file already exists
if exist(analysis_fname,'file')
	fprintf('Loading analyses from file: %s\n', analysis_fname);
  load(analysis_fname);
else
  an_data = metric_readout('init_outdata');
end
ac = set_var_col_const(an_data.vars);

prevEventWindowDuration = 400;

params.load_wav = false;

% create parameter check list - pcl
pcl = {};
nmetrics = length(params.calcMetrics);
nreadout = length(params.readoutWindows);
nspaces  = length(params.useSpaces);
for m=1:nmetrics
  for r=1:nreadout
    for s=1:nspaces
      pcl(end+1,:) = {params.useSpaces{s},...
          params.readoutWindows(r).name,params.calcMetrics{m}};
    end % for s=1:nspaces
  end % for r=1:nreadout
end % for m=1:nmetrics

% initialize filepath structure
fpath_st = ensemble_init_data_struct;
fpath_st.vars = {'stimpath'};

% iterate over datasets, load paths
ndataset = length(params.datasets);
ntconst  = length(params.timeConstantList);
for k=1:ndataset
  fprintf(1,'analyzing %s\n',params.datasets(k).id);
  
  % Specify which readout windows we want to use
  emptyOnsetMask = cellfun('isempty', {params.readoutWindows.event_onset_ms});

  % Assign target event onsets
  targOnset = params.datasets(k).event_onsets(end);
  [params.readoutWindows(emptyOnsetMask).event_onset_ms] = deal(targOnset);

  % For previous event comparisons, assign pre_ms windows
  emptyPreMask = cellfun('isempty', {params.readoutWindows.pre_ms});
  prepostMask = ismember({params.readoutWindows.type},'prepost');
  pre_ms = params.datasets(k).event_onsets(end-1)-targOnset;
  pre_ms = [pre_ms pre_ms+prevEventWindowDuration];
  [params.readoutWindows(emptyPreMask&prepostMask).pre_ms] = deal(pre_ms);

  % update params for this dataset
  params.experiment_name = params.datasets(k).id;
  
  % create file path input structure for this dataset
  if ~isfield(params.datasets(k), 'subsets') || isempty(params.datasets(k).subsets)
	  fpath_st.data = {build_path_list(fullfile(params.paths.data_root,params.datasets(k).id))};
  else
  	  fpath_st.data = {build_path_list(fullfile(params.paths.data_root,params.datasets(k).id), ...
 	    'targ_dirs', {params.datasets(k).subsets.id})};
  end

  % check an_data for the proposed analyses
  dsmask = ismember(an_data.data{ac.experiment},params.datasets(k).id);
  if ~isempty(an_data.data{ac.tc_local})
    an_tcpairs = [an_data.data{[ac.tc_local ac.tc_global]}];
  else
    an_tcpairs = [];
  end

  % iterate over requested time constant pairs
  for t=1:ntconst
	tcpairMask = ismember(an_tcpairs,params.timeConstantList{t},'rows');
    
	if ~isempty(tcpairMask)
		compositeMask = tcpairMask & dsmask;
	else
		compositeMask = false(size(dsmask));
    end
    
    if ~any(compositeMask)
      fprintf(1,['\nno analysis results found for %s, (%0.2f/%0.2f time '...
          'constant pair), running analysis\n'],params.datasets(k).id,...
          params.timeConstantList{t}(1),params.timeConstantList{t}(2));
      
      lparams = params;
      lparams.timeConstantList = params.timeConstantList(t);
      result_st = metric_readout(fpath_st,lparams);
      an_data = ensemble_merge_data({an_data,result_st},struct());
    else
      % check spaces, readout_windows and metrics
      mdata = [an_data.data{ac.space}(compositeMask) ...
          an_data.data{ac.readout_window}(compositeMask) ...
          an_data.data{ac.metric}(compositeMask)];
      notfound = find(~all(ismember(pcl,mdata),2));
      if ~isempty(notfound)
        % set missing params in lparams
        lparams = params;
        lparams.timeConstantList = params.timeConstantList(t);
				for n=1:length(notfound)
					
					% THIS PRINT STATEMENT NEEDS FIXING
					fprintf(1,['\nno analysis results found for %s, '...
						'(%0.2f/%0.2f time constant pair, %s space, '...
						'%s readout window, %s metric), running analysis\n'],...
						params.datasets(k).id,params.timeConstantList{t}(1),...
						params.timeConstantList{t}(2),pcl{notfound(n),1},pcl{notfound(n),2},pcl{notfound(n),3});
				end % for n=notfound
				
				% set spaces, readout window, metrics
				roidx = find(ismember({params.readoutWindows.name},pcl(notfound,2)));
				lparams.readoutWindows = params.readoutWindows(roidx);
				lparams.useSpaces = unique(pcl(notfound,1));
				lparams.calcMetrics = unique(pcl(notfound,3));
				
				% analyze, merge
				result_st = metric_readout(fpath_st,lparams);
				an_data = ensemble_merge_data({an_data,result_st},struct());
      end % if ~all(notfound
    end % if ~any(mask
  end % for t=1:ntconst
	
  % Save the data now so that we don't have to recompute if we crash on a
  % subsequent experiment
  save_vars = {'an_data','params'};
  fprintf(1,'Saving data: %s\n', analysis_fname);
  save(analysis_fname,save_vars{:});
	
end % for k=1:ndataset

fprintf(1,'CollinsEtAl_metric_compile: done processing. data saved to %s\n',...
    analysis_fname);
