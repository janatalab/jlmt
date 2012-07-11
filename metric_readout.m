function out_st = metric_readout(data_st, params)
% Given a list of stimulus paths specified in data_st, this function will
% extract metrics based on several fields in the params structure.
% 
% Copyright (c) 2012 The Regents of the University of California
% All Rights Reserved.
%

% 01Nov2011 Petr Janata
% 07Nov2011 Fred Barrett - added init_outdata, check each stimpath to
% verify that it points to a directory

try 
	verbose = params.verbose;
catch ME
	verbose = 0;
end

%% Initialize an output data structure
out_st = ensemble_init_data_struct;
outdims = {'stimulus','stimulus_path','wav','wavname','space','tc_local',...
    'tc_global','readout_window','timevect','metric','value','experiment'};
out_st.vars = outdims;
out_st.type = 'metric_readout';
outcols = set_var_col_const(out_st.vars);
out_st.data = cell(1,length(out_st.vars));

if ischar(data_st) && strmatch(data_st,'init_outdata'), return, end

%% Initialize readout parameters
% Get the list of stimulus projection spaces that we want to calculate the
% metrics on
useSpaces = params.useSpaces;
numSpaces = length(useSpaces);

% Get the list of time constant pairs (or single time constants) that we
% want to return metrics for
timeConstantList = params.timeConstantList;
numTimeConstantPairs = length(timeConstantList);

% Get the list of readout windows. This is a structure array
readoutWindowList = params.readoutWindows;
nWind = length(readoutWindowList);

% Get list of metrics
calcMetrics = params.calcMetrics;
numMetrics = length(calcMetrics);

% Get list of stimuli from the data structure
dsCols = set_var_col_const(data_st.vars);

stimPaths = data_st.data{dsCols.stimpath};
nstims = length(stimPaths);

%% Loop over stims, spaces, time constants, metrics
for istim = 1:nstims
	currStimPath = stimPaths{istim};
	[~, stimName] = fileparts(currStimPath);
    if ~exist(currStimPath,'dir')
      fprintf(['Not processing stimulus (%d/%d: %s) - this is not a '...
          'stimulus directory'],istim,nstims,stimName);
      continue
    else
	  fprintf('Processing stimulus (%d/%d): %s\n', istim, nstims, stimName);
    end
	
	% Locate stimulus file
	if params.load_wav
		out_st.data{outcols.wav}{end+1,1} = wavread(fullfile(currStimPath,'audio',sprintf('%s.wav', stimName)));
		out_st.data{outcols.wavname}{end+1,1} = stimName;
	end
	
	% Loop over spaces (ci, toract)
	for ispace = 1:numSpaces
		currSpace = useSpaces{ispace};
		spaceDir = fullfile(currStimPath, currSpace);
		if ~exist(spaceDir)
			error('Desired analysis space (%s) does not exist for stimulus (%s)', currSpace, stimName);
		end
		
		% Loop over time-constant pairs.
		
		% Maintain a list of data we have already loaded (cached) and use it.
		haveTCData = [];
		
		for iTCPair = 1:numTimeConstantPairs
			% Look through mat-files for the current space to find one containing
			% all the necessary time-constants
			currTCPair = timeConstantList{iTCPair};
			
			% Find existing analysis
			clear matchParams;
			switch currSpace
				case 'ci'
					varName = 'cntxt_img';
					tcName = 'HalfDecayTimes';
					
				case 'pc_ci'
					varName = 'pc_cntxt_img';
					tcName = 'HalfDecayTimes';

                case 'toract'
					varName = 'toract';
					tcName = 'HalfDecayTimes';
					
				otherwise
					varName = currSpace;
			end
			matchParams.varName = varName;
			matchParams.paramFind.ani = params.jlmt.ani;
			matchParams.paramFind.pp = params.jlmt.pp;
            if ~strcmp(currSpace, 'pc_ci')
                matchParams.paramFind.(currSpace) = params.jlmt.(currSpace);
                matchParams.paramFind.(currSpace).(tcName) = currTCPair;
            end
			matchParams.ignore = {'ani.aniPath','ci','Fs','pc'};
			matchParams.subsetIsOK = {tcName};
			
			[previousCalcFname, violation, violationReason] = check_anal_exist(spaceDir,matchParams);
			
			if isempty(previousCalcFname)
				error(['Could not locate MAT-file with matching parameters in %s\n' ...
				'Violation reason: %s'], spaceDir, violationReason)
					% This is where we can add calculation of the analysis via
					% jlmt_proc_series
			else
				calcInfo = load(fullfile(spaceDir,previousCalcFname{1}));
			end
		
			% Find the datastruct
            if strcmp(currSpace, 'pc_ci')
                currSpaceLabel = 'ci';
            else
                currSpaceLabel = currSpace;
            end
			idx = ensemble_find_analysis_struct(calcInfo.(varName), struct('type',currSpaceLabel));
			if isempty(idx)
				error('Could not find analysis of type: %s', currSpace)
			end
			
			% Extract the data
			an_st = calcInfo.(varName);
			anCols = set_var_col_const(an_st.vars);
			
			% Get appropriate time-constant indices
			halfDecayTimes = an_st.data{anCols.params}.(currSpaceLabel).HalfDecayTimes;
			[tcMask, tcIdxs] = ismember(currTCPair, halfDecayTimes);
			if ~all(tcMask)
				error('Not all desired time-constants were found')
			end
			
			% Refer to the shorter time-constant as the target data. Handle case
			% in which we are doing a comparison using a single time-constant
			[~,tmpIdx] = min(halfDecayTimes(tcIdxs));
			targIdx = tcIdxs(tmpIdx);
			
			if ~any(diff(tcIdxs))
				refIdx = targIdx;
			else
				refIdx = setdiff(tcIdxs, targIdx);
			end
			
			switch currSpace
				case 'ci'
					dataField = 'signals';
				case 'toract'
					dataField = 'activations';
				case 'pc_ci'
					dataField = 'signals';
				otherwise
					dataField = '';
			end
			
			% Make copies of the data sources
			targData = an_st.data{anCols.(dataField)}{targIdx};
			refData = an_st.data{anCols.(dataField)}{refIdx};
			
			switch currSpaceLabel
				case 'ci'
					Fs = an_st.data{anCols.Fs};
				case 'toract'
					Fs = an_st.data{anCols.params}.(currSpace).Fs;
					if ndims(targData) == 3
						dims = size(targData);
						targData = reshape(targData, prod(dims(1:2)), dims(3));
						refData = reshape(refData, prod(dims(1:2)), dims(3));
					end
			end
			
			timevect = (0:(size(targData,2)-1))/Fs*1000; % get time vector in milliseconds
			
			% Loop over desired readout windows, and handle each type
			for iwind = 1:nWind
				currWind = readoutWindowList(iwind);
				readoutType = currWind.type;				
						
				% Loop over metrics
				for imetric = 1:numMetrics
					currMetric = calcMetrics{imetric};
			
					fh = str2func(currMetric);
					if ~isfield(params,currMetric)
						params.(currMetric) = struct();
					end
					
					% Make sure this metric is compatible with this readoutType
					if ismember(currMetric, {'jlPre_target_baseline'}) && ismember(readoutType,{'post','prepost'})
						if verbose
							fprintf('Metric %s is incompatible with readout window type %s\n', currMetric, readoutType);
						end
						continue
					end
					
					% Copy window information to the params struct for this metric
					params.(currMetric).wininfo = currWind;
					
					% Get the appropriate chunks of data
					switch readoutType
						case 'pre'
							startSampRef = find(timevect >= currWind.event_onset_ms+currWind.pre_ms(1), 1, 'first' );
							stopSampRef = find(timevect <= currWind.event_onset_ms+currWind.pre_ms(2), 1, 'last');
							
							startSampTarg = startSampRef;
							stopSampTarg = stopSampRef;
							
						case 'post'
							startSampRef = find(timevect >= currWind.event_onset_ms+currWind.post_ms(1), 1, 'first' );
							stopSampRef = find(timevect <= currWind.event_onset_ms+currWind.post_ms(2), 1, 'last');
							
							startSampTarg = startSampRef;
							stopSampTarg = stopSampRef;
							
							% Make a copy of the current time vector available to the
							% function that's called
							currTimevectRelOnset = timevect(startSampTarg:stopSampTarg)-currWind.event_onset_ms;
							params.(currMetric).timevect = currTimevectRelOnset;
							
							
						case 'prepost'
							
							% Specify start and stop samples depending on whether the metric relies on target and
							% reference data to be registered in time
							switch currMetric
								% both chunks span from earliest pre-stimulus to latest
								% post-stimulus
								case {'jlPost_target_max','jlPost_target_sum','jlPost_target_mean','map_max'}
									startSampRef = find(timevect >= currWind.event_onset_ms+currWind.pre_ms(1), 1, 'first' );
									stopSampRef = find(timevect <= currWind.event_onset_ms+currWind.post_ms(2), 1, 'last');
									
									startSampTarg = startSampRef;
									stopSampTarg = stopSampRef;
									
									% Make a copy of the current time vector
									currTimevectRelOnset = timevect(startSampTarg:stopSampTarg)-currWind.event_onset_ms;
									params.(currMetric).timevect = currTimevectRelOnset;
									
								otherwise
									% reference chunk is pre-stimulus; target chunk is
									% post-stimulus
									startSampRef = find(timevect >= currWind.event_onset_ms+currWind.pre_ms(1), 1, 'first' );
									stopSampRef = find(timevect <= currWind.event_onset_ms+currWind.pre_ms(2), 1, 'last');
									
									startSampTarg = find(timevect >= currWind.event_onset_ms+currWind.post_ms(1), 1, 'first' );
									stopSampTarg = find(timevect <= currWind.event_onset_ms+currWind.post_ms(2), 1, 'last');
							end
							
						otherwise
							error('Unknown readout type: %s', readoutType)
					end % switch readoutType
					
					% Call the function
					val = fh(refData(:,startSampRef:stopSampRef), targData(:,startSampTarg:stopSampTarg), params.(currMetric));
					
					% Write the result to the output structure
					% {'stimulus','space','tc_pair','readoutWindows','metric','value'};
					out_st.data{outcols.stimulus}{end+1,1} = stimName;
                    out_st.data{outcols.stimulus_path}{end+1,1} = currStimPath;
					out_st.data{outcols.space}{end+1,1} = currSpace;
					out_st.data{outcols.tc_local}(end+1,1) = currTCPair(1);
					out_st.data{outcols.tc_global}(end+1,1) = currTCPair(2);
					out_st.data{outcols.readout_window}{end+1,1} = currWind.name;
					out_st.data{outcols.metric}{end+1,1} = currMetric;
					out_st.data{outcols.timevect}{end+1,1} = timevect(startSampTarg:stopSampTarg);
					out_st.data{outcols.value}{end+1,1} = val;
					
					% Add optional writing/appending of this result to the
					% metric_readout var in the analysis space matfile
					
				end % for imetric
			end % for iwind = 
		end % for iTCPair=
	end % ispace
end % for istim=

if isfield(params,'experiment_name')
  out_st.data{outcols.experiment} = cell(length(out_st.data{outcols.stimulus}),1);
  [out_st.data{outcols.experiment}{:}] = deal(params.experiment_name);
end

return