function outData = jlmt_proc_series(inData,params)
% Performs a series of analyses for the Janata Lab Music Toolbox 
% 
%   outData = jlmt_proc_series(inData,params)
%
% INPUTS
%  inData: this could be one of many different input types:
% 
%           1) a string specifying a single path and filename (for
%                   processing just one input)
%           2) a cell array of filenames, with paths included
%           3) an Ensemble data struct with filename and path variables
%           4) a cell array of 'aud' structs (see init_aud.m and new_aud.m)
% 
%           (the following are only applicable if you have installed the
%           "Ensemble" database)
% 
%           5) an Ensemble data struct with a "stimulus_id" variable, which
%                   specifies a list of stimulus IDs in Ensemble
%           6) a vector of stimulus IDs in Ensemble
%
%       if 'getDefaultParams' tag os passed to jlmt_proc_series as inData,
%       the default parameter structure for jlmt_proc_series (defined in a
%       sub-function at the end of this script) is returned, along with any
%       specific parameters that are set in the 'params' argument (the
%       second argument to jlmt_proc_series). For example, if I execute:
% 
%           jlmt_proc_series('getDefaultParams')
% 
%       it will return the default parameters for the script. If I execute:
% 
%           jlmt_proc_series('getDefaultParams',initial_parameters)
% 
%       it will return the default parameters for the script, with specific
%       values set for parameters that are provided in 'initial_parameters'
%
%  params.glob.lfid: Specifies a log file ID (obtained from fopen). All
%       standard output will be sent to this log file. If this parameter is
%       not set, the output will be sent to the Matlab console.
%
%  params.glob.process: a cell array of cell arrays of strings indicating
%       series of analysis steps to run on the input stimuli. string values
%       will be evaluated to function handles in the form 'calc_%s'. The
%       input signal from a given stimulus will be used as the input for
%       the first analysis step, and subsequent analysis steps will be fed
%       the output from the previous analysis step.
% 
%       possible values: 'ani','pp','li','toract','pc','rp'
% 
%       example:
%           params.glob.process = {{'ani','pp','li','toract'},...
%               {'ani','pp','pc','toract'}};
%           With the above-given process array, two job series will be
%           executed for each stimulus:
%               1) ani->pp->li->toract, and
%               2) ani->pp->pc->toract
%           In each case, the first step ('ani') will receive as input the
%           signal from the stimuli provided to jlmt_proc_series. In the
%           first analysis series, input to 'pp' will be output from 'ani',
%           input to 'li' will be the output from 'pp', and the input from
%           'toract' will be the output from 'li'. In the second analysis
%           series, input to 'pp' will be the output from 'ani', input to
%           'pc' will be the output from 'pp', and input to 'toract' will
%           be the output from 'pc'.
%
%  params.glob.force_recalc: same possible values as params.glob.process.
%       cell array of strings that specify which jobs will force a
%       recalculation, even if a saved calculation with matching parameters
%       is found.
%
%  params.glob.ignore: a cell array of param fields to ignore when finding
%       a match for an identical saved calculation. This allows you to
%       ignore params that are not relevant to the desired calculation.
%       (e.g. params.glob.ignore = {'ani.PlotFlag','rp.temporalExpect'})
% 
%  params.(analysis_step): same possible values for 'analysis_step' as
%       params.glob.process. Each analysis step requires parameters that
%       govern the execution of that particular analysis step. You can
%       choose to pass set these parameters and pass them in as a single
%       structure or a structure array into params.(analysis_step) for any
%       given step that you wish to execute in an analysis series. If you
%       do not provide any parameters for a given analysis step,
%       jlmt_proc_series will load the default parameters from that step
%       (see 'getDefaultParams' for the given step's calc_%s function) and
%       attempt to execute the job using those parameters. If you provide a
%       single parameter structure for a given job type in
%       params.(analysis_step), that parameter structure will be passed on
%       to any instance of a particular job type that is attempted. If you
%       provide an array of structures for a given job type,
%       jlmt_proc_series will try to find a structure in the array that has
%       a 'prev_steps' field that matches the previous steps in the
%       particular analysis series in question. If it finds a matching
%       parameter structure, it will pass that into the job at analysis
%       time. Otherwise, it will pass in the default parameters for the
%       given job type (see 'getDefaultParams' for the given calc function)
% 
%       specific (optional) sub-fields for each params.(analysis_step) that
%       are referenced in and used by jlmt_proc_series include:
% 
%       .prev_steps - a cell array of strings specifying jobs executed
%           before the current analysis step. If you provide an array of
%           parameter structures in params.(analysis_step), .prev_steps
%           will be used to find the particular parameter structure in the 
%           array to pass on in an instance of a job in a job series. For
%           more clarification, see the example below.
% 
%       .inDataType - if you want to specify a different input for a given
%           analysis step than the output of the previous analysis step,
%           specify the variable name of the desired input data (the
%           analysis step from a previous analysis) in
%       	params.(analysis_step).inDataType
% 
%           example: if you have the following proc series
% 
%               params.glob.process = {{'ani','pp','pc','li','toract'}}
%           
%           by default, the input to the pc step will be pp output, the
%           input to the li step will be pc output, and the input to the
%           toract step will be li output. If, however, you wanted to
%           calculate pc data but then feed the pp input to the li step,
%           you could set the following parameter:
% 
%               params.li.inDataType = 'pp'
% 
%           and instead of using the output of the previous step in the
%           series (pc) as the input to li, the output from the step
%           specified in 'params.li.inDataType' (pp) would be used.
% 
%       .metrics - if you want to calculate metrics based upon the output
%           of a particular analysis step, specify those in
%           params.(analysis_step).metrics. For each metric you would like
%           to calculate for a given analysis step, add a fieldname to
%           params.(analysis_step).metrics that matches the name of the
%           metric [calculated by metric_readout()] that you want to
%           generate. The value of that fieldname will be used as the
%           parameter input to the 'metric_readout()' function for that
%           metric calculation. See metric_readout() for more information.
% 
%           The output of metrics steps are saved to the output files
%           created after an analysis step, and they are not checked when
%           comparing current proposed analysis steps to previously run
%           analyses. Therefore, if you've run a previous analysis, haven't
%           changed the parameters of that analysis, and you want to
%           calculate new metrics based on that previously saved
%           calculation, you can enter those metrics into
%           params.(analysis_step).metrics and re-run jlmt_proc_series, and
%           your new metrics (as well as any previously calculated metrics)
%           will be returned, without having to re-run the original
%           analysis step.
% 
%       .plotfun - if you have a custom function that you want to use to
%           generate plots based on the output of an analysis step and/or
%           the output of metric calculations, you can pass in a string
%           reference to that function here. The function will be called in
%           the following manner:
% 
%           plotfun({analysis_output},params.(analysis_step))
%               - or, if metrics have been calculated -
%           plotfun({analysis_output,metric_output},params.(analysis_step))
% 
%       example:
% 
%           params.glob.process = {{'ani','pp','li','toract'},...
%               {'ani','pp','pc','toract'}};
%           params.ani = calc_ani('getDefaultParams');
%           params.pc = calc_pc('getDefaultParams');
%           params.toract(1) = calc_toract('getDefaultParams');
%           params.toract(2) = one_set_of_predefined_toract_params;
%           params.toract(2).prev_steps = {'ani','pp','pc'};
%           params.toract(3) = another_set_of_predefined_toract_params;
%           params.toract(3).prev_steps = {'ani','pp','li'};
% 
%           In the above case, when running the 'ani' and the 'pc' steps,
%           params.ani and params.pc, respectfully, will be used to conduct
%           those analyses. When 'pp' and 'li' are run in both analysis
%           series, jlmt_proc_series will dynamically execute
%           calc_pp('getDefaultParams') and calc_li('getDefaultParams'),
%           respectively, and use the parameters returned by those function
%           calls when executing calc_pp. In both analysis series, when
%           attempting to run the 'toract' analysis, jlmt_proc_series will
%           compare the previous jobs run in the series with the
%           '.prev_steps' in each params.toract case, and use the
%           parameters that match the previous jobs run in each series.
%           This will result in params.toract(3) being used in the first
%           analysis series ({'ani','pp','li','toract'}), params.toract(2)
%           being used in the second analysis series
%           ({'ani','pp','pc','toract'}), and params.toract(1) not being
%           used at all.
% 
%  params.ensemble.conn_id: a connection ID to the ensemble
%                           database. Only useful when retrieving
%                           stimuli from ensemble. If
%                           no conn_id is provided, a new database 
%                           connection will be made.
%
% OUTPUTS
%
% An output structure is returned for each job specified in
% params.glob.process. Each row of each output structure will contain the
% parameters, including a prev_steps field, that were used to generate the
% given data.
%
% Copyright (c) 2006-2011 The Regents of the University of California
% All Rights Reserved.
%
% Authors:
% December 2006, Petr Janata
% April 26, 2007, Stefan Tomic - added the ability to input different data
%                                types, including ensemble data
%                                structs. 
%                              - added ability to match parameters
%                                to load a previous calculation.
%                              - added population of missing params
%                                with default values
%                              - changed the way force_recalc is
%                                specified (now a cell array of strings)
%                              -
%
% 17, December, 2007 Stefan Tomic - added handling of mat files for
%                                   input signal
% 27, August, 2008 Stefan Tomic - Now supports a cell array of aud
%                                 structs. Also supports a cell array of mat filenames
% November 2, 2009 Stefan Tomic - Added support for ignoring parameters when
%                                param matching through params.glob.ignore
% May 6, 2011 Fred Barrett - added calc_pitchclass.m to the proc series,
%                            also cleaned up deprecated functionality
%                            handling "correlation" ... this has been moved
%                            into metrics calculated in calc_toract and
%                            calc_pitchclass, added comments, did a little
%                            bit of cleaning
% March 29, 2012 FB - added pc_ci handling (leaky integrate pc vectors)
% May 24, 2012 FB - cleaned up code, expanded to allow separate processing
% series, abstracted each job block so they can be run in the order
% specified by the processing series

%%
% Make sure IPEM setup has been run
if isempty(getpref('IPEMToolbox'))
  switch computer
    case {'PCWIN','WINPC'}
      IPEMSetup;
    otherwise
      IPEMSetupLinuxOSX;
  end
end

% get a list of all current jlmt calc functions
rootdir = fileparts(which('jlmt_proc_series'));
calcs = dir(fullfile(rootdir,'calc*m'));
calcfuncs = cell(1,length(calcs));
for k=1:length(calcs)
  lname = regexp(calcs(k).name,'^calc_(.*)\.m$','tokens');
  calcfuncs{k} = lname{1};
end % for k=1:length(calcs

% get/return default jlmt parameter structure?
if(ischar(inData) && strcmp(inData,'getDefaultParams'))
  outData = getDefaultParams('',calcfuncs);
  return
end

% This parameter sets the fieldnames to ignore when comparing the
% parameters of previously run and saved analyses with current requested
% analyses. For more information, see check_anal_exist()
if(~isfield(params,'glob') || ~isfield(params.glob,'ignore'))
  params.glob.ignore = {};
end

% If Ensemble is installed and available in the path, ensemble_globals will
% set the variables stimulus_root and stimulus_jlmt_analysis_root. These
% variables are used when processing stimuli from the Ensemble stimulus
% database. If Ensemble is not installed, these variables will not be used.
try
  ensemble_globals;
catch
  stimulus_root = '';
  stimulus_jlmt_analysis_root = '';
end

% % % % FIXME: change ipem -> jlmt

% initialize the output data structure
outData = ensemble_init_data_struct;
outData.type = 'jlmt_proc_series';
outData.name = 'jlmt_proc_series';
outData.vars = calcfuncs;
outDataCols = set_var_col_const(outData.vars);
outData.data = repmat({{[]}},size(outData.vars));

% if params were passed in, fill in any missing parameters with
% defaults. Otherwise, if no params passed in, just set them all to defaults.
if(exist('params','var'))
  params = getDefaultParams(params,calcfuncs);
else
  params = getDefaultParams('',calcfuncs);
end

%option to set log file ID in params
%if not set, then output is sent to stdout (id 1)
try
  lfid = params.glob.lfid;
catch
  lfid = 1;
end

%force_recalc global param is a cell array of strings that
%specifies which calculations should be forced to recalculate, even
%if a previous matching saved calculation was found
try
  params.glob.force_recalc;
catch
  params.glob.force_recalc = {};
end

% outputType dictates whether the output will contain the actual data, or
% a filepath
try
  params.glob.outputType;
catch
  params.glob.outputType = 'data';
end

%%  identify input format
if(ischar(inData)) 

  % a single file path as a string
  [~,~,ext] = fileparts(inData);
  switch ext
   case {'.wav','.aif','.mp3'}
    inDataType = 'aud_file_string';
   case {'.mat'}
    inDataType = 'mat_file_string';
  end
  nfiles = 1;

elseif(iscell(inData) && ischar(inData{1}))

  % cell array of filepaths as strings
  [~,~,ext] = fileparts(inData{1});
  switch lower(ext)
   case {'.wav','.aif','.mp3'}
    inDataType = 'aud_file_string_cell_array';
   case {'.mat'}
    inDataType = 'mat_file_string_cell_array';
  end
  nfiles = length(inData);

elseif(iscell(inData) && all(isfield(inData{1},{'vars','type','data'})) ...
       && strcmp(inData{1}.type,'aud'))

  % aud cell array as input
  inDataType = 'aud_cell_array';
  nfiles = length(inData);

elseif isstruct(inData) && all(isfield(inData{1},{'vars','type','data'}))

  % deal with ensemble data struct types here
  if(ismember('path',inData.vars) && ismember('filename',inData.vars))
    
    % ensemble structure with path info
    inDataType = 'file_data_struct';
    cols = set_var_col_const(inData.vars);
    nfiles = length(inData.data{cols.filename});

  elseif(strcmp(inData.type,'aud'))
    
    % ensemble structure with aud info
    inDataType = 'aud';
    cols = set_var_col_const(inData.vars);
    nfiles = length(inData);
    
  elseif(ismember('stimulus_id',inData.vars))
    
    % ensemble structure with stimulus ids : must consult Ensemble database
    inDataType = 'stimulus_id';
    cols = set_var_col_const(inData.vars);
    stimIDList = inData.data{cols.stimulus_id};
    nfiles = length(stimIDList);
    
    % check database connection
    try 
      params.ensemble.conn_id;
    catch
      params.ensemble.conn_id = 0;
      mysql_make_conn('','',params.ensemble.conn_id);
      tmp_conn_id = 1;
    end
    
    % get location field rows for the stimuli
    stimLocs = mysql_extract_data('table','stimulus','stimulus_id', ...
				  stimIDList,'extract_vars',{'location'},'conn_id',params.ensemble.conn_id);
    stimLocs = stimLocs{1};
    destStimLocs = check_stim_dirs(stimLocs,'srcroot',stimulus_root,'destroot',stimulus_ipem_analysis_root,'verbose',false);
    
  end
  
elseif(isnumeric(inData)) || (iscell(inData) && isnumeric(inData{1}))

  % if inData is numeric, assuming a list of Ensemble stim IDs
  % must consult the Ensemble database
  inDataType = 'stimulus_id';
  if iscell(inData)
    stimIDList = inData{1};
  else
    stimIDList = inData;
  end
  nfiles = length(stimIDList);
  
  % get location field rows for the stimuli
  try 
    params.ensemble.conn_id;
  catch
    params.ensemble.conn_id = 0;
    mysql_make_conn('','',params.ensemble.conn_id);
    tmp_conn_id = 1;
  end
  stimLocs = mysql_extract_data('table','stimulus','stimulus_id', ...
				stimIDList,'extract_vars',{'location'},'conn_id',params.ensemble.conn_id);
  stimLocs = stimLocs{1};
  destStimLocs = check_stim_dirs(stimLocs,'srcroot',stimulus_root,'destroot',stimulus_ipem_analysis_root,'verbose',false);

elseif(ischar(inData) && strcmp(inData,'getDefaultParams'))
    
  % if inData is a string equal to 'getDefaultParams', return the default params
  outData = getDefaultParams('',calcfuncs);
  return
  
end

%% process the list of audio files
for ifile = 1:nfiles

  switch inDataType
   case 'aud_file_string'
    [path,name,ext] = fileparts(inData);
    filename = [name ext];
    paramsSig = params_aud('path',path,'filename',filename);
    sig_st = new_aud(paramsSig);
    descript = '(single audio file)';
    
   case 'aud_file_string_cell_array'
    [path,name,ext] = fileparts(inData{ifile});
    filename = [name ext];
    paramsSig = params_aud('path',path,'filename',filename);
    sig_st = new_aud(paramsSig);
    descript = fullfile(path,filename);
   
   case {'mat_file_string','mat_file_string_cell_array'}
    % currently, mat input is only supported for rhythm profiles
    
    % doing check_stim_dirs here, since the path needs to be
    % separated from the filenames. Need to extract name and ext first
    if(strcmp(inDataType,'mat_file_string_cell_array'))
      [path,name,ext] = fileparts(inData{ifile});
    else
      [path,name,ext] = fileparts(inData);
    end
    destStimLoc = check_stim_dirs([name ext],'srcroot',path,'destroot',path);

    % obtaining the path of destStimLoc
    [path,name,ext] = fileparts(destStimLoc{1});
      
    filename = [name ext];
    paramsSig = params_aud('path',path,'filename',filename,'Fs',params.glob.insigFs,'var_name',params.glob.matVarName);
    sig_st = new_aud(paramsSig);
    descript = '(single mat file)';
    
   case 'aud'
    sig_st = inData;
    descript = 'aud structure';
    path = inData.data{cols.path};
    filename = inData.data{cols.filename};
    
   case 'aud_cell_array'
    sig_st = inData{ifile};
    cols = set_var_col_const(inData{ifile}.vars);
    descript = ['aud_struct ' num2str(ifile)];
    path = inData{ifile}.data{cols.path};
    filename = inData{ifile}.data{cols.path};
    
   case 'file_data_struct'
    path = inData.data{cols.path}{ifile};
    filename = inData.data{cols.filename}{ifile};
    paramsSig = params_aud('path',path,'filename',filename);
    sig_st = new_aud(paramsSig);
    descript = fullfile(path,filename);
    
   case 'stimulus_id'
    [path,fstub,ext] = fileparts(destStimLocs{ifile});
    filename = [fstub ext];
    paramsSig = params_aud('path',path,'filename',filename);
    sig_st = new_aud(paramsSig);
    descript = sprintf(' Stimulus ID %d ',stimIDList(ifile));
    
  end
  
  % if empty waveform, go to next audio file (or aud struct) in list
  if(exist('sig_st','var'))
    sig_st_cols = set_var_col_const(sig_st.vars);
    if isempty(sig_st.data{sig_st_cols.wvf})
      continue
    end
  end % if(exist('sig_st
  
  fprintf(lfid,'jlmt_proc_series: \n\nWorking on file %d/%d: %s\n',...
      ifile, nfiles, descript);

  % iterate over job series
  for iseries = 1:nproc
    series = params.glob.process{iseries};
    if isempty(series), continue, end
    if length(unique(series)) < length(series)
      error('duplicate steps in a series are not permitted');
    end
    series_params = struct();

    fprintf(lfid,'jlmt_proc_series: running job series %d/%d (%s)\n',...
        iseries,nproc,cell2str(series,'-'));
    
    % iterate over steps in this series
    for istep = 1:length(series)
      proc = series{istep};
      proc_fh = parse_fh(sprintf('calc_%s',proc));
      
      % get parameters for this calc step
      nprocparam = length(params.(proc));
      if nprocparam > 1
        %%% if there is more than one parameter structure provided for this
        %%% step, we will attempt to find the parameter structure for this
        %%% step that has a matching cell array of strings in the
        %%% 'prev_steps' field indicating the previous jobs used in this
        %%% series
        pmatch = cellfun(@ismember,{params.(proc).prev_steps},...
            repmat(series(1:istep),1,nprocparam),'UniformOutput',false);
        pidx = find(cellfun(@all,pmatch,1,'first'));
        if isempty(pidx)
          wstr = sprintf(['multiple param structs specified for %s, '...
              'however none matched the previous jobs for the current '...
              'series (%s). Default parameters are being used\n'],proc,...
              cell2str(series(1:istep-1),','));
          warning(wstr);
          fprintf(lfid,wstr);
          lparams = proc_fh('getDefaultParams');
        else
          if length(pidx) > 1
            warning('multiple parameter structs matched for %s, using first',...
                proc);
            pidx = pidx(1);
          end
          
          fprintf(lid,'using params index %d for %s\n',pidx,proc);
          lparams = params.(proc)(pidx);
        end
      else
        % use the only set of parameters available for this step
        lparams = params.(proc);
      end % if isfield(params,proc
      
      % save the previous steps in this proc series
      if istep > 1, lparams.prev_steps = series(1:istep); end
      
      % save this step's parameters to the series params
      series_params.(proc) = lparams;

      % do the parameters specify an input data type?
      if istep == 1
        % the first job in a series must take the stimulus signal as its
        % input
        input_data = 'sig_st';
      elseif isfield(lparams,'inDataType') && ~isempty(lparams.inDataType)
        input_data = lparams.inDataType;
      else
        input_data = series{istep-1};
      end % if isfield(lparams,'inDataType
      
      % look for previously run jobs that match this parameter set
      lfname = construct_analfname(fullfile(path,filename),proc);
      lPath = fileparts(lfname);

      % find existing analysis
      clear matchParams;
      matchParams.varName = proc;
      matchParams.paramFind = series_params;
      matchParams.ignore = [params.glob.ignore 'metrics' 'ani.aniPath'];
      previousCalcFname = check_anal_exist(lPath,matchParams);
    
      % run this step?
      if ~ismember(proc,params.glob.force_recalc) && ~isempty(previousCalcFname)
        % previous job with matching parameters exists, and this job hasn't
        % been listed for a force recalc, so load the previous results
        fprintf(lfid,['jlmt_proc_series: Loading a previously calculated '...
            '%s with matching parameters...\n'],proc);
        lfname = fullfile(lPath,previousCalcFname{1});
        load(lfname,proc)
      else
        % run this step!
        fprintf(lfid,'jlmt_proc_series: Calculating %s ...\n\n',proc);
        eval(sprintf('%s = proc_fh(%s,lparams)',proc,input_data));
     
        % save the output?
        if ismember(proc,params.glob.save_calc)
          % if a previous calculation exists and we are recalculating,
          % use previous calculation filename and overwrite it
          if ~isempty(previousCalcFname)
            fprintf(lfid,['FORCE RECALC WAS SET. Overwriting a previous' ...
                ' %s calculation\n\n'],proc);
      	    lfname = fullfile(lPath,previousCalcFname{1});
          end
        
          fprintf(lfid,'jlmt_proc_series: Saving %s to %s\n\n',proc,lfname);
          save(lfname,proc);
        end % if ismember(proc,params.glob.save_calc
      end % if ~ismember(proc,params.glob.force_recalc) && ~isempty(...
    
      % calculate metrics?
      if isfield(lparams,'metrics')
        mets_to_calc = fieldnames(lparams.metrics);
        nmet_to_calc = length(mets_to_calc);
        for imet=1:nmet_to_calc
          currmet = mets_to_calc{imet};
          
          %FIXME: check if requested metrics have already been calculated
          
          % get current metric parameters
          lmetparam = lparams.metrics.(currmet);
          
          % calculate metrics
          eval(sprintf('metric_output = metric_readout(%s,lmetparam)',proc));
          
          % FIXME: update 'meta' variable for current analysis output, to
          % reflect metrics that have been calculated
        end % for imet=1:nmet_to_calc
      end
      
      % plot function?
      if isfield(lparmas,'plotfun') && ~isempty(lparams.plotfun)
        plotfun = parse_fh(lparams.plotfun);
        eval(sprintf('plotfun(%s,lparams)',proc));
        %% FIXME: add functionality to support multiple inputs to plotfun
      end % if isfield(lparams,'plotfun...
      
      % output
      if strcmp(params.glob.outputType,'filepath') && ...
              ismember(proc,params.glob.save_calc) && ...
              exist(lfname,'file')
        outData.data{outDataCols.(proc)}{ifile,1} = lfname;
      else
        eval(sprintf('outData.data{outDataCols.%s}{ifile,1} = %s',proc,proc));
      end
    end % for istep = 1:length(series
  end % for iseries = 1:nproc

%%%% plot code from previous versions for rhythm profiler: must confirm
%%%% that it integrates well with the new coding scheme
%     if(ismember('rp',params.plot.process))
%       plotRPInput{1} = rp;
%       plotRPInput{2} = sig_st;
%       rpPlotParams = params.plot.rp;
%       rpPlotParams.fig.title = sig_st.data{sig_st_cols.tag};
%       if(isfield(params.plot,'addInputFname') && params.plot.addInputFname == 1)
% 	[dummyvar,plotFstub] = fileparts(params.plot.rp.plotFname);
% 	[dummyvar,sigFstub] = fileparts(sig_st.data{sig_st_cols.filename});
% 	rpPlotParams.plotFname = [plotFstub '_' sigFstub];
%       end
%       plot_rhythmProfile(plotRPInput,rpPlotParams);
%     end

end % for ifile=

if(exist('tmp_conn_id','var'))
  mysql(params.ensemble.conn_id,'close');
end



% % % % % % 
% % % % % %         subfunctions
% % % % % % 

function params = getDefaultParams(params,calc_names)

% getDefaultParams returns default params for each calc step where params
% were not otherwise specified

def.glob.force_recalc = {};
def.glob.ignore = {};
def.glob.save_calc = {'ani','pp','li','toract'};
def.glob.process = {{'ani','pp','li','toract'},...
    {'ani','pp','pc','li','toract'}};
def.glob.outputType = 'data';

for k=1:length(calc_names)
  fh = parse_fh(calc_names{k});
  if nargin && isstruct(params) && isfield(params,calc_names{k})
    def.(calc_names{k}) = fh('getDefaultParams',params.(calc_names{k}));
  else
    def.(calc_names{k}) = fh('getDefaultParams');
  end
end % for k=1:length(calc_names

if (nargin)
  %only assign defaults to the processes that we are performing
  if(isfield(params,'glob') && isfield(params.glob,'process'))
    defTypes = [params.glob.process{:} {'glob'}];
  else
    defTypes = fieldnames(def);    
  end
    
  for iType = 1:length(defTypes)

    type = defTypes{iType};  

    if ~isfield(params,type)
      params.(type) = struct();
    end
    
    %see if there are any fields in params that are not in
    %defs. This should be reported since it may have been caused by
    %a misspelling. All possible fields should also have defaults
    %(even if the default is simply empty)
    
    % should we be checking sub-structs of the main analysis parameters? If
    % we do this, then we need to make sure that for every toract metric
    % function we create, we update the defaults for params.metrics in
    % params_toract.m to include the defaults for the given metric. If we
    % do not want to do this, we will have to either not descend into
    % parameters for each ipem_proc_series function, or we will have to
    % find another solution to checking the defaults. - FB 01/10/2011
%     [validParams,badFields,reason] = compare_structs(params.(type),def.(type),'values',0,'substruct',0,'types',0);
    [validParams,badFields,~] = compare_structs(params.(type),def.(type),'values',0,'substruct',1,'types',0);
    if(~validParams)
      if iscell(badFields), badFields = cell2str(badFields,', '); end
      warning(['Field(s) %s for ''%s'' calculation is not specified in '...
			 'default parameters. Check spelling or add'...
			 ' the field(s) to the list of defaults\n\n'],badFields,type);
    end
    
    %populate missing param fields with defaults
    params.(type) = struct_union(params.(type),def.(type));

  end
  
else
  params = def;
end
