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
%       .plotfun - if you have a custom function that you want to use to
%           generate plots based on the output of an analysis step, you can
%           pass in a string reference to that function here. The function
%           will be called in the following manner:
% 
%           plotfun({analysis_output},params.(analysis_step))
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
% Copyright (c) 2006-2012 The Regents of the University of California
% All Rights Reserved.

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
%                            added comments, did a little bit of cleaning
% March 29, 2012 FB - added pc_ci handling (leaky integrate pc vectors)
% May 24, 2012 FB - cleaned up code, expanded to allow separate processing
% series, abstracted each job block so they can be run in the order
% specified by the processing series, calc_ci/pc_ci is now calc_li
% 26Oct2012 PJ - updated to generate default behavior in the event that
%                only a list of file names is passed in. Also changed path
%                to fpath.
% 29Oct2012 PJ - fixed handling of input as directory or single file. Now
%                calls jlmt_prepare_dirs()
% 30Oct2012 PJ - Cleaned up handling of default parameter specification
%                with regard to context specificity. Could still use some
%                cleanup of the logic associated with returning default
%                parameters.

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
calcs = dir(fullfile(rootdir,'proc','calc*m'));
calcfuncs = cell(1,length(calcs));
for k=1:length(calcs)
  lname = regexp(calcs(k).name,'^calc_(.*)\.m$','tokens');
  calcfuncs{k} = lname{1}{1};
end % for k=1:length(calcs

% get/return default jlmt parameter structure?
if ischar(inData) && strcmp(inData,'getDefaultParams')
  if exist('params','var')
    if isfield(params,'glob') && isfield(params.glob,'process') && ~isempty(params.glob.process)
      outData = getDefaultParams(params);
    else
      outData = getDefaultParams(params, calcfuncs);
    end
  else
    outData = getDefaultParams('',calcfuncs);
  end
  return
end

if ~exist('params','var')
  fprintf('No param structure specified. Using defaults ...\n')
  params = getDefaultParams;
end

% This parameter sets the fieldnames to ignore when comparing the
% parameters of previously run and saved analyses with current requested
% analyses. For more information, see check_anal_exist()
if(~isfield(params,'glob') || ~isfield(params.glob,'ignore'))
  params.glob.ignore = {};
end

% If Ensemble is installed and available in the path, ensemble_globals will
% set the variables stimulus_root and stimulus_ipem_analysis_root. These
% variables are used when processing stimuli from the Ensemble stimulus
% database. If Ensemble is not installed, these variables will not be used.
if exist('ensemble_globals','file')
  ensemble_globals;
end

% initialize the output data structure
outData = ensemble_init_data_struct;
outData.type = 'jlmt_proc_series';
outData.name = 'jlmt_proc_series';
outData.vars = ['stimulus_path' 'proc_series' calcfuncs];
outDataCols = set_var_col_const(outData.vars);
outData.data = repmat({{}},size(outData.vars));

% if params were passed in, fill in any missing parameters with
% defaults. Otherwise, if no params passed in, just set them all to defaults.
if(exist('params','var'))
  params = getDefaultParams(params);
else
  params = getDefaultParams;
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

if ~iscell(params.glob.process)
  params.glob.process = {params.glob.process};
end
if ~iscell(params.glob.process{1})
  params.glob.process = {params.glob.process};
end
nproc = length(params.glob.process);

%%  identify input format
if ischar(inData) && exist(inData)
  % Check to see if this is a directory, in which case we would try to
  % process all .wav, .aif, .mp3 and .aiff files
  if isdir(inData)
    inData = jlmt_prepare_dirs(inData);
    inDataType = 'aud_file_string_cell_array';
    nfiles = size(inData,1);
  else
    [fpath,fstub,ext] = fileparts(inData);
    switch ext
      case {'.wav','.aif','.mp3'}
        inDataType = 'aud_file_string';
      case {'.mat'}
        inDataType = 'mat_file_string';
    end
    
    % a single file path as a string
    inData = jlmt_prepare_dirs(inData);
    inData = inData{1};
    
    nfiles = 1;
  end
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

elseif isstruct(inData) && all(isfield(inData,{'vars','type','data'}))

  % deal with ensemble data struct types here
  if(ismember('stimulus_id',inData.vars))
    
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
  elseif(strcmp(inData.type,'aud'))
    
    % ensemble structure with aud info
    inDataType = 'aud';
    cols = set_var_col_const(inData.vars);
    nfiles = length(inData);
    
  elseif(ismember('path',inData.vars) && ismember('filename',inData.vars))
    
    % ensemble structure with path info
    inDataType = 'file_data_struct';
    cols = set_var_col_const(inData.vars);
    nfiles = length(inData.data{cols.filename});

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
if ~nfiles
  fprintf('No files found ...\n');
  return
end

for ifile = 1:nfiles

  switch inDataType
   case 'aud_file_string'
    [fpath,name,ext] = fileparts(inData);
    filename = [name ext];
    paramsSig = params_aud('path',fpath,'filename',filename);
    sig_st = new_aud(paramsSig);
    descript = '(single audio file)';
    
   case 'aud_file_string_cell_array'
    [fpath,name,ext] = fileparts(inData{ifile});
    filename = [name ext];
    paramsSig = params_aud('path',fpath,'filename',filename);
    sig_st = new_aud(paramsSig);
    descript = fullfile(fpath,filename);
   
   case {'mat_file_string','mat_file_string_cell_array'}
    % currently, mat input is only supported for rhythm profiles
    
    % doing check_stim_dirs here, since the path needs to be
    % separated from the filenames. Need to extract name and ext first
    if(strcmp(inDataType,'mat_file_string_cell_array'))
      [fpath,name,ext] = fileparts(inData{ifile});
    else
      [fpath,name,ext] = fileparts(inData);
    end
    destStimLoc = check_stim_dirs([name ext],'srcroot',fpath,'destroot',fpath);

    % obtaining the path of destStimLoc
    [fpath,name,ext] = fileparts(destStimLoc{1});
      
    filename = [name ext];
    paramsSig = params_aud('path',fpath,'filename',filename,'Fs',params.glob.insigFs,'var_name',params.glob.matVarName);
    sig_st = new_aud(paramsSig);
    descript = '(single mat file)';
    
   case 'aud'
    sig_st = inData;
    descript = 'aud structure';
    fpath = inData.data{cols.path};
    filename = inData.data{cols.filename};
    
   case 'aud_cell_array'
    sig_st = inData{ifile};
    cols = set_var_col_const(inData{ifile}.vars);
    descript = ['aud_struct ' num2str(ifile)];
    fpath = inData{ifile}.data{cols.path};
    filename = inData{ifile}.data{cols.path};
    
   case 'file_data_struct'
    fpath = inData.data{cols.path}{ifile};
    filename = inData.data{cols.filename}{ifile};
    paramsSig = params_aud('path',fpath,'filename',filename);
    sig_st = new_aud(paramsSig);
    descript = fullfile(fpath,filename);
    
   case 'stimulus_id'
    [fpath,fstub,ext] = fileparts(destStimLocs{ifile});
    filename = [fstub ext];
    paramsSig = params_aud('path',fpath,'filename',filename);
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
    pseries = params.glob.process{iseries};
    if isempty(pseries), continue, end
    series_params = struct();
    nseries = length(pseries);

    fprintf(lfid,'jlmt_proc_series: running job series %d/%d (%s)\n',...
        iseries,nproc,cell2str(pseries,'-'));
    
    sidx = length(outData.data{1}) + 1;
    outData.data{outDataCols.stimulus_path}{sidx,1} = fullfile(fpath,filename);
    outData.data{outDataCols.proc_series}{sidx,1} = cell2str(pseries,'->');
    
    % iterate over steps in this series
    for istep = 1:nseries
      proc = pseries{istep};
      proc_fh = parse_fh(sprintf('calc_%s',proc));
      
      % get parameters for this calc step
      nprocparam = length(params.(proc));
      if nprocparam > 1

        lparams = '';
        for iproc = 1:nprocparam
          if (istep > 1 && (isfield(params.(proc)(iproc),'prev_steps') && ...
                  compare_cells(pseries(1:istep-1),params.(proc)(iproc).prev_steps))) || ...
                  (isfield(params.(proc)(iproc),'future_steps') && ...
                  compare_cells(pseries(istep+1:nseries),params.(proc)(iproc).future_steps))
            fprintf(lfid,'using params index %d for %s\n',iproc,proc);
            lparams = params.(proc)(iproc);
            break
          elseif istep == 1 && (~isfield(params.(proc)(iproc),'prev_steps') ...
                  || isempty(params.(proc)(iproc).prev_steps))
            lparams = params.(proc)(iproc);
          end
        end % for iproc = 1:nprocparam

        if isempty(lparams)
          wstr = sprintf(['multiple param structs specified for %s, '...
              'however none matched the previous jobs for the current '...
              'series (%s). Default parameters are being used\n'],proc,...
              cell2str(pseries(1:istep-1),','));
          warning(wstr);
          fprintf(lfid,wstr);
          lparams = proc_fh('getDefaultParams','prev_steps',pseries(1:istep-1));
        end
      else
        % use the only set of parameters available for this step
        lparams = params.(proc);
      end % if isfield(params,proc
      
      % save the previous steps in this proc series
      if istep > 1, lparams.prev_steps = pseries(1:istep); end
      
      % save this step's parameters to the series params
      if isfield(series_params,proc)
        series_params.(proc)(end+1) = lparams;
      else
        series_params.(proc) = lparams;
      end

      % do the parameters specify an input data type?
      if istep == 1
        % the first job in a series must take the stimulus signal as its
        % input
        input_data = 'sig_st';
      elseif isfield(lparams,'inDataType') && ~isempty(lparams.inDataType)
        input_data = lparams.inDataType;
      else
        input_data = pseries{istep-1};
      end % if isfield(lparams,'inDataType
      
      % look for previously run jobs that match this parameter set
      lfname = construct_analfname(fullfile(fpath,filename),proc);
      lPath = fileparts(lfname);

      % find existing analysis
      clear matchParams;
      matchParams.varName = proc;
      matchParams.paramFind = series_params;
      matchParams.ignore = [params.glob.ignore 'ani.aniPath' 'future_steps'];
      previousCalcFname = check_anal_exist(lPath,matchParams);
    
      % run this step?
      if ~ismember(proc,params.glob.force_recalc) && ~isempty(previousCalcFname)
        % previous job with matching parameters exists, and this job hasn't
        % been listed for a force recalc, so load the previous results
        fprintf(lfid,['jlmt_proc_series: Loading a previously calculated '...
            '%s with matching parameters...\n'],proc);
        lfname = fullfile(lPath,previousCalcFname{1});
        load(lfname,proc)
      elseif eval(sprintf('isempty(%s)',input_data))
        fprintf('Input data (%s) for %s is empty! Skipping ...\n', input_data, proc);
        continue
      else
        % run this step!
        fprintf(lfid,'jlmt_proc_series: Calculating %s ...\n\n',proc);
        eval(sprintf('%s = proc_fh(%s,lparams);',proc,input_data));
     
        % save the output?
        if ismember(proc,params.glob.save_calc{iseries})
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
      
      % plot function?
      if isfield(lparams,'plotfun') && ~isempty(lparams.plotfun)
        plotfun = parse_fh(lparams.plotfun);
        eval(sprintf('plotfun(%s,lparams);',proc));
        %% FIXME: add functionality to support multiple inputs to plotfun
      end % if isfield(lparams,'plotfun...
      
      % output
      if strcmp(params.glob.outputType,'filepath') && ...
              ismember(proc,params.glob.save_calc) && ...
              exist(lfname,'file')
        outData.data{outDataCols.(proc)}{sidx,1} = lfname;
      else
        eval(sprintf('outData.data{outDataCols.%s}{sidx,1} = %s;',proc,proc));
      end
    end % for istep = 1:nseries

%%% FIXME: plot code from previous versions for rhythm profiler: must
%%% confirm that it integrates well with the new coding scheme
    if(ismember('rp',pseries)) && isfield(params.rp,'plot')
      rpPlotParams = params.rp.plot;
      rpPlotParams.fig.title = sig_st.data{sig_st_cols.tag};
      if(isfield(params.rp.plot,'addInputFname') && params.rp.plot.addInputFname == 1)
    	[dummyvar,plotFstub] = fileparts(params.rp.plot.plotFname);
        [dummyvar,sigFstub] = fileparts(sig_st.data{sig_st_cols.filename});
    	rpPlotParams.plotFname = [plotFstub '_' sigFstub];
      end
      plot_rhythmProfile({rp,sig_st},rpPlotParams);
    end

  end % for iseries = 1:nproc
end % for ifile=

% fill out short columns
for k=1:length(outData.vars)
  if length(outData.data{k}) < sidx
    outData.data{k}{sidx,1} = [];
  end
end

if(exist('tmp_conn_id','var'))
  mysql(params.ensemble.conn_id,'close');
end



% % 
%           subfunctions


function params = getDefaultParams(params,calc_names)
% getDefaultParams returns default params for each calc step where params
% were not otherwise specified.
%
% We need to handle two situations here. One situation is that in which
% no parameters have been passed in whatsoever and we have to populate
% everything with default parameters.  The other is the situation in which
% a partial set of parameters has been passed in and we need to complete
% the set.

VERBOSE = 0;  % flag for debugging

def.glob.force_recalc = {};
def.glob.ignore = {};
def.glob.save_calc = {'ani','pp','li','toract','pc','rp'};
def.glob.process = {{'ani','pp','li','toract'},...
    {'ani','pp','pc','li','toract'}, ...
    {'ani','rp'}};
def.glob.outputType = 'data';

% Check to see whether we have a params struct and initialize any defaults
% as necessary
if ~nargin || ~isstruct(params)
  params = struct;
end

% Check to see if we have the various globals fields and if not, copy the
% default
checkGlobTypes = {'process','save_calc','outputType','force_recalc','ignore'};
for iglob = 1:length(checkGlobTypes)
  currType = checkGlobTypes{iglob};
  if ~isfield(params,'glob') || ~isfield(params.glob, currType)
    params.glob.(currType)= def.glob.(currType);
  end
end

% Clear the defaults
clear def

if exist('calc_names','var')
  if iscell(calc_names{1})
    params.glob.process = calc_names;
  else
    params.glob.process = {calc_names};
  end
end

if ~iscell(params.glob.process{1})
  params.glob.process = {params.glob.process};
end

numSeries = length(params.glob.process);
for iseries = 1:numSeries
  procCalc = params.glob.process{iseries};
  for k=1:length(procCalc)
    currCalc = procCalc{k};
    fh = parse_fh(['calc_' currCalc]);
    
    % Determine whether any parameters have been specified for this
    % calculation whatsoever, and if not, get the defaults
    if ~isfield(params, currCalc)
      params.(currCalc) = fh('getDefaultParams','prev_steps',procCalc(1:k-1));
    else
      % If a field for the current calculation already exists, determine
      % whether any instances of it match the current context

      nspec = length(params.(currCalc));
      contextIdx = [];
      nullContextIdx = [];
      for ispec=1:nspec
        % Copy the current spec
        currSpec = params.(currCalc)(ispec);
        
        if ~isfield(currSpec, 'prev_steps') || isempty(currSpec.prev_steps)
          currSpec.prev_steps = {};
          params.(currCalc)(ispec).prev_steps = {};
          nullContextIdx = ispec;
        end
        
        % If the context matches, then assume this is the one we want;
        % otherwise get a new spec for this context. If the context was
        % empty, assume that we want to use those parameters throughout
        if length(currSpec.prev_steps) == length(procCalc(1:k-1)) ...
            && k > 1 && all(strcmp(currSpec.prev_steps, procCalc(1:k-1)))
          contextIdx = ispec;
          break     
        end
      end % for l=1:nspec
      
      if isempty(contextIdx)
        if ~isempty(nullContextIdx)
          % Copy this as a default spec)
          defSpec = params.(currCalc)(nullContextIdx);
          defSpec.prev_steps = procCalc(1:k-1);
          extraParams = [fieldnames(defSpec) struct2cell(defSpec)]';
        else
          extraParams = {'prev_steps',procCalc(1:k-1)};
        end
          
        tmpParams = fh('getDefaultParams',extraParams{:});
        params.(currCalc)(end+1) = tmpParams;
      end
      
    end % if ~isfield(params, currCalc)
    
    if 0
      tmpdef = fh('getDefaultParams','prev_steps',procCalc(1:k-1));
      if ~isfield(params, currCalc)
        params.(currCalc) = tmpdef;
      else
        % compare structures and see if we need to append this spec
        
        % Set fields to ignore
        switch currCalc
          case 'ani'
            ignoreFlds = 'aniPath';
          otherwise
            ignoreFlds = {};
        end
        
        [structsMatch, firstViolation, violationReason] = ...
          compare_structs(params.(currCalc), tmpdef, 'ignore_fieldnames',ignoreFlds, 'verbose', VERBOSE);
        if structsMatch
          % Perform any cleanup necessary due to redundancy
          switch currCalc
            case 'ani'
              if exist(tmpdef.aniPath,'dir'), rmdir(tmpdef.aniPath), end
          end
        else
          
          % Append this new definition
          params.(currCalc)(end+1) = tmpdef;
        end % if compare_structs
      end % if ~isfield(def, currCalc)
    end % if 0
  end % for k=1:length(calc_names
end % for iseries


if 0

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
  
  [validParams,badFields,~] = compare_structs(params.(type),def.(type),'values',0,'substruct',1,'types',0);
  if(~validParams)
    if iscell(badFields), badFields = cell2str(badFields,', '); end
    warning(['Field(s) %s for ''%s'' calculation is not specified in '...
      'default parameters. Check spelling or add'...
      ' the field(s) to the list of defaults\n\n'],badFields,type);
  end
  
  %populate missing param fields with defaults
  nspec = length(def.(type));
  for l=1:nspec
    if l == 1
      tmpstruct = struct_union(params.(type)(l),def.(type)(l));
    else
      tmpstruct(l) = struct_union(params.(type)(l),def.(type)(l));
    end
  end % for l=1:nparam
  params.(type) = tmpstruct;
end

end % if 0

return