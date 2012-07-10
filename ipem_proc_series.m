function outData = ipem_proc_series(inData,params)
% Performs a series of calculations dependent on the IPEM toolbox
% 
% outData = ipem_proc_series(inData,params)
%
%
% INPUTS
%  inData: this could be one of many different input types:
%           1) an ensemble data struct with a "stimulus_id"
%              variable, which specifies a list of stimulus IDs in ensemble
%           2) a vector of stimulus IDs in ensemble
%           3) a cell array of 'aud' structs (see init_aud.m  and
%              new_aud.m )
%           4) an ensemble data struct with filename and path variables
%           5) a cell array of filenames, with paths included
%           6) a string specifying a single path and filename (for
%              processing just one input)
%
%  params.glob.lfid: Specifies a log file ID (obtained from
%                    fopen). All standard output will be sent to this log file. If
%                    this parameter is not set, the output will be sent to the Matlab console.
%
%  params.glob.process: possible values: 'ani','pp','ci','toract','pc','rp'
%                       a cell array of strings specifying which
%                       calculations to perform. If any dependent calculations are not
%                       specified, they will run anyway.
%
%  params.glob.force_recalc: same possible values as params.glob.process.
%  			     cell array of strings that specify
%                            which jobs will force a recalculation, even if a saved
%                            calculation with matching parameters is found. The matching
%                            saved calculation will be overwritten if params.glob.save_calc
%                            is set.
%
%  params.glob.ignore: a cell array of param fields to ignore when finding a match for
%                      an identical saved calculation 
%                      (e.g. params.glob.ignore = {'ani.PlotFlag','rp.temporalExpect'})
%                      This allows ignoring params that are not relevant to the desired
%                      calculation
%
%  params.ani:    a structure containing ani params (see params_ani.m)
%  params.ci:     structure containing contextuality index
%                 parameters (see params_ci.m)
%  params.pp:     structure containing periodicity pitch analysis
%                 parameters (see params_pp.m)
%  params.toract: structure containing torus activation parameters
%                 (see calc_toract.m)
%  params.pc:     structure containing pitch class projection parameters
%                 (see params_pc.m)
%  params.pc_ci:  structure containing contextuality index parameters for
%                 application to pitch class vectors (see params_ci.m)
%  params.rp      parameter structure for rhythm_profiler (Beyond
%                 the Beat paper code).See
%                 (params_rhythm_profiler.m)
%  params.ensemble.conn_id: a connection ID to the ensemble
%                           database. Only useful when retrieving
%                           stimuli from ensemble. If
%                           no conn_id is provided, a new database 
%                           connection will be made.
%
% if 'getDefaultParams' tag was passed as inData, the default
% parameter structure is returned
%
%
% OUTPUTS
%
%
% The calculations to be performed are specified in the params structure
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

if(ischar(inData) && strcmp(inData,'getDefaultParams'))
  outData = getDefaultParams;
  return
end

if(~isfield(params,'glob') || ~isfield(params.glob,'ignore'))
  params.glob.ignore = {};
end

%need to catch an error in case ensemble_globals is not in the path
%this provides for public use of ipem_proc_series (without Ensemble)
if exist('ensemble_globals','file')
  ensemble_globals
end

outData = ensemble_init_data_struct;
outData.type = 'ipem_proc_series';
outData.name = 'ipem_proc_series';
outData.vars = {'ani','pp','ci','toract','pc','pc_ci','rp'};
outDataCols = set_var_col_const(outData.vars);
outData.data = repmat({{[]}},size(outData.vars));

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

%outputType dictates whether the output will contain the actual data, or
%a filepath
try
  params.glob.outputType;
catch
  params.glob.outputType = 'data';
end

try
  params.plot.process;
catch
  params.plot.process = {};
end

%Make sure all prerequisite analyses are performed.

if(ismember('pp',params.glob.process))
  params.glob.process = union(params.glob.process,'ani');
end
  
if(ismember('ci',params.glob.process))
  params.glob.process = union(params.glob.process,{'ani','pp'});
end

if(ismember('toract',params.glob.process))
  params.glob.process = union(params.glob.process,{'ani','pp','ci'});
end

if(ismember('pc',params.glob.process))
  params.glob.process = union(params.glob.process,{'ani','pp'});
end

if(ismember('pc_ci',params.glob.process))
  params.glob.process = union(params.glob.process,{'ani','pp','pc'});
end

%% find out what type of inData was passed
if isstruct(inData) && isfield(inData,'vars') && isfield(inData,'type')
  %deal with ensemble data struct types here
  
  if(ismember('stimulus_id',inData.vars))
    inDataType = 'stimulus_id';
    cols = set_var_col_const(inData.vars);
    stimIDList = inData.data{cols.stimulus_id};
    if isfield(cols,'location')
      stimLocs = inData.data{cols.location};
    else
      % check database settings
      try 
        params.ensemble.conn_id;
      catch
        params.ensemble = mysql_login();
        mysql_make_conn(params.ensemble);
        tmp_conn_id = 1;
      end
      
      %get location field rows for the stimuli
      stimLocs = mysql_extract_data('table','stimulus','stimulus_id', ...
	  			    stimIDList,'extract_vars',{'location'},'conn_id',params.ensemble.conn_id);
      stimLocs = stimLocs{1};
    end
    destStimLocs = check_stim_dirs(stimLocs,'srcroot',stimulus_root,'destroot',stimulus_ipem_analysis_root,'verbose',false);
    nfiles = length(destStimLocs);
    
  elseif(strcmp(inData.type,'aud'))
    inDataType = 'aud';
    cols = set_var_col_const(inData.vars);
    nfiles = length(inData);
    
  elseif(ismember('path',inData.vars) && ismember('filename',inData.vars))
    inDataType = 'file_data_struct';
    cols = set_var_col_const(inData.vars);
    nfiles = length(inData.data{cols.filename});
  end
elseif(iscell(inData) && all(isfield(inData{1},{'vars','type','data'})) ...
       && strcmp(inData{1}.type,'aud'))
  inDataType = 'aud_cell_array';
  nfiles = length(inData);
elseif(iscell(inData) && ischar(inData{1}))
  %cell array of filepaths as strings
  [~,~,ext] = fileparts(inData{1});
  switch lower(ext)
   case {'.wav','.aif','.mp3'}
    inDataType = 'aud_file_string_cell_array';
   case {'.mat'}
    inDataType = 'mat_file_string_cell_array';
  otherwise
    error('unknown indata type for extension %s',ext);
  end
  nfiles = length(inData);
    
elseif(ischar(inData) && strcmp(inData,'getDefaultParams'))
  %if inData is a string equal to 'getDefaultParams', return
  %the default params
  outData = getDefaultParams;
  return
elseif(ischar(inData)) 
  %a single file path as a string
  [~,~,ext] = fileparts(inData);
  switch ext
   case {'.wav','.aif','.mp3'}
    inDataType = 'aud_file_string';
   case {'.mat'}
    inDataType = 'mat_file_string';
  otherwise
    error('unknown indata type for extension %s',ext);
  end
  nfiles = 1;
  
elseif(isnumeric(inData)) || (iscell(inData) && isnumeric(inData{1}))
    
  %if inData is numeric, assuming a list of stim IDs
  inDataType = 'stimulus_id';
  if iscell(inData)
    stimIDList = inData{1};
  else
    stimIDList = inData;
  end
  nfiles = length(stimIDList);
  
  %get location field rows for the stimuli
  try 
    params.ensemble.conn_id;
  catch
    params.ensemble = mysql_login();
    mysql_make_conn(params.ensemble);
    tmp_conn_id = 1;
  end
  stimLocs = mysql_extract_data('table','stimulus','stimulus_id', ...
				stimIDList,'extract_vars',{'location'},'conn_id',params.ensemble.conn_id);
  stimLocs = stimLocs{1};
  destStimLocs = check_stim_dirs(stimLocs,'srcroot',stimulus_root,'destroot',stimulus_ipem_analysis_root,'verbose',false);

end

%% process the list of audio files
for ifile = 1:nfiles

  switch inDataType
   case 'stimulus_id'
    [path,fstub,ext] = fileparts(destStimLocs{ifile});
    filename = [fstub ext];
    paramsSig = params_aud('path',path,'filename',filename);
    sig_st = new_aud(paramsSig);
    descript = sprintf(' Stimulus ID %d ',stimIDList(ifile));
    
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
    descript = sprintf('%s/%s',path,filename);
    
   case 'aud_file_string_cell_array'
    [path,name,ext] = fileparts(inData{ifile});
    filename = [name ext];
    paramsSig = params_aud('path',path,'filename',filename);
    sig_st = new_aud(paramsSig);
    descript = sprintf('%s/%s',path,filename);
   
   case 'aud_file_string'
    [path,name,ext] = fileparts(inData);
    filename = [name ext];
    paramsSig = params_aud('path',path,'filename',filename);
    sig_st = new_aud(paramsSig);
    descript = '(single audio file)';
    
   case {'mat_file_string','mat_file_string_cell_array'}
    %currently, mat input is only supported for rhythm profiles
    
    %doing check_stim_dirs here, since the path needs to be
    %separated from the filenames. Need to extract name and ext first
    if(strcmp(inDataType,'mat_file_string_cell_array'))
      [path,name,ext] = fileparts(inData{ifile});
    else
      [path,name,ext] = fileparts(inData);
    end
    destStimLoc = check_stim_dirs([name ext],'srcroot',path,'destroot',path);
    %obtaining the path of destStimLoc
    [path,name,ext] = fileparts(destStimLoc{1});
      
    filename = [name ext];
    paramsSig = params_aud('path',path,'filename',filename,'Fs',params.glob.insigFs,'var_name',params.glob.matVarName);
    sig_st = new_aud(paramsSig);
    descript = '(single mat file)';
    
  end
  
  %if empty waveform, go to next audio file (or aud struct) in list
  if(exist('sig_st','var'))
    sig_st_cols = set_var_col_const(sig_st.vars);
    if isempty(sig_st.data{sig_st_cols.wvf})
      continue
    end
  end % if(exist('sig_st
  
  fprintf(lfid,'ipem_proc_series: \n\nWorking on file %d/%d: %s\n', ifile, nfiles, descript);

% %     calculate auditory nerve image
  if( ismember('ani',params.glob.process))
    
    %if path and filename are not empty, we can check for previous
    %analyses, otherwise, this step is skipped and we can also not
    %save any analyses
    if(~isempty(path) && ~isempty(filename))
      anifname = construct_analfname(fullfile(path,filename),'ani');
      aniPath = fileparts(anifname);

      %find existing analysis
      clear matchParams;
      matchParams.varName = 'ani';
      matchParams.paramFind.ani = params.ani;
      matchParams.ignore = [params.glob.ignore {'ani.aniPath'}];
      previousCalcFname = check_anal_exist(aniPath,matchParams);
    else
      previousCalcFname = {};
    end
      
    force_recalc_ani = ismember('ani',params.glob.force_recalc);
    
    if ~force_recalc_ani && ~isempty(previousCalcFname)
      fprintf(lfid,'ipem_proc_series: Loading a previously calculated ANI with matching parameters...\n');
      anifname = fullfile(aniPath,previousCalcFname{1});
      load(anifname,'ani')
    else
      % if a previous calculation exists and we are recalculating,
      % use previous calculation filename and overwrite it
      if ~isempty(previousCalcFname)
	fprintf(lfid,['FORCE RECALC WAS SET. Overwriting a previous' ...
		      ' ANI calculation\n\n']);
	anifname = fullfile(aniPath,previousCalcFname{1});
      end
      fprintf(lfid,'ipem_proc_series: Calculating ANI ...\n\n');
      ani = calc_ani(sig_st,params.ani);
     
      if ismember('ani',params.glob.save_calc)
	
	if(exist('anifname','var'))
	  fprintf(lfid,'ipem_proc_series: Saving ANI to %s\n\n',anifname);
	  save(anifname,'ani');
	else
	  fprintf(lfid,'ipem_proc_series: Cannot save ANI. Filename not specified\n\n');
	end
	
      end
    end
 
    switch (params.glob.outputType)
     case 'data'
      outData.data{outDataCols.ani}{ifile,1} = ani;
     case 'filepath'
      outData.data{outDataCols.ani}{ifile,1} = anifname;
    end
      
  end %if proc ani

% %   calculate periodicity pitch image
  if( ismember('pp',params.glob.process) )  
    if(~isempty(path) && ~isempty(filename))
      ppfname = construct_analfname(fullfile(path,filename),'pp');
      ppPath = fileparts(ppfname);
    
      %find existing analysis
      clear matchParams;
      matchParams.varName = 'pp';
      matchParams.paramFind.ani = params.ani;
      matchParams.paramFind.pp  = params.pp;
      matchParams.ignore = [params.glob.ignore {'ani.aniPath'}];
      previousCalcFname = check_anal_exist(ppPath,matchParams);
    else      
      previousCalcFname = {};
    end
    
    force_recalc_pp = ismember('pp',params.glob.force_recalc);
    
    if ~force_recalc_pp && ~isempty(previousCalcFname)
      ppfname = fullfile(ppPath,previousCalcFname{1});
      fprintf(lfid,['Loading a previously calculated Periodicity Pitch Image' ...
	    ' with matching parameters...\n\n']);
      load(ppfname,'pp')
    else
      
      % if a previous calculation exists and we are recalculating,
      % use previous calculation filename and overwrite it
      if ~isempty(previousCalcFname)
	fprintf(lfid,['FORCE RECALC WAS SET. Overwriting a previous' ...
		      ' periodicity pitch calculation\n\n']);
	ppfname = fullfile(ppPath,previousCalcFname{1});
      end
      
      fprintf(lfid,'ipem_proc_series: Calculating periodicity pitch image ...\n\n');
      pp = calc_pp(ani,params.pp);
      if ismember('pp',params.glob.save_calc)
	
	if(exist('ppfname','var'))
	  fprintf(lfid,'ipem_proc_series: Saving periodicity pitch image to %s\n\n',ppfname);
	  save(ppfname,'pp')
	else
	  fprintf(lfid,'ipem_proc_series: Cannot save pitch image. Filename not specified.\n\n');
	end
	
      end
    end
    
    switch (params.glob.outputType)
     case 'data'
      outData.data{outDataCols.pp}{ifile,1} = pp;
     case 'filepath'
      outData.data{outDataCols.pp}{ifile,1} = ppfname; 
    end
    
  end %if proc pp
    
% %     calculate context images
  if( ismember('ci',params.glob.process) )
    
    if(~isempty(path) && ~isempty(filename))
      cifname = construct_analfname(fullfile(path,filename),'ci');
      ciPath = fileparts(cifname);
    
      %find existing analysis
      clear matchParams;
      matchParams.varName = 'cntxt_img';
      matchParams.paramFind.ani = params.ani;
      matchParams.paramFind.pp  = params.pp;
      matchParams.paramFind.ci  = params.ci;
      matchParams.ignore = [params.glob.ignore {'ani.aniPath'}];
      previousCalcFname = check_anal_exist(ciPath,matchParams);
    else
      previousCalcFname = {};
    end
    
    force_recalc_ci = ismember('ci',params.glob.force_recalc);
    
    if ~force_recalc_ci && ~isempty(previousCalcFname)
      cifname = fullfile(ciPath,previousCalcFname{1});
      fprintf(lfid,['Loading a previously calculated Context Image' ...
	    ' with matching parameters...\n\n']);
      load(cifname,'cntxt_img');
       
    else
      
      % if a previous calculation exists and we are recalculating,
      % use previous calculation filename and overwrite it
      if ~isempty(previousCalcFname)
	fprintf(lfid,['FORCE RECALC WAS SET. Overwriting a previous' ...
		      ' context image calculation\n\n']);
	cifname = fullfile(ciPath,previousCalcFname{1});
      end

      fprintf(lfid,'ipem_proc_series: Calculating context images ...\n\n');
      cntxt_img = calc_context(pp,params.ci);
      if ismember('ci',params.glob.save_calc)
	
	if(exist('cifname','var'))
	  fprintf(lfid,'ipem_proc_series: Saving context image to %s\n\n',cifname);
	  save(cifname,'cntxt_img')
	else
	  fprintf(lfid,'ipem_proc_series: Cannot save context image. Filename not specified.\n\n');
	end
	
      end
    end
    
    switch (params.glob.outputType)
     case 'data'
      outData.data{outDataCols.ci}{ifile,1} = cntxt_img;
     case 'filepath'
      outData.data{outDataCols.ci}{ifile,1} = cifname;
    end
    
  end %if proc ci
  
% %     project context images to the torus
  if( ismember('toract',params.glob.process) )  
   
    if(~isempty(path) && ~isempty(filename))
      toract_fname = construct_analfname(fullfile(path,filename),'toract');
      toractPath = fileparts(toract_fname);
    
      %find existing analysis
      clear matchParams;
      matchParams.varName = 'toract';
      matchParams.paramFind.ani = params.ani;
      matchParams.paramFind.pp  = params.pp;
      matchParams.paramFind.ci  = params.ci;
      matchParams.paramFind.toract  = params.toract;
      matchParams.ignore = [params.glob.ignore {'ani.aniPath' 'toract.ci' 'toract.Fs'}];
      previousCalcFname = check_anal_exist(toractPath,matchParams);
    else
      previousCalcFname = {};
    end
    
    force_recalc_toract = ismember('toract',params.glob.force_recalc);
    
    if (~force_recalc_toract && ~isempty(previousCalcFname))
      toract_fname = fullfile(toractPath,previousCalcFname{1});
      fprintf(lfid,'ipem_proc_series: Loading a previously calculated torus activation with matching parameters...\n');
      load(toract_fname,'toract')
    else

      % if a previous calculation exists and we are recalculating,
      % use previous calculation filename and overwrite it
      if ~isempty(previousCalcFname)
	fprintf(lfid,['FORCE RECALC WAS SET. Overwriting a previous' ...
		      ' torus activation calculation\n\n']);
	toract_fname = fullfile(toractPath,previousCalcFname{1});
      end

      fprintf(lfid,'ipem_proc_series: Calculating projection to torus ...\n\n');
      toract = calc_toract(cntxt_img, params.toract);
      if ismember('toract',params.glob.save_calc)
	
	if(exist('toract_fname','var'))
	  fprintf(lfid,'ipem_proc_series: Saving torus activation image to %s\n\n',toract_fname);
	  save(toract_fname,'toract','params')
	else
	  fprintf(lfid,['Cannot save torus activation image. Filename' ...
			' not specified.\n\n.']);
	end
      end
    end
    
    switch (params.glob.outputType)
     case 'data'
      outData.data{outDataCols.toract}{ifile,1} = toract;
     case 'filepath'
      outData.data{outDataCols.toract}{ifile,1} = toract_fname;
    end
    
  end %if proc toract
  
% %     project periodicity pitch or context images to pitch class space  
  if( ismember('pc',params.glob.process) )  
   
    if(~isempty(path) && ~isempty(filename))
      pc_fname = construct_analfname(fullfile(path,filename),'pc');
      pcPath = fileparts(pc_fname);
    
      %find existing analysis
      clear matchParams;
      matchParams.varName = 'pc';
      matchParams.paramFind.ani = params.ani;
      matchParams.paramFind.pp  = params.pp;
      if isfield(params.pc,'ci_siglist') && ~isempty(params.pc.ci_siglist)
        matchParams.paramFind.ci  = params.ci;
      end
      matchParams.paramFind.pc  = params.pc;
      matchParams.ignore = [params.glob.ignore {'ani.aniPath' 'pc.ci' 'pc.Fs' }];
      previousCalcFname = check_anal_exist(pcPath,matchParams);
    else
      previousCalcFname = {};
    end
    
    force_recalc_pc = ismember('pc',params.glob.force_recalc);
    
    if (~force_recalc_pc && ~isempty(previousCalcFname))
      pc_fname = fullfile(pcPath,previousCalcFname{1});
      fprintf(lfid,'ipem_proc_series: Loading a previously calculated pitch class projection with matching parameters...\n');
      load(pc_fname,'pc')
    else

      % if a previous calculation exists and we are recalculating,
      % use previous calculation filename and overwrite it
      if ~isempty(previousCalcFname)
	fprintf(lfid,['FORCE RECALC WAS SET. Overwriting a previous' ...
		      ' pitch class projection\n\n']);
	pc_fname = fullfile(pcPath,previousCalcFname{1});
      end

      fprintf(lfid,'ipem_proc_series: Projecting to pitch class space ...\n\n');
      
      if isfield(params.pc,'ci_siglist') && ~isempty(params.pc.ci_siglist)
        pc = calc_pitchclass(cntxt_img,params.pc);
      else
        pc = calc_pitchclass(pp,params.pc);
      end
      
      if ismember('pc',params.glob.save_calc)
	
	if(exist('pc_fname','var'))
	  fprintf(lfid,'ipem_proc_series: Saving pitch class projections to %s\n\n',pc_fname);
	  save(pc_fname,'pc','params')
	else
	  fprintf(lfid,['Cannot save pitch class projections. Filename' ...
			' not specified.\n\n.']);
	end
      end
    end
    
    switch (params.glob.outputType)
     case 'data'
      outData.data{outDataCols.pc}{ifile,1} = pc;
     case 'filepath'
      outData.data{outDataCols.pc}{ifile,1} = pc_fname;
    end
    
  end %if proc pc

% %     calculate context images from pitch class vectors
  if( ismember('pc_ci',params.glob.process) )
    
    if(~isempty(path) && ~isempty(filename))
      pccifname = construct_analfname(fullfile(path,filename),'pc_ci');
      pcciPath = fileparts(pccifname);
    
      %find existing analysis
      clear matchParams;
      matchParams.varName = 'pc_cntxt_img';
      matchParams.paramFind.ani = params.ani;
      matchParams.paramFind.pp  = params.pp;
      matchParams.ignore = [params.glob.ignore {'ani.aniPath'}];
      previousCalcFname = check_anal_exist(pcciPath,matchParams);
    else
      previousCalcFname = {};
    end
    
    force_recalc_pc_ci = ismember('pc_ci',params.glob.force_recalc);
    
    if ~force_recalc_pc_ci && ~isempty(previousCalcFname)
      pccifname = fullfile(pcciPath,previousCalcFname{1});
      fprintf(lfid,['Loading a previously calculated PC Context Image' ...
	    ' with matching parameters...\n\n']);
      load(pccifname,'pc_cntxt_img');
    else
      
      % if a previous calculation exists and we are recalculating,
      % use previous calculation filename and overwrite it
      if ~isempty(previousCalcFname)
	fprintf(lfid,['FORCE RECALC WAS SET. Overwriting a previous' ...
		      ' pc context image calculation\n\n']);
	pccifname = fullfile(pcciPath,previousCalcFname{1});
      end

      fprintf(lfid,'ipem_proc_series: Calculating pc context images ...\n\n');
      pc_cntxt_img = calc_context(pc,params.pc_ci);
      if ismember('pc_ci',params.glob.save_calc)
	
	if(exist('pccifname','var'))
	  fprintf(lfid,'ipem_proc_series: Saving pc context image to %s\n\n',pccifname);
	  save(pccifname,'pc_cntxt_img')
	else
	  fprintf(lfid,'ipem_proc_series: Cannot save pc context image. Filename not specified.\n\n');
	end
	
      end
    end
    
    switch (params.glob.outputType)
     case 'data'
      outData.data{outDataCols.pc_ci}{ifile,1} = pc_cntxt_img;
     case 'filepath'
      outData.data{outDataCols.pc_ci}{ifile,1} = pccifname;
    end
    
  end %if proc ci
  
% %     calculate the rhythm profile    
  if( ismember('rp',params.glob.process) )  
    
    if(~isempty(path) && ~isempty(filename))
      rp_fname = construct_analfname(fullfile(path,filename),'rp');
      rpPath = fileparts(rp_fname);
      
      %find existing analysis
      clear matchParams;
      matchParams.varName = 'rp';
      
      if(strcmp(params.rp.inDataType,'ani'))
	matchParams.paramFind.ani = params.ani;
	matchParams.ignore = [params.glob.ignore {'ani.aniPath'}];
      else
	matchParams.ignore = params.glob.ignore;
      end
      
      matchParams.paramFind.rp  = params.rp;
      previousCalcFname = check_anal_exist(rpPath,matchParams);
    else
      previousCalcFname = {};
    end
    
    force_recalc_rp = ismember('rp',params.glob.force_recalc);
    
    if (~force_recalc_rp && ~isempty(previousCalcFname))
      rp_fname = fullfile(rpPath,previousCalcFname{1});
      fprintf(lfid,'ipem_proc_series: Loading a previously calculated rhythm profile with matching parameters...\n');
      load(rp_fname,'rp')
    else
      
      % if a previous calculation exists and we are recalculating,
      % use previous calculation filename and overwrite it
      if ~isempty(previousCalcFname)
	fprintf(lfid,['ipem_proc_series: FORCE RECALC WAS SET. Overwriting a previous' ...
		      ' rhythm profile calculation\n\n']);
	rp_fname = fullfile(rpPath,previousCalcFname{1});
      end
      
      fprintf(lfid,'ipem_proc_series: Calculating rhythm profile ...\n\n');

      
      %if ani was not processed, passing the sig_st directly
      %through the rhythm profiler. Otherwise, pass the ani

      if(ismember('ani',params.glob.process))
	fprintf(lfid,'ipem_proc_series: passing ANI through rhythm profiler.\n\n');
	rp = rhythm_profiler_v2(ani,params.rp);
      else
	fprintf(lfid,'ipem_proc_series: passing signal directly through rhythm profiler.\n\n');
	rp = rhythm_profiler_v2(sig_st,params.rp);
      end
      
      if ismember('rp',params.glob.save_calc)
	
	if(exist('rp_fname','var'))
	  fprintf(lfid,'ipem_proc_series: Saving rhythm profile to %s\n\n',rp_fname);
	  save(rp_fname,'rp','params')
	else
	  fprintf(lfid,'ipem_proc_series: Cannot save rhythm profile. Filename not specified.\n\n');
	end
      end
      
      
    end

    if(ismember('rp',params.plot.process))
      plotRPInput{1} = rp;
      plotRPInput{2} = sig_st;
      rpPlotParams = params.plot.rp;
      rpPlotParams.fig.title = sig_st.data{sig_st_cols.tag};
      if(isfield(params.plot,'addInputFname') && params.plot.addInputFname == 1)
	[dummyvar,plotFstub] = fileparts(params.plot.rp.plotFname);
	[dummyvar,sigFstub] = fileparts(sig_st.data{sig_st_cols.filename});
	rpPlotParams.plotFname = [plotFstub '_' sigFstub];
      end
      plot_rhythmProfile(plotRPInput,rpPlotParams);
    end
    
    switch (params.glob.outputType)
     case 'data'
      outData.data{outDataCols.rp}{ifile,1} = rp;
     case 'filepath'
      outData.data{outDataCols.rp}{ifile,1} = rp_fname;
    end
    
  end %if proc rp
  
end % for ifile=

if(exist('tmp_conn_id','var'))
  mysql(params.ensemble.conn_id,'close');
end


end % ipem_proc_series

%populates empty params with defaults. This needs to be done before
%searching directories for equivalent calcs
function params = getDefaultParams(params)

def.ani = calc_ani('getDefaultParams');

def.pp = params_pp;

def.ci = params_ci;

def.toract = params_toract;

def.pc = params_pc;

def.pc_ci = params_ci;

if exist('params','var') && isfield(params,'rp')
  def.rp = rhythm_profiler_v2('getDefaultParams',params.rp);
else
  def.rp = rhythm_profiler_v2('getDefaultParams');
end

def.glob.save_calc = {'ani','pp','ci','toract','pc','pc_ci'};
def.glob.force_recalc = {};
def.glob.ignore = {};
def.glob.process = {'ani','pp','ci','toract','pc','pc_ci'};
def.glob.outputType = 'data';

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

end %getDefaultParams
