function outData = ipem_proc_series(inData,params)
% Performs a series of IPEM calculations on a list of input wav-files
%
% The calculations to be performed are specified in the params structure
%
% Copyright (c) 2007 The Regents of the University of California
% All Rights Reserved.
%
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
% if 'getDefaultParams' tag was passed, then just return the
% default params



if(ischar(inData) & strcmp(inData,'getDefaultParams'))
  outData = getDefaultParams;
  return
end

%need to catch an error in case ensemble_globals is not in the path
%this provides for public use of ipem_proc_series (without Ensemble)
try
  ensemble_globals;
catch
  stimulus_root = [];
  audio_stim_root = [];
end

outData = ensemble_init_data_struct;
outData.vars = {'ani','pp','ci','toract','rp'};
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


%find out what type of inData was passed
if(isfield(inData,'vars') & isfield(inData,'type')) 
  %deal with ensemble data struct types here
  
  if(ismember('stimulus_id',inData.vars))
    inDataType = 'stimulus_id';
    cols = set_var_col_const(inData.vars);
    stimIDList = inData.data{cols.stimulus_id};
    nfiles = length(stimIDList);
    try 
      params.ensemble.conn_id;
    catch
      params.ensemble.conn_id = 0;
      mysql_make_conn('','',params.ensemble.conn_id);
      tmp_conn_id = 1;
    end
    %get location field rows for the stimuli
    stimLocs = mysql_extract_data('table','stimulus','stimulus_id', ...
				  stimIDList,'extract_vars',{'location'},'conn_id',params.ensemble.conn_id);
    stimLocs = stimLocs{1};
    destStimLocs = check_stim_dirs(stimLocs,'srcroot',stimulus_root,'destroot',stimulus_ipem_analysis_root,'verbose',false);
    
  elseif(strcmp(inData.type,'aud'))
    inDataType = 'aud';
    cols = set_var_col_const(inData.vars);
    nfiles = length(inData);
    
  elseif(ismember('path',inData.vars) & ismember('filename',inData.vars))
    inDataType = 'file_data_struct';
    cols = set_var_col_const(inData.vars);
    nfiles = length(indata.data{cols.filename});
  end
  
elseif(iscell(inData) & ischar(inData{1}))
  %cell array of filepaths as strings
  [p,f,ext] = fileparts(inData{1});
  switch ext
   case {'.wav','.aif','.mp3'}
    inDataType = 'aud_file_string_cell_array';
   case {'.mat'}
    inDataType = 'mat_file_string_cell_array';
  end
  nfiles = length(inData);
    
elseif(ischar(inData) & strcmp(inData,'getDefaultParams'))
  %if inData is a string equal to 'getDefaultParams', return
  %the default params
  outData = getDefaultParams;
  return
elseif(ischar(inData)) 
  %a single file path as a string
  [p,f,ext] = fileparts(inData);
  switch ext
   case {'.wav','.aif','.mp3'}
    inDataType = 'aud_file_string';
   case {'.mat'}
    inDataType = 'mat_file_string';
  end
  nfiles = 1;
  
elseif(isnumeric(inData))
  %if inData is numeric, assuming a list of stim IDs
  inDataType = 'stimulus_id';
  stimIDList = inData;
  nfiles = length(stimIDList);
  %get location field rows for the stimuli
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
    
end
  
%process the list of audio files
for ifile = 1:nfiles

  switch inDataType
   case 'stimulus_id'
    [path,fstub,ext] = fileparts(destStimLocs{ifile});
    filename = [fstub ext];
    paramsSig = params_sig('path',path,'filename',filename);
    sig_st = new_sig(paramsSig);
    descript = sprintf(' Stimulus ID %d ',stimIDList(ifile));
    
   case 'aud'
    sig_st = inData;
    descript = 'aud structure';
    path = inData.data{cols.path};
    filename = inData.data{cols.filename};
    
   case 'file_data_struct'
    path = inData.data{cols.path}{ifile};
    filename = inData.data{cols.filename}{ifile};
    paramsSig = params_sig('path',path,'filename',filename);
    sig_st = new_sig(paramsSig);
    descript = sprintf('%s/%s',path,filename);
    
   case 'aud_file_string_cell_array'
    [path,name,ext] = fileparts(inData{ifile});
    filename = [name ext];
    paramsSig = params_sig('path',path,'filename',filename);
    sig_st = new_sig(paramsSig);
    descript = sprintf('%s/%s',path,filename);
   
   case 'aud_file_string'
    [path,name,ext] = fileparts(inData);
    filename = [name ext];
    paramsSig = params_sig('path',path,'filename',filename);
    sig_st = new_sig(paramsSig);
    descript = '(single audio file)';
    
   case 'mat_file_string'
    %currently, mat input is only supported for rhythm profiles
    %this will need to be expanded
    [path,name,ext] = fileparts(inData);
    filename = [name ext];
    paramsSig = params_sig('path',path,'filename',filename,'Fs',params.insigFs);
    sig_st = new_sig(paramsSig);
    descript = '(single mat file)';
    
  end
  
  %if empty waveform, go to next audio file (or aud struct) in list
  if(exist('sig_st','var'))
    sig_st_cols = set_var_col_const(sig_st.vars);
    if(length(sig_st.data{sig_st_cols.wvf}) == 0)
      continue
    end
  end
  
  fprintf(lfid,'\n\nWorking on file %d/%d: %s\n', ifile, nfiles, descript);
  
  if( ismember('ani',params.glob.process))
    %calculating auditory nerve image
    
    %if path and filename are not empty, we can check for previous
    %analyses, otherwise, this step is skipped and we can also not
    %save any analyses
    if(~isempty(path) & ~isempty(filename))
      anifname = construct_analfname(fullfile(path,filename),'ani');
      aniPath = fileparts(anifname);

      %find existing analysis
      clear matchParams;
      matchParams.varName = 'ani';
      matchParams.paramFind.ani = params.ani;
      matchParams.ignore = {'ani.aniPath'};
      previousCalcFname = check_anal_exist(aniPath,matchParams);
    else
      previousCalcFname = 0;
    end
      
    force_recalc_ani = ismember('ani',params.glob.force_recalc);
    
    if ~force_recalc_ani & previousCalcFname
      fprintf(lfid,'Loading a previously calculated ANI with matching parameters...\n');
      load(fullfile(aniPath,previousCalcFname),'ani')
    else
      % if a previous calculation exists and we are recalculating,
      % use previous calculation filename and overwrite it
      if previousCalcFname
	fprintf(lfid,['FORCE RECALC WAS SET. Overwriting a previous' ...
		      ' ANI calculation\n\n']);
	anifname = fullfile(aniPath,previousCalcFname);
      end
      fprintf(lfid,'Calculating ANI ...\n\n');
      ani = calc_ani(sig_st,params.ani);
     
      if ismember('ani',params.glob.save_calc)
	
	if(exist('anifname','var'))
	  fprintf(lfid,'Saving ANI to %s\n\n',anifname);
	  save(anifname,'ani');
	else
	  fprintf(lfid,'Cannot save ANI. Filename not specified\n\n');
	end
	
      end
    end
 
    outData.data{outDataCols.ani}{ifile} = ani;
    
  end %if proc ani
    
  if( ismember('pp',params.glob.process) )  
    %Calculating periodicity pitch image
    if(~isempty(path) & ~isempty(filename))
      ppfname = construct_analfname(fullfile(path,filename),'pp');
      ppPath = fileparts(ppfname);
    
      %find existing analysis
      clear matchParams;
      matchParams.varName = 'pp';
      matchParams.paramFind.ani = params.ani;
      matchParams.paramFind.pp  = params.pp;
      matchParams.ignore = {'ani.aniPath'}; 
      previousCalcFname = check_anal_exist(ppPath,matchParams);
    else      
      previousCalcFname = 0;
    end
    
    force_recalc_pp = ismember('pp',params.glob.force_recalc);
    
    if ~force_recalc_pp & previousCalcFname
      fprintf(lfid,['Loading a previously calculated Periodicity Pitch Image' ...
	    ' with matching parameters...\n\n']);
      load(fullfile(ppPath,previousCalcFname),'pp')
    else
      
      % if a previous calculation exists and we are recalculating,
      % use previous calculation filename and overwrite it
      if previousCalcFname
	fprintf(lfid,['FORCE RECALC WAS SET. Overwriting a previous' ...
		      ' periodicity pitch calculation\n\n']);
	ppfname = fullfile(ppPath,previousCalcFname);
      end
      
      fprintf(lfid,'Calculating periodicity pitch image ...\n\n');
      pp = calc_pp(ani,params.pp);
      if ismember('pp',params.glob.save_calc)
	
	if(exist('ppfname','var'))
	  fprintf(lfid,'Saving periodicity pitch image to %s\n\n',ppfname);
	  save(ppfname,'pp')
	else
	  fprintf(lfid,'Cannot save pitch image. Filename not specified.\n\n');
	end
	
      end
    end
    
    outData.data{outDataCols.pp}{ifile} = pp;
    
  end %if proc pp
    
    
  if( ismember('ci',params.glob.process) )
    % Calculate a context image
    
    
    if(~isempty(path) & ~isempty(filename))
      cifname = construct_analfname(fullfile(path,filename),'ci');
      ciPath = fileparts(cifname);
    
      %find existing analysis
      clear matchParams;
      matchParams.varName = 'cntxt_img';
      matchParams.paramFind.ani = params.ani;
      matchParams.paramFind.pp  = params.pp;
      matchParams.paramFind.ci  = params.ci;
      matchParams.ignore = {'ani.aniPath'}; 
      previousCalcFname = check_anal_exist(ciPath,matchParams);
    else
      previousCalcFname = 0;
    end
    
    force_recalc_ci = ismember('ci',params.glob.force_recalc);
    
    if ~force_recalc_ci & previousCalcFname
      cifname = previousCalcFname;
      fprintf(lfid,['Loading a previously calculated Context Image' ...
	    ' with matching parameters...\n\n']);
      load(fullfile(ciPath,cifname),'cntxt_img');
       
    else
      
      % if a previous calculation exists and we are recalculating,
      % use previous calculation filename and overwrite it
      if previousCalcFname
	fprintf(lfid,['FORCE RECALC WAS SET. Overwriting a previous' ...
		      ' context image calculation\n\n']);
	cifname = fullfile(ciPath,previousCalcFname);
      end

      fprintf(lfid,'Calculating context images ...\n\n');
      cntxt_img = calc_context(pp,params.ci);
      if ismember('ci',params.glob.save_calc)
	
	if(exist('cifname','var'))
	  fprintf(lfid,'Saving context image to %s\n\n',cifname);
	  save(cifname,'cntxt_img')
	else
	  fprintf(lfid,'Cannot save context image. Filename not specified.\n\n');
	end
	
      end
    end
    
    outData.data{outDataCols.ci}{ifile} = cntxt_img;
    
  end %if proc ci
  
    
  if( ismember('toract',params.glob.process) )  
    % Project the context image to the torus
   
    if(~isempty(path) & ~isempty(filename))
      toract_fname = construct_analfname(fullfile(path,filename),'toract');
      toractPath = fileparts(toract_fname);
    
      %find existing analysis
      clear matchParams;
      matchParams.varName = 'toract';
      matchParams.paramFind.ani = params.ani;
      matchParams.paramFind.pp  = params.pp;
      matchParams.paramFind.ci  = params.ci;
      matchParams.paramFind.toract  = params.toract;
      matchParams.ignore = {'ani.aniPath'};
      previousCalcFname = check_anal_exist(toractPath,matchParams);
    else
      previousCalcFname = 0;
    end
    
    force_recalc_toract = ismember('toract',params.glob.force_recalc);
    
    if (~force_recalc_toract & previousCalcFname)
      toract_fname = previousCalcFname;
      fprintf(lfid,'Loading a previously calculated torus activation with matching parameters...\n');
      load(fullfile(toractPath,toract_fname),'toract')
    else

      % if a previous calculation exists and we are recalculating,
      % use previous calculation filename and overwrite it
      if previousCalcFname
	fprintf(lfid,['FORCE RECALC WAS SET. Overwriting a previous' ...
		      ' torus activation calculation\n\n']);
	toract_fname = fullfile(toractPath,previousCalcFname);
      end

      fprintf(lfid,'Calculating projection to torus ...\n\n');
      toract = calc_toract(cntxt_img, params.toract);
      if ismember('toract',params.glob.save_calc)
	
	if(exist('toract_fname','var'))
	  fprintf(lfid,'Saving torus activation image to %s\n\n',toract_fname);
	  save(toract_fname,'toract','params')
	else
	  fprintf(lfid,['Cannot save torus activation image. Filename' ...
			' not specified.\n\n.']);
	end
      end
    end
    
    outData.data{outDataCols.toract}{ifile} = toract;
    
  end %if proc toract
    
  
  if( ismember('rp',params.glob.process) )  
    % Calculate the rhythm profile
    
    if(~isempty(path) & ~isempty(filename))
      rp_fname = construct_analfname(fullfile(path,filename),'rp');
      rpPath = fileparts(rp_fname);
      
      %find existing analysis
      clear matchParams;
      matchParams.varName = 'rp';
      
      if(strcmp(params.rp.inDataType,'ani'))
	matchParams.paramFind.ani = params.ani;
	matchParams.ignore = {'ani.aniPath'};
      end
      matchParams.paramFind.rp  = params.rp;
      previousCalcFname = check_anal_exist(rpPath,matchParams);
    else
      previousCalcFname = 0;
    end
    
    force_recalc_rp = ismember('rp',params.glob.force_recalc);
    
    if (~force_recalc_rp & previousCalcFname)
      rp_fname = previousCalcFname;
      fprintf(lfid,'Loading a previously calculated rhythm profile with matching parameters...\n');
      load(fullfile(rpPath,rp_fname),'rp')
    else
      
      % if a previous calculation exists and we are recalculating,
      % use previous calculation filename and overwrite it
      if previousCalcFname
	fprintf(lfid,['FORCE RECALC WAS SET. Overwriting a previous' ...
		      ' rhythm profile calculation\n\n']);
	rp_fname = fullfile(rpPath,previousCalcFname);
      end
      
      fprintf(lfid,'Calculating rhythm profile ...\n\n');
      
      %allow for passing ani to the rp model, or passing an aud
      %struct directly to the model, bypassing the ani filters
      switch(params.rp.inDataType)
       case 'ani'
	rp = rhythm_profiler(ani,params.rp);
       case 'aud'
	rp = rhythm_profiler(sig_st,params.rp);
       case 'mat'
	rp = rhythm_profiler(matFilePath,params.rp);
      end
      
      if ismember('rp',params.glob.save_calc)
	
	if(exist('rp_fname','var'))
	  fprintf(lfid,'Saving rhythm profile to %s\n\n',rp_fname);
	  save(rp_fname,'rp','params')
	else
	  fprintf(lfid,'Cannot save rhythm profile. Filename not specified.\n\n');
	end
      end
    end
    
    outData.data{outDataCols.rp}{ifile} = rp;
    
  end %if proc rp
  
end % for ifile=

if(exist('tmp_conn_id','var'))
  mysql(params.ensemble.conn_id,'close');
end


end % ipem_proc_series

%populates empty params with defaults. This needs to be done before
%searching directories for equivalent calcs
function params = getDefaultParams(params)

%this list of defs is not currently complete or correct

def.ani = calc_ani('getDefaultParams');

def.pp = params_pp;

def.ci = params_ci;

def.toract.ci.siglist = {'sig1','sig2'};
def.toract.calc_spher_harm = 1;
def.toract.spher_harm.nharm_theta = 3;
def.toract.spher_harm.nharm_phi = 4;
def.toract.spher_harm.min_rsqr = 0.95;
def.toract.norm = [];
def.toract.som.fname = [];
def.toract.SampleFreq = [];

if(exist('params') & isfield(params,'rp'))
  def.rp = rhythm_profiler('getDefaultParams',params.rp);
else
  def.rp = rhythm_profiler('getDefaultParams');
end

def.glob.save_calc = {'ani','pp','ci','toract'};
def.glob.force_recalc = {};
def.glob.process = {'ani','pp','ci','toract'};

if (nargin)
    
  %only assign defaults to the processes that we are performing
  if(isfield(params,'glob') & isfield(params.glob,'process'))
    defTypes = {params.glob.process{:},'glob'};
  else
    defTypes = fieldnames(def);    
  end
    
  for iType = 1:length(defTypes)

    type = defTypes{iType};  
    fnames = fieldnames(def.(type));

    if ~isfield(params,type)
      params.(type) = [];
    end
    
    %populate newParams with the params values, only if the default was
    %not previously defined or the params is not empty (i.e. don't
    %overwrite a default with an empty value, but do allow empty values
    %to be written in if no default param exists)
    for ifld = 1:length(fnames)
      if(~isfield(params.(type),fnames{ifld}) | isempty(params.(type).(fnames{ifld})) )
	params.(type).(fnames{ifld}) = def.(type).(fnames{ifld});
      end
    end
  end
  
else
  params = def;
end

end %getDefaultParams
