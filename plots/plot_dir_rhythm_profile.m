function outData = plot_dir_rhythm_profile(inData,params)
% Runs the BTB algorithm and plots the resulting periodicity surfaces
%
% outData = plot_dir_rhythm_profile(inData,params)
%
% Cycles through a list of wav or midi files and plots resonator
% outputs, average periodicity surfaces, band periodicity surfaces,
% and mean periodicity profiles. Any subset of these will be plotted
% depending on how the params are set. This function calls
% jlmt_proc_series to either find an existing saved calculation or
% to calculate new ANI and rp, so this function may be used to
% perform the calculations as well as plotting.
%
% INPUT/PARAMS 
%  inData: cell array of files to cycle through. Expects the full
%          paths of the files
%  params.plot.savePath:     directory where the plot(s) should be saved
%  params.plot.plotFname:    filename where the plot(s) should be saved
%  params.separatePlotFiles: flag (0 or 1) that determines whether
%                            or not an extra filename should be 
%                            tagged to the plot name (for
%                            separating plots for each input wav 
%                            or midi file)
%  params.jlmt:              parameter structure for
%                            jlmt_proc_series. The data structures
%                            for ani and rp calculations will
%                            either be searched for and found, or
%                            will be calculated. params.jlmt.glob,
%                            params.jlmt.ani, params.jlmt.rp each
%                            should be populated with the desired
%                            jlmt parameters (see jlmt_proc_series)
%
% Copyright (c) 2007-2013 The Regents of the University of California
% All Rights Reserved
%
% Author(s):
% Stefan Tomic 8/07

if(~isfield(params,'separatePlotFiles'))
  params.separatePlotFiles = 0;
end

plotFile = fullfile(params.plot.savePath,params.plot.plotFname);
[p,plotFStub,plotext] = fileparts(plotFile);
if(params.plot.savePlot  & exist(plotFile))
    
  confirmOverwrite = '';
  while(~ismember(confirmOverwrite,{'y','n'}))
    confirmOverwrite = input(sprintf('File %s exists, overwrite? (y/n): ',plotFile),'s');
  end
  
  if(strcmp(confirmOverwrite,'y'))
    delete(plotFile);
  else
    disp('OK, not saving the plot. It will be displayed on the screen.');
    params.plot.savePlot = 0;
  end
  
end
  
%support processing mat filenames in a cell array of strings, or
%stored in an ensemble data struct
if(iscell(inData) && ischar(inData{1}))
  fList = inData;
  nFiles = length(fList);
elseif(isstruct(inData) && isfield(inData,'type'))
  inDataCols = set_var_col_const(inData.vars);
  paths = inData.data{inDataCols.path};
  fnames = inData.data{inDataCols.filename};
  for iFile = 1:length(fnames)
    fList{iFile} = fullfile(paths{iFile},fnames{iFile});
  end
  nFiles = length(fList);
end


for iFile = 1:nFiles;
 
  thisLocation = fList{iFile};
  [path,fstub,ext] = fileparts(thisLocation);
  
  rpPath = fullfile(path,fstub,'rp');
  
  params.plot.fig.title = sprintf('%s',[fstub ext]);
  params.plot.inDataPath = fullfile(path,fstub);

  clear inputPlotData;
  clear audStruct;

  switch lower(ext)
   case {'.wav','.mat','.csv','.mp3'}

    if(isfield(params,'srcMatVarName'))
      pa =  params_aud('path',path,'filename',[fstub ext],'var_name',params.srcMatVarName);
    elseif(isfield(params,'csvDelim'))
      pa =  params_aud('path',path,'filename',[fstub ext],'delim',params.csvDelim,'Fs',params.Fs);
    else
      pa = params_aud('path',path,'filename',[fstub ext]);
    end
    
    audStruct = new_aud(pa);
    inputPlotData = {audStruct};
   case {'.mid','*.midi'}
    
    params.n2a_params.Fs = params.jlmt.rp.Fs;
    nmat_cols = nmat_columns;
    tr = javamidi2tr(thisLocation);
    %unique midi notes are used to differentiate left and right hands in
    %our current tapping trials
    unique_midi_notes = unique(tr.nmat(:,nmat_cols.midi_pitch));
    trialSets = [];
    %separate left and right hands for plotting the input
    for inote = 1:length(unique_midi_notes)
      params.n2a_params.extractNotes = unique_midi_notes(inote);
      inputPlotData{inote} = nmat2aud(tr.nmat,params.n2a_params);
      trialSets = [trialSets;inote mod(inote-1,2)+1];
      params.plot.trialSets = trialSets;
    end
    %create an aud struct for inputting to a resonator bank (both hands)
    params.n2a_params.extractNotes = unique_midi_notes;
    audStruct = nmat2aud(tr.nmat,params.n2a_params);
    audStructCols = set_var_col_const(audStruct.vars);
    audStruct.data{audStructCols.path} = path;
    audStruct.data{audStructCols.filename} = [fstub ext];
  end

  %params.jlmt.glob.force_recalc = {};
  %params.jlmt.glob.save_calc = params.jlmt.glob.process;
  disp('Running jlmt_proc_series...');
  jlmtOut = jlmt_proc_series(audStruct,params.jlmt);
  disp('Finished jlmt_proc_series...');
  jlmtOutCols = set_var_col_const(jlmtOut.vars);
  
  rp = jlmtOut.data{jlmtOutCols.rp}{1};
  
  plotData{1} = rp;
  plotData(2:length(inputPlotData)+1) = inputPlotData;
  
  if(params.separatePlotFiles)
    params.plot.plotFname = [plotFStub '_' fstub plotext];
  else
    params.plot.appendToPlot = 1;
  end
  
  plot_rhythmProfile(plotData,params.plot);
  
end

outData = plotFile;

return

