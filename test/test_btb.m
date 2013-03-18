% BASIC BTB test script
%
% This script tests jlmt_proc_series and JLMT rhythm profiler analyses on
% .wav and .mp3 files to verify that everything is working. 
%
% Copyright (c) 2007-2013 The Regents of the University of California
% All Rights Reserved
%
% Author:
% Stefan Tomic, Fred Barrett

% SET THE FOLLOWING LINE to the path where you would like to store
% data and output from the test scripts (eg. '/home/username/data' or
% 'C:\Documents and Settings\Username\data')
dest_dir = '';

if isempty(dest_dir)
  error('Please specify "dest_dir" in test/test_btb.m');
end

rp_proc_dir = fullfile(fileparts(which('jlmt_proc_series')),...
    'data','test_jlmt');
input_file_pattern = {'*.wav','*.mp3'};
% input_file_pattern = {'*.mid'};

for k=1:length(input_file_pattern)
  [pth,fn,ext] = fileparts(input_file_pattern{k});

  dirStruct = dir(fullfile(rp_proc_dir,input_file_pattern{k}));
  fnames = {dirStruct.name};
  mapped_audioFileList = check_stim_dirs(fnames,'srcroot',rp_proc_dir,'destroot',dest_dir);

  params.jlmt.ani  = ani_paramGroups;

  switch lower(ext)
    case {'.wav','.mp3'}
      params.jlmt.rp   = rp_paramGroups_v2('param_group',...
				       'reson_filterQSpacing_periodBasedDecay',...
				       'input_type','ani','gain_type','beta_distribution','prev_steps',{'ani'});
  
      params.jlmt.glob.process = {{'ani','rp'}};
      %note that the peak ratio bargraph plot is disabled for audio
      %files since the automated procedure for finding peak ratios was mainly
      %tested and used with MIDI files. In order to use the peak ratio
      %plotter, simply add 'peakRatios' to params.plot.perform.
      params.plot.perform = {'plotBands','plotBandOutputs','plotBandEnergies', 'resonEnergy','meanEnergy','stdEnergy','findPeaks'};
      params.plot.inputSigType = 'wvf'; 

    case {'.mid','.midi'}
      params.jlmt.rp = rp_paramGroups_v2('param_group',...
                           'reson_filterQSpacing_periodBasedDecay',...
                           'input_type','aud','gain_type','beta_distribution');

      params.jlmt.glob.process = {'rp'};

      params.plot.perform = {'resonOutput','resonEnergy','meanEnergy','stdEnergy','findPeaks','peakRatios','showInput'};
      params.plot.inputSigType = 'impulse';

    case {'.csv'}
      params.Fs = 30; %make sure to set this to the appropriate Fs for your csv file
      params.jlmt.rp = rp_paramGroups_v2('param_group',...
                         'reson_filterQSpacing_periodBasedDecay',...
                         'input_type','aud','gain_type','beta_distribution',...
                         'Fs',params.Fs);
      params.jlmt.rp.invertInputSig = 0;
      params.jlmt.glob.process = {'rp'};
      params.plot.perform = {'resonOutput','resonEnergy','meanEnergy','stdEnergy','findPeaks','peakRatios','showInput'};
      params.plot.inputSigType = 'wvf';
      params.plot.inputPlot.normalize = 1; %only normalizes the plot,not the input signal
      params.csvDelim = '\t';
      
  end % switch lower(ext

  %set this to a list of calculations you wish to force recalculation
  %(saved calculations with identical parameters will be overwritten)
  params.jlmt.glob.force_recalc = {};

  %list of calculations that will be saved. Currently set to all
  %processed calculations
  params.jlmt.glob.save_calc = params.jlmt.glob.process;

  % the following parameters govern plotting
  params.n2a_params.timingColumn = 'onset_beats';
  params.n2a_params.beatsPerSec = 2;
  params.n2a_params.ioiThresh = 0.050; %IOI threshold below which to remove double strikes
  params.n2a_params.scaleByVelocity = 1;
  params.separatePlotFiles = 1;
  params.plot.freqLegend = 1;
  params.plot.legendParams.minHz = 1.5;
  params.plot.legendParams.maxHz = 5;
  params.plot.legendParams.numKeys = 7;
  params.plot.colorbar = 1;
  params.plot.savePath = dest_dir;
  params.plot.timeTickRes = 30.0;
  %params.plot.timeLim = [0 15];
  params.plot.numResonTicks = 15;
  params.plot.maxResonFreqPlot = 10;
  params.plot.plotMeans = 1;
  params.plot.color = 1;
  params.plot.separatePages = 0;
  params.plot.savePlot = 1;
  params.plot.peakFinder.thresh = 0.25;
  params.plot.peakFinder.perform = {};
  params.plot.peakFinder.ratioLabels = params.jlmt.rp.resonatorBank.peakFinder.ratioLabels;
  params.plot.peakFinder.ratioVals = params.jlmt.rp.resonatorBank.peakFinder.ratioVals;
  params.plot.meanEnergy.numticks = 4;
  params.plot.plotFname = 'btb.ps';

  plot_dir_rhythm_profile(mapped_audioFileList,params);

end % for k=1:length(input_file_pattern
