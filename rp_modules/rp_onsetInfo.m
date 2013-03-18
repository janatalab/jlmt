function outData = rp_onsetInfo(varargin)
%
% Calculates and returns various onset information (times, heights) of onsets.
%
% outData = rp_onsetInfo(varargin)
%
% 
% INPUTS:
% tag, value pairs:
%      'rp': the rhythm profile to process
%
%      relevant params:
%            params.peakfinder.thresh: threshold for the peak
%                                      finder for finding onsets
%            params.Fs: the sampling rate of the reson filters
%            params.predictCurrentOutput: params struct for predictCurrentOutput.m
%            
% Copyright (c) 2009-2013 The Regents of the University of California
% All Rights Reserved.
%
% Author: Stefan Tomic 03/19/2009


for iArg = 1:2:nargin
  switch(varargin{iArg})
   case('rp')
    rp = varargin{iArg+1};
   case('inputSig')
    inputSig = varargin{iArg+1};
  end
end

%if this is an ensemble data struct, convert to tree
if(isfield(rp,'vars') && isfield(rp,'data'))
  rp = ensemble_datastruct2tree(rp);
end

durSamps = size(rp.resonatorOut,3);

params = rp.params.rp;

%onsetDetect.thresh contains the threshold value for detecting
%onset times
try
  params.onsetInfo.thresh;
catch
  warning('onset detector threshold not set. Setting to 0.05');
  params.onsetInfo.thresh = 0.05;
end

Fs = params.Fs;

nResons = length(rp.resonatorFreqs);

nBands = size(rp.resonatorOut,1);
nTimeSamps = size(rp.resonatorOut,3);

if(~isempty(params.timeLimSecs))
  timeStartSamps = round(params.timeLimSecs(1)*params.Fs) + 1;
  timeStopSamps = min(round(params.timeLimSecs(2)*params.Fs),nTimeSamps);
else
  timeStartSamps = 1;
  timeStopSamps = nTimeSamps;
end

nFrames = size(rp.resonatorEnergy,2);

for iBand = 1:nBands
  warning off
  onsetTimesSamps = find_peaks(inputSig(iBand,timeStartSamps:timeStopSamps),params.onsetInfo);
    
  %Some ANI signals have exhibited large spikes at the beginning of the signal
  %where there probably shouldn't have been one. So, for now, params should be
  %set to throw out the first detected onset for audio files until this problem is solved.
  if(isfield(params.onsetInfo,'throwOutFirstOnset') && params.onsetInfo.throwOutFirstOnset)
    onsetTimesSamps(1) = [];
  end
  
  onsetHeights = peak_heights('signal',inputSig(iBand,:),'peakIdxs',onsetTimesSamps);
  pwParams.calcAt = 'base_high';
  
  try
    pwParams.gtThreshold = params.onsetInfo.peakBaseThreshold;
  catch
    pwParams.gtThreshold = 0.001;
  end
  
  
  onsetBorders =  peak_widths('signal',inputSig(iBand,:),'peakIdxs',onsetTimesSamps,'params',pwParams);
  
  if(size(onsetHeights,1) > size(onsetHeights,2))
    onsetHeights = onsetHeights';
  end
  
  onsetTimesSampsByBand{iBand} = onsetTimesSamps;
  onsetHeightsByBand{iBand} = onsetHeights;
  onsetBordersByBand{iBand} = onsetBorders;
  


end %for iBand

outData = ensemble_init_data_struct;
outData.vars = {'onsetTimesSampsByBand','onsetHeightByBand','onsetBordersByBand'};
outDataCols = set_var_col_const(outData.vars);
outData.data{outDataCols.onsetTimesSampsByBand} = onsetTimesSampsByBand;
outData.data{outDataCols.onsetHeightByBand} = onsetHeightsByBand;
outData.data{outDataCols.onsetBordersByBand} = onsetBordersByBand;