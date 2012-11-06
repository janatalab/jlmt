% Obtains generic rhythm profile (rp) struct info for plotting
%
% script that is called at the beginning of many of the
% rhythm_profiler functions. Obtains generic rhythm profile plot info
% useful for multiple plots (rp output, rp energy, mean energy, etc)
%
%
%
% Copyright (c) 2007-2012 The Regents of the University of California
% All Rights Reserved
%
% Author(s):
% Stefan Tomic 09/2007

thisRPCols = set_var_col_const(thisRP.vars);
Fs = thisRP.data{thisRPCols.Fs};
resonFreqs = thisRP.data{thisRPCols.resonatorFreqs};

try
  timeTickRes = params.timeTickRes;
catch
  timeTickRes = 2.0;
end

%make sure resonFreqs is a column vector to construct the
%proper label string matrix
if(size(resonFreqs,2) > 1)
  resonFreqs = resonFreqs';
end
  
nResons = length(resonFreqs);
freqMin = 1;
if(isfield(params,'maxResonFreqPlot'))
  freqMax = find(resonFreqs < params.maxResonFreqPlot,1,'last');
else
  freqMax = length(resonFreqs);
end

resIdxVector = [1:nResons];
resonFreqTick = [freqMin:floor(freqMax./params.numResonTicks):freqMax];
resonFreqLabel = num2str(resonFreqs(resonFreqTick),'%.2f');
resEnergy = thisRP.data{thisRPCols.resonatorEnergy};

%if timelim is not set, then we can use the stored stdEnergy and
%meanEnergy. Otherwise, we need to recalculate them to correspond
%to the time limits being displayed. Same goes for our peak picker.


if(~isfield(params,'timeLim'))
  stdEnergy = thisRP.data{thisRPCols.stdResonatorEnergy};
  meanEnergy = thisRP.data{thisRPCols.meanResonatorEnergy};
  
  try
    peakInfo = thisRP.data{thisRPCols.peakInfo};
  catch
    peakInfo = [];
  end
  
  if isempty(peakInfo)
    peakInfo = get_resonPeakInfo(thisRP.data{thisRPCols.meanResonatorEnergy},...
      resonFreqs,params.peakFinder);
    peakInfo = ensemble_tree2datastruct(peakInfo);
    %temporary solution until ensemble_tree2datastruct is fixed
    if(iscell(peakInfo.data{1}))
      for iTmp = 1:length(peakInfo.data)
        peakInfo.data{iTmp} = peakInfo.data{iTmp}{1};
      end
    end
  end
  
  peakInfoCols = set_var_col_const(peakInfo.vars);
  peakIdxs = peakInfo.data{peakInfoCols.peakIdxs};

else
  maxTime = size(thisRP.data{thisRPCols.resonatorOut},3)./Fs - 1./Fs;
  timeLim = min(params.timeLim,maxTime);
  tempTimeV = [timeLim(1)*Fs:timeLim(2)*Fs]+1;%add one to allow
                                               %for starting at
                                               %time t=0
					       

					       
  stdEnergy = std(thisRP.data{thisRPCols.resonatorEnergy}(:,tempTimeV),[],2);
  meanEnergy =  mean(thisRP.data{thisRPCols.resonatorEnergy}(:, ...
						  tempTimeV),2);
  
  if(~isfield(params,'peakFinder'))
    params.peakFinder = struct();
    params.peakFinder.minPeakHeight = 0.1;
  end
  
  peakInfo = get_resonPeakInfo(thisRP.data{thisRPCols.meanResonatorEnergy},...
		    resonFreqs,params.peakFinder);
  peakInfo = ensemble_tree2datastruct(peakInfo);
  %temporary solution until ensemble_tree2datastruct is fixed
  if(iscell(peakInfo.data{1}))
    for iTmp = 1:length(peakInfo.data)
      peakInfo.data{iTmp} = peakInfo.data{iTmp}{1};
    end
  end
  
  thisRP.data{thisRPCols.peakInfo} = peakInfo;
  
end

minMeanEnergy = min(meanEnergy);
maxMeanEnergy = max(meanEnergy);
energyMaxVal = max(resEnergy(:));

nTimeSamps = size(resEnergy,2);

%some rp structs only have mpp info (no time info)
%this occurs generally in cases where mpps were calculated from
%multiple mpps (e.g. average across multiple mpps)
if(~isempty(Fs) && (nTimeSamps ~= 0))
  timeVector = [0:nTimeSamps-1]./Fs;
else
  timeVector = [];
end

if(~isempty(timeVector))
  timeTick = [0:timeTickRes:timeVector(end)];
else
  timeTick = 0;
end



