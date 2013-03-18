function [approxRatios,ratioBinIdx,errorVec] = approx_peak_ratios(varargin)
% Approximates the peak frequency ratios relative to the highest peak.
%
% approxRatios = approx_peak_ratios(varargin)
%
% input arguments are tag/value pairs that can be passed in any
% order.

% peakFreqs: A vector of frequency values for each peak.
% maxPeakFreq: The peakFreq that has the highest amplitude.
% params: a parameter structure with the following fields.
%              ratioVals: vector of ratio values to try and fit to the
%                         peak frequency values.
%            ratioLabels: cell array of ratio labels (as strings)
%                         that are index matched to ratioVals. 
%                maxMult: This parameter is only used if
%                         params.relativeTo='minFreq'
%                         The highest multiple used for approximation. Any
%                         peaks that are at higher multiples will not be
%                         assigned a value. 
%         errorTolerance: error percentage within which to regard a
%                         peak as fitting to a certain
%                         ratio. e.g. .05 means that a frequency
%                         value can be fit to a ratio or multiple
%                         within a 5% error. If that frequency
%                         value doesn't fit to any ratio within
%                         that margin of error, it is regarded as
%                         not falling within any of the given ratios.
%             relativeTo: 'maxPeakFreq' or 'minFreq'.
%                         if set to 'maxPeakFreq', then ratios will
%                         be relative to params.maxPeakFreq.
%                         If set to 'minFreq', ratios will be
%                         relative to the smallest frequency value
%                         in params.peakFreqs.
%    
%
% Copyright (c) 2008-2013 The Regents of the University of California
% All Rights Reserved
%
% Author(s):
% Stefan Tomic 9/08
%

for iarg = 1:2:nargin

  switch varargin{iarg}
   case 'peakFreqs'
    peakFreqs = varargin{iarg+1};
   case 'params'
    params = varargin{iarg+1};
   case 'maxPeakFreq'
    maxPeakFreq = varargin{iarg+1};
  end

end

try
  maxMult = params.maxMult;
catch
  maxMult = 16;
end

try
  errorTolerance = params.errorTolerance;
catch
  errorTolerance = .05;
end

try
  params.relativeTo;
catch
  params.relativeTo = 'maxPeakFreq';
end

try
  doReport = params.doReport;
catch
  doReport = 1;
end

switch(params.relativeTo)
  
 case 'maxPeakFreq'
  
  labels = params.ratioLabels;
  ratioVals = params.ratioVals;
  
  peakFreqs = peakFreqs./maxPeakFreq;
  
  for iPeakFreq = 1:length(peakFreqs)
    diffFreqs = abs(peakFreqs(iPeakFreq) - ratioVals);
    %if the peak frequency was equidistant from two ratios
    %just taking the higher ratio. This is a very rare occurrence,
    %and will probably be thrown out anyway because the ratios are beyond
    %the error tolerance
    ratioIdx = find(diffFreqs == min(diffFreqs),1,'last');
    errorPct = abs(ratioVals(ratioIdx) - peakFreqs(iPeakFreq))./ratioVals(ratioIdx);
    errorVec(iPeakFreq) = errorPct;
    if(errorPct <= errorTolerance)
      approxRatios{iPeakFreq,1} = labels{ratioIdx}; 
      ratioBinIdx(iPeakFreq,1) = ratioIdx;
    else
      approxRatios{iPeakFreq,1} = 'N/A';
      ratioBinIdx(iPeakFreq,1) = NaN;
      
      if(doReport)
	fprintf('\napprox_peak_ratios: unable to approximate peak ratios. Value was %.3f\n\n',peakFreqs(iPeakFreq));
      end
    end
    
  end
  
 case 'minFreq'

  %express peakFreqs in relation to min frequency.
  peakFreqs = peakFreqs./min(peakFreqs);
  
  %we're going to try and find the
  %least multiplier that will nicely
  %characterize all the peak frequencies as ratios
  for iMult = 1:maxMult
    
    checkRatios = peakFreqs.*iMult;
    roundCheckRatios = round(checkRatios);
    
    %if any rounded ratios are 0, we know this can't be right, so go
    %to the next loop iteration
    if(any(roundCheckRatios == 0))
      continue;
    else
      errorPct = abs(roundCheckRatios - checkRatios)./roundCheckRatios;
    end
  
    if(all(errorPct <= repmat(errorTolerance,size(checkRatios))))
      approxRatios = roundCheckRatios;
      return
    end
  
  end

  %if we didn't find good ratios in the loop, just return a vector of
  %NaNs for now. Eventually, we may want to isolate those peaks that
  %can't be expressed as simple ratios and still return those that we can
  approxRatios = repmat(NaN,size(peakFreqs));

  
    
end







