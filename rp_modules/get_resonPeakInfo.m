function outData = get_resonPeakInfo(meanEnergy,resonFreqs,params)
% Retrieves peak frequencies and peak widths of an RP structure
%
% outData = get_resonPeakInfo(thisRP,params)
%
% This script calls various other functions for plotting resonator
% bank outputs, periodicity surfaces, average periodicity surfaces
% (APS), mean periodicity profiles (MPP), and input signals.
%
% This function passes the params struct to find_peaks.m
% Please see the discription for that function for details on params.
%
% Copyright (c) 2007-2013 The Regents of the University of California
% All Rights Reserved
%
% Author(s):
% Stefan Tomic 3/07
%


if(size(resonFreqs,2) > size(resonFreqs,1))
  resonFreqs = resonFreqs';
end

peakIdxs = find_peaks(meanEnergy,params);

if(~isempty(peakIdxs))
  phParams.transformation = 'proportion';
  normHeights = peak_heights('signal',meanEnergy,'peakIdxs',peakIdxs,'xVals',resonFreqs,'params',phParams);
  if(isfield(params,'minPeakHeight'))
    goodPeakMask = normHeights >= params.minPeakHeight;
    peakIdxs = peakIdxs(goodPeakMask,1);
    normHeights = normHeights(goodPeakMask);
  end
end

if(~isempty(peakIdxs))
  peakFreqs = resonFreqs(peakIdxs);
  maxPeakEnergyFreq = resonFreqs(find(meanEnergy == max(meanEnergy(peakIdxs)),1,'first'));
  normPeakFreqs = peakFreqs ./ maxPeakEnergyFreq;
  widths = get_rp_peakWidths('rpEnergy',meanEnergy,'peakFreqs',peakFreqs,'resonFreqs',resonFreqs);
  [approxRatios,ratioBinIdx,errorVec] = approx_peak_ratios('peakFreqs',peakFreqs,'maxPeakFreq',maxPeakEnergyFreq,'params',params);

  %get peak areas
  peakAreas = peak_areas('signal',meanEnergy,'xVals',resonFreqs,'peakIdxs',peakIdxs);
  %retrieve actual peak heights (not normalized)
  phParams.transformation = 'none';
  peakHeights = peak_heights('signal',meanEnergy,'peakIdxs', ...
			     peakIdxs,'xVals',resonFreqs,'params',phParams);
  
  [dummyVar,peakIdxs] = intersect(resonFreqs,peakFreqs);

  %get base peak widths
  pwParams.calcAt = 'base_high';
  pwParams.interpolate.method = 'spline';
  [baseFreqs] = peak_widths('signal',meanEnergy,'peakIdxs',peakIdxs,'xVals',resonFreqs,'params',pwParams);

  baseWidths = baseFreqs(2:2:end) - baseFreqs(1:2:end-1);

  
else
  peakFreqs = [];
  normPeakFreqs = [];
  widths = [];
  baseWidths = [];
  peakHeights = [];
  normHeights = [];
  peakAreas = [];
  approxRatios = [];
  errorVec = [];
  ratioBinIdx = [];
end

outData.peakFreq = peakFreqs;
outData.peakIdxs = peakIdxs;
outData.normPeakFreq = normPeakFreqs;
outData.width = widths;
outData.baseWidths = baseWidths;
outData.peakHeight = peakHeights;
outData.normHeight = normHeights;
outData.peakAreas = peakAreas;
outData.approxRatio = approxRatios;
outData.ratioBinIdx = ratioBinIdx;
outData.errorVec = errorVec;
return
  

  
