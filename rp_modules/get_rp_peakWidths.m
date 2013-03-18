function [outData] = get_rp_peakWidths(varargin)
% Calculates the widths of peaks in a signal
%
% outData = get_rp_peakWidths('rpEnergy',periodicity_profile,'peakFreqs',peakFreqs,'resonFreqs',resonFreqs)
%
% uses peak_widths function to return values that are meaningful
% for periodicity profiles. The returned widths are at a point
% somewhat equivalent to a full-width half max (FWHM). The half max
% point is calculated from the higher of the left and right bases of
% the peak.
%
% INPUT
%  tag/value pairs. Expected tag/values are:
%   'rpEnergy'  : corresponds to value of periodicity profile
%   'peakFreqs' : corresponds to vector of peak frequencies
%                 (returned from find_peaks.m)
%   'resonFreqs': corresponds to vector of all reson filter
%                 frequencies
%
% OUTPUT
%  outData: vector of peak widths (in Hz) for the peak values
%  corresponding to peakFreqs.
%
% Copyright (c) 2007-2013 The Regents of the University of California
% All Rights Reserved
%
% Author:
% Stefan Tomic 8/07

for iarg = 1:2:nargin
  
  switch(varargin{iarg})
   case 'rpEnergy'
    rpEnergy = varargin{iarg+1};
   case 'peakFreqs'
    peakFreqs = varargin{iarg+1};
   case 'resonFreqs'
    resonFreqs = varargin{iarg+1};
  end
 
end


[dummyVar,peakIdxs] = intersect(resonFreqs,peakFreqs);

pwParams.calcAt = 'half_height';
pwParams.interpolate.perform = 'spline';
pwParams.interpolate.strength = 20;
[hhFreqs] = peak_widths('signal',rpEnergy,'peakIdxs',peakIdxs,'xVals',resonFreqs,'params',pwParams);

hhFreqWidths = hhFreqs(2:2:end) - hhFreqs(1:2:end-1);

outData = hhFreqWidths;
