function outData = plot_findPeaks(ax,inData,params)
% Finds and displays the peaks in a mean periodicity profile 
%
% outData = plot_findPeaks(ax,inData,params)
%
% Description:
% Finds the peaks in a mean periodicity profile and displays the
% peak frequencies on the specified axes
%
% INPUTS
%  ax: the axes handle to display the values of peak frequencies
%  inData: the rhythm profile (rp) struct containing the data to
%          analyze
%  params.peakFinder.perform: cell array. If it contains "lowPass",
%                             then a low pass filter will be
%                             applied to the signal before peak detection
%  params.peakFinder.lpFreq: low pass cutoff frequency if
%                            performing low pass filtering for peak detection
%  params.peakFinder.thresh: float between 0 and 1. fraction of max
%                            value of the signal to regard as a
%                            threshold below which peak
%                            detection will not be performed
%
%
% Copyright (c) 2007-2012 The Regents of the University of California
% All Rights Reserved
%
% Author(s):
% Stefan Tomic 3/07

try
  fontSize = params.fontSize;
catch
  fontSize = 10;
end

thisRP = inData;
thisRPCols = set_var_col_const(thisRP.vars);
get_rpPlotInfo;

peakInfo = thisRP.data{thisRPCols.peakInfo};
peakInfoCols = set_var_col_const(peakInfo.vars);
peakFreqs = peakInfo.data{peakInfoCols.peakFreq};
normPeakFreqs = peakInfo.data{peakInfoCols.normPeakFreq};
widths = peakInfo.data{peakInfoCols.width};
heights = peakInfo.data{peakInfoCols.normHeight};
approxRatios = peakInfo.data{peakInfoCols.approxRatio};

if(isfield(params,'fpTitle'))
  peakFindText = sprintf([params.fpTitle '\n\n']);
else
  peakFindText = '';
end

peakFindHeader = sprintf('Peak   Norm   Width   normHeight   Ratio\n');

peakFindText = [peakFindText peakFindHeader];

for idx = 1:length(peakFreqs)
  if(iscell(approxRatios))  
    if(isempty(approxRatios))
      r = 'N/A';
    else
      r = approxRatios{idx};
    end
    
    if(iscell(r))
      r = r{1};
    end
    
    newRow = sprintf('%0.2f%8.2f%9.2f%8.2f%14s\n', ...
		     peakFreqs(idx),normPeakFreqs(idx),widths(idx),heights(idx),r);
  else
    newRow = sprintf('%0.2f%8.2f%9.2f%8.2f%14d\n', ...
		     peakFreqs(idx),normPeakFreqs(idx),widths(idx),heights(idx),approxRatios(idx));
  end
  peakFindText = [peakFindText newRow];
end
axes(ax);
set(ax,'visible','off');
t = text(0,1,peakFindText,'VerticalAlignment','top','FontSize',fontSize);
  

  
