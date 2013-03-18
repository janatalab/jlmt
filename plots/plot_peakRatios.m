function outData = plot_peakRatios(ax,inData,params)
% Plots a bargraph that illustrates relative heights of peaks,
% expressed as frequency ratios relative to the highest peak
%
% outData = plot_peakRatios(ax,inData,params)
%
% Description:
%
% INPUTS
%  ax: the axes handle to display the values of peak frequencies
%  inData: the rhythm profile (rp) struct containing the data to
%          analyze
%
%
% Copyright (c) 2008-2013 The Regents of the University of California
% All Rights Reserved
%
% Author(s):
% Stefan Tomic 8/20/2008


try
  ratioLabels = params.peakFinder.ratioLabels;
  numVals = length(ratioLabels);
catch
  try
    numVals = params.numVals;
  catch
    numVals = 16;
  end
end

rp = inData;
rpCols = set_var_col_const(rp.vars);
peakInfo = rp.data{rpCols.peakInfo};
peakInfoCols = set_var_col_const(peakInfo.vars);
ratios = peakInfo.data{peakInfoCols.approxRatio};
heights = peakInfo.data{peakInfoCols.normHeight};

vals = zeros(1,numVals);
if(isnumeric(ratios))
  vals(ratios) = heights;
  bar(ax,vals,0.4,'k');
else
  [dummy,idxs1,idxs2] = intersect(ratioLabels,ratios);
  vals(idxs1) = heights(idxs2);
  bar(ax,vals,0.4,'k');
  set(ax,'Xlim',[0 numVals+1]);
  set(ax,'XTick',[1:numVals]);
  set(ax,'XTickLabel',ratioLabels);
  set(ax,'FontSize',8);
end
  

set(ax,'xlim',[0 numVals + 1]);
set(ax,'ylim',[0 1]);
set(ax,'xtick',1:numVals);

return
