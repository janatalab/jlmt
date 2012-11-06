function plot_rpMeanResonEnergy(ax,inData,params)
% Plots a mean periodicity profile (MPP)
% 
% plot_rpMeanResonEnergy(ax,thisRP,params)
%
% plots a mean periodicity profile (MPP) on a specified axes given
% a rp struct. Standard deviation bars can optionally be plotted.
%
% INPUT
% ax: the axes handle to plot to
% thisRP: rhythm profile (rp) data struct
% params.meanEnergy.numTicks: number of y ticks
% params.meanEnergy.yTickLabelFmt: C String style format for y tick labels
% params.meanEnergy.yLim: ylim specification for MPP
% params.perform: if this contains a cell array with 'stdEnergy',
%                 then standard deviation bars will be plotted
% params.color: 0 or 1. Whether or not to plot in color
%
% Copyright (c) 2007-2012 The Regents of the University of California
% All Rights Reserved
%
% Author(s):
% Stefan Tomic 4/07

try
  plotStyle = params.plotStyle;
catch
  plotStyle = 'line';
end

try
  showResonTickLabels = params.showResonTickLabels;
catch
  showResonTickLabels = 0;
end

%multiple RPs are passed in as an array of cell structs
if(iscell(inData))
  nData = length(inData);
else
  nData = 1;
end
  
try(params.plotStd)
  params.plotStd;
catch
  params.plotStd = ones(nData,1);
end

setYLim = [NaN NaN];

for idx = 1:nData

  if(iscell(inData))
    thisRP = inData{idx};
  else
    thisRP = inData;
  end
  
  if(isfield(params,'lineSpec'))
    lspec = params.lineSpec{idx};
  elseif(strcmp(plotStyle,'fill'))
    lspec = [0 0 1];
  elseif(strcmp(plotStyle,'line'))
    lspec = 'k-';
  end
 
  get_rpPlotInfo;
  axes(ax);
  
  if(~isfield(params,'meanEnergy') || (~isfield(params.meanEnergy,'numticks')))
    params.meanEnergy.numticks = 4;
  end
  
  if(~isfield(params,'meanEnergy') || (~isfield(params.meanEnergy,'yTickLabelFmt')))
    params.meanEnergy.yTickLabelFmt = '%.2f';
  end
  
  if(params.plotStd(idx))
    
    
    if(params.color)
      colorSpec = [0 0.8 0];
    else
      colorSpec = [0.5 0.5 0.5];
    end
    
     switch(plotStyle)
     case 'line'      
      errorbar(ax,meanEnergy,stdEnergy,lspec,'color',colorSpec);
      lineH(idx) = plot(ax,meanEnergy,lspec,'lineWidth',.5);
     case 'fill'
      axes(ax);
      lineH(idx) = jbfill([1:length(meanEnergy)],meanEnergy'+stdEnergy',...
	     meanEnergy'-stdEnergy',lspec,'k',1,0.5);
     hold on
     plot(ax,meanEnergy,'color',lspec,'lineWidth', 1);
    end
  else
    lineH(idx) = plot(ax,meanEnergy,'color',lspec,'lineWidth',1);
  end

  hold on;
  
  if(ismember('stdEnergy',params.perform))
    lowerBound = min(meanEnergy - stdEnergy);
    if sign(lowerBound) < 1
      setYLim(1) = min(setYLim(1),1.2*lowerBound);      
    else
      setYLim(1) = min(setYLim(1),0.8.*lowerBound);
    end
    setYLim(2) = max(setYLim(2),1.2.*max(meanEnergy + stdEnergy));
  else
    setYLim(1) = min(setYLim(1),0.8.*min(meanEnergy));
    setYLim(2) = max(setYLim(2),1.2.*max(meanEnergy));
  end
  
  if(ismember('plotPeaks',params.perform))
    plot(peakIdxs,meanEnergy(peakIdxs),'r*');
    nPeaks = length(peakIdxs);
    for iPeak = 1:nPeaks
      peakLabel = sprintf('%.2f',resonFreqs(peakIdxs(iPeak)));
      text(peakIdxs(iPeak),meanEnergy(peakIdxs(iPeak)),peakLabel,...
	   'horizontalalignment','center','verticalalignment','bottom');
    end
  end
  
  
  
end %for idx
  
if(isfield(params,'meanEnergy') && isfield(params.meanEnergy,'ylim'))
  ylim = params.meanEnergy.ylim;
else
  ylim = setYLim;
end


set(ax,'ylim',ylim);
set(ax,'xlim',[freqMin freqMax]);
  
set(ax,'Xtick',resonFreqTick);

if(showResonTickLabels)
  set(ax,'Xticklabel',resonFreqLabel);
else
  set(ax,'XtickLabel',[]);
end

ytick = [ylim(1):(ylim(2)-ylim(1))/(params.meanEnergy.numticks-1):ylim(2)]';

set(ax,'ytick',ytick);
set(ax,'yticklabel',num2str(ytick,params.meanEnergy.yTickLabelFmt));
    
if(isfield(params,'legendLabels'))
  legend(lineH,strrep(params.legendLabels,'_','\_'),'Location','NorthWest');
end

return
