function plot_rpResonEnergy(ax,thisRP,params)
% Plots an average periodicity surface (APS)
%
% plot_rpResonEnergy(ax,thisRP,params)
%
%
% INPUTS
% ax:                    axes handle to the axes for plotting the APS
% thisRP:                the rhythm profile data struct to plot the APS from
% params.colorbar:       0 or 1. Whether or not to plot a "colorbar"
%                        which indicates the values that correspond 
%                        to color or gray gradations
% params.color:          0 or 1.  Whether or not to plot the APS in color.
% params.resEnergy.clim: the clim values to use for the plot. If
%                        not specified, the clim minimum is set 
%                        to the minimum of the mean periodicity
%                        profile (MPP) and the maximum clim is set to the
%                        maximum value in the APS.
% params.timeLim:        The min and max time (in seconds) to plot.
% params.xminortick:     0 or 1. Whether or not to plot minor ticks.   
%
%
% Copyright (c) 2007-2012 The Regents of the University of California
% All Rights Reserved
%
% Author(s):
% Stefan Tomic 8/07

if(~isfield(params,'colorbar'))
  params.colorbar = 0;
end

get_rpPlotInfo;
resEnergy = thisRP.data{thisRPCols.resonatorEnergy};

axes(ax);
tag = get(ax,'tag');
imagesc(timeVector,resIdxVector,resEnergy);
if(params.color)
  colormap(hot);
else
  colormap(flipud(gray));
end

set(ax,'tag',tag);  
if(isfield(params,'resEnergy') && isfield(params.resEnergy,'clim'))
  setClim = params.resEnergy.clim;
else
  setClim = [minMeanEnergy energyMaxVal];
end

set(ax,'FontSize',12);
set(ax,'clim',setClim);
set(ax,'Xtick',timeTick);
set(ax,'Ytick',resonFreqTick);
set(ax,'YtickLabel',resonFreqLabel);
set(ax,'YLim',[freqMin freqMax]);
set(ax,'tickdir','out');
if(isfield(params,'timeLim'))
  set(ax,'Xlim',params.timeLim);
end

if(isfield(params,'xminortick'))
  set(ax,'xminortick',params.xminortick);
end

xlabel('time (s)');
ylabel('Resonator Freq (Hz)');
if(params.colorbar)
  pos = get(ax,'position');
  cb = colorbar('position',[pos(1)+pos(3)+0.17 pos(2) 0.0347 0.3405]);
end

return
