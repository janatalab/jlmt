function plot_rpResonOutput(ax,thisRP,params)
% Produces a surface plot of a resonator bank output for a single band implementation.
%
% Produces a surface plot of resonator outputs. This function is
% only used for cases in which only one band is employed (e.g. MIDI
% inputs). For plotting resonator outputs of individual bands in
% cases that employ multiple bands, use plot_rpBandSigs.m instead.
%
% 
% INPUT
% ax:                   axes handle of axes to plot to
% thisRP:               rhythm profile (rp) data struct that contains the
%                       resonator output data.
% params.colorbar:      0 or 1. Whether or not to plot a color bar that
%                       indicates values for gradations of color or gray in the plot
% params.color:         0 or 1. Whether or not to plot in color. 
%                       Grayscale is used if 0.
% params.timeTickRes:   Time tick resolution in seconds. Determines
%                       the spacing of time ticks.
% params.resonOut.clim: Vector with 2 elements. clim values to use
%                       for plot.
% params.timeLim:       Vector with 2 elements. Time limit(xlim)
%                       values in seconds to use for plot.
% 
% OUTPUT
% No output variables. Plots to an axes.
%
% Copyright (c) 2007-2012 The Regents of the University of California
% All Rights Reserved
%
% Author:
% Stefan Tomic 8/07



if(~isfield(params,'colorbar'))
  params.colorbar = 0;
end

get_rpPlotInfo;

rOut = squeeze(thisRP.data{thisRPCols.resonatorOut});

nTimeSamps = size(rOut,2);
timeVector = [0:nTimeSamps-1]./Fs;
timeTick = [0:params.timeTickRes:timeVector(end)];
tag = get(ax,'tag');
axes(ax);
imagesc(timeVector,resIdxVector,rOut);

if(params.color)
  load new_seismic;
  colormap(new_seismic);
else
  colormap(gray);
end
set(ax,'tag',tag);
maxRoutVal = max(max(abs(rOut))); 
if(isfield(params,'resonOut') && isfield(params.resonOut,'clim'))
  setClim = params.resonOut.clim;
else
  setClim = [-maxRoutVal maxRoutVal];
end
if(isfield(params,'timeLim'))
  set(ax,'Xlim',params.timeLim);
else
  set(ax,'xlim',[min(timeVector) max(timeVector)]);
end

set(ax,'FontSize',12);
set(ax,'clim',setClim);
set(ax,'Xtick',timeTick);
set(ax,'Ytick',resonFreqTick);
set(ax,'YtickLabel',resonFreqLabel);
set(ax,'YLim',[freqMin freqMax]);
xlabel('time (s)');
ylabel('Resonator Freq (Hz)');
if(params.colorbar)
  pos = get(ax,'position');
  cb = colorbar('position',[pos(1)+pos(3)+0.02 pos(2) 0.0347 0.3405]);
end


return
