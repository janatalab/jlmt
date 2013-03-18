function plot_rpBandSigs(thisRP,params)
% Plots the resonator bank outputs or periodicity surfaces of each band
%
%
% plot_rpBandSigs(thisRP,params)
%
% Plots band periodicity surfaces or band resonator outputs on a
% 5x1 subplot. params.plotRPBandSigs.plotSig determines whether to
% plot the resonator output or the periodicity surface.
%
% INPUT
%  thisRP: a rhythm profile (rp) data structure to plot
%  params.plotRPBandSigs.plotSig: either 'resonatorOut' to plot
%                                 resonator output or
%                                 'resonatorBandEnergy' for
%                                 periodicity surface
%  params.fig: figure handle to plot to
%  params.bands.climMode: 'scaledAcrossBands' scales clim according
%                         to the min and max value across all bands
%                         'scaledPerBand' scales clim according to
%                         min and max for each band (clims for each
%                         band are independent)
%                         'meanAcrossBands' uses the mean min value
%                         and mean max value to scale the clim of
%                         all bands (makes various peaks clearer in
%                         some instances)
%  params.color: flag (0 or 1) 1: plots in color, 
%                              0: plots in grayscale
%  params.timeLim: vector with two elements that specify the
%                  time range to plot ([min max])
%  params.savePlot: flag (0 or 1) whether to save the plot or just
%                   display to screen
%  params.plotFname: file name to save plot to
%  params.separatePages: (0 or 1) whether or not to plot reson
%                        outputs and periodicity surfaces to separate files
%  params.savePath: file path to save plots (if separatePages is
%                   set, an extra tag will be added to separate the files)
%
% Copyright (c) 2007-2013 The Regents of the University of California
% All Rights Reserved
%
% Author:
% Sept, 2007 Stefan Tomic




get_rpPlotInfo;

if(~isfield(params,'plotRPBandSigs') || ~isfield(params.plotRPBandSigs,'plotSig'))
  params.plotRPBandSigs.plotSig = 'resonatorOut';
end

sigToPlot = thisRP.data{thisRPCols.(params.plotRPBandSigs.plotSig)};

meanChannelCBU = thisRP.data{thisRPCols.meanChannelCBU};
nBands = size(sigToPlot,1);
maxRBEnergy = max(abs(sigToPlot(:)));
  

for iBand = 1:nBands
  
  if(mod(iBand-1,5) == 0)
    [figH,axHList] = rp_figure(params.fig);
  end
  
  ax = get_axhandle(axHList,sprintf('subplot_%d',mod(iBand-1,5)+1));
  axes(ax);
      
  imagesc(timeVector,resIdxVector,squeeze(sigToPlot(iBand,:,:)));

  switch(params.plotRPBandSigs.plotSig)
   case 'resonatorOut'
    plotName = '_rpBandOut_';
    if(strcmp(params.bands.climMode,'scaledAcrossBands'))
      climMin = min(mean(min(sigToPlot,[],3)));
      climMax = max(mean(max(sigToPlot,[],3)));
    elseif(strcmp(params.bands.climMode,'scaledPerBand'))
      climMin = min(min(sigToPlot(iBand,:,:)));
      climMax = max(max(sigToPlot(iBand,:,:)));
    elseif(strcmp(params.bands.climMode,'meanAcrossBands'))
      climMin = mean(mean(min(sigToPlot,[],3)));
      climMax = mean(mean(max(sigToPlot,[],3)));
    end
    if(params.color)
      load new_seismic;
      colormap(new_seismic);
    else
      colormap(gray);
    end
    
    
   case 'resonatorBandEnergy'
    plotName = '_rpBandEnergy_';
    if(strcmp(params.bands.climMode,'scaledAcrossBands'))
      climMin = min(min(mean(sigToPlot,3)));
      climMax = max(max(mean(sigToPlot,3)));
    elseif(strcmp(params.bands.climMode,'scaledPerBand'))
      climMin = min(min(sigToPlot(iBand,:,:)));
      climMax = max(max(sigToPlot(iBand,:,:)));
    elseif(strcmp(params.bands.climMode,'meanAcrossBands'))
      climMin = mean(min(mean(sigToPlot,3),[],2));
      climMax = mean(max(mean(sigToPlot,3),[],2));
    end
    if(params.color)
      colormap(hot);
    else
      colormap(flipud(gray));
    end
    
  end
  
  
  fontSize = .1;
  
  freqTickSpacer = freqMax./7;
  resonFreqTick = floor([freqMin+freqTickSpacer:freqTickSpacer:freqMax-2]);
  
  
  set(ax,'clim',[climMin climMax]);
  set(ax,'Xtick',timeTick);
  set(ax,'Ytick',resonFreqTick);
  set(ax,'YLim',[freqMin freqMax]);
  if(isfield(params,'timeLim'))
    set(ax,'Xlim',params.timeLim);
  else
    set(ax,'XLim',[min(timeVector) max(timeVector)]);
  end

  if(iBand == 1)
    ylabel('Resonator Frequency (Hz)','FontUnits','normalized','FontSize',fontSize);
    xlabel('Time (s)','FontUnits','normalized','FontSize',fontSize);
  end
  
  set(ax,'XtickLabel',sprintf('%d|',timeTick));
  set(ax,'YtickLabel',sprintf('%.1f|',resonFreqs(resonFreqTick)));
  set(ax,'tickdir','out');
  if(mod(iBand,5) == 0 || iBand == nBands)
    if(isfield(params,'savePlot') & params.savePlot)
        [pth,fstub,ext] = fileparts(params.plotFname);
	if(params.separatePages)
	  pfname = [fstub plotName ext];
	  print(params.format,fullfile(params.savePath,pfname));
	else
	  pfname = [fstub ext];
	  print(params.format,'-append',fullfile(params.savePath,pfname));
	end
	disp(sprintf('Writing new plots to %s\n\n',fullfile(params.savePath,pfname)));  
	close(figH);
    end
  
  end
  
end
