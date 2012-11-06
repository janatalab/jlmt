function [fig_h,ax_h] = rp_figure(params)
% Produces a figure with axes locations commonly used with BTB
%
% [fig_h,ax_h] = rp_figure(params)
% 
% Creates a figure for plotting resonator outputs, periodicity
% surfaces, and mean periodicity profiles. Axes are added to the
% figure at locations specified by params.
%
% INPUT
%  params.figh:    figure handle to add axes to. Defaults to next
%                  available figure number if not specified
%  params.figSize: if 'letter', position and paperposition WxH
%                  properties will be 8.5"x11"
%                  if 'ipod': position WxH is set to 640x480
%                  pixels.
%  params.title:   string specifying the figure title
%  params.axLocs: cell array of strings specifying the axes you want.
%  params.rowLabels: cell array of strings for subplot row labels
%  params.colLabels: cell array of strings for subplot column labels
%
% Copyright (c) 2007-2012 The Regents of the University of California
% All Rights Reserved
%
% Author(s):
% Stefan Tomic 6/07


if(isfield(params,'figh'))
  fig_h = figure(figh);
else
  fig_h = figure;
end

ax_h = [];

if(~isfield(params,'rectPlotWidth'))
  rectPlotWidth = 0.65;
else
  rectPlotWidth = params.rectPlotWidth;
end

if(~isfield(params,'rectPlotHeight'))
  rectPlotHeight = 0.34;
else
  rectPlotHeight = params.rectPlotHeight;
end

squarePlotDim = 0.5;
leftMargin = 0.1;
rightMargin = 0.1;
bottomMargin = 0.05;
spacer = 0.03;

topMargin = 0.08;
middle = bottomMargin./2 + (1 - topMargin)./2;

if(~isfield(params,'figSize'))
  params.figSize = 'letter';
end


switch(params.figSize)
 case 'letter'
  set(fig_h,'units','inches');
  set(fig_h,'position',[1 1 8 10.5]); 
  set(fig_h,'paperposition',[0.25 0.25 8 10.5]); 
 case 'ipod'
  set(fig_h,'units','pixels');
  set(fig_h,'position',[100 100 640 480]);
end

%assign a title to the figure
if(isfield(params,'title'))
  titleStr = strrep(params.title,'_','\_');
  axes('position',[leftMargin 0.95 rectPlotWidth 0.02],'visible','off');
  txt_h = text(0.5,0.5,titleStr);
  set(txt_h,'HorizontalAlignment','center','FontWeight','bold', ...
	    'FontSize',14);
end

if(~isfield(params,'aspect'))
  params.aspect = 'rectangle';
end
  
if(ismember('rectPlot_top',params.axLocs))
  ax_h(end+1) = axes('position',[leftMargin middle rectPlotWidth rectPlotHeight]);
  set(ax_h(end),'tag','rectPlot_top');
end

if(ismember('smallRectPlot_top',params.axLocs))
  ax_h(end+1) = axes('position',[leftMargin middle+(rectPlotHeight/2)+spacer*3 rectPlotWidth rectPlotHeight/2]);
  set(ax_h(end),'tag','smallRectPlot_top');
end

if(ismember('smallRectPlot_middle',params.axLocs))
  ax_h(end+1) = axes('position',[leftMargin middle+spacer*2 rectPlotWidth+0.20 rectPlotHeight/2]);
  set(ax_h(end),'tag','smallRectPlot_middle');
end

if(ismember('rectPlot_bottom',params.axLocs))
  ax_h(end+1) = axes('position',[leftMargin bottomMargin rectPlotWidth rectPlotHeight]);
  set(ax_h(end),'tag','rectPlot_bottom');
end

if(ismember('rectPlot_bottomRightVert',params.axLocs))
  ax_h(end+1) = axes('position',[leftMargin+rectPlotWidth bottomMargin 0.15 rectPlotHeight]);
  set(ax_h(end),'tag','rectPlot_bottomRightVert');  
end

if(ismember('infoBox_top',params.axLocs))
  ax_h(end+1) = axes('position',[leftMargin 0.81 0.775 0.1]);
  set(ax_h(end),'tag','infoBox_top');
  set(ax_h(end),'visible','off');
end

if(ismember('smallLegend_topRight',params.axLocs))
  ax_h(end+1) = axes('position',[0.8 0.87 0.15 0.05]);
  set(ax_h(end),'tag','smallLegend_topRight');
  set(ax_h(end),'visible','off');
end

if(ismember('smallLegend_bottomRight',params.axLocs))
  ax_h(end+1) = axes('position',[0.8 (bottomMargin+squarePlotDim-0.1) 0.15 0.05]);
  set(ax_h(end),'tag','smallLegend_bottomRight');
  set(ax_h(end),'visible','off');
end

if(ismember('smallBox_top',params.axLocs))
  ax_h(end+1) = axes('position',[leftMargin middle+rectPlotHeight rectPlotWidth 0.1]);
  set(ax_h(end),'tag','smallBox_top');
  set(ax_h(end),'xtick',[],'ytick',[]);
end

if(ismember('smallBox_bottom',params.axLocs))
  ax_h(end+1) = axes('position',[leftMargin bottomMargin+rectPlotHeight rectPlotWidth 0.1]);
  set(ax_h(end),'tag','smallBox_bottom');
  set(ax_h(end),'xtick',[],'ytick',[]);
end

if(ismember('squarePlot_top',params.axLocs))
  ax_h(end+1) = axes('position',[ (1-squarePlotDim)/2 middle squarePlotDim squarePlotDim]);
  axis square;
  set(ax_h(end),'tag','squarePlot_top');
end

if(ismember('squarePlot_bottom',params.axLocs))
  dims = [ (1-squarePlotDim)/2 bottomMargin squarePlotDim squarePlotDim];
  ax_h(end+1) = axes('position',dims);
  axis square
  set(ax_h(end),'tag','squarePlot_bottom');
end
    
if(ismember('movie_img',params.axLocs))
  ax_h(end+1) = axes('position',[0 0 1 1]);
  set(ax_h(end),'tag','movie_img');
end

if(ismember('movie_tapdata',params.axLocs))
  img_h = get_axhandle(ax_h,'movie_img');
  set(img_h,'position',[0 0 1 .8]);
  ax_h(end+1) = axes('position',[0 .8 1 .2]);
  set(ax_h(end),'tag','movie_tapdata');
  set(ax_h(end),'xtick',[],'ytick',[],'ydir','reverse');
end

%note that this subplot provides different subplot spacing than the
%"subplot" function. This is done to provide more room for labels.
if(ismember('subplot',params.axLocs))
  nRows = params.nRows;
  nCols = params.nCols;
  plotSpacing = (1 - topMargin - bottomMargin)/nRows;
  
  %only supports one column for now
  for iRow = 1:nRows
    
    for iCol = 1:nCols

      
      plotBottom = 1 - topMargin - iRow*plotSpacing;
      plotLeft = leftMargin + 0.07;
      
      ax_h(end+1) = axes('position',[plotLeft plotBottom rectPlotWidth plotSpacing*0.7]);
      set(ax_h(end),'tag',sprintf('subplot_%d',(iRow-1)*nCols + iCol));
      
      if(isfield(params,'colLabels') && (iRow == 1))
	axPos = get(ax_h(end),'position');
	axPos(1) = axPos(1)+axPos(3)/2 -0.04;
	axPos(2) = axPos(2) + axPos(4) + .02;
	axPos(3) = 0.08;
	axPos(4) = 0.02;
	ax_lbl_h = axes('position',axPos,'visible','off');
	thisLabel = params.colLabels{iCol};
	if(isnumeric(thisLabel))
	  thisLabel = num2str(thisLabel);
	end
	txt_h = text(0.5,0.5,thisLabel);
	set(txt_h,'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold');
	
      end
      
      if(iRow == 1 && iCol == 1)
	
	axPos = get(ax_h(end),'position');
	fontSize = .1;
	
	if(isfield(params,'xLabel'))
	  xlabel(ax_h(end),params.xLabel,'FontUnits','normalized','FontSize',fontSize);
	end
	if(isfield(params,'yLabel'))
	  ylabel(ax_h(end),params.yLabel,'FontUnits','normalized','FontSize',fontSize);
	end
      end
      
    end
    
    %display row labels
    if(isfield(params,'rowLabels'))
      axPos = get(ax_h(end),'position');
      axPos(1) = 0.01;
      axPos(2) = axPos(2) + axPos(4)/2;
      axPos(3) = 0.08;
      axPos(4) = 0.02;
      ax_lbl_h = axes('position',axPos,'visible','off');
      thisLabel = params.rowLabels{iRow};
      if(isnumeric(thisLabel))
	thisLabel = num2str(thisLabel);
      end
      txt_h = text(0.5,0.5,thisLabel);
      set(txt_h,'HorizontalAlignment','center','VerticalAlignment','middle','FontWeight','bold');
    end
    
  end
  
end   


return
