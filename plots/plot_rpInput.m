function plot_rpInput(ax,audArray,params)
% Plots the input signal of a BTB calculation
% 
% INPUTS
% ax:                         axes handle to plot the input signal
% audArray:                   array of aud structs used as inputs. This provides the
%                             ability to either plot multiple subject tapping or simply left and right
%                             hand tapping (if using midi data converted to aud using Dirac impulses)
% params.xlim:                xlim for the axes. Usually, an xlim is obtained from
%                             the periodicity surface so that the xlims match.
% params.inputPlot.normalize: 0 or 1. Whether or not to normalize
%                             the inputs
% params.trialSets:           separates trial types by grouping them and
%                             plotting them in alternating colors by type.
%                             This is a 2 column matrix where the first
%                             column indicates input index and the second
%                             column indicates a trial index.
%
%
% Copyright (c) 2007-2012 The Regents of the University of California
% All Rights Reserved
%
% Author:
% Stefan Tomic 4/07


xlimAx = params.xlim;

if(~isfield(params,'inputPlot') || ~isfield(params.inputPlot,'normalize'))
  params.inputPlot.normalize = 0;
end

axes(ax);
hold on

wvfPlotNum = 1;

audArrayCols = set_var_col_const(audArray{1}.vars);

wvfFs = audArray{1}.data{audArrayCols.Fs}(1);
    
for iAudGroup = 1:length(audArray)
  audStruct = audArray{iAudGroup};
  thisWvf = audStruct.data{audArrayCols.wvf};
  if(iscell(thisWvf))
    thisWvf = thisWvf{1};
  end
  
  thisWvfFs = audStruct.data{audArrayCols.Fs};
  
  if(any(diff(thisWvfFs)))
    error('Wvf overlays on RP output plot do not have matching sampling rates.');
  else      
    thisWvfFs = thisWvfFs(1);
  end
  
  if(isfield(params,'trialSets'))
    findTrialSet = params.trialSets(find(iAudGroup == params.trialSets(:,1)),2);
    thisColor = [1 0 0]*(mod(findTrialSet-1,2));
  else
    thisColor = [0 0 0];
  end
      
  if(~isempty(thisWvf))
    
    if(thisWvfFs ~= wvfFs)
	  error('Wvf overlays on RP output plot do not have matching sampling rates.');
    end
	  
    switch (params.inputSigType)
     case 'impulse'
      if(isfield(params,'overlayInput') && params.overlayInput)
	yMin = 0;
	yMax = 1;
      else
	yMin = wvfPlotNum-1;
	yMax = wvfPlotNum;
      end
      wvfPointsX = repmat(find(thisWvf)',2,1);
      %wvfPointsY = [yMin ; yMax] * ones(1,size(wvfPointsX,2));
      wvfPointsY = [yMax - (yMax - yMin) * thisWvf(find(thisWvf))';
		    yMax .* ones(1,size(wvfPointsX,2))];
      plot(wvfPointsX,wvfPointsY,'Color',thisColor);
      
     case 'wvf'
      yMin = 0;
      yMax = 1;
      %convert wvf to mono if necessary
      thisWvf = mean(thisWvf,2);
      %normalize
      if(params.inputPlot.normalize)
	thisWvf = thisWvf./max(abs(thisWvf)).*0.98;
      end

      %scale and offset wvf to be between 0 and 1
      %thisWvf = (thisWvf+1)./2;
      thisWvf = 1 - thisWvf;
      plot(thisWvf,'k');
	  
     case 'mid'
	  
      
    end	  
	  
   
    wvfPlotNum = wvfPlotNum + 1;
  else
    
    yMax = 1;
  end
end
 
set(ax,'xlim',xlimAx*wvfFs);
set(ax,'ylim',[0 yMax]);
set(ax,'xtick',[]);
set(ax,'xticklabel',[]);
set(ax,'ytick',[]);
set(ax,'yticklabel',[]);
set(ax,'ydir','reverse');
%set(ax,'visible','off');

return
