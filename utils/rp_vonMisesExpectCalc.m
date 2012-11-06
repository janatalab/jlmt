function outData = rp_vonMisesExpectCalc(varargin)
% Calculates an expectation error based on von Mises distributions across metric levels
% 
% outData = rp_vonMisesExpectCalc(inData,params)
% 
% INPUTS (tag/value pairs):
%         rp - a rhythm profiler struct
%   inputSig - the input signal to the resonator model (denv, sumDenv, or any
%              other signal)
%
%     RELEVANT PARAMS:
%              params.vmExpect.distResolution: The bin size of the
%                 von Mises distribution
%
% Copyright (c) 2009-2012 The Regents of the University of California
% All Rights Reserved.
%
% Author: Stefan Tomic 04/23/2009
% ~9/2009 S.T. - Considerably cleaned up code 
% 2/22/10 S.T. - Fixed handling of multiple bands after cleanup
  
PLOT_METER_PREDICTION = 1; 
  
for iarg = 1:2:nargin
  switch varargin{iarg}
   case 'rp'
    rp = varargin{iarg+1};
   case 'inputSig'
    inputSig = varargin{iarg+1};
  end
end

rpCols = set_var_col_const(rp.vars);
params = rp.data{rpCols.params}.rp;

%for better results, use a distResolution that is odd (results in an even
%number of bins)
try
  params.vmExpect.distResolution = 100;
catch
  warning('von Mises distribution resolution not set. Setting to 100.');
  params.vmExpect.distResolution = 100;
end

distResolution = params.vmExpect.distResolution;
angleDistVector = [-pi:2*pi/distResolution:pi];


%obtain idx
if(iscell(rp))
  rp = rp{1};
end

rp = ensemble_datastruct2tree(rp);

onsetInfo = rp.onsetInfo;

peakFinderParams =  params.resonatorBank.peakFinder;
resonFreqs = rp.resonatorFreqs;
nResons = length(resonFreqs);
nBands = length(onsetInfo.onsetTimesSampsByBand);

PLOT_METER_PREDICTION = PLOT_METER_PREDICTION && (nBands == 1); % a little kludge so that I can
                                            % plot on one band in one
                                            % environment while I run real
                                            % data in another environment



if(params.perform.vmExpect.useInterpPeak)
  bandPeakInfo = rp.interpMPPStruct.bandPeakInfo;
else
  bandPeakInfo = rp.peakInfo.bandPeakInfo;
end

for iBand = 1:nBands

  onsetTimesThisBand = onsetInfo.onsetTimesSampsByBand{iBand};
  
  if(length(onsetTimesThisBand) <=1)
    vmPredictionPerOnset{iBand} = [];
    vmFulfillPerOnset{iBand} = [];
    vmLogFulfillPerOnset{iBand} = [];
    vmAreaPerPeakList{iBand} = [];
    vmProportionsPerPeakList{iBand} = [];
    vmPeakFreqs{iBand} = [];
    vmK{iBand} = [];
    vmTimingDev{iBand} = [];
    vmDistPerPeakPerOnset = [];
    continue;
  end
     
  nOnsets = length(onsetTimesThisBand);
  
  %obtain the denv magnitudes at onset times
  denvAtOnsets = squeeze(inputSig(iBand,onsetTimesThisBand));
  denvThisBand = inputSig(iBand,:);
  
  peakFreqs = bandPeakInfo.peakFreq{iBand};
  peakWidths = bandPeakInfo.baseWidths{iBand};
  peakHeights = bandPeakInfo.normHeight{iBand};
    
  %normalizing by RMS here
  peakIdxs = rp.peakInfo.peakIdxs; %get peak idxs of non-interpolated peaks
  rmsResonThisBand = sqrt(mean(squeeze(rp.resonatorOut(iBand,peakIdxs,:)).^2,2));
  meanRMSThisBand = mean(rmsResonThisBand);

  onsetExpectByBand = rp.complexExpectAtOnsets{iBand};
  
  %removing this normalization for now
  %onsetExpectByBand = onsetExpectByBand ./ meanRMSThisBand;

  %denvEnergy = sqrt(mean(denvThisBand.^2));
  %adjDenvAtOnsets = denvAtOnsets ./ denvEnergy;
    
  
  clear normDenvPerOnset    maxProbPerOnset normProbPerOnset normProbPerOnset2;
  
  for iOnset = 1:nOnsets
    
    thisComplexResonOut = onsetExpectByBand(:,iOnset);
    thisComplexResonMags = abs(thisComplexResonOut);
    thisComplexResonPhases = angle(thisComplexResonOut);
  
    clear peakPhases peakMags
    if(params.perform.vmExpect.useInterpPeak)
      
      nInterpPeaks = length(peakFreqs);
      for iInterpPeak = 1:nInterpPeaks
	thisPeak = peakFreqs(iInterpPeak);
	lBorderReson = find(resonFreqs <= thisPeak,1,'last');
	rBorderReson = find(resonFreqs >= thisPeak,1,'first');
	pct(iInterpPeak) = (thisPeak - resonFreqs(lBorderReson)) ./ ...
	    (resonFreqs(rBorderReson) - resonFreqs(lBorderReson));
	if(~isfinite(pct(iInterpPeak)))
	  pct(iInterpPeak) = 0;
	end
		
	lAngle = angle(thisComplexResonOut(lBorderReson));
	rAngle = angle(thisComplexResonOut(rBorderReson));
	
	lMag = thisComplexResonMags(lBorderReson);
	rMag = thisComplexResonMags(rBorderReson);
	diffMag = diff([lMag rMag]);
	
	peakPhases(iInterpPeak) = interpAngle(lAngle,rAngle,pct(iInterpPeak));
	peakMags(iInterpPeak) = lMag + pct(iInterpPeak) .* diffMag;
	
      end      
   
    else
      peakPhases = angle(thisComplexResonOut(peakIdxs))';
      peakMags = thisComplexResonMags(peakIdxs);
    end
    
    nPeaks = length(peakFreqs);
    
    
    clear k weight variance vmProbPerPeak vmDistPerPeak vmAreaPerPeak vmAdjProbPerPeak ...
	vmAdjProbPerPeakTest metricLevel pwProportion evalAngleIdxPerPeak;
  
    for iPeak = 1:nPeaks
	
      pwProportion(iPeak) = peakWidths(iPeak)./ peakFreqs(iPeak)./peakHeights(iPeak);
      variance(iPeak) = pwProportion(iPeak);

      k(iPeak)  = 1./(variance(iPeak));
      vmDistPerPeak(iPeak,:) =  von_mises_pdf(angleDistVector,0,k(iPeak));

      diffAngle = peakPhases(iPeak) - angleDistVector;
      evalAngleIdx = find(abs(diffAngle) == min(abs(diffAngle)));
      evalAngleIdxPerPeak(iPeak) = evalAngleIdx;
     
      %USING PEAK HEIGHTS HERE INSTEAD OF RESONMAGS
      %RESONMAGS ARE GOING TO SHIFT AROUND A LOT
      sumResonMags = sum(peakMags);
      sumPeakHeights = sum(peakHeights);
      
      %if(sumPeakHeights ~= 0)
      if(sumResonMags ~= 0)
	weight(iPeak) = peakMags(iPeak) ./ sumResonMags;
	%weight(iPeak) = peakHeights(iPeak) ./ sumPeakHeights;
      else
	weight(iPeak) = 0;
      end

      vmProbPerPeak(iPeak) = vmDistPerPeak(iPeak,evalAngleIdx).*weight(iPeak);

      
      %find the indices to evaluate for the area of the von mises
      %dist (timing dev calc)
      %not including index corresponding to angle zero. This could result
      %in a value greater than .5 (since zero bin is directly in the middle
      %when using odd number of bins).
      if(peakPhases(iPeak) > 0)
	areaIdxs = find((angleDistVector > 0) & (angleDistVector <= peakPhases(iPeak)));
      else
	areaIdxs = find((angleDistVector >= peakPhases(iPeak)) & (angleDistVector < 0));
      end
     
      vmAreaPerPeak(iPeak) = sum(vmDistPerPeak(iPeak,areaIdxs).*...
				 weight(iPeak).*(pi/distResolution));    
      
     
      if(PLOT_METER_PREDICTION)
	
	nPhasesStartOffset = - 0.1;
	startPhasePredict = 2*pi*nPhasesStartOffset;
	nPhasesPredict = 2;
	endPhasePredict = 2*pi*nPhasesPredict;
	pPlotDistRes = 500;
	
	
      	%constructing an overall distribution matrix
	if(iPeak == 1)
	  clear metricLevel;
	  basePeriod = 1./peakFreqs(iPeak);
	  nPredictPlotSec = basePeriod * nPhasesPredict;
	  startPredictSecPlot = nPhasesStartOffset * basePeriod;
	  evalDist = [startPhasePredict:(endPhasePredict-startPhasePredict)./pPlotDistRes:endPhasePredict];
	  distThisPeak = von_mises_pdf(evalDist,peakPhases(iPeak),k(iPeak));
	  metricLevel(1,:) = distThisPeak.*weight(iPeak);	  
	  numSamps = size(metricLevel,2);
	else
	  thisPeriod = 1./peakFreqs(iPeak);
	  periodRatio = basePeriod/thisPeriod;
	  angleBegin = startPhasePredict * periodRatio;
	  angleEnd = endPhasePredict*periodRatio;
	  evalDist = [angleBegin:(angleEnd-angleBegin)/pPlotDistRes:angleEnd];
	  distThisPeak = von_mises_pdf(evalDist,peakPhases(iPeak),k(iPeak));
	  %distThisPeak = (distThisPeak - min(distThisPeak)) ./ ...
	   %   (max(distThisPeak) - min(distThisPeak));
	  %resampDist = angleBegin:(angleEnd-angleBegin)/(numSamps-1):angleEnd;
	  %distThisPeakResamp = spline(evalDist,distThisPeak,resampDist);
	  metricLevel(iPeak,:) = distThisPeak .* weight(iPeak);
	end
      end %PLOT_METER_PREDICTION
      
    end %for iPeak
    
    if(PLOT_METER_PREDICTION)
      Fs = rp.Fs;
      thisOnsetTimeSamps = onsetTimesThisBand(iOnset);
      startPredictSampsPlot = startPredictSecPlot * Fs;
      startTimeSamps = onsetTimesThisBand(iOnset) + startPredictSampsPlot;
      endTimeSamps = onsetTimesThisBand(iOnset) + nPredictPlotSec*Fs;
      figure(9);
      set(gcf,'position',[100 100 800 1000],'paperposition',[0.25 0.5 8 10]);
      
      resAx = axes('position',[0.05 .5 .9 .2]);
      hold off;
      plot([startPhasePredict:(endPhasePredict-startPhasePredict)./(numSamps-1):endPhasePredict],sum(metricLevel,1));
      set(gca,'xlim',[startPhasePredict endPhasePredict]);
      set(gca,'xtick',[0:pi/2:nPhasesPredict*2*pi],'xticklabel',{'0','pi/2','pi',...
		    '3/2 pi','2 pi','2 1/2 pi','3 pi','3 1/2 pi','4 pi'});

      inputAx = axes('position',[0.05 .7 .9 .2]);
      bar([startPredictSecPlot:1/Fs:nPredictPlotSec],inputSig(startTimeSamps:endTimeSamps));
      hold on;
      %plot the current onset in red
      bar([-1:1]./Fs,inputSig(thisOnsetTimeSamps-1:thisOnsetTimeSamps+1),'r');
      set(gca,'xlim',[startPredictSecPlot nPredictPlotSec]);
      set(gca,'xtick',[0:.5:endTimeSamps./Fs]);
      set(gca,'ytick',[0:.2:1]);

    
      clf;
    
    end
    
    %normalize denvAtOnsets and vmProbPerPeak by a specified window size
    %for fulfillment calculation
    normWindow = params.vmExpect.normalizeWindow;
    normWindowSamps = round(normWindow .* params.Fs);
    windowStart = max(1,(onsetTimesThisBand(iOnset) - normWindowSamps));
    onsetIdxsThisWindow = find(onsetTimesThisBand >= windowStart & ...
			       onsetTimesThisBand <= ...
			       onsetTimesThisBand(iOnset));
    
    normDenvPerOnset(iOnset) = denvAtOnsets(iOnset) ./ max(denvAtOnsets(onsetIdxsThisWindow));
    vmPredictionPerOnset{iBand}(iOnset) = sum(vmProbPerPeak);
    
    %normalizing Prob by maximum calculated prediction value over a window
    %normProbPerOnset(iOnset) = vmPredictionPerOnset{iBand}(iOnset) ./ max(vmPredictionPerOnset{iBand}(onsetIdxsThisWindow));
    
    %normalizing Prob by maximum possible distribution value in window
    %maxProbPerOnset(iOnset) = sum(max(vmDistPerPeak,[],2) .* weight');
    %maxProbThisWindow = max(maxProbPerOnset(onsetIdxsThisWindow));
    %normProbPerOnset(iOnset) = (sum(vmProbPerPeak) ./ maxProbThisWindow);  

    %normalizing Prob by maximum Distribution only for that onset
    maxProb = sum(max(vmDistPerPeak,[],2) .* weight');
    normProbPerOnset(iOnset) = (sum(vmProbPerPeak) ./ maxProb);
        
    %fulfillThisOnset = normDenvPerOnset(iOnset) ./ normProbPerOnset(iOnset);
    fulfillThisOnset = denvAtOnsets(iOnset) ./ normProbPerOnset(iOnset);

    %trying simple inverse of probability
    % fulfillThisOnset = 1 ./ normProbPerOnset(iOnset);
    
    
    %vmPredictionPerOnset{iBand}(iOnset) = sum(vmProbPerPeak);
    vmFulfillPerOnset{iBand}(iOnset) = fulfillThisOnset;
    vmLogFulfillPerOnset{iBand}(iOnset) = log(fulfillThisOnset);
    vmAreaPerPeakList{iBand,iOnset} = vmAreaPerPeak;
    vmProportionsPerPeakList{iBand,iOnset} = variance;
    vmK{iBand,iOnset} = k;
    vmTimingDev{iBand}(iOnset) = sum(vmAreaPerPeak);
    vmDistPerPeakPerOnset{iOnset} = vmDistPerPeak;
    
  end %for iOnset 
  
  vmPeakFreqs{iBand} = peakFreqs;
  
    
end % for iBand

expectationData = ensemble_init_data_struct;
expectationData.vars{end+1} = 'vmAngles';
expectationData.vars{end+1} = 'vmPeakFreqs';
expectationData.vars{end+1} = 'vmPredictionPerOnset';
expectationData.vars{end+1} = 'vmAreaPerPeak';
expectationData.vars{end+1} = 'vmFulfillPerOnset';
expectationData.vars{end+1} = 'vmLogFulfillPerOnset';
expectationData.vars{end+1} = 'vmTimingDev';
expectationData.vars{end+1} = 'vmK';
expectationData.vars{end+1} = 'proportions';
expectationData.vars{end+1} = 'vmDistPerOnset';
expectationData.vars{end+1} = 'vmAngles';
expDataCols = set_var_col_const(expectationData.vars);
expectationData.data{expDataCols.vmAngles} = angleDistVector';
expectationData.data{expDataCols.vmPeakFreqs} = vmPeakFreqs;
expectationData.data{expDataCols.vmPredictionPerOnset} = vmPredictionPerOnset;
expectationData.data{expDataCols.vmFulfillPerOnset} = vmFulfillPerOnset;
expectationData.data{expDataCols.vmLogFulfillPerOnset} = vmLogFulfillPerOnset;
expectationData.data{expDataCols.vmAreaPerPeak} = vmAreaPerPeakList;
expectationData.data{expDataCols.proportions} = vmProportionsPerPeakList;
expectationData.data{expDataCols.vmTimingDev} = vmTimingDev;
expectationData.data{expDataCols.vmDistPerOnset} = vmDistPerPeakPerOnset;
expectationData.data{expDataCols.vmK} = vmK;
expectationData.data{expDataCols.vmAngles} = angleDistVector;
outData = expectationData;
