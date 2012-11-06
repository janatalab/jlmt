function [mppByFrame,stdmpByFrame,ratioBinsByFrame,ratioEntropyByFrame,ratioTally] = calc_mppByFrame(resonEnergy,resonFreqs,params)
% calculates a mean metric profile (MMP) on a frame-by-frame basis
%
% [mppByFrame,stdmpByFrame,ratioBinsByFrame,ratioEntropyByFrame,ratioTally] = calc_mppByFrame(resonEnergy,resonFreqs,params)function outData = rp_findResonPhases(inData,params)
%
%  INPUT
% 
%
%  OUTPUT
%
% Copyright (c) 2009-2012 The Regents of the University of California
% All Rights Reserved 
%
% Author:
% Jan, 2009 - Stefan Tomic 


try
  fid = params.logFileFid;
catch
  fid = 1;
end


Fs = params.Fs;
frameSizeSamps = params.mppByFrame.frameSize*Fs;
hopSizeSamps = params.mppByFrame.hopSize*Fs;
nSamps = size(resonEnergy,2);
nFrames = floor((nSamps - frameSizeSamps - 1)/hopSizeSamps);
ratioEntropyByFrame = nan(nFrames,1);

fprintf(fid,'Calculating mean metric profiles by frame...\n');
fprintf(fid,'Frame ');

nRatios = length(params.resonatorBank.peakFinder.ratioLabels);
ratioBinsByFrame = zeros(nRatios,nFrames);

ratioTally = [];


for iFrame = 1:nFrames
  
  if(mod(iFrame,100) == 0)
    fprintf(fid,[num2str(iFrame) ' ']);
  end
  
  startSamp = (iFrame-1)*hopSizeSamps+1;
  stopSamp = startSamp + frameSizeSamps;
  
  mmpThisFrame = mean(resonEnergy(:,startSamp:stopSamp),2);
  stdmpThisFrame = std(resonEnergy(:,startSamp:stopSamp),[],2);
  
  mppByFrame(:,iFrame) = mmpThisFrame;
  stdmpByFrame(:,iFrame) = stdmpThisFrame;
  
  params.resonatorBank.peakFinder.doReport = 0;
  peakInfo = get_resonPeakInfo(mppByFrame(:,iFrame),resonFreqs,params.resonatorBank.peakFinder);
  
  %use only peaks that successfully binned  
  idxsBinSuccesses = peakInfo.ratioBinIdx(~isnan(peakInfo.ratioBinIdx));
  binSuccessHeights = peakInfo.normHeight(~isnan(peakInfo.ratioBinIdx));
  
  ratioBinsByFrame(idxsBinSuccesses,iFrame) = binSuccessHeights;
  ratioEntropyByFrame(iFrame) = sigEntropy(ratioBinsByFrame(:,iFrame));

  
  heightRatios = peakInfo.normHeight./max(peakInfo.normHeight);
  nNewRatios = length(heightRatios);
  nTally = size(ratioTally,1);
  ratioTally(nTally+1:nTally+nNewRatios,1) = peakInfo.normPeakFreq;
  ratioTally(nTally+1:nTally+nNewRatios,2) = heightRatios;
  ratioTally(nTally+1:nTally+nNewRatios,3) = peakInfo.ratioBinIdx;
  ratioTally(nTally+1:nTally+nNewRatios,4) = peakInfo.errorVec;
  
end

