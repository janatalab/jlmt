function hfunc = resonatorBanks(select)
% Filters each channel of a multi channel signal through a resonator bank
%
% hfunc = resonatorBanks(select)
%
% returns a function handle to one of the two subfunctions for
% implementing comb filters or reson filters.
% h = resonatorBanks('reson') returns a handle to the reson subfunction.
% h = resonatorBands('comb') returns a handle to the comb subfunction.
%
% SUBFUNCTION RESON:
%  [wvfOut,resonatorFreqs,BList,AList,QList] = reson(wvfIn,params)
%
% INPUTS
%  wvfIn:                       waveform used as input to the filters. Expects a matrix
%                               where the rows are individual band signals and time is the x-axis.
%  params.freqLimits:           Two element vector. Frequency limits of reson filters
%  params.resonFreqSpacing:     'log2'    uses log base 2 spacing. 
%                                         params.numFilters determines the
%                                         number of filters
%                               'linear'  uses linear spacing. The
%                                         number of filters are determined by params.spacingHz.
%                               'filterQ' spaces filters apart by half
%                                         a bandwidth. Bandwidth is determined by
%                                         the filters' Q factor.
%  params.numFilters:           Specifies the number of filters to
%                               use. This parameter is only used for 'log2' spacing.
%  params.spacingHz:            Used for spacing the filters when using
%                               linear spacing. 
%  params.numPeriodDecay:       Filter Q-factor (also the number of
%                               periods the filters decay by approx. -27 dB, see BTB paper)   
%  params.gainType:             'constant'
%                               'weighted'
%                               'normalized_response'
%  params.gain:                 gain value if using constant gain.
%  params.gainFunction.method:  determines shape of rise of gain function
%                               across resonators. Currently 
%  params.gainFunction.minGain: minimum gain value.
%  params.gainFunction.maxGain: maximum gain value.
%  params.gainFunction.maxFreq: frequency at which gain rise
%                               reaches the maximum gain value.
%  params.Fs:                   sampling rate of input signal 
%                               (also sampling rate of output signal)
%  params.bw_measurement:       'constant'           - uses a constant
%                                                      bandwidth (in Hz) across all filters
%                               'halfPowerBandwidth' - determines filter bandwidth so that filter magnitude
%                                                      meet at their half power points.
%                               'periodBasedDecay'   - A constant Q factor across all the filters is
%                                                      used to determine their bandwidths.
%  params.spacingFactor:        when params.resonFreqSpacing is set
%                               to 'filterQ', this is multiplied by the bandwidth of the
%                               filters to determine the spacing of filters.
%  params.bandwidth:            when using a constant bandwidth in  Hz (params.bw_measurement = 'constant'),
%                               this is the bandwidth used in Hz.
%
% Copyright (c) 2006 The Regents of the University of California
% All Rights Reserved
%
% Author(s):
% Stefan Tomic 12/06

  switch select
   case {'reson','reson_r','reson_z'}
    hfunc=@reson;
   case 'comb'
    hfunc=@comb;
   case 'partialHanning'
    hfunc=@partialHanning;
   case 'betaDistribution'
    hfunc=@betaDistribution;
   otherwise
    hfunc=[];
    error('Unknown Input');
  end
return

function [wvfOut,resonatorFreqs,BList,AList,QList,RList] = reson(wvfIn,params)
  disp('Filtering signal(s) through reson filterbank(s)...');

  freqLimits = params.freqLimits;
 
  switch(params.resonFreqSpacing)
   case 'log2'
    numFilters = params.numFilters; 
    resonatorFreqs = logspace2(log2(freqLimits(1)),log2(freqLimits(2)),numFilters);

   case 'linear'
    spacingHz = params.spacingHz; 
    resonatorFreqs = [freqLimits(1):spacingHz:freqLimits(2)];
    numFilters = length(resonatorFreqs);
    
   case 'filterQ'
    Q = params.numPeriodDecay;
    freqIdx = 1;
    resonatorFreqs(freqIdx) = freqLimits(1);
    while(resonatorFreqs(freqIdx) <= freqLimits(2))
      freqIdx = freqIdx + 1;
      lastBW = resonatorFreqs(freqIdx-1)./Q;
      resonatorFreqs(freqIdx) = resonatorFreqs(freqIdx-1) + params.spacingFactor*lastBW;
    end
    numFilters = length(resonatorFreqs);
    
  end
  
  if(ismember(params.gainType,{'weighted','prenorm_weighted'}))
    gainFunc = str2func(params.gainFunction.method);
    params.gainFunction.resonatorFreqs = resonatorFreqs;
    gainVector = gainFunc(params.gainFunction);
  end
    
  if(~isfield(params,'bw_measurement'))
    params.bw_measurement = 'constant';
  end
   
  for ibank = 1:size(wvfIn,1)
    disp(['channel ' num2str(ibank) ' ']);
    for iReson = 1:numFilters
      
      switch(params.bw_measurement)
       case 'constant'
	bw = params.bandwidth;

       case 'halfPowerBandwidth'
	if(iReson == numFilters)
	  percentBw = (resonatorFreqs(iReson) - resonatorFreqs(iReson-1))./resonatorFreqs(iReson-1);
	  bw = percentBw*resonatorFreqs(iReson);
	else
	  bw = (resonatorFreqs(iReson+1) - resonatorFreqs(iReson));
	end
       
       case 'periodBasedDecay'
	Q = params.numPeriodDecay;
	bw = resonatorFreqs(iReson)./Q;
      end

      R = 1 - bw./params.Fs*pi;
      theta = resonatorFreqs(iReson)./params.Fs*(2*pi);    

      switch(params.gainType)
       case 'constant'
	gainFactor = params.gain;
       case 'weighted'
	gainFactor = gainVector(iReson);
       case 'prenorm_weighted'
	%normalizing before multiplying by beta gain
	gainFactor = (1-R^2)/2.*gainVector(iReson);	
       case 'normalized_response'
	switch(params.method)
	  %normalization method depends on the type of reson filter used
	 case 'reson_r'
	  gainFactor = 1-R; 
	 case 'reson_z'
	  gainFactor = (1-R^2)/2;
	 case 'reson'
	  gainFactor = (1 - R^2)*sin(theta);
	end
      end

      % approximating actual peak magnitude response from theta
      % (Steiglitz, 1994)

      switch(params.method)
       case 'reson_z'
	B_unscaled = [1 0 -1];
	psi = acos((2*R)./(1+R^2).*cos(theta));
       case 'reson'
	B_unscaled = [1];
	psi = acos((1+R^2)./(2*R)*cos(theta));
       case 'reson_r'
	psi = theta; %we have no peak response correction equation for reson_r
	B_unscaled = [1 0 -R];
      end

      realFreq = psi./(2*pi).*params.Fs;
      resonatorFreqs(iReson) = realFreq;
      
      B = B_unscaled.*gainFactor;
      A = [1 -2*R*cos(theta) R^2];
      
      wvfOut(ibank,iReson,:) = filter(B,A,wvfIn(ibank,:));
    
      %the filter poles, zeros, and Q-factors can be returned for testing/tuning
      BList{iReson} = B;
      AList{iReson} = A;
      QList(iReson) = resonatorFreqs(iReson)./bw;
      RList(iReson) =  R;
    end
  end
return



function [wvfOut,resonatorFreqs,BList,AList,QList] = comb(wvfIn,params)
    disp('Filtering signal(s) through comb filterbank(s)...');
   
    freqLimits = params.freqLimits;
    numFilters = params.numFilters;
    resonatorFreqs = logspace2(log2(freqLimits(1)),log2(freqLimits(2)),numFilters);
  
    tSamps = params.halfEnergyTimeSecs*params.Fs;
 
    resonatorPeriodSecs = 1./resonatorFreqs;
    resonatorPeriodSamps = floor(resonatorPeriodSecs*params.Fs);
    
    %calculation for alpha (gain)
    alpha = 0.5.^(resonatorPeriodSamps./tSamps);

    if(strcmp(params.gainType,'weighted'))
      gainFunc = str2func(params.gainFunction.method);
      params.gainFunction.resonatorFreqs = resonatorFreqs;
      gainVector = gainFunc(params.gainFunction);
    end
  
    for ibank = 1:size(wvfIn,1)
      disp(['channel ' num2str(ibank) ' ']);
      for iComb = 1:numFilters
	A = [1 zeros(1,resonatorPeriodSamps(iComb)-1) -1*alpha(iComb)];
	B = 1-alpha(iComb);

	if(strcmp(params.gainType,'weighted'))
	  B = B.*gainVector(iComb);
	end
	
	wvfOut(ibank,iComb,:) = filter(B,A,wvfIn(ibank,:));
	
	%the filter poles and zeros can be returned for testing/tuning
	BList{iComb} = B;
	AList{iComb} = A;
	QList(iComb) = NaN;
      end
    end

return


function gainVector = partialHanning(params)
%weighting function for resonator gains
gainVector = zeros(size(params.resonatorFreqs));
gainVector(find(params.resonatorFreqs >= params.maxFreq)) = params.maxGain;  
gainCurveLength = sum(params.resonatorFreqs < params.maxFreq);
hanningCurve = hanning(gainCurveLength*2);
gainVector(find(params.resonatorFreqs < params.maxFreq)) = params.minGain + (params.maxGain-params.minGain).*hanningCurve(1:gainCurveLength);
  
return

function gainVector = betaDistribution(params)

freqRange = params.resonatorFreqs - min(params.resonatorFreqs);
xvals = freqRange./max(freqRange);
betaDist = betapdf(xvals,params.betaDistribution.alpha,params.betaDistribution.beta);

if(isfield(params.betaDistribution,'flattenAt') && ~isempty(params.betaDistribution.flattenAt))
  [smallestDiff,closestFreqIdx] = min(abs(params.resonatorFreqs - params.betaDistribution.flattenAt));
  scaleFactor = (params.maxGain - params.minGain)./betaDist(closestFreqIdx);
  gainVector = betaDist.*scaleFactor + params.minGain;
  gainVector(find(gainVector > params.maxGain)) = params.maxGain;
else
  gainVector = betaDist./max(betaDist).*(params.maxGain - params.minGain) + params.minGain;
end

return
