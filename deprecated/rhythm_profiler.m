function outData = rhythm_profiler(inData,params)
%  This is the entry point for the BTB algorithm
%
%  outData = rhythm_profiler(inData,params)
%
%
%  if inData is the string 'getDefaultParams', the default
%  parameter structure will be returned.
%
%  Otherwise, inData can either be an ani struct or an aud
%  struct. The model will process both data structs the same way.
%  Passing in an aud struct allows for bypassing the ANI process,
%  which may not be applicable for inputs such as impulse waveforms
%  that mark onset patterns.
%  
%  return the default parameter struct.
%
%                        bypassOnsetDetect: logical 1 or 0. If 0,
%                                           then the onset detection section of the
%                                           model is skipped (for processing data when
%                                           onset detection isn't appropriate, such as
%                                           Dirac impulses)
%                       convertInputToMono: logical 1 or 0. Only
%                                           applicable if bypassOnsetDetect is
%                                           true. Will convert the input signal to mono
%                                           by taking the mean signal of all the channels.
%                            ipemRMSMethod: logical 1 or 0 (determines whether to use IPEM
%                                           RMS method or not)
%                               inDataType: 'ani' or 'aud', specifies the type of inData used 
%                                           (see comment above)
%                    downSampleInputFactor: used to downsample 'aud' inputs before processing
%                                           (set to empty vector [] for no downsampling)
%                                       Fs: output sampling rate for RMS calculation (or
%                                           converted sample rate if onset detection is bypassed)
%                            rmsFrameWidth: framewidth in seconds for calculating RMS
%                            sumBandsAfter: 'RMS','diff', or
%                                           'HWR'. Whether to sum the ANI channels
%                                           after calculating RMS, difference, or
%                                           half-wave rectification. Normal value
%                                           is 'HWR'. Others were only tested for comparison.
%                         sumAdjacentBands: number of adjacent bands to sum.
%                 mulawCompression.perform: logical 0 or 1. Whether or not to perform mu-law
%                                           compression, as per Klapuri, et al. (2006). 
%                                           Normally, this is  turned off
%                      mulawCompression.mu: value to use for mu if performing mulawcompression
%                                           since it did not improve performance.
%                         params.hwrWeight: If this parameter is present, will perform a weighted average
%                                           between the signal and its result after
%                                           HWR and differentiation (see Klapuri et
%                                           al. 2006, eq 3). This did not improve results
%                                           so it is normally turned off.
%                     resonatorBank.method: 'reson' or 'comb', specifies
%                                           type of resonator filters to use
%                 resonatorBank.freqLimits: a 2 element vector that descibes
%                                           the frequency limits in Hz of the resonatorbank (e.g. [0.5 20])
%                         resonatorBank.Fs: sampling rate of the resonatorbank (nominally the same as Fs)
%                       resonatorBank.gain: gain of the resonator filters
%                  resonatorBank.bandwidth: bandwidth of the resonator filters
%         resonatorBank.halfEnergyTimeSecs: half time decay of resonator (comb filter only)
%           resonatorBank.energyWindowSecs: Window size of energy calculation in seconds.
%             resonatorBank.bw_measurement: string that describes
%                                           the method of bandwidth
%                                           calculation (either
%                                           'constant' or
%                                           'halfPowerBandWidth'
%                                           are currently
%                                           supported)
%                 resonatorBank.peakFinder: param struct for peak
%                                           finder function (see get_resonPeakInfo.m and find_peaks.m)
%
% Copyright (c) 2006-2013 The Regents of the University of California
% All Rights Reserved.
%
% Author: Stefan Tomic 12/2006

if(ischar(inData) & strcmp(inData,'getDefaultParams'))
  if(exist('params','var'))
    outData = getDefaultParams(params);
  else
    outData = getDefaultParams;
  end
  return
end

if(params.bypassOnsetDetect)
  sigStruct = inData;
  sigStructCols = set_var_col_const(sigStruct.vars);
  inSig.wvf = sigStruct.data{sigStructCols.wvf}';
  inSig.Fs = sigStruct.data{sigStructCols.Fs};

  if(size(inSig.wvf,1) > 1 && params.convertInputToMono)
    inSig.wvf = mean(inSig.wvf,1);
  end
  
  origFs = sigStruct.data{sigStructCols.Fs};
  
  if(params.Fs ~= inSig.Fs)
    %since resampling the input signal can drastically change some
    %signals (such as impulses, throw a warning here)
    warning(['RESONATOR MODEL IS RESAMPLING INPUT SIGNAL. This can' ...
	     ' drastically change inputs, especially impulses. Make' ...
	      ' sure that you check the resulting resampled signal.']);
    
    inSig.wvf = resample(inSig.wvf,params.Fs,inSig.Fs); 
    inSig.Fs = params.Fs;
  end
 
else
  aniStruct = inData;
  aniStructCols = set_var_col_const(aniStruct.vars);
  inSig.wvf = aniStruct.data{aniStructCols.sig};
  inSig.Fs  = aniStruct.data{aniStructCols.Fs};
end

if(params.invertInputSig)
  inSig.wvf = invert_sig(inSig.wvf);
end

%initialize rp struct

rp = init_rp_struct;
rpCols = set_var_col_const(rp.vars);

%save ani and rp params structs in the rp struct
if(strcmp(params.inDataType,'ani'))
  saveParams.ani = aniStruct.data{aniStructCols.params}.ani;
end
saveParams.rp = params;
rp.data{rpCols.params} = saveParams;


% do not perform RMS and denv if aud structure is passed in
% assumes that aud struct contains impulses
if(~params.bypassOnsetDetect)
  %%%%%% Calc RMS %%%%%%%%%
  
  disp('Calculating RMS...');
  
  if(params.ipemRMSMethod)
    
    [outRms,outRmsFs] = IPEMCalcRMS(inSig.wvf,inSig.Fs,params.rmsFrameWidth,1/params.Fs);
    
    rp.data{rpCols.rms} = outRms;
    rp.data{rpCols.Fs} = outRmsFs;
    
    disp('IPEM RMS method');
    
  else
  
    disp('RMS coded by S.T.');
    framewidthsamps = round(params.rmsFrameWidth*inSig.Fs);
    rms_rows = size(inSig.wvf,1);
    rms_cols = size(inSig.wvf,2)-framewidthsamps;
    
    rp.data{rpCols.rms} = zeros(rms_rows,rms_cols);
    for(irms = 1:rms_cols)
      frame_idxs = [irms:irms+framewidthsamps];
      rms_frame = inSig.wvf(:,frame_idxs);
      rp.data{rpCols.rms}(:,irms) = sqrt(sum(rms_frame.^2,2)./size(rms_frame,2));    
    end

    rp.data{rpCols.rms} = resample(rp.data{rpCols.rms}',params.Fs,inSig.Fs)';
    rp.data{rpCols.Fs} = params.Fs;
    
    Flp = 10; %cutoff freq in Hz
    Wn = Flp*2/rp.data{rpCols.Fs};
    [b,a] = butter(12,Wn);
    nbands = size(rp.data{rpCols.rms},1);
    for iband = 1:nbands
      rp.data{rpCols.rms}(iband,:) = filtfilt(b,a,rp.data{rpCols.rms}(iband,:));
    end  
    
  end
  
  if(strcmp(params.sumBandsAfter,'RMS'))
    rp.data{rpCols.rms} = sumBands(rp.data{rpCols.rms},params.sumAdjacentBands);
    rp.data{rpCols.meanChannelCBU} = meanBandCBU(rp.data{rpCols.params});
  end
  
  if(params.mulawCompression.perform)
    mu = params.mulawCompression.mu;
    rp.data{rpCols.rmsCompressed} = log(1+mu.*rp.data{rpCols.rms})./log(1+mu);
    rmsOut = rp.data{rpCols.rmsCompressed};
  else
    rmsOut = rp.data{rpCols.rms};
  end
  
  rp.data{rpCols.rms} = rmsOut;
  
  %%%%%%%% Calc Diff %%%%%%%%
  
  disp('Calculating difference of the RMS...');
    
  denv = diff(rmsOut,1,2); %simple difference		
  denv = [zeros(size(denv,1),1) denv];
  
  if(strcmp(params.sumBandsAfter,'diff'))
    denv = sumBands(denv,params.sumAdjacentBands);
    rp.data{rpCols.meanChannelCBU} = meanBandCBU(rp.data{rpCols.params});
  end
  
  denv(find(denv < 0)) = 0; %half-wave rectify the difference 
			      
  if(isfield(params,'hwrWeight'))			      
    %performing weighting of HWR(denv) Klapuri 2006
    lambda = params.hwrWeight;
    denv = (1-lambda).*rmsOut + (lambda.*(params.Fs)./Flp.*denv);
  end
      
  rp.data{rpCols.denv} = denv;
  
  if(strcmp(params.sumBandsAfter,'HWR'))
    denv = sumBands(denv,params.sumAdjacentBands);
    rp.data{rpCols.meanChannelCBU} = meanBandCBU(rp.data{rpCols.params});
  end
  
  rp.data{rpCols.sumDenv} = denv;
  inputToResonators = denv;
    
  
else
  
  denv = inSig.wvf;  
  rp.data{rpCols.Fs} = inSig.Fs;
  rp.data{rpCols.denv} = denv;
  rp.data{rpCols.sumDenv} = denv;
  inputToResonators = denv;
  
end %if(~params.bypassOnsetDetect)
%%%%%%%% Calc Reson Banks %%%%%%%%%%%%%%

disp('Calculating Reson Filterbanks for each frequency band...');

%resonator bank Fs = rms signal Fs
params.resonatorBank.Fs = params.Fs;

hfunc_resonatorBanks = resonatorBanks(params.resonatorBank.method);
[rp.data{rpCols.resonatorOut},rp.data{rpCols.resonatorFreqs},BList,AList,QList,RList] = ...
    hfunc_resonatorBanks(inputToResonators,params.resonatorBank);


%%%%%%%% Calculate reson bank energy %%%%%%%%%%%%%%%%%%%%%
  
disp('Calculating Resonator RMS...');

%calc reson energy for each frequency band and hilbert transform of
%each band
nBands = size(rp.data{rpCols.resonatorOut},1);
for iBand = 1:nBands
  thisBandResonatorOut = rp.data{rpCols.resonatorOut}(iBand,:,:);
  rpThisBand = init_rp_struct('resonatorOut',thisBandResonatorOut,'Fs',params.Fs,'resonatorFreqs',rp.data{rpCols.resonatorFreqs});
  rpChannelEnergy  = resonatorEnergy(rpThisBand,params.resonatorBank);
  rp.data{rpCols.resonatorBandEnergy}(iBand,:,:) = rpChannelEnergy;

  %rp.data{rpCols.complexResonatorOut}(iBand,:,:) = (hilbert(squeeze(thisBandResonatorOut)'))';
  mppThisBand = mean(rpChannelEnergy,2);
  rp.data{rpCols.mppByBand} =mppThisBand;
  bandPeakInfo{iBand} = get_resonPeakInfo(mppThisBand,rp.data{rpCols.resonatorFreqs},params.resonatorBank.peakFinder);
  
end

rp.data{rpCols.resonatorEnergy} = squeeze(mean(rp.data{rpCols.resonatorBandEnergy},1));
rp.data{rpCols.meanResonatorEnergy} = mean(rp.data{rpCols.resonatorEnergy},2);
rp.data{rpCols.stdResonatorEnergy} = std(rp.data{rpCols.resonatorEnergy},[],2);

peakInfo = get_resonPeakInfo(rp.data{rpCols.meanResonatorEnergy},rp.data{rpCols.resonatorFreqs},...
					   params.resonatorBank.peakFinder);
%peakInfo.bandPeakInfo = bandPeakInfo;
peakInfo = ensemble_tree2datastruct(peakInfo);
%temporary solution until ensemble_tree2datastruct is fixed
if(iscell(peakInfo.data{1}))
  for iTmp = 1:length(peakInfo.data)
    peakInfo.data{iTmp} = peakInfo.data{iTmp}{1};
  end
end
peakInfo.vars{end+1} = 'bandPeakInfo';
peakInfoCols = set_var_col_const(peakInfo.vars);
peakInfo.data{peakInfoCols.bandPeakInfo} = bandPeakInfo;
rp.data{rpCols.peakInfo} = peakInfo;
  
%%%%%%%% Post BTB Paper Modules and Calculations %%%%%%%%%%%%%%

%[rp.data{rpCols.mmpByFrame},rp.data{rpCols.stdmpByFrame},rp.data{rpCols.ratioBinsByFrame}, ...
%rp.data{rpCols.ratioEntropyByFrame},rp.data{rpCols.ratioTally}] = calc_mmpByFrame(rp.data{rpCols.resonatorEnergy},rp.data{rpCols.resonatorFreqs},params);

%[rp.data{rpCols.expectationError}] = rp_temporalExpectationError('rp',rp,'BList',BList,'AList',AList);

%params.vmExpect.scaleByDenv = 0;
%[rp.data{rpCols.vmExpect}] = rp_vonMisesExpectCalc(rp,params);
%[rp.data{rpCols.coherence}] = rp_phaseCoherenceByPeak(rp,params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outData = rp;

return %rhythm_profiler


function averageCBU = meanBandCBU(params)

aniParams = params.ani;
rpParams = params.rp;
firstCBU = aniParams.FirstCBU;
CBUStep  = aniParams.CBUStep;
numChannels = aniParams.NumOfChannels;
numChannelsSummed = rpParams.sumAdjacentBands;
origChannelsCBU = [firstCBU:CBUStep:(firstCBU+(numChannels-1)*CBUStep)];

for iSummed = 1:ceil(numChannels/numChannelsSummed)
  firstIdx = (iSummed-1)*numChannelsSummed + 1;
  lastIdx = min(iSummed*numChannelsSummed,length(origChannelsCBU));
  averageCBU(iSummed) = mean(origChannelsCBU(firstIdx:lastIdx));  
end

return


function outWvf = sumBands(inWvf,numToSum)

m0 = numToSum;
c0 = size(inWvf,1)/m0;

for c = 1:c0
  b = [ (c - 1)*m0+1 : c*m0 ]; 
  outWvf(c,:) = sum(inWvf(b,:),1);
end

return




function defParams = getDefaultParams(params)

if(exist('params','var') & isfield(params,'inDataType'))
  switch(params.inDataType)
    
   case 'ani'
    defParams = rp_paramGroups('param_group','reson_filterQSpacing_periodBasedDecay',...
			       'input_type','ani','gain_type','beta_distribution','perform_mulawCompression',0);
    
   case 'aud'
    defParams = rp_paramGroups('param_group','reson_filterQSpacing_periodBasedDecay','input_type','aud','gain_type','beta_distribution','perform_mulawCompression',0);
   case 'mat'
    defParams = rp_paramGroups('param_group','reson_filterQSpacing_periodBasedDecay','input_type','mat','gain_type','beta_distribution','perform_mulawCompression',0);
  
  end

else
  
   defParams = rp_paramGroups('param_group','reson_filterQSpacing_periodBasedDecay',...
			       'input_type','ani','gain_type','beta_distribution','perform_mulawCompression',0);
  
end
  
return
