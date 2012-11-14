function outData = calc_rp(inData, varargin)
% This is the entry point for the BTB algorithm
%
%   outData = calc_rp(inData,params)
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
%                           params.perform: 'convertInputToMono'
%                                              Will convert the input signal to mono
%                                              by taking the mean signal of all the channels.
%                                           'calcRMS'
%                                           'calcDiffHWR'
%                                           'calcResonBanks'
%                                           'calcResonRMS'
%                                           'calcAPS'
%                                           'calcMPP'
%                                           'calcPeakInfo'
%                                           
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
% Copyright (c) 2006-2012 The Regents of the University of California
% All Rights Reserved.
%
% Author: Stefan Tomic 12/2006
% 6/16/09 - Stefan Tomic. restructured params to use a cell array of strings to
%           determine which calculations to perform. Cleaned up handling of
%           post JASA paper calculations, and set up handling of rp structs
%           as inputs.
% 30Oct2012 PJ - added support for varargin and handling of default
%                parameters


if nargin > 1 && isstruct(varargin{1})
  params = varargin{1};
end

%params for converting tree to ensemble data struct
etParams.encapsulate_in_cells = 1;
  
%if string 'getDefaultParams' is passed in, then just return the
%default params struct
if ischar(inData) && strcmp(inData,'getDefaultParams')
  if exist('params','var')
    outData = getDefaultParams('params',params);
  else
    outData = getDefaultParams(varargin{:});
  end
  return
end

% if perform substruct wasn't specified, throw an error.
if(~isfield(params,'perform'))
  error('perform substructure was not specified in params. Perhaps this is an outdated params struct.');
end

%see if inData is an ensemble data struct. If not, throw an error
if(~isfield(inData,'type'))
  error(['input data (inData) must be an ensemble data struct of' ...
	 ' type ''aud'',''ani'', or ''rp''']);
end

%deal with different input types here.
switch(inData.type)
  
 case 'aud'
  %if the input was actually a mat file, this should have been
  %converted to a struct of type 'aud' at this point
  
  sigStruct = inData;
  sigStructCols = set_var_col_const(sigStruct.vars);
  inSig.wvf = sigStruct.data{sigStructCols.wvf}';
  inSig.Fs = sigStruct.data{sigStructCols.Fs};

  if(size(inSig.wvf,1) > 1 && params.perform.convertInputToMono)
    inSig.wvf = mean(inSig.wvf,1);
  end
  
  origFs = sigStruct.data{sigStructCols.Fs};
  
  if(params.Fs ~= inSig.Fs)
    warning('Resampling input signal. This may result in undesirable artifacts');
    inSig.wvf = resample(inSig.wvf,params.Fs,inSig.Fs);
    inSig.Fs = params.Fs;
  end

  %initialize rp struct
  rp = init_rp_struct;
  rpCols = set_var_col_const(rp.vars);
  saveParams.rp = params;
  rp.data{rpCols.params} = saveParams;
  rp.data{rpCols.Fs} = params.Fs;
  
 case 'ani'
  
  aniStruct = inData;
  aniStructCols = set_var_col_const(aniStruct.vars);
  inSig.wvf = aniStruct.data{aniStructCols.sig};
  inSig.Fs  = aniStruct.data{aniStructCols.Fs};
  
  %initialize rp struct
  rp = init_rp_struct;
  rpCols = set_var_col_const(rp.vars);
  saveParams.ani = aniStruct.data{aniStructCols.params}.ani;
  saveParams.rp = params;
  rp.data{rpCols.params} = saveParams;
  rp.data{rpCols.Fs} = params.Fs;
  
 case 'rp'
  
  rp = inData;
  
  %make sure that the params of the rp are contained within the
  %passed in params
  rpParams = rp.data{rpCols.params};
  [paramsAreContained] = compare_structs(rpParams,params,'substruct',1,'values',1);
  if(~paramsAreContained)
    error(['The passed in params does not comprise a superset of the params' ...
	   ' of the input RP']);
  end
  
  %populate rp struct with missing vars
  empty_rp = init_rp_struct;
  missingVars = setdiff(empty_rp.vars,rp.vars);
  rp.vars = {rp.vars{:} missingVars{:}};
  rpCols = set_var_col_const(rp.vars);
  
  
end

%invertInputSig may be useful in some cases where the signal is not
%symmetrical about DC, and we are interested in the periodicities
%of troughs the signal
if(params.perform.invertInputSig)
  inSig.wvf = invert_sig(inSig.wvf);
end

if(params.perform.calcRMS)
  %%%%%% Calc RMS %%%%%%%%%
  
  disp('Calculating RMS...');
  
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

    
  Flp = 10; %cutoff freq in Hz
  Wn = Flp*2/rp.data{rpCols.Fs};
  [b,a] = butter(12,Wn);
  nbands = size(rp.data{rpCols.rms},1);
  for iband = 1:nbands
    rp.data{rpCols.rms}(iband,:) = filtfilt(b,a,rp.data{rpCols.rms}(iband,:));
  end  
    
  if(params.perform.sumChannels && strcmp(params.sumChannels.after,'RMS'))
    rp.data{rpCols.rms} = sumBands(rp.data{rpCols.rms},params.sumChannels.num);
    rp.data{rpCols.meanChannelCBU} = meanBandCBU(rp.data{rpCols.params});
  end
  
end % perform calcRMS


if(params.perform.calcDiff)

  %%%%%%%% Calc Diff %%%%%%%%
  
  if((length(rp.data) < rpCols.rms) || isempty(rp.data{rpCols.rms}))
    error('Missing required rms data for calcDiff');
  end
  
  disp('Calculating difference of the RMS...');
  
  denv = diff(rp.data{rpCols.rms},1,2); %simple difference		
  denv = [zeros(size(denv,1),1) denv];
  
end % perform calcDiff
    
if(params.perform.calcHWR)
  denv(find(denv < 0)) = 0; %half-wave rectify the difference 

end % perform calcHWR

%if calcDiff or calcHWR were performed, assign the result to denv
%column here. Note that denv variable is BEFORE normalizing or summing across bands
if(params.perform.calcDiff || params.perform.calcHWR)
  rp.data{rpCols.denv} = denv;
end

%summing across adjacent bands
%result is stored in sumDenv
if(params.perform.sumChannels && strcmp(params.sumChannels.after,'HWR'))
  denv = sumBands(denv,params.sumChannels.num);
  rp.data{rpCols.meanChannelCBU} = meanBandCBU(rp.data{rpCols.params});
  rp.data{rpCols.sumDenv} = denv;
end

if(params.perform.normalize.denv)
  nBands = size(denv,1);
  nSamps = size(denv,2);
  %set the normalization window size to the period of the lowest reson filter freq.
  %we need to determine this by accessing the params, since the resonators have
  %not been calculated yet.
  minFreq = min(rp.data{rpCols.params}.rp.resonatorBank.freqLimits);
  winSize = (1 ./ minFreq) .* rp.data{rpCols.Fs};
  for iBand = 1:nBands
    for iSamp = 1:nSamps
      winStart = max(1,iSamp-winSize);
      window = [winStart:iSamp];
      maxInWindow = max(abs(denv(iBand,window)));
      if(maxInWindow ~= 0)
	normDenv(iBand,iSamp) = denv(iBand,iSamp) ./ maxInWindow;
      else
	normDenv(iBand,iSamp) = 0;
      end
    end
  end

  denv = normDenv;
  clear normDenv;
  
  %assign denv back to either sumDenv or denv variable in the rp structure
  %the choice of which one will depend on whether summation was performed or
  %not.
  if(params.perform.sumChannels && strcmp(params.sumChannels.after,'HWR'))
    rp.data{rpCols.sumDenv} = denv;
  else
    rp.data{rpCols.denv} = denv;
  end
end

%now determine what variable to use as input to resonators
if(params.perform.calcDiff || params.perform.calcHWR)
  inputToResonators = denv;
elseif(params.perform.calcRMS)
  inputToResonators = rp.data{rpCols.rms};
else
  %if none of RMS, diff, HWR were performed, the input signal is
  %passed directly to the resonators
  inputToResonators = inSig.wvf;
end
  
if(params.perform.calcResonBanks)
  %%%%%%%% Calc Reson Banks %%%%%%%%%%%%%%
  disp('Calculating Reson Filterbanks for each frequency band...');

  %resonator bank Fs = rms signal Fs
  params.resonatorBank.Fs = params.Fs;

  hfunc_resonatorBanks = resonatorBanks(params.resonatorBank.method);
  [rp.data{rpCols.resonatorOut},rp.data{rpCols.resonatorFreqs}] = ...
      hfunc_resonatorBanks(inputToResonators,params.resonatorBank);

end

bandPeakInfo = {};
if(params.perform.calcResonRMS)
  %%%%%%%% Calculate reson bank RMS %%%%%%%%%%%%%%%%%%%%%
  
  disp('Calculating Resonator RMS...');

  %calc reson energy for each frequency band
  nBands = size(rp.data{rpCols.resonatorOut},1);
  for iBand = 1:nBands
    rpThisBand = init_rp_struct('resonatorOut',rp.data{rpCols.resonatorOut}(iBand,:,:),'Fs',params.Fs,'resonatorFreqs',rp.data{rpCols.resonatorFreqs});
    rpChannelEnergy  = resonatorEnergy(rpThisBand,params.resonatorBank);
    rp.data{rpCols.resonatorBandEnergy}(iBand,:,:) = rpChannelEnergy;
    mppThisBand = mean(rpChannelEnergy,2);
    rp.data{rpCols.mppByBand}{iBand,1} = mppThisBand;
    bpInfo(iBand,1) = get_resonPeakInfo(mppThisBand,rp.data{rpCols.resonatorFreqs},params.resonatorBank.peakFinder);
    
  end  

  bandPeakInfo = ensemble_tree2datastruct(bpInfo,etParams);
end

if(params.perform.calcAPS)
  if(length(rp.data) < rpCols.resonatorBandEnergy)
    error('required data ''resonatorBandEnergy'' must be calcuated for calcAPS');
  end
  rp.data{rpCols.resonatorEnergy} = squeeze(mean(rp.data{rpCols.resonatorBandEnergy},1));  
end

if(params.perform.calcMPP)
  if(length(rp.data) < rpCols.resonatorEnergy)
    error('required data ''resonatorEnergy'' must be calculated for calcMPP');
  end
  rp.data{rpCols.meanResonatorEnergy} = mean(rp.data{rpCols.resonatorEnergy},2);
  rp.data{rpCols.stdResonatorEnergy} = std(rp.data{rpCols.resonatorEnergy},[],2);
end
  
if(params.perform.calcPeakInfo)  
  if(length(rp.data) < rpCols.meanResonatorEnergy)
    error('required data ''meanResonatorEnergy'' must be calculated for calcPeakInfo');
  end
  
  peakInfo = get_resonPeakInfo(rp.data{rpCols.meanResonatorEnergy},rp.data{rpCols.resonatorFreqs},...
					     params.resonatorBank.peakFinder);
  peakInfo.bandPeakInfo = bandPeakInfo;
  peakInfo = ensemble_tree2datastruct(peakInfo);
 
  rp.data{rpCols.peakInfo} = peakInfo;
end

if(params.perform.interpMPP)
  if( (length(rp.data) < rpCols.meanResonatorEnergy) || ...
      isempty(rp.data{rpCols.meanResonatorEnergy})) 
     error('required data ''meanResonatorEnergy'' must be calculated for calcPeakInfo');
  end

  resonFreqs = rp.data{rpCols.resonatorFreqs};
  interpMPPInterval = params.resonatorBank.interpMPPInterval;
  interpFreqs = resonFreqs(1):interpMPPInterval:resonFreqs(end);
  
  for iBand = 1:nBands
    band_pp(iBand,1) = spline(rp.data{rpCols.resonatorFreqs},rp.data{rpCols.mppByBand}{iBand});
    bandInterpMPP{iBand,1} = ppval(band_pp(iBand),interpFreqs);  
    bandInterpPeakInfo(iBand,1) = get_resonPeakInfo(bandInterpMPP{iBand},interpFreqs,params.resonatorBank.peakFinder);
  
  end
    
  bandInterpPeakInfo = ensemble_tree2datastruct(bandInterpPeakInfo,etParams);
  
  interpMPP_pp = spline(rp.data{rpCols.resonatorFreqs},rp.data{rpCols.meanResonatorEnergy});

  interpMPPStruct = ensemble_init_data_struct;
  interpMPPStruct.vars = {'pp','freqs','interpMPP','peakInfo','band_pp','bandInterpMPP','bandPeakInfo'};
  imsCols = set_var_col_const(interpMPPStruct.vars);

  interpMPPStruct.data{imsCols.freqs} = interpFreqs;
  interpMPPStruct.data{imsCols.pp} = interpMPP_pp;
  interpMPP = ppval(interpMPP_pp,interpFreqs);
 
  interpMPPStruct.data{imsCols.interpMPP} = interpMPP;
  
  interpPeakInfo = get_resonPeakInfo(interpMPP,interpFreqs,params.resonatorBank.peakFinder);
  interpPeakInfo = ensemble_tree2datastruct(interpPeakInfo);

  interpMPPStruct.data{imsCols.peakInfo} = interpPeakInfo;
  interpMPPStruct.data{imsCols.band_pp} = band_pp;
  interpMPPStruct.data{imsCols.bandInterpMPP} = bandInterpMPP;
  interpMPPStruct.data{imsCols.bandPeakInfo} = bandInterpPeakInfo;
  rp.data{rpCols.interpMPPStruct} = interpMPPStruct;
  
end
	  
%%%%%%%% Post JASA 2008 Modules and Calculations %%%%%%%%%%%%%%

if(params.perform.calcMPPByFrame)
  [rp.data{rpCols.mppByFrame},rp.data{rpCols.stdmpByFrame},rp.data{rpCols.ratioBinsByFrame}, ...
   rp.data{rpCols.ratioEntropyByFrame},rp.data{rpCols.ratioTally}] = rp_mppByFrame(rp.data{rpCols.resonatorEnergy},rp.data{rpCols.resonatorFreqs},params);
end

if(params.perform.calcOnsetInfo)
  [rp.data{rpCols.onsetInfo}] = rp_onsetInfo('inputSig',inputToResonators,'rp',rp);
end

if(params.perform.calcComplexOutput)
  [rp.data{rpCols.complexExpectAtOnsets}] = rp_complexExpectAtOnsets(rp,params);
end

if(params.perform.calcVonMises)
  [rp.data{rpCols.vmExpect}] = rp_vonMisesExpectCalc('rp',rp,'inputSig',inputToResonators);
end

if(params.perform.calcPhaseCoherence)
  [rp.data{rpCols.coherence}] = rp_phaseCoherenceByPeak(rp);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outData = rp;

return %rhythm_profiler


function averageCBU = meanBandCBU(params)
% obtains the mean CBU of a series of adjacent channels
% (useful when wanting to label mean CBU of summed channels)
aniParams = params.ani;
rpParams = params.rp;
firstCBU = aniParams.FirstCBU;
CBUStep  = aniParams.CBUStep;
numChannels = aniParams.NumOfChannels;
numChannelsSummed = rpParams.sumChannels.num;
origChannelsCBU = [firstCBU:CBUStep:(firstCBU+(numChannels-1)*CBUStep)];

for iSummed = 1:ceil(numChannels/numChannelsSummed)
  firstIdx = (iSummed-1)*numChannelsSummed + 1;
  lastIdx = min(iSummed*numChannelsSummed,length(origChannelsCBU));
  averageCBU(iSummed) = mean(origChannelsCBU(firstIdx:lastIdx));  
end

return


function outWvf = sumBands(inWvf,numToSum)
%sums adjacent channels

m0 = numToSum;
c0 = size(inWvf,1)/m0;

for c = 1:c0
  b = [ (c - 1)*m0+1 : c*m0 ]; 
  outWvf(c,:) = sum(inWvf(b,:),1);
end

return




function defParams = getDefaultParams(varargin)
% returns a default parameter struct
% this can be called from external functions by passing 
% 'getDefaultParams' to rhythm_profiler

inDataType = '';
for iarg = 1:2:nargin
  switch varargin{iarg}
      case 'inDataType'
          inDataType = varargin{iarg+1};
    case 'params'
      params = varargin{iarg+1};
    case 'prev_steps'
      prev_steps = varargin{iarg+1};
  end
end


if ~exist('inDataType') && exist('params','var') && isfield(params,'inDataType')
    inDataType = params.inDataType;
end

if ~isempty(inDataType)
  switch inDataType
    
   case 'ani'
    defParams = rp_paramGroups_v2('param_group','reson_filterQSpacing_periodBasedDecay',...
			       'input_type','ani','gain_type','beta_distribution','prev_steps',prev_steps);
    
   case 'aud'
    defParams = rp_paramGroups_v2('param_group','reson_filterQSpacing_periodBasedDecay', ...
        'input_type','aud', ...
        'gain_type','beta_distribution', ...
        'prev_steps', {});
   case 'mat'
    defParams = rp_paramGroups_v2('param_group','reson_filterQSpacing_periodBasedDecay',...
        'input_type','mat',...
        'gain_type','beta_distribution', ...
        'prev_steps', {});
  
  end

else
  
   defParams = rp_paramGroups_v2('param_group','reson_filterQSpacing_periodBasedDecay',...
			       'input_type','ani','gain_type','beta_distribution','prev_steps',prev_steps);
  
end
  
return
