function rp = rp_paramGroups(varargin)
% Returns an rp parameteter structure
%
% Returns an rp parameter structure. The parameters are grouped so
% that commonly used parameters for various tasks do not need to
% be explicitly specified.
%
% INPUTS
% Accepts tag/value pairs, where possible tags are:
%   'input_type':                'ani','aud',or 'mat'. 'aud' and 'mat' set the
%                                appropriate parameters to bypass ANI calculation
% 
%   'gain_type':                 'normalized_response', 'constant',
%                                'weighted_partialHanning', or 'beta_distribution'
%                                'normalized_response' sets the reson filter gain
%                                so that the peak response of all filters is 0 dB.
%                                'beta_distribution' is the gain specified
%                                in the BTB paper.
%                                'constant' specifies that a constant gain is used
%                                across all resons
%    'perform_mulawCompression': 0 or 1. Whether or not to perform ...
%	                         muLawCompression (see Klapuri,2006). 
%                                This did not improve our results 
%                                so it is normally turned off.
%
%    'Fs':                       sampling rate used for RP model.
%    'param_group':              mechanism for grouping parameters
%                                to simplify switching between sets of parameters.
%                                The BTB paper uses 'reson_filterQSpacing_periodBasedDecay'.
%                                Descriptions:
%                                 'reson_filterQSpacing_periodBasedDecay'
%                                    Sets all reson filters with
%                                    constant Q-factor. Spaces
%                                    filters so the magnitude
%                                    responses meet at their
%                                    half-power bandwidth.
%                                    Sets all decay rates based on
%                                    the number of periods (13 periods currently used).
%                                 'reson_linearSpacing_halfPowerBandwidth'
%                                    Sets all reson filters with
%                                    the same bandwidth. The reson
%                                    filters are spaced so that
%                                    their  magnitude responses
%                                    meet at their half power bandwidth.
%                                 'reson_linearSpacing_constantBandwidth'
%                                    Sets all reson filters with
%                                    the same bandwidth. The
%                                    filters are spaced so that
%                                    their center frequencies are a
%                                    constant distance apart
%                                 'reson_linearSpacing_periodBasedDecay'
%                                    Sets all reson filters with
%                                    constant Q-factor (therefore
%                                    bandwidth increases with
%                                    center frequency). Spaces
%                                    filters by a constant value.
%                                 'reson_log2spacing_constBW'
%                                    Sets reson filter bandwidths
%                                    to a constant value. Uses a
%                                    logarithmic spacing
%                                 'comb'
%                                    Uses comb filters instead of
%                                    reson filters. Uses linear
%                                    spacing of filters and
%                                    constant half-energy time.
%
%
% Copyright (c) 2007-2013 The Regents of the University of California
% All Rights Reserved.
%
% Author: Stefan Tomic 10/18/2007


for iarg = 1:2:nargin
  
  switch(varargin{iarg})  
   case 'param_group'
    param_group = varargin{iarg+1};
    
   case 'input_type'
    input_type = varargin{iarg+1};
    
   case 'gain_type'
    gain_type = varargin{iarg+1};
    
   case 'perform_mulawCompression'
    if(varargin{iarg+1} == 1)
      rp.mulawCompression.perform = 1;
      rp.mulawCompression.mu = 100;
    else
      rp.mulawCompression.perform = 0;
    end
    
   case 'Fs'
    rp.Fs = varargin{iarg+1};
    
  end
  
end

%SETTINGS THAT ARE IN COMMON WITH ALL PARAM SETS
rp.perform = {};
rp.calcDiff = 1;
rp.invertInputSig = 0;
rp.sumBandsAfter = 'HWR';
rp.resonatorBank.freqLimits = [0.25 10];
rp.resonatorBank.energyWindowSecs = 2.0; 
rp.resonatorBank.periodBasedEnergy = 0;
rp.resonatorBank.peakFinder.thresh = 0;
rp.resonatorBank.peakFinder.perform = {};
rp.resonatorBank.peakFinder.minPeakHeight = 0.05;
rp.resonatorBank.peakFinder.errorTolerance = 0.05;
rp.resonatorBank.peakFinder.ratioVals = ...
    [ 1/8 1/6 1/5 1/4 1/3 3/8 2/5 1/2 3/5 5/8 2/3 3/4 4/5 5/6 7/8 1 5/4 4/3 ...
      3/2 5/3 2 5/2 3:8];
rp.resonatorBank.peakFinder.ratioLabels = ...
    { '1/8' '1/6' '1/5' '1/4' '1/3' '3/8' '2/5' '1/2' '3/5' '5/8' '2/3' '3/4' ...
      '4/5' '5/6' '7/8' '1' '5/4' '4/3' '3/2' '5/3' '2' '5/2' '3' '4' '5' '6' '7' '8' };
rp.mmpByFrame.frameSize = 4.0;
rp.mmpByFrame.hopSize = 0.1;

%INPUT SIGNAL TYPE (ANI OR AUD)
switch(input_type)
  
 case 'ani'
  rp.inDataType = 'ani';
  rp.sumAdjacentBands = 8;
  rp.ipemRMSMethod = 0;
  try rp.Fs; catch  rp.Fs = 100; end
  rp.rmsFrameWidth = 0.050;
  rp.bypassOnsetDetect = 0;
  rp.convertInputToMono = 0;
 case 'aud'
  rp.inDataType = 'aud';
  rp.sumAdjacentBands = 1;
  rp.ipemRMSMethod = 0;
  try rp.Fs; catch rp.Fs = 100; end
  rp.bypassOnsetDetect = 1;
  rp.convertInputToMono = 1;
 case 'mat'
  rp.inDataType = 'mat';
  rp.sumAdjacentBands = 1;
  rp.ipemRMSMethod = 0;
  try rp.Fs; catch rp.Fs = 100; end
  rp.bypassOnsetDetect = 1;
  rp.convertInputToMono = 0;
  
end

%GAIN TYPE
switch(gain_type)
 case 'constant'
  rp.resonatorBank.gainType = 'constant';
  rp.resonatorBank.gain = 0.4;
 case 'weighted_partialHanning'
  rp.resonatorBank.gainType = 'weighted';
  rp.resonatorBank.gainFunction.method = 'partialHanning';
  rp.resonatorBank.gainFunction.maxGain = 0.4;
  rp.resonatorBank.gainFunction.minGain = 0.1;
  rp.resonatorBank.gainFunction.maxFreq = 1.66;
 case 'normalized_response'
  rp.resonatorBank.gainType = 'normalized_response';
 case 'beta_distribution'
  rp.resonatorBank.gainType = 'weighted';
  rp.resonatorBank.gainFunction.method = 'betaDistribution';
  rp.resonatorBank.gainFunction.maxGain = 0.4;
  rp.resonatorBank.gainFunction.minGain = 0.1;
  rp.resonatorBank.gainFunction.betaDistribution.alpha = 2;
  rp.resonatorBank.gainFunction.betaDistribution.beta = 5;
 case 'prenorm_beta_distribution'
  rp.resonatorBank.gainType = 'prenorm_weighted';
  rp.resonatorBank.gainFunction.method = 'betaDistribution';
  rp.resonatorBank.gainFunction.maxGain = 0.4;
  rp.resonatorBank.gainFunction.minGain = 0.1;
  rp.resonatorBank.gainFunction.betaDistribution.alpha = 2;
  rp.resonatorBank.gainFunction.betaDistribution.beta = 5;

end

rp.resonatorBank.Fs = rp.Fs;

%PROCESSING PARAM GROUPS
switch(param_group)
  
 case 'reson_filterQSpacing_periodBasedDecay';
  rp.resonatorBank.method = 'reson_z';
  rp.resonatorBank.resonFreqSpacing = 'filterQ';
  rp.resonatorBank.numPeriodDecay = 13;
  rp.resonatorBank.bw_measurement = 'periodBasedDecay';
  rp.resonatorBank.spacingFactor = 0.5;
  
 case 'reson_linearSpacing_halfPowerBandwidth'
  rp.resonatorBank.method = 'reson_z';
  rp.resonatorBank.resonFreqSpacing = 'linear';
  rp.resonatorBank.spacingHz = .06;
  rp.resonatorBank.bw_measurement = 'halfPowerBandwidth';
  
 case 'reson_linearSpacing_constantBandwidth'
  rp.resonatorBank.method = 'reson_z';
  rp.resonatorBank.resonFreqSpacing = 'linear';
  rp.resonatorBank.spacingHz = .06;
  rp.resonatorBank.bw_measurement = 'constant';
  rp.resonatorBank.bandwidth = 0.12;
  
 case 'reson_linearSpacing_periodBasedDecay';
  rp.resonatorBank.method = 'reson_z';
  rp.resonatorBank.resonFreqSpacing = 'linear';
  rp.resonatorBank.spacingHz = .06;
  rp.resonatorBank.bw_measurement = 'periodBasedDecay';

 case 'reson_log2spacing_constBW'  
  rp.resonatorBank.method = 'reson_z';
  rp.resonatorBank.resonFreqSpacing = 'log2';
  rp.resonatorBank.numFilters = 90;
  rp.resonatorBank.bw_measurement = 'constant';
  rp.resonatorBank.bandwidth = 0.12;
 
 case 'comb'
  rp.resonatorBank.method = 'comb';
  rp.resonatorBank.resonFreqSpacing = 'linear';
  rp.resonatorBank.halfEnergyTimeSecs = 1.5;
  rp.resonatorBank.numFilters = 99;
end
