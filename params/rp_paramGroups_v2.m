function rp = rp_paramGroups_v2(varargin)
% Returns an rp parameteter structure
%
% Returns an rp parameter structure. The parameters are grouped so
% that commonly used parameters for various tasks do not need to
% be explicitly specified.
%
% INPUTS
% Accepts tag/value pairs, where possible tags are:
%   'prev_steps':   (optional)   a cell array of strings specifying the
%                                analysis steps previous to calc_rp that
%                                should be encountered when using the given
%                                parameters
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
% Changes:
%  TC, 11.13.2012. Added some default parameters for onset detection, line
%   142 onwards.
%  BH 22Dec2013 - Added rp.perform.calcComplexExpectAtOnsets as a default 'perform' flag

% Initialize an empty RP structure
rp = struct('Fs',[], ...
    'perform', struct, ...
    'timeLimSecs', [], ...
    'resonatorBank', struct, ...
    'onsetInfo', struct, ...
    'resonatorBand', [], ...
    'temporalExpect', struct, ...
    'vmExpect', struct, ...
    'mppByFrame', struct, ...
    'inDataType', '', ...
    'prev_steps', {});

for iarg = 1:2:nargin
  
  switch(varargin{iarg})  
   case 'param_group'
    param_group = varargin{iarg+1};
    
   case 'input_type'
    input_type = varargin{iarg+1};
    
   case 'gain_type'
    gain_type = varargin{iarg+1};
    
   case 'Fs'
    rp(1).Fs = varargin{iarg+1};
    
   case 'perform'
    rp(1).perform = varargin{iarg+1};
    
   case 'prev_steps'
    rp(1).prev_steps = varargin{iarg+1};
  end
  
end

if ~exist('rp') || isempty(rp) || ~isfield(rp,'Fs') || isempty(rp.Fs)
    rp(1).Fs = 100; 
end

%assign default 'perform' flags.
%these will be overridden based 
%on param_group and input_type
rp.perform.convertInputToMono = 0;
rp.perform.invertInputSig = 0;
rp.perform.calcRMS = 0;
rp.perform.calcDiff = 0;
rp.perform.calcHWR = 0;
rp.perform.sumChannels = 0;
rp.perform.calcResonBanks = 1;
rp.perform.calcResonRMS = 1;
rp.perform.calcAPS = 1;
rp.perform.calcMPP = 1;
rp.perform.calcPeakInfo = 1;
rp.perform.calcMPPByFrame = 1;
rp.perform.calcOnsetInfo = 0;
rp.perform.calcComplexOutput = 0;
rp.perform.calcComplexExpectAtOnsets = 0;
rp.perform.calcVonMises = 0;
rp.perform.calcPhaseCoherence = 0;
rp.perform.interpMPP = 0;
rp.perform.vmExpect.useInterpPeak = 1;
rp.perform.phaseCoherence.useInterpPeak = 1;
rp.perform.normalize.denv = 0;
rp.perform.complexOut.useOnsetInfo = 0;

%SETTINGS THAT ARE IN COMMON WITH ALL PARAM SETS
rp.timeLimSecs = [];
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
%rp.resonatorBank.interpMPPInterval = 0.0005;
rp.resonatorBank.interpMPPInterval = 0.01;
rp.onsetInfo.thresh = 0.05;
% Adding some default parameters for onset detection.
rp.onsetInfo.peak_height_thresh = 0.05;
rp.onsetInfo.peak_height_window = 0.05;
rp.onsetInfo.throwOutFirstOnset = 1;
rp.onsetInfo.Fs = rp.Fs;
rp.resonatorBand = 1;
% End of additions.
rp.onsetInfo.peakBaseThreshold = 0;
rp.temporalExpect.predictSec = 5.0;
rp.vmExpect.normalizeWindow = 4;

if(rp.perform.calcMPPByFrame)
  rp.mppByFrame.frameSize = 4.0;
  rp.mppByFrame.hopSize = 0.1;
end



%INPUT SIGNAL TYPE (ANI OR AUD)
switch(input_type)
  
 case 'ani'
  rp.inDataType = 'ani';
  rp.sumChannels.num = 8;
  rp.sumChannels.after = 'HWR';
  %throw out first onset for ANI (see rp_temporalExpectationError code)
  rp.temporalExpect.throwOutFirstOnset = 1; 
  rp.rmsFrameWidth = 0.050;
  rp.perform.calcRMS = 1;
  rp.perform.calcDiff = 1;
  rp.perform.calcHWR = 1;
  rp.perform.sumChannels = 1;
 case 'aud'
  rp.inDataType = 'aud';
  rp.perform.convertInputToMono = 1;
   rp.temporalExpect.throwOutFirstOnset = 0; 
 case 'mat'
  rp.inDataType = 'mat';
   rp.temporalExpect.throwOutFirstOnset = 0; 
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
  rp.resonatorBank.numPeriodDecay = 13;
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

% Add previous steps info
