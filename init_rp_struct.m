function rp = init_rp_struct(varargin)
% Initializes a rhythm profiler (rp) data struct
% 
% rp = init_rp_struct(varargin)
% 
% INPUTS:
% tag, value pairs that can be used to initialize variables in the
% struct
%
% valid variables are: params,meanChannelCBU,rms,rmsCompressed,Fs,denv,sumDenv,resonatorFreqs,
%	               resonatorOut,resonatorEnergy,resonatorBandEnergy,
%	               meanResonatorEnergy,stdResonatorEnergy,peakFreqs,maxFilt
%
%
% Copyright (c) 2007 The Regents of the University of California
% All Rights Reserved.
%
% Author: Stefan Tomic 10/18/2007

rp = ensemble_init_data_struct;
rp.type = 'rhythm_profile';
rp.vars = {'params','meanChannelCBU','rms','rmsCompressed','Fs','denv','sumDenv','resonatorFreqs', ...
	  'resonatorOut','resonatorEnergy','resonatorBandEnergy','mppByBand', ...
	   'meanResonatorEnergy','stdResonatorEnergy',...
	   'mppByFrame','stdmpByFrame','ratioBinsByFrame','ratioEntropyByFrame',...
	   'peakInfo','interpMPPStruct','ratioTally','onsetInfo','complexOutput','complexExpectAtOnsets','vmExpect','coherence'};
%peakFreqs is depreceated but left in the struct for
%compatibility. Use peakInfo (which contains another data struct
%with peak widths and  norm freq instead
rpCols = set_var_col_const(rp.vars);

for ivar = 1:2:length(varargin)
  
  tag = varargin{ivar};
  if(~ismember(tag,rp.vars))
    error(sprintf('%s is not a valid variable for type rp',tag));
  end
  rp.data{rpCols.(tag)} = varargin{ivar+1};
  
end
