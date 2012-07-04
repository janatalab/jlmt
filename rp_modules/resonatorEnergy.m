function rpEnergy = resonatorEnergy(rp,params)
% Calculates the RMS of a resonator bank output (periodicity surface)
%
% rpEnergy = resonatorEnergy(rp,params)
%
%
% Accepts an rp struct, calculates the RMS of the resonator
% outputs, and returns the result in rpEnergy. The RMS is
% calculated for each frame with a hop size of one sample.
%
% INPUT
% rp - an rp struct that contains the resonator output (see init_rp.m)
% params.energyWindowSecs  - the frame size for the RMS
% calculation. 
% params.periodBasedEnergy - 0 or 1. Whether or not to calculate
% the RMS with a framesize equal to the period of the reson
% filter. If true, params.energyWindowSecs is ignored. 
%
% OUTPUT
% rpEnergy - result of RMS calculation
%
% Copyright (c) 2006 The Regents of the University of California
% All Rights Reserved
%
% Author(s):
% Stefan Tomic 12/06


rpCols = set_var_col_const(rp.vars);
Fs = rp.data{rpCols.Fs};
energyWindowSamps = round(params.energyWindowSecs*Fs);

if(ndims(rp.data{rpCols.resonatorOut}) == 3)
  summedFilters = squeeze(sum(rp.data{rpCols.resonatorOut},1));
else
  summedFilters = rp.data{rpCols.resonatorOut};
end

numResonators = size(summedFilters,1);
numTimeSamps = size(summedFilters,2);


rpEnergy = zeros(numResonators,numTimeSamps);    

if(~params.periodBasedEnergy)
  for(iSamp = 1:numTimeSamps)
    %start samp is the current sample - the energy window size
    %if we haven't progressed a full energy window then take sample
    %1 as the start samp
    startSamp =  max(iSamp - energyWindowSamps + 1, 1);
    endSamp = iSamp;
    windowEnergy = sqrt(sum(summedFilters(:,startSamp:endSamp).^2,2)./energyWindowSamps);
    rpEnergy(:,iSamp) = windowEnergy;
  end
else
  for(iReson = 1:numResonators)
    %energyWindowSamps = the period duration in samples
    energyWindowSamps = round(1./rp.data{rpCols.resonatorFreqs}(iReson).*Fs);
    for(iSamp = 1:numTimeSamps)
      startSamp = max(iSamp - energyWindowSamps + 1, 1);
      endSamp = iSamp;
      windowEnergy = sqrt(sum(summedFilters(iReson,startSamp:endSamp).^2,2)./energyWindowSamps);
      rpEnergy(iReson,iSamp) = windowEnergy;
    end
  end
end

