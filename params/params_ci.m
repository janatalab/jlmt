function params = params_ci(varargin)
% Returns a parameter structure for IPEM Toolbox Contextuality Index calculations
%
% params = params_ci(varargin);
% 
% Input parameters for IPEM Toolbox Contextuality Index calculations
%
%        PeriodicityPitch: []
%              SampleFreq: []
%                 Periods: []
%                SnapShot: []
%             Enlargement: []
%                PlotFlag: []
%          HalfDecayTimes: []
% % % %         HalfDecayChords: []
% % % %    HalfDecayToneCenters: []
%
% Copyright (c) 2006-2013 The Regents of the University of California
% All Rights Reserved.
%
% 12/4/06 Petr Janata
% 2010.02.18 FB - adapted to use HalfDecayTimes
% 2011.01.28 FB - added field for metrics

fields = {...
      'PeriodicityPitch', ...
      'SnapShot', ...
      'HalfDecayTimes',...
      'Enlargement', ...
      'PlotFlag',...
      'metrics'};

params = mkstruct(fields,varargin);

if isempty(params.metrics)
  params.metrics = struct();
end

end % params_ci
