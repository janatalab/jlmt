function params = params_li(varargin)
% Returns a parameter structure for IPEM Toolbox Leaky Integration calculations
%
% params = params_li(varargin);
% 
% Input parameters for IPEM Toolbox Leaky Integration calculations
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
% Copyright (c) 2006 The Regents of the University of California
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
      'PlotFlag'};

params = mkstruct(fields,varargin);

end % params_li
