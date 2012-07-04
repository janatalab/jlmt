function params = params_pp(varargin)
% Returns a parameter structure for IPEM Toolbox Periodicity Pitch calculations
% 
% params = params_pp(varargin);
% 
% input arguments are optional and are in the form of tag/value
% pairs.
%
% Input parameters for IPEM Toolbox Periodicity Pitch calculations
%
%           Matrix: []
%       SampleFreq: []
%     LowFrequency: []
%       FrameWidth: []
%    FrameStepSize: []
%         PlotFlag: []
%            Atten: [] - handle to weighting function to apply over periodicity dimension
%
% Copyright (c) 2006 The Regents of the University of California
% All Rights Reserved.
%
%
% Author:
% 12/4/06 Petr Janata

fields = {...
    'LowFrequency', ...
    'FrameWidth', ...
    'FrameStepSize', ...
    'PlotFlag', ...
    'Atten'};

params = mkstruct(fields,varargin);

end % params_pp
