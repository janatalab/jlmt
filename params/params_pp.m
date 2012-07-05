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
%       inDataType: specify the type of input for calc_pp (optional)
%       prev_steps: specify the analysis steps previous to calc_pp that
%                   should be encountered when using the given
%                   parameters (optional)
%
% Copyright (c) 2006-2012 The Regents of the University of California
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
    'Atten', ...
    'prev_steps', ...
    'inDataType'};

params = mkstruct(fields,varargin);

end % params_pp
