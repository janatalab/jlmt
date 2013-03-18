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
%       FrameWidth: 0.0381
%    FrameStepSize: 0.0381
%         PlotFlag: 0
%            Atten: 'ipem_squash_hf' - handle to weighting function to apply over periodicity dimension
%       inDataType: specify the type of input for calc_pp (optional)
%       prev_steps: specify the analysis steps previous to calc_pp that
%                   should be encountered when using the given
%                   parameters (optional)
%
% Copyright (c) 2006-2013 The Regents of the University of California
% All Rights Reserved.
%
%
% Author:
% 12/4/06 Petr Janata
% Changes.
%  2012.11.13 TC. Altered prev_steps to empty cell array, not double.

fields = {...
    'LowFrequency', ...
    'FrameWidth', ...
    'FrameStepSize', ...
    'PlotFlag', ...
    'Atten', ...
    'prev_steps', ...
    'inDataType'};

params = mkstruct(fields,varargin);

def.LowFrequency = [];
def.FrameWidth = 0.0381;
def.FrameStepSize = 0.0381;
def.PlotFlag = 0;
def.Atten = 'ipem_squash_hf';
def.prev_steps = {};
def.inDataType = '';

for ifld = 1:length(fields)
  if isempty(params.(fields{ifld}))
    params.(fields{ifld}) = def.(fields{ifld});
  end
end
