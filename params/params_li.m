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
%          HalfDecayTimes: 2
% % % %         HalfDecayChords: []
% % % %    HalfDecayToneCenters: []
%              inDataType: specify the type of input for calc_li (optional)
%              prev_steps: specify the analysis steps previous to calc_li
%                          that should be encountered when using the given
%                          parameters (optional)

% Copyright (c) 2006-2013 The Regents of the University of California
% All Rights Reserved.
%
% 12/4/06 Petr Janata
% 2010.02.18 FB - adapted to use HalfDecayTimes
% 2011.01.28 FB - added field for metrics
% 2012.07.05 FB - added inDataType, prev_steps
% 2012.11.13 TC. Altered prev_steps to empty cell array, not double.
% 2019.07.29 PJ - Fixed handling of varargin 

fields = {...
  'PeriodicityPitch', ...
  'SnapShot', ...
  'HalfDecayTimes',...
  'Enlargement', ...
  'PlotFlag', ...
  'prev_steps', ...
  'inDataType'};

% Convert params structure to argument list if necessary
if strcmp(varargin{1},'params')
  args = [fieldnames(varargin{2}) struct2cell(varargin{2})]';
else
  args = varargin;
end
params = mkstruct(fields,args(:));

% set defaults if not otherwise specified
def.PeriodicityPitch = [];
def.SnapShot = [];
def.HalfDecayTimes = [0.2 2];
def.Enlargement = [];
def.PlotFlag = [];
def.prev_steps = {};
def.inDataType = [];

for ifld = 1:length(fields)
  if isempty(params.(fields{ifld}))
    params.(fields{ifld}) = def.(fields{ifld});
  end
end
