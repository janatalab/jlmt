function params = params_tonreg_fmri(varargin)
% Returns a paramneter struct for use with calc_tonreg_fmri
%
%   params = params_tonreg_fmri(varargin);
%
% This function accepts a variable number of tag/value pairs. If no
% argument is passed in, then a parameter structure with empty values
% is returned.
%
% Parameter structure for controlling fmri tonality tracking regressor generation
%      scanner: TR (time to collect one volume) and dt (micro-timing slices
%               within each TR)
%      use_sig: string name of integration time constant [created by
%               calc_li('calc_names',<integration time constant vector>)]
%               to use when extracting torus activation timecourses
%   inDataType: specify the type of input for calc_pp (optional)
%   prev_steps: specify the analysis steps previous to calc_pp that should
%               be encountered when using the given parameters (optional)
% 
% Initialization values must be passed in as parameter/value pairs
%
% Copyright (c) 2012 The Regents of the University of California
% All Rights Reserved.
%
% Author(s):
% Fred Barrett - 2012.07.05

%this is the parameter list
%if you need to add more parameters, add them here
fields = {'scanner','use_sig','inDataType','prev_steps'};

params = mkstruct(fields,varargin);

% set defaults if not otherwise specified
def.scanner = struct('TR',2,'dt',16); % Janata Lab defaults: 2s TR, 1/16 dt
def.use_sig = calc_li('calc_names',[2]); % timescale used in Janata 2009 Cerebral Cortex
def.inDataType = '';
def.prev_steps = [];

for ifld = 1:length(fields)
  if isempty(params.(fields{ifld}))
    params.(fields{ifld}) = def.(fields{ifld});
  end
end
