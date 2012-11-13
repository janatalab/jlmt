function params = params_pc(varargin)
% Returns a parameter structure for projections of pp or ci images to pitch class space
% 
% params = params_pc(varargin);
% 
% input arguments are optional and are in the form of tag/value
% pairs.
%
% Input parameters for pitch class projection:
%
%       wmtx.fname: []
%   HalfDecayTimes: []
%       ci_siglist: []
%             norm: []
%       inDataType: specify the type of input for calc_pc (optional)
%       prev_steps: specify the analysis steps previous to calc_pc that
%                   should be encountered when using the given
%                   parameters (optional)
%
% Copyright (c) 2011-2012 The Regents of the University of California
% All Rights Reserved.
%
%
% Author:
% 05/06/2011 - FB
% 2012.07.05 FB - added inDataType and prev_steps.
% 2012.11.13 TC - if empty, set prev_steps to an empty cell array rather
%  than an empty double.

fields = {...
    'wmtx',...
    'HalfDecayTimes',...
    'ci_siglist',...
    'norm', ...
    'prev_steps', ...
    'inDataType'};
params = mkstruct(fields,varargin);
if isempty(params.prev_steps)
  params.prev_steps = {};
end

if isempty(params.wmtx)
  root_path = fileparts(which('jlmt_proc_series'));
  params.wmtx = struct('fname',fullfile(root_path,'maps','pp2pitchclass_W_20121105.mat'));
end

end % params_pc
