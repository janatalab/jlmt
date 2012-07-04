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
%       nnet.fname: []
%   HalfDecayTimes: []
%       ci_siglist: []
%             norm: []
%          metrics: []
%
% Copyright (c) 2011 The Regents of the University of California
% All Rights Reserved.
%
%
% Author:
% 05/06/2011 - FB

fields = {...
    'nnet',...
    'HalfDecayTimes',...
    'ci_siglist',...
    'norm',...
    'metrics'};
params = mkstruct(fields,varargin);

if isempty(params.nnet)
  params.nnet = struct('fname',[]);
end
if isempty(params.metrics)
  params.metrics = struct();
end

end % params_pc
