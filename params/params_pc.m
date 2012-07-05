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
% 2012.07.05 FB - added inDataType and prev_steps

fields = {...
    'nnet',...
    'HalfDecayTimes',...
    'ci_siglist',...
    'norm', ...
    'prev_steps', ...
    'inDataType'};
params = mkstruct(fields,varargin);

if isempty(params.nnet)
  params.nnet = struct('fname',[]);
end

end % params_pc
