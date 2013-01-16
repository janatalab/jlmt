function params = params_mode_estimate(varargin)
% Returns a parameter structure for mode estimate from torus projection in jlmt_proc_series.m
%
%   params = params_mode_estimate(varargin);
% 
% Input parameters for mode estimate from torus projection within jlmt_proc_series.m
%
%     torus_mode_map.fname = 'jlmt/maps/toract_mode_map_20130116.mat'
%
% Copyright (c) 2006-2013 The Regents of the University of California
% All Rights Reserved.
%
% 2013.01.16 FB - adapted from params_toract.m

fields = {...
    'toract_mode_map',...
    };
params = mkstruct(fields,varargin);

jlmt_root = fileparts(which('jlmt_proc_series'));
def.toract_mode_map.fname = fullfile(jlmt_root,'maps','toract_mode_map_20130116.mat');

def.inDataType = '';
def.prev_steps = {};

for ifld = 1:length(fields)
  if isempty(params.(fields{ifld}))
    params.(fields{ifld}) = def.(fields{ifld});
  end
end
