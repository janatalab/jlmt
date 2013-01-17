function params = params_mode_estimate(varargin)
% Returns a parameter structure for mode estimate from torus projection in jlmt_proc_series.m
%
%   params = params_mode_estimate(varargin);
% 
% Input parameters for mode estimate from torus projection within jlmt_proc_series.m
%
%     torus_mode_map.fname = 'jlmt/maps/toract_mode_map_20130116.mat'
% 	            li_siglist = []; % tc_2 if not specified
% 	        HalfDecayTimes = []; % 2 if not specified
%
% Copyright (c) 2006-2013 The Regents of the University of California
% All Rights Reserved.
%
% 2013.01.16 FB - adapted from params_toract.m

fields = {...
    'toract_mode_map',...
    'li_siglist',...
    'ci_siglist',... % for backwards compability
    'HalfDecayTimes',...
    'prev_steps',...
    'inDataType',...
    };
params = mkstruct(fields,varargin);

jlmt_root = fileparts(which('jlmt_proc_series'));
def.toract_mode_map.fname = fullfile(jlmt_root,'maps','toract_mode_map_20130116.mat');
def.li_siglist = [];
def.ci_siglist = [];
def.HalfDecayTimes = [0.1 2];

def.inDataType = '';
def.prev_steps = {};

for ifld = 1:length(fields)
  if isempty(params.(fields{ifld}))
    params.(fields{ifld}) = def.(fields{ifld});
  end
end

if ~isempty(params.ci_siglist)
  % for backwards compatibility; ci_siglist is deprecated
  params.li_siglist = params.ci_siglist;
end
if isempty(params.HalfDecayTimes) && isempty(params.li_siglist)
  params.HalfDecayTimes = [2];
end
if isempty(params.li_siglist)
  params.li_siglist = calc_li('calc_names',params.HalfDecayTimes);
end