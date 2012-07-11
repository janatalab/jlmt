function params = params_toract(varargin)
% Returns a parameter structure for torus projection in jlmt_proc_series.m
%
%   params = params_toract(varargin);
% 
% Input parameters for torus projection within jlmt_proc_series.m
%
% 	            li_siglist = []; % tc_2 if not specified
% 	        HalfDecayTimes = []; % 2 if not specified
% 	       calc_spher_harm = 1;
%   spher_harm.nharm_theta = 3;
%     spher_harm.nharm_phi = 4;
% 	   spher_harm.min_rsqr = 0.95;
% 	                  norm = 'x./repmat(sum(x),size(x,1),1)';
% 	             som.fname = []; % modulating-melody-trained som if not specified
% 	                    Fs = [];
%               inDataType = [];
%               prev_steps = [];
%
% Copyright (c) 2006-2012 The Regents of the University of California
% All Rights Reserved.
%
% 2011.01.10 FB - adapted from params_ci.m
% 2012.07.03 FB - move from ci_siglist to li_siglist, add inDataType and
%   prev_steps

fields = {...
    'li_siglist',...
    'ci_siglist',... % for backwards compability
    'HalfDecayTimes',...
    'calc_spher_harm',...
    'spher_harm',...
    'norm',...
    'som',...
    'Fs', ...
    'inDataType', ...
    'prev_steps'};
params = mkstruct(fields,varargin);

root_path = fileparts(fileparts(which('jlmt_proc_series')));

def.li_siglist = [];
def.ci_siglist = [];
def.HalfDecayTimes = [];
def.calc_spher_harm = 1;
def.spher_harm = struct('nharm_theta',3,'nharm_phi',4,'min_rsqr',0.95);
def.norm = 'x./repmat(sum(x),size(x,1),1)';
def.som = struct('fname',fullfile(root_path,'data',...
    'map_10-Dec-2006_16:18.mat'));
def.Fs = [];
def.inDataType = '';
def.prev_steps = [];

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
