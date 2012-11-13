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
% 	                    Fs = [];
%               inDataType = [];
%               prev_steps = [];
% 
% 	             som.fname = []; % mus specify either:
%                   the pp->li->toract at:
%                       jlmt/data/maps/map_10-Dec-2006_16:18.mat
%                   OR the pp->pc->li->toract map at:
%                       jlmt/data/maps/pc_ci2toract_map_12-Jun-2012_15:47.mat
%
% Copyright (c) 2006-2012 The Regents of the University of California
% All Rights Reserved.
%
% 2011.01.10 FB - adapted from params_ci.m
% 2012.07.03 FB - move from ci_siglist to li_siglist, add inDataType and
%   prev_steps
% 2012.10.30 PJ - added context dependent assignment of weight matrices
% 2012.11.13 TC - Altered prev_steps to empty cell array, not double.

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

def.li_siglist = [];
def.ci_siglist = [];
def.HalfDecayTimes = [0.1 2];
def.calc_spher_harm = 1;
def.spher_harm = struct('nharm_theta',3,'nharm_phi',4,'min_rsqr',0.95);
def.norm = 'x./repmat(sum(x),size(x,1),1)';
def.som = [];
% Assign the default weight matrix based on context
switch cell2str(params.prev_steps,'_')
  case 'ani_pp_li'
    def.som.fname = 'map_10-Dec-2006_16:18.mat';
  case 'ani_pp_pc_li'
    def.som.fname = 'pc_ci2toract_map_12-Jun-2012_15:47.mat';
  otherwise
    def.som.fname = '';
end
def.Fs = [];
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
