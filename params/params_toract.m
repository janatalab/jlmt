function params = params_toract(varargin)
% Returns a parameter structure for torus projection in jlmt_proc_series.m
%
%   params = params_toract(varargin);
% 
% Input parameters for torus projection within jlmt_proc_series.m
%
% 	            li_siglist = [];
% 	        HalfDecayTimes = [];
% 	       calc_spher_harm = 1;
%   spher_harm.nharm_theta = 3;
%     spher_harm.nharm_phi = 4;
% 	   spher_harm.min_rsqr = 0.95;
% 	                  norm = 'x./repmat(sum(x),size(x,1),1)';
% 	             som.fname = [];
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

if ~isempty(params.ci_siglist)
  % for backwards compatibility; ci_siglist is deprecated
  params.li_siglist = params.ci_siglist;
end
if isempty(params.som)
  params.som = struct('fname',[]);
end
if isempty(params.spher_harm)
  params.spher_harm = struct('nharm_theta',[],'nharm_phi',[],'min_rsqr',[]);
end
if isempty(params.norm)
  params.norm = 'x./repmat(sum(x),size(x,1),1)';
end
