function params = params_toract(varargin)
% Returns a parameter structure for torus projection in jlmt_proc_series.m
%
%   params = params_toract(varargin);
% 
% Input parameters for torus projection within jlmt_proc_series.m
%
% 	                  li_siglist = [];
% 	              HalfDecayTimes = [];
% 	             calc_spher_harm = 1;
%     	  spher_harm.nharm_theta = 3;
%       	spher_harm.nharm_phi = 4;
% 	         spher_harm.min_rsqr = 0.95;
% 	                        norm = [];
% 	                   som.fname = [];
% 	                          Fs = [];
%
% Copyright (c) 2006 The Regents of the University of California
% All Rights Reserved.
%
% 2011.01.10 FB - adapted from params_ci.m
% 2012.07.03 FB - move from ci_siglist to li_siglist

fields = {...
    'li_siglist',...
    'ci_siglist',... % for backwards compability
    'HalfDecayTimes',...
    'calc_spher_harm',...
    'spher_harm',...
    'norm',...
    'som',...
    'Fs'};
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
%   params.metrics = struct('toract_corr',struct('tc_pair',[]));
%   params.metrics = struct('toract_kldist',struct('tc_pair',[]));
% note: if we stick with verifying sub-structs of main analysis params in
% ipem_proc_series.getDefaultParameters, we will have to add a default
% definition for every toract metric function that we use
