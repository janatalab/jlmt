function [mode_tc,avg,pct_major,pct_minor] = mode_est_max_ratio(planes,grp_mtx)

% estimate a mode timecourse from the max torus projection across time
% 
%   [mode_tc,avg,pct_major,pct_minor] = mode_est_max(planes,grp_mtx)
% 
% This function identifies the maximum value in major regions and the
% maximum value in minor regions of a given torus projection timecourse at
% each sample, and calculates estimated mode as the ratio of the max major
% to max minor value. It returns a timecourse of this ratio, the average of
% this ratio for the entire piece, and the percentage of time that the
% piece is estimated to have been in major and minor mode.
% 
% This is a sub-function, to be used with calc_mode_estimate.
% 
% REQUIRES
%   planes
%   grp_mtx
% 
% RETURNS
%   mode_tc - time course of ratio of max major to max minor toract value
%   avg - average ratio max major to max minor value
%   pct_major - percentage of samples estimated to be in major mode
%   pct_minor - percentage of samples estiamted to be in minor mode
% 
% Copyright (c) 2013 The Regents of the University of California
% All Rights Reserved.

% FB 2013.02.05

mode_tc = nan(1,size(planes,3));
lgrp_mtx = logical(grp_mtx);
for k=1:size(planes,3)
  new_plane = squeeze(planes(:,:,k));
  new_plane_maj = new_plane(lgrp_mtx);
  new_plane_min = new_plane(~lgrp_mtx);
  
  % get max value for this timepoint
  majidx = new_plane_maj == max(new_plane_maj);
  minidx = new_plane_min == max(new_plane_min);
  
  mode_tc(k) = new_plane_maj(majidx)/new_plane_min(minidx);
  
  plot_plane;
end % for k=1:size(planes,3

avg = mean(mode_tc);
pct_major = sum(mode_tc > 1)/length(mode_tc);
pct_minor = sum(mode_tc < 1)/length(mode_tc);
