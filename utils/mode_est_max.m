function [mode_tc,avg,pct_major,pct_minor] = mode_est_max(planes,grp_mtx)

% estimate a mode timecourse from the max torus projection across time
% 
%   [mode_tc,avg,pct_major,pct_minor] = mode_est_max(planes,grp_mtx)
% 
% This function identifies the region in which the maximum value for the
% torus projection is found in each sample. It returns a timecourse
% indicating the region (major = 1, minor = 0) containing the max value at
% each sample. It also returns the average of this timecourse for the
% entire piece and the percentage of time that the piece is estimated to
% have been in major and minor mode.
% 
% This is a sub-function, to be used with calc_mode_estimate.
% 
% REQUIRES
%   planes
%   grp_mtx
% 
% RETURNS
%   mode_tc
%   avg - average mode of the piece (major = 1, minor = 0)
%   pct_major - percentage of samples in major mode
%   pct_minor - percentage of samples in minor mode
% 
% Copyright (c) 2013 The Regents of the University of California
% All Rights Reserved.

% FB 2013.02.05

mode_tc = nan(1,size(planes,3));

for k=1:size(planes,3)
  new_plane = squeeze(planes(:,:,k));

  % get max value for this timepoint
  [myidx,mxidx] = ...
      find(new_plane == repmat(max(new_plane(:)),size(new_plane)));
  mode_tc(k) = grp_mtx(myidx,mxidx);
  plot_plane;
end % for k=1:size(planes,3

avg = mean(mode_tc);
pct_major = sum(mode_tc == 1)/length(mode_tc);
pct_minor = sum(mode_tc == 0)/length(mode_tc);
