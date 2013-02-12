function [mode_tc,avg,pct_major,pct_minor] = mode_est_ratio(planes,grp_mtx)

% estimate a mode timecourse from the ratio of mean major and minor activation
% 
%   [mode_tc,avg,pct_major,pct_minor] = mode_est_ratio(planes,grp_mtx)
% 
% This function calculates the ratio of the average pixel value in major
% areas compared to the average pixel value in minor areas of a torus
% projection timecourse, at each sample. It returns a timecourse of this
% ratio, the average of this ratio for the entire piece, and the percentage
% of time that the piece is estimated to have been in major and minor mode.
% 
% This is a sub-function, to be used with calc_mode_estimate.
% 
% REQUIRES
%   planes
%   grp_mtx
% 
% RETURNS
%   mode_tc
%   avg - average ratio of strength of major to minor areas
%   pct_major - percentage of samples estimated to be in major mode
%   pct_minor - percentage of samples estiamted to be in minor mode
% 
% FB 2013.02.05

mode_tc = nan(1,size(planes,3));

for k=1:size(planes,3)
  new_plane = squeeze(planes(:,:,k));

  % get ratio of major to minor
  mode_tc(k) = (mean(new_plane(logical(grp_mtx)))+eps)/...
      (mean(new_plane(logical(~grp_mtx)))+eps);
  plot_plane;
end % for k=1:size(planes,3

avg = mean(mode_tc);
pct_major = sum(mode_tc > 1)/length(mode_tc);
pct_minor = sum(mode_tc < 1)/length(mode_tc);
