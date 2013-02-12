function [mode_tc,avg,pct_major,pct_minor] = mode_est_hr_ratio_2(planes,grp_mtx)

% estimate a mode timecourse from ratio of mean major and minor activation on a separate half-rectified torus projection
% 
%   [mode_tc,avg,pct_major,pct_minor] = mode_est_hr_ratio_2(planes,grp_mtx)
% 
% This function calculates the ratio of the average pixel value in major
% areas compared to the average pixel value in minor areas of a torus
% projection timecourse, at each sample. Before doing so, each sample is
% half-rectified separately for major and minor mode areas, meaning that
% values below the mean for the major areas, and separately for the minor
% mode areas, are set to 0 before the mode-region average is calculated. 
% This function returns a timecourse of this ratio, the average of this
% ratio for the entire piece, and the percentage of time that the piece is
% estimated to have been in major and minor mode.
% 
% This is a sub-function, to be used with calc_mode_estimate.
% 
% REQUIRES
%   planes
%   grp_mtx
% 
% RETURNS
%   mode_tc
%   avg - average ratio strength of major to minor areas after rectification
%   pct_major - percentage of samples estimated to be in major mode
%   pct_minor - percentage of samples estiamted to be in minor mode
% 
% FB 2013.02.05

mode_tc = nan(1,size(planes,3));

for k=1:size(planes,3)
  new_plane = squeeze(planes(:,:,k));
  maj_plane = new_plane(logical(grp_mtx));
  min_plane = new_plane(logical(~grp_mtx));

  % get half-rectified ratio
  maj_plane(maj_plane<mean(maj_plane(:))) = NaN;
  min_plane(min_plane<mean(min_plane(:))) = NaN;
  mode_tc(k) = (nanmean(maj_plane)+eps)/(nanmean(min_plane)+eps);
  plot_plane;
end % for k=1:size(planes,3

avg = mean(mode_tc);
pct_major = sum(mode_tc > 1)/length(mode_tc);
pct_minor = sum(mode_tc < 1)/length(mode_tc);
