function mode_tc = mode_est_hr_ratio(planes,grp_mtx)

% estimate a mode timecourse from ratio of mean major and minor activation on a half-rectified torus projection
% 
%   mode_tc = mode_est_hr_ratio(planes,grp_mtx)
% 
% This is a sub-function, to be used with calc_mode_estimate.
% 
% REQUIRES
%   planes
%   grp_mtx
% 
% RETURNS
%   mode_tc
% 
% FB 2013.02.05

mode_tc = nan(1,size(planes,3));

for k=1:size(planes,3)
  new_plane = squeeze(planes(:,:,k));

  % get half-rectified ratio
  new_plane(new_plane<0) = 0;
  mode_tc(k) = (mean(new_plane(logical(grp_mtx)))+eps)/...
      (mean(new_plane(logical(~grp_mtx)))+eps);
end % for k=1:size(planes,3
