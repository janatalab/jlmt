function mode_tc = mode_est_max(planes,grp_mtx)

% estimate a mode timecourse from the max torus projection across time
% 
%   mode_tc = mode_est_max(planes,grp_mtx)
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

  % get max value for this timepoint
  [myidx,mxidx] = ...
      find(new_plane == repmat(max(new_plane(:)),size(new_plane)));
  mode_tc(k) = grp_mtx(myidx,mxidx);
end % for k=1:size(planes,3
