% load/init variables
if ~exist('sM','var')
  try mode_map_fname = params.toract_mode_map.fname;
  catch mode_map_fname = fullfile(fileparts(which('jlmt_proc_series')),...
          'maps','toract_mode_map_20130118.mat');
  end

  if exist('grp_mtx','var')
    grp_mtx_bak = grp_mtx;
  end
  load(mode_map_fname);
  if exist('grp_mtx_bak','var')
    grp_mtx = grp_mtx_bak;
  end
  map_coords = som_vis_coords(sM.topol.lattice,sM.topol.msize);
  try labelSize = params.labelSize; catch labelSize=14; end
  try pause_length = params.pause_length; catch pause_length = 0.25; end
end % if exist('sM

figure(1);
clf;
subplot(1,2,1)
imagesc(new_plane);

coord_translation = size(new_plane)./sM.topol.msize;
label_coords = (map_coords(label_idxs,:)-0.5).*repmat(coord_translation,length(label_idxs),1);
t = text(label_coords(:,1),label_coords(:,2),new_label_array, ...
  'fontweight','bold','fontsize', labelSize, ...
  'horizontalalignment', 'center', ...
  'verticalalign', 'middle');

for l=1:size(bound_params,1)
  avg_slope = mean(bound_params(l,1));
  avg_int = mean(bound_params(l,2));
  x = [0 -avg_int/avg_slope];
  y = [avg_int 0];
  line(x,y,'color','w');
end % for l=1:size(bound_params,1

subplot(1,2,2)
imagesc(grp_mtx);

if exist('mxidx','var')
  subplot(1,2,1);
  text(mxidx,myidx,'X','color','w','fontsize',15,'fontweight','bold');
  subplot(1,2,2)
  text(mxidx,myidx,'X','color','w','fontsize',15,'fontweight','bold');
else
  title(mode_tc(k));
end % if exist('mxidx

pause(0.25);
  
  