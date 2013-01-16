% this script finds approximate demarcation lines between major and minor
% strata in the flattened torus projection. It estimates the slope and
% intercept of each path through major and minor chords in the circle of
% fifths, then approximates demarcation lines between the regions following
% these paths as the average intercept and slope of adjacent path lines.
% 
% this script generates a static file that is referenced by
% 'get_toract_mode_estimate.m' in the 'utils' directory of the jlmt distro.
% 
% FB 2013.01.16

% initialize some variables
SAVE_NEW_FILE = 0;
if SAVE_NEW_FILE
  out_root = fileparts(which('jlmt_proc_series'));
  out_fname = fullfile(out_root,'maps',sprintf('toract_mode_map_%s.mat',...
      datestr(now(),'yyyymmdd')));
  if exist(out_fname,'file'), delete(outfname); end
end % if SAVE_NEW_FILE

% load toract data, for visualization purposes only
load('/home/fbarrett/data/partest/data.mat','outData'); %from jlmt analysis
oc = set_var_col_const(outData.vars);
toract = outData.data{oc.toract}{1};
tc = set_var_col_const(toract.vars);
map_fname = toract.data{tc.params}.toract.som.fname;
load(map_fname);
planes = toract.data{tc.activations}{1};
new_plane = planes(:,:,1);

% initialize labels
label_idxs = [];
new_label_array = {};
for(idx = 1:length(new_labels))
  if ~isempty(new_labels{idx})
    new_label_array{end+1} = new_labels{idx};
    label_idxs = [label_idxs idx];
  end
end

% plot torus activation
figure();
imagesc(new_plane);
title('torus projection with key centers labeled')

% plot key centers
try labelSize = params.labelSize; catch labelSize=14; end
map_coords = som_vis_coords(sM.topol.lattice,sM.topol.msize);

% filched from plot_torusMovieFrames.m
coord_translation = size(new_plane)./sM.topol.msize;
label_coords = (map_coords(label_idxs,:)-0.5).*repmat(coord_translation,length(label_idxs),1);
t = text(label_coords(:,1),label_coords(:,2),new_label_array, ...
     'fontweight','bold','fontsize', labelSize, ...
     'horizontalalignment', 'center', ...
     'verticalalign', 'middle');

% 
line_order = {...
    {'c','f'},...
    {'Bb','Eb','Ab'},...
    {'bb','d#','g#','c#'},...
    {'Db','Gb','B','E'},...
    {'f#','b','e'},...
    {'A','D','G','C'},...
    {'a','d','g'},...
    {'F'},...
    };
nlines = length(line_order);
line_class = [-1 1 -1 1 -1 1 -1 1]; % -1 = minor, 1 = major

% solve slope and intercept for n-1 lines, as the final line contains only
% one point
line_params = zeros(nlines-1,2);
for k=1:nlines-1
  idxs = find(ismember(new_label_array,line_order{k}));
  line_params(k,:) = polyfit(label_coords(idxs,1),label_coords(idxs,2),1);
end % for k=1:length(line_order

% draw up to n-1 lines, as the final line only contains one point
% each line is the average of two adjacent lines, excluding the final line
bound_params = nan(size(line_params));
for k=1:nlines-2
  avg_slope = mean(line_params(k:k+1,1));
  avg_int = mean(line_params(k:k+1,2));
  x = [0 -avg_int/avg_slope];
  y = [avg_int 0];
  line(x,y,'color','w');
  bound_params(k,:) = [avg_slope avg_int];
end % for k=1:size(line_params,1

% assume average slope, solve for intercept with last point, draw final line
final_slope = mean(line_params(:,1));
final_int = -(label_coords(end,1)*final_slope-label_coords(end,2));
avg_slope = mean([line_params(end,1) final_slope]);
avg_int = mean([line_params(end,2) final_int]);
x = [0 -avg_int/avg_slope];
y = [avg_int 0];
line(x,y,'color','w');
bound_params(end,:) = [avg_slope avg_int];

% get the average value for all cells within a given region
areas = cell(size(line_params,1)+1,1);
grp_mtx = nan(size(new_plane));
for k=1:size(planes,1)
  for l=1:size(planes,2)
    % iterate over lines, find first line that cell (k,l) is above
    for m=1:nlines-1
      ylim = bound_params(m,1)*l+bound_params(m,2);
      xlim = (k-bound_params(m,2))/bound_params(m,1);
      
      if l <= round(xlim) && k <= round(ylim)
        areas{m}(end+1,:) = planes(k,l,:);
        grp_mtx(k,l) = m;
        break
      end
    end % for m=1:nlines-1
    
    if isnan(grp_mtx(k,l))
      areas{m+1}(end+1,:) = planes(k,l,:);
      grp_mtx(k,l) = m+1;
    end % if isnan(grp_mtx(k,l
  end % for l=1:size(new_plane,2
end % for k=1:size(new_plane,1

figure();
imagesc(grp_mtx);
title('torus pixel group membership');

tcs = cell2mat(cellfun(@nanmean,areas,'UniformOutput',false));
mode_tc = line_class*tcs;
figure();
plot(mode_tc);
title('timecourse of mode estimate (maj > 0, min < 0)');

if SAVE_NEW_FILE
  save(out_fname,'line_params','bound_params','grp_mtx',...
      'line_order','line_class','label_coords','line_class','sM');
end % if SAVE_NEW_FILE
