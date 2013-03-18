function plot_toract_w_mode_map(inData,params)

% plot a torus projection with an overlaid mode map
% 
%   plot_toract_w_mode_map(inData,params)
% 
% Copyright (c) 2006-2013 The Regents of the University of California
% All Rights Reserved.
%
% REQUIRES
%   inData - output from calc_toract
%   params
%       .toract_mode_map.fname - path to file containing the mode map for
%           the given toract model. This must contain the variables
%           'bound_params' and 'grp_mtx', which can be generated using the
%           utility file 'make_toract_mode_map.m'
%       .li_siglist - names of images in inData to use
%       .pause_length (default: 0.25s) - length in seconds to pause between
%           displaying successive toract frames

% FB 2013.01.17

% load/init variables
try mode_map_fname = params.toract_mode_map.fname;
catch mode_map_fname = fullfile(fileparts(which('jlmt_proc_series')),...
        'maps','toract_mode_map_20130118.mat');
end
load(mode_map_fname);

map_coords = som_vis_coords(sM.topol.lattice,sM.topol.msize);
try labelSize = params.labelSize; catch labelSize=14; end
try pause_length = params.pause_length; catch pause_length = 0.25; end
tc = set_var_col_const(inData.vars);
nsig = length(params.li_siglist);

% iterate over signals, display them
figure();
for isig=1:nsig
    signame = params.li_siglist{isig};
    lidx = find(ismember(inData.data{tc.labels},signame));
    if isempty(lidx)
      warning('no context image found for %s, SKIPPING',signame);
      continue
    end
    
    % iterate over frames of the signal
    planes = inData.data{tc.activations}{lidx};
    for k=1:size(planes,3)
      tic
      new_plane = planes(:,:,k);

      % plot torus activation for the current frame
      clf();
      imagesc(new_plane);
      title(sprintf('%s, sample %d/%d',...
          regexprep(signame,'_',' '),k,size(planes,3)));

      % label keys
      coord_translation = size(new_plane)./sM.topol.msize;
      label_coords = (map_coords(label_idxs,:)-0.5).*...
          repmat(coord_translation,length(label_idxs),1);
      text(label_coords(:,1),label_coords(:,2),new_label_array, ...
          'fontweight','bold','fontsize', labelSize, ...
          'horizontalalignment', 'center', ...
          'verticalalign', 'middle');

      % draw mode boundaries
      for l=1:size(bound_params,1)
        avg_slope = bound_params(l,1);
        avg_int = bound_params(l,2);
        x = [0 -avg_int/avg_slope];
        y = [avg_int 0];
        line(x,y,'color','w');
      end % for l=1:size(bound_params,1

      pause(pause_length);
      toc
    end % for k=1:size(planes,3
    
end % for isig=1:nsig
