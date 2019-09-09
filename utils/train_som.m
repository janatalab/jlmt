function [sM, sD, sT, outfname] = train_som(train_mtx, params)
% Trains a self-organizing map (SOM) using the Finnish somtoolbox
%
% Input:
% training_mtx - a matrix in which each row is an input representation
% params -
%
% If the string 'getDefaultParams' is passed in as the single argument, a
% structure containing the default parameter values is returned

% 08Sep2019 PJ - adapted initial script from train_som_w_melody.m

if nargin == 1 && strcmp(train_mtx, 'getDefaultParams')
  sM = getDefaultParams;
  return
end  

% Normalize the training mtx
if ~isempty(params.norm_func)
  train_mtx = params.norm_func(train_mtx);
end

% Initialize the self-organizing map data structure
disp('Initializing SOM structures');
sD = som_data_struct(train_mtx);

% Set up the surface of the SOM
input_space_dim = size(sD.data,2);
sM = som_map_struct(input_space_dim,'msize', params.map_dims);
sM.topol = som_topol_struct('msize', params.map_dims, 'shape', params.map_shape, 'lattice', ...
  params.map_lattice);

% Initialize the weights
sM = som_randinit(sD.data,sM);

% Do the actual training
disp('Training network');
[sM, sT] = som_batchtrain(sM,sD, ...
  'tracking', params.tracking, ...
  'neigh', params.neighborhood, ...
  'radius_ini', params.start_radius, ...
  'radius_fin', params.final_radius, ...
  'trainlen', params.trainlen);

% Optional handling of the results
if params.results.save
  % Construct a filename if the filestub is empty
  fname = params.results.fname;
  if isempty(fname)
    fname = sprintf('map-%s.mat', sM.trainhist(end).time);
    fname = regexprep(fname,'\s','-');
  end
    
  outfname = fullfile(params.results.path, fname);
    
  save_vars = {'sM','sT'};
  fprintf('Saving map data to %s\n', outfname);
  save(outfname, save_vars{:})
end

return % train_som()
end

function params = getDefaultParams()

params = struct();
params.norm_func = @(mtx) (mtx ./ repmat(sum(mtx,2),1,size(mtx,2)));
params.map_dims = [12 16]; % [12 16], [15 20], [18, 24], [24 32]
params.map_shape = 'toroid';
params.map_lattice = 'hexa';
params.neighborhood = 'gaussian';
params.start_radius = 3; % 6
params.final_radius = 1; % 2
params.trainlen = 200;
params.tracking = 2;
params.results.save = 0;
params.results.path = '/tmp';
params.results.fname = 'trained_som.mat';

return % getDefaultParams()
end