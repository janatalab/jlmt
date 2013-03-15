% Copyright (c) 2012 The Regents of the University of California
% All Rights Reserved.

% This script is intended to accompany:
%
%  Tom Collins, Barbara Tillmann, Frederick S. Barrett, Charles Delbé, and
%   Petr Janata. A combined model of sensory and cognitive representations
%   underlying tonal expectation: from audio signals to behavior.
%   Submitted, 2012.
%
% Please cite the paper if you use/adapt this code in your own work. This
% script sets parameters and then iterates over the datasets defined in
% the function CollinsEtAl_datasets, creating path and variable names of
% audio files for subsequent analysis. You can place your own audio files
% in the jlmt data folder and have these analyzed as well/instead.
% Conforming to our existing folder structure when placing audio will help
% to avoid errors, and it will also be necessary to extend/modify the
% function CollinsEtAl_datasets to include your audio files. Event onsets
% (or at least the time of the target event) need to be specified there in
% order to run the function CollinsEtAl_calc_attributes. If you do not have
% this information available, it is possible to modify our code to estimate
% onsets instead.
%
% The function jlmt_proc_series creates the representations of audio that
% are discussed in the paper (periodicity pitch, chroma vector, and tonal
% space), as well as some others. The function CollinsEtAl_calc_attributes
% iterates over the output of jlmt_proc_series, and calculates attributes
% of the representations in a temporal region of target events. For
% example, one attribute calculates the correlation between short- and
% long-term tonal space activations, averaged over a 0-200 ms post-target
% window. Please see the above paper and the code in
% CollinsEtAl_calc_attributes for more details. In particular, at the top
% of this function it is possible to experiment with the default parameter
% values.

% Tom Collins, 2012.07.21.

%% Initialize parameters
params = CollinsEtAl_globals;

% Initialize structure that will be passed to jlmt_proc_series
in_data = ensemble_init_data_struct;
in_data.vars = {'path', 'filename', 'path_no_ext', 'name_no_ext'};
ic = set_var_col_const(in_data.vars);
in_data.data{ic.path} = {};
in_data.data{ic.filename} = {};
in_data.data{ic.path_no_ext} = {};
in_data.data{ic.name_no_ext} = {};

%% Iterate over datasets, collect stimuli to analyze.
datasets = params.datasets;
nsets = length(datasets);
for iset = 1:nsets
  fprintf('\n\nDataset %d\n', iset);
  fprintf('Identifier: %s\n', datasets(iset).id);
  % Determine file path
  location_list = {};
  location_root = fullfile(params.paths.data_root,...
    params.datasets(iset).id);
  if isfield(datasets(iset), 'subsets')
    nsub = length(datasets(iset).subsets);
    for isub = 1:nsub
      fprintf('\tSubset %d: %s\n', isub, datasets(iset).subsets(isub).id);
      location_list{isub} = fullfile(location_root,...
        params.datasets(iset).subsets(isub).id);
    end
  else
    location_list = {location_root};
  end
  for iloc = 1:length(location_list)
    location = location_list{iloc};
    params.paths.stim_root = location;
    % Now create a data structure that is compatible with jlmt_proc_series,
    % and run jlmt_proc_series.
    contents = dir(location);
    ndir = size(contents, 1);
    for idir = 1:ndir
      % If the name folder is one of the following, do nothing.
      if strcmp(contents(idir).name, '.') || ...
          strcmp(contents(idir).name, '..') || ...
          strcmp(contents(idir).name, '.DS_Store') || ...
          ~isempty(regexp(contents(idir).name, '.txt', 'once'));
        continue
      end
      % Get parts of path and filename.
      currPath = fullfile(location, contents(idir).name, 'audio');
      flist = listFilesOfType(currPath, {'wav','mp3'});
      if length(flist) > 1
        error('Found %d files in %s\n', length(flist), currPath);
      end
      in_data.data{ic.path}(end+1) = {currPath};
      in_data.data{ic.filename}(end+1) = flist;
      in_data.data{ic.path_no_ext}(end+1) = {fullfile(location,...
        contents(idir).name)};
      in_data.data{ic.name_no_ext}(end+1) = {contents(idir).name};
    end % for idir=
  end % for iloc
end % for iset=

%% Run jlmt_proc_series.
jlmt_out = jlmt_proc_series(in_data, params.jlmt);

%% Calculate metrics.
% Parameters within this script can be experimented with if you wish.
results = CollinsEtAl_calc_attributes(in_data, params, jlmt_out);
