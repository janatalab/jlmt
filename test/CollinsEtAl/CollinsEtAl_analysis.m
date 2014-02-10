% Copyright (c) 2013 The Regents of the University of California
% All Rights Reserved.

% This script is intended to accompany:
%
%  Tom Collins, Barbara Tillmann, Frederick S. Barrett, Charles Delbe, and
%   Petr Janata. A combined model of sensory and cognitive representations
%   underlying tonal expectation: from audio signals to behavior.
%   Psychological Review, 121(1):33-65, 2014.
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
zeroMeanRTobs = [];
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
          ~isempty(regexp(contents(idir).name, '.txt', 'once')) || ...
          strcmp(contents(idir).name,'.svn')
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
  zeroMeanRTobs = [zeroMeanRTobs; datasets(iset).zeroMeanRTobs];
end % for iset=

%% Run jlmt_proc_series.
jlmt_out = jlmt_proc_series(in_data, params.jlmt);

%% Calculate metrics.
% Parameters within this script can be experimented with if you wish.
results = CollinsEtAl_calc_attributes(in_data, params, jlmt_out);

%% Predict zero-mean response times using our regression equation (2) from
% the paper. The variables in this model are:
exp_vars2 = {...
  'TS' [.1 4] 'MC' 'rel' [0 200] 7;...            % x_TS
  'TS' [.1 4] {'CL' 'hypo'} 'NA' [100 300] 25;... % p_clos
  'PP' [.1 4] 'MC' 'rel' [0 200] 3;...            % x_PP
  'PP' [.1 4] 'MC' 'abs' [201 600] 2;...          % z_PP
  'PP' 4 'MV' 'abs' [201 600] 14;...              % y_PP
  'CV' 4 'MV' 'abs' [0 200] 21};                  % z_PP
% The commented column above is the label used in the paper. The last
% column in the cell is for checking that the current results match a
% previous analysis.
nresults = size(results, 1);
stim_names2 = in_data.data{4}';
% Coefficient values from equation (2).
alphaHat = -57.370601599275098;
betaHat = [-1.117975643676257e+02 -1.415175120627715e+04...
  -3.570081944978055e+02 2.258324648389623e+02 0.141210000049002...
  -3.444465129538324];
nexp_vars = size(exp_vars2, 1);
nstim = size(stim_names2, 1);
zeroMeanRTfit = zeros(nstim, 1);
X2 = zeros(nstim, 6);
% For each stimulus, get the regression variable values and form a linear
% combintation with the coefficients to estimate zero-mean response time.
for istim = 1:nstim
  x = zeros(nexp_vars, 1);
  for ivar = 1:nexp_vars
    % Find the relevant row for this stimulus and variable combination.
    rel_row = [];
    irow = 2;
    while irow <= nresults
      % Test match on stimulus string.
      if strcmp(results{irow, 2}, stim_names2{istim}) &&...
          ... % Test match on representational space.
          strcmp(results{irow, 3}, exp_vars2{ivar, 1}) &&...
          ... % Test match on time constants.
          sum(results{irow, 4} == exp_vars2{ivar, 2})...
          == size(exp_vars2{ivar, 2}, 2) &&...
          ... % Test match on calculation type.
          sum(strcmp(results{irow, 5}, exp_vars2{ivar, 3}))...
          == size(exp_vars2{ivar, 2}, 3) &&...
          ... % Test match on window comparison.
          strcmp(results{irow, 6}, exp_vars2{ivar, 4}) &&...
          ... % Test match on post-target window.
          sum(results{irow, 7} == exp_vars2{ivar, 5})...
          == size(exp_vars2{ivar, 5}, 2)
        rel_row = irow;
        irow = nresults;
      end
      irow=irow+1;
    end
    % Found the relevant row. Now assign the value to the variable.
    x(ivar) = results{rel_row, 8};
  end
  % Form a linear combintation with the coefficients to estimate zero-mean
  % response time.
  zeroMeanRTfit(istim) = alphaHat + betaHat*x;
end

% Now plot the fitted and observed RTs and add a line of best fit.
plot(zeroMeanRTobs, zeroMeanRTfit, 'ko')
hold on
pfit = polyfit(zeroMeanRTobs, zeroMeanRTfit, 1);
l = pfit(1)*(-250:1:250) + pfit(2);
plot(-250:1:250, l, '-');
hold off
title('Plot of Fitted v Observed RT Using Eq. 2 from Paper')
xlabel('Observed Response Time (Zero-Mean)')
ylabel('Fitted Response Time (Zero-Mean)')
