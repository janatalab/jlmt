function results = CollinsEtAl_calc_attributes(in_data, params, jlmt_out)

% Copyright (c) 2012 The Regents of the University of California
% All Rights Reserved.

% This function calculates attributes of audio representations in a
% temporal region of target events. The types of audio representation used
% and attributes calcuated may be altered by the user.

% INPUT
%  in_data is a structure consisting of path and file names of audio that
%   has been analysed by jlmt_proc_series.
%  params is a structure containig path names and parameter values for
%   jlmt_proc_series.
%  jlmt_out is the result of running jlmt_proc_series. Because there are
%   three procssing routes for this project, there are three cells per
%   audio stimulus in the output. We are interested in the first (direct
%   route to torus), second (via chroma vector), and third (rhythm profile)
%   of each of these. NB, there is only one cell per audio stimulus in
%   in_data.

% Tom Collins, 2012.07.21.
% Petr Janata, 2011.11.08.

% Default parameters.
integration_constants = [.1 2 4];
representational_space = {'PP' 'CV' 'TS'};
calculation_type = {'MC' 'MV' 'CL'};
window_comparison = {'abs' 'rel'};
post_target_window = {[0 200] [201 600]};
pre_target_window = [-100 0];
params.closure.post_event_window = [100 300];
closure_distributions = {'hypo' 'empr' 'smth'};

% Make sure integration constants are in ascending order.
integration_constants = sort(integration_constants, 'ascend');
nstim = size(in_data.data{1}, 2);
Nstim = size(jlmt_out.data{1}, 1);
nconstant = size(integration_constants, 2);
nspace = size(representational_space, 2);
ntype = size(calculation_type, 2);
ncomparison = size(window_comparison, 2);
nwindow = size(post_target_window, 2);
ndb = size(closure_distributions, 2);
% Preallocate results. There will be fewer results due to some
% combinations not being defined.
nrow = nstim*nconstant*nspace*ntype*ncomparison*nwindow +...
  nstim*ndb*nconstant*(nconstant-1)/2;
% Define output table.
results = cell(nrow, 8);

% Load closure distributions.
if ismember('CL', calculation_type)
  ndb = size(closure_distributions, 2);
  distrbn = cell(1, ndb);
  for idb = 1:ndb
  closure_struct = load(fullfile(params.paths.closure_struct,...
      ['closure_struct_' closure_distributions{idb} '.mat']));
    distrbn{idb} = closure_struct.closure_struct;
  end
end

jres = 2; % Iterate over results table (leave row 1 blank for titles).
ic = set_var_col_const(in_data.vars);
jc = set_var_col_const(jlmt_out.vars);
exp_names = {params.datasets.id};
nexp = size(exp_names, 2);

%% Iterate over stimuli.
for Istim = 1:3:Nstim
  istim = (Istim + 2)/3; % Counter for in_data (one row per stimulus).
  fprintf('Calculating attributes for stimulus %d of %d.\n', istim, nstim)
  % Get ontime of target note/chord. Begin by finding its experiment.
  iexp = 1; % Increment over exp_names.
  while iexp <= nexp
    if ~isempty(regexp(jlmt_out.data{jc.stimulus_path}{Istim},...
        exp_names{iexp}, 'once'))
      exp_idx = iexp;
      iexp = nexp;
    end
    iexp = iexp + 1;
  end
  if exp_idx <= nexp
    targ_on = params.datasets(exp_idx).event_onsets(end);
  else
    error('calc_attributes:event_onsets',...
      'No event onsets found for this stimulus.')
  end
  %% Iterate over pairs of integration constants.
  for iconstant = 1:nconstant
    for jconstant = iconstant+1:nconstant
      %% Iterate over representational space.
      for ispace = 1:nspace
        % Get local and global variables representing this space,
        % referred to as svar_loc and svar_glo.
        switch representational_space{ispace}
          case 'PP'
            % Variables in pp folder.
            svar = jlmt_out.data{jc.li}{Istim};
            sc = set_var_col_const(svar.vars);
            % Half decay times and sampling rate.
            hdt = svar.data{sc.params}.li.HalfDecayTimes;
            Fs = svar.data{sc.Fs};
            % Find local time constant, define space variable.
            [~, idx_loc] = ismember(integration_constants(iconstant),hdt);
            svar_loc = svar.data{sc.signals}{idx_loc};
            % Find global time constant, define space variable.
            [~, idx_glo] = ismember(integration_constants(jconstant),hdt);
            svar_glo = svar.data{sc.signals}{idx_glo};
            nsamp = size(svar_loc, 2); % Samples.
          case 'CV'
            % Variables in cv folder.
            svar = jlmt_out.data{jc.li}{Istim+1};
            sc = set_var_col_const(svar.vars);
            % Half decay times and sampling rate.
            hdt = svar.data{sc.params}.li.HalfDecayTimes;
            Fs = svar.data{sc.Fs};
            % Find local time constant, define space variable.
            [~, idx_loc] = ismember(integration_constants(iconstant),hdt);
            svar_loc = svar.data{sc.signals}{idx_loc};
            % Find global time constant, define space variable.
            [~, idx_glo] = ismember(integration_constants(jconstant),hdt);
            svar_glo = svar.data{sc.signals}{idx_glo};
            nsamp = size(svar_loc, 2); % Samples.
          case 'TS'
            % Variables in toract folder.
            svar = jlmt_out.data{jc.toract}{Istim};
            sc = set_var_col_const(svar.vars);
            % Half decay times and sampling rate.
            hdt = svar.data{sc.params}.toract.HalfDecayTimes;
            Fs = svar.data{sc.Fs};
            % Find local time constant, define space variable.
            [~, idx_loc] = ismember(integration_constants(iconstant),hdt);
            sarr_loc = svar.data{sc.activations}{idx_loc};
            % Find global time constant, define space variable.
            [~, idx_glo] = ismember(integration_constants(jconstant),hdt);
            sarr_glo = svar.data{sc.activations}{idx_glo};
            % Reshape these arrays to matrices.
            sizarr = size(sarr_loc);
            svar_loc = zeros(sizarr(1)*sizarr(2), sizarr(3));
            svar_glo = zeros(sizarr(1)*sizarr(2), sizarr(3));
            for isamp = 1:sizarr(3)
              sarr_temp = sarr_loc(:, :, isamp);
              svar_loc(:, isamp) = sarr_temp(:);
              sarr_temp = sarr_glo(:, :, isamp);
              svar_glo(:, isamp) = sarr_temp(:);
            end
            nsamp = size(svar_loc, 2); % Samples.
        end
        %% Iterate over calculation type.
        for itype = 1:ntype
          scal = zeros(nsamp, 1);
          switch calculation_type{itype}
            case 'MC'
              for isamp = 1:nsamp
                corr_temp = corrcoef(svar_loc(:, isamp),...
                  svar_glo(:, isamp));
                scal(isamp) = corr_temp(2);
              end
            case 'MV'
              scal = svar_glo;
            case 'CL'
              % Closure probabilities only calculated for TS. Get the
              % correlation time course and rhythm profile.
              if strcmp(representational_space{ispace}, 'TS')
                for isamp = 1:nsamp
                  corr_temp = corrcoef(svar_loc(:, isamp),...
                    svar_glo(:, isamp));
                  scal(isamp) = corr_temp(2);
                end
                rvar = jlmt_out.data{jc.rp}{Istim+2};
                ecm = CollinsEtAl_closure_event_corr_mean(...
                  scal, Fs, rvar, params);
                for idb = 1:ndb
                  % Need to check this function.
                  p = CollinsEtAl_closure_probability( ecm, distrbn{idb});
                  results{jres, 1} = in_data.data{ic.path_no_ext}{istim};
                  results{jres, 2} = in_data.data{ic.name_no_ext}{istim};
                  results{jres, 3} = representational_space{ispace};
                  results{jres, 4} = [integration_constants(iconstant)...
                    integration_constants(jconstant)];
                  results{jres, 5} = {calculation_type{itype},...
                    closure_distributions(idb)};
                  results{jres, 6} = 'NA';
                  results{jres, 7} = params.closure.post_event_window;
                  results{jres, 8} = p;
                  jres = jres + 1;
                end
              end % if strcmp(rep_space{ispace}, 'TS')
          end
          if ismember(calculation_type{itype}, {'MC', 'MV'})
            %% Iterate over post-target window.
            tpts = (0:nsamp-1)/Fs*1000; % Times in msec.
            for iwindow = 1:nwindow
              % Calculate the post-target attribute value.
              post_start_idx = find(tpts >= targ_on +...
                post_target_window{iwindow}(1), 1, 'first');
              post_end_idx = find(tpts <= targ_on +...
                post_target_window{iwindow}(2), 1, 'last');
              switch calculation_type{itype}
                case 'MC'
                  post_val = mean(scal(...
                    post_start_idx:...
                    post_end_idx));
                case 'MV'
                  [post_val, post_val_row] = max(mean(scal(:,...
                    post_start_idx:post_end_idx), 2));
              end
              %% Iterate over window comparison.
              for icomparison = 1:ncomparison
                % Updating results.
                results{jres, 1} = in_data.data{ic.path_no_ext}{istim};
                results{jres, 2} = in_data.data{ic.name_no_ext}{istim};
                results{jres, 3} = representational_space{ispace};
                switch calculation_type{itype}
                  case 'MC'
                    results{jres, 4} = [integration_constants(iconstant)...
                      integration_constants(jconstant)];
                  case 'MV'
                    results{jres, 4} = integration_constants(jconstant);
                end
                results{jres, 5} = calculation_type{itype};
                results{jres, 7} = post_target_window{iwindow};
                switch window_comparison{icomparison}
                  case 'abs'
                    results{jres, 6} = 'abs';
                    results{jres, 8} = post_val;
                  case 'rel'
                    % The pre-target attribute value.
                    pre_start_idx = find(tpts >= targ_on +...
                      pre_target_window(1), 1, 'first' );
                    pre_end_idx = find(tpts <= targ_on +...
                      pre_target_window(2), 1, 'last');
                    switch calculation_type{itype}
                      case 'MC'
                        pre_val = mean(scal(pre_start_idx:pre_end_idx));
                      case 'MV'
                        pre_val = mean(...
                          scal(post_val_row,pre_start_idx:pre_end_idx), 2);
                    end
                    results{jres, 6} = 'rel';
                    results{jres, 8} = post_val - pre_val;
                end
                jres = jres + 1;
              end % for iwindow
            end % for icomparison
          end % if ismember(calc_type{itype}, {'MC', 'MV'})
        end % for itype
      end % for ispace
    end % for jconstant
  end % for iconstant
end % for Istim

% Tidy up results table and give column titles.
results{1, 1} = 'Stimulus Path';
results{1, 2} = 'Stimulus Name';
results{1, 3} = 'Represetational Space';
results{1, 4} = 'Time Constants';
results{1, 5} = 'Calculation Type';
results{1, 6} = 'Window Comparison';
results{1, 7} = 'Post-Target Window';
results{1, 8} = 'Attribute Value';
results = results(1:jres-1, :);

end
