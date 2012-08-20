function [in_data2, p] = CollinsEtAl_closure_probability(...
    tonmodcomp_params, varargin)

% This function retrieves the event correlations for a time-constant pair
% and calculates the probabilities.

% Tom Collins, 2012.06.26

tonmodcomp_params = tonmodcomp_globals;
restrict = [1:5 8 11];
% restrict = 12;
tonmodcomp_params.datasets = tonmodcomp_params.datasets(restrict);
% tc_restrict = {'tc_0p1000' 'tc_2'};
tc_restrict = {'tc_0p1000' 'tc_4'};


%% Get event correlations.
% if nargin >= 2
%     restrict = varargin{1};
%     tonmodcomp_params.datasets = tonmodcomp_params.datasets(restrict);
%     if nargin == 4
%         tc_restrict = varargin{2};
%     else
%         tc_restrict = [];
%     end
% end

in_data = ensemble_init_data_struct;
in_data.vars = {'path', 'filename'};
ic = set_var_col_const(in_data.vars);
in_data.data{ic.path} = {};
in_data.data{ic.filename} = {};
% Loop over datasets.
datasets = tonmodcomp_params.datasets;
nsets = length(datasets);
fprintf('Found %d datasets:\n', nsets);
for iset = 1:nsets
	fprintf('\n\nDataset %d\n', iset);
	fprintf('Identifier: %s\n', datasets(iset).id);
    % Determine file path
	location_list = {};
	location_root = fullfile(tonmodcomp_params.paths.project_root, tonmodcomp_params.datasets(iset).id);
	if isfield(datasets(iset), 'subsets')
        nsub = length(datasets(iset).subsets);
		for isub = 1:nsub
			fprintf('\tSubset %d: %s\n', isub, datasets(iset).subsets(isub).id);
			location_list{isub} = fullfile(location_root, tonmodcomp_params.datasets(iset).subsets(isub).id);
		end
	else
		location_list = {location_root};
	end
	for iloc = 1:length(location_list)
		location = location_list{iloc};
		tonmodcomp_params.paths.stim_root = location;
		% Now create a data structure that is compatible with ipem_proc_series, and
		% run ipem_proc_series.
		contents = dir(location);
		ndir = size(contents, 1);
		for idir = 1:ndir
			% If the name folder is one of the following, do nothing.
			if strcmp(contents(idir).name, '.') || ...
               strcmp(contents(idir).name, '..') || ...
               strcmp(contents(idir).name, '.DS_Store') || ...
               ~isempty(regexp(contents(idir).name, '.txt'));
				continue
			end
			% Get parts of path and filename.
			currPath = fullfile(location, contents(idir).name, 'audio');
			flist = listFilesOfType(currPath, {'wav','mp3'});
			if length(flist) > 1
				error('Found %d files in %s\n', length(flist), currPath);
			end
			in_data.data{ic.path}(end+1) = {currPath};
			strcat(location, '/', contents(idir).name, '/audio');
			in_data.data{ic.filename}(end+1) = flist;
		end % for idir=
	end % for iloc
end % for iset=

% Now restrict ipem_out to the stimulus categories of interest.
% Stimulus categories of interest:
cats_int = {{'Mediant' 'LeadingTone'} ...
    {'PianoTonicDull' 'PianoSubdominantDull' 'PureTonicA' ...
    'PureSubdominantA'} ...
    {'RelatedConsonantInB' 'RelatedConsonantInC' ...
    'UnrelatedConsonantInB' 'UnrelatedConsonantInC'} ...
    {'Tonic' 'Dominant' 'Subdominant' 'TonicBaseline' ...
    'DominantBaseline' 'SubdominantBaseline'} ...
    {'TonicConsonant' 'TonicBaselineConsonant' 'SubdominantConsonant' ...
    'SubdominantBaselineConsonant'} ...
    {'TonicNTC' 'SubdominantNTC'} ...
    {'Tonic' 'Subdominant'}};
% cats_int = {{'TonicA' 'TonicB' 'SubdominantA' 'SubdominantB'}};
nexp = size(cats_int, 2);
dc = set_var_col_const(in_data.vars);
nstim = size(in_data.data{1}, 2);
fNames = cell(1, nstim);
for istim = 1:nstim
    [~, stim_name, ~] = fileparts(in_data.data{dc.filename}{istim});
    fNames{istim} = stim_name;
end
rel_idx = zeros(nstim, 1);
idxi = 1; % Increment to create rel_idx.
for istim = 1:nstim
    iexp = 1;
    while iexp <= nexp
        ncats_int = size(cats_int{iexp}, 2);
        icat = 1;
        while icat <= ncats_int
            % Get contents of this category.
            cat_conts =...
                tonmodcomp_params.datasets(iexp).stimCategoryMembership.(cats_int{iexp}{icat});
            % See if current stimulus is a member of this category.
            [tf, loc] = ismember(fNames(istim), cat_conts);
            if tf
                % If its a member, include the index in rel_idx and move on
                % to next stimulus.
                rel_idx(idxi) = istim;
                idxi = idxi + 1;
                icat = ncats_int + 1;
                iexp = nexp + 1;
            else
                % Carry on looking in the next category.
                icat = icat + 1;
            end
        end
        iexp = iexp + 1;
    end
end
rel_idx = rel_idx(1:idxi-1);

in_data2 = in_data;
in_data2.data{1} = in_data.data{1}(rel_idx);
in_data2.data{2} = in_data.data{2}(rel_idx);
nstim = size(in_data2.data{1}, 2);
% Update fNames and get fPaths.
dc = set_var_col_const(in_data2.vars);
fNames = cell(1, nstim);
fPaths = cell(1, nstim);
for istim = 1:nstim
    [~, stim_name, ~] = fileparts(in_data2.data{dc.filename}{istim});
    fNames{istim} = stim_name;
    [path_str, ~, ~] = fileparts(in_data2.data{dc.path}{istim});
    fPaths{istim} = path_str;
end

% Iterate over stims and retrieve event correlations, placing them in an
% nstim x 2 matrix EC.
EC = zeros(nstim, 2);
for istim = 1:nstim
    contents = dir(fullfile(fPaths{istim}, 'ecm', 'ecm_20120702T16*'));
    % contents = dir(fullfile(fPaths{istim}, 'ecm', 'ecm_20120806T13*'));
    ecm = load(fullfile(fPaths{istim}, 'ecm', contents.name));
    ecm = ecm.ecmi;
    ec = set_var_col_const(ecm.vars);
    event_corr_mtx = ecm.data{ec.event_corr_mtx};
    t_const = ecm.data{ec.time_constants};
    tc_pair = ecm.data{ec.tc_pair};
    rel_idx = find(strcmp(tc_restrict(1), tc_pair(:,1)).*...
        strcmp(tc_restrict(2), tc_pair(:,2)));
    if isempty(rel_idx)
        error('Could not find time-constant pair amongst existing analysis')
    end
    itc = find(strcmp(tc_restrict(1), t_const));
    jtc = find(strcmp(tc_restrict(2), t_const));
    for kev = 1:2
        EC(istim, kev) = event_corr_mtx(itc, jtc, kev);
    end
end

%% Load the closure matrix and calculate probabilities of event
% correlations. Theoretical:
closure_struct = load(fullfile(tonmodcomp_params.paths.closure_struct,...
    'closure_struct-hypoth_tom_20120703T111311.mat'));
% closure_struct = load(fullfile(tonmodcomp_params.paths.closure_struct,...
%     'closure_struct-theor_tom_20120703T111311.mat'));
% Empirical, 0.1 and 4 s time constants:
% closure_struct = load(fullfile(tonmodcomp_params.paths.closure_struct,...
%     'closure_struct-0p1_4_empir_20120710T153125.mat'));
% Empirical, 0.1 and 2 s time constants:
% closure_struct = load(fullfile(tonmodcomp_params.paths.closure_struct,...
%     'closure_struct-0p1_2_empir_20120711T145528.mat'));
closure_struct = closure_struct.closure_struct;
closure_mtx = closure_struct.closure_mtx;
% imagesc(closure_mtx);
% axis square
bin_size = closure_struct.bin_size;
vals = -1:bin_size:1;
m = size(vals, 2);
p = zeros(nstim, 1);
for istim = 1:nstim
    r_idx = zeros(1, 2);
    for kev = 1:2
        r_raw = EC(istim, kev);
        [~, idx] = min(abs(r_raw - vals));
        r_idx(kev) = idx;
        % May need to correct bin index:
        if idx > 1 && r_raw <= vals(idx)
            r_idx(kev) = idx - 1;
        end
    end
    p(istim) = closure_mtx(m - r_idx(1), r_idx(2));
end

% Save p.
closure_vector = struct;
closure_vector.in_data = in_data2;
closure_vector.fNames = fNames;
closure_vector.fPaths = fPaths;
closure_vector.EC = EC;
closure_vector.p = p;
% Save closure vector.
f_name = ['closure_vector-tillmann_etal_ejocp_2006_0p1_4_', datestr(now,30)];
save(fullfile(tonmodcomp_params.paths.project_root, 'closure_work',...
    'closure_vector', f_name), 'closure_vector');


% Some plots for tillmann_etal_ejocp_2006
start_idx = 1;
end_idx = 48;
ncat = 4;
color_spec = {'b' 'c' 'g' 'y'};
title_str = ['tillman_etal_ejocp_2006, I_A = blue, I_B = cyan,'...
    ' IV_A = green, IV_B = yellow'];
figure
hold on
for icat = 1:ncat
    % plot(EC(12*(icat-1)+1:12*icat, 2), ...
    %     EC(12*(icat-1)+1:12*icat, 1), ['.' color_spec{icat}],...
    %     'MarkerSize', 12);
    plot(p(12*(icat-1)+1:12*icat), color_spec{icat});
end
hold off
grid on
title(title_str, 'Interpreter', 'None')
xlim([0.5 1])
ylim([0.5 1])
axis square
xlabel('r_n')
ylabel('r_{n-1}')


% Some plots for tillmann_etal_jep_2008
start_idx = 136;
end_idx = 207;
ncat = 6;
color_spec = {'r' 'k' 'g' 'k' 'b' 'c'};
title_str = ['tillman_etal_jep_2008, I = blue, V = red, IV = green,'...
    ' base = black'];

figure
hold on
for icat = 1:ncat
    plot(EC(start_idx+icat-1:ncat:end_idx, 2), ...
        EC(start_idx+icat-1:ncat:end_idx, 1), ['.' color_spec{icat}],...
        'MarkerSize', 12);
    % plot(p(start_idx+icat-1:ncat:end_idx), color_spec{icat});
end
hold off
grid on
title(title_str, 'Interpreter', 'None')
xlim([0.5 1])
ylim([0.5 1])
axis square
xlabel('r_n')
ylabel('r_{n-1}')


% 
% 
% % Some plots for tillmann_etal_jep_2003
start_idx = 208;
end_idx = 255;
ncat = 4;
color_spec = {'k' 'g' 'k' 'b'};
title_str = 'tillman_etal_jep_2003, I = blue, IV = green, base = black';

% 
% Some plots for bigand_etal_jep_2003
start_idx = 256;
end_idx = 279;
ncat = 2;
color_spec = {'g' 'b'};
title_str = 'bigand_etal_jep_2003, I = blue, IV = green';



print('-dpsc2', '-append', fullfile('/data2', 'tonmodcomp',...
    'closure_work', 'closure_struct',...
    'event_correlations_tc_0p1_2_win_100_300.ps'));

end