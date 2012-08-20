% Analysis of example files from Collins et al 201X tonality model comparisons
% 
% 

%% initialize parameters
params = CollinsEtAl_globals;

% Initialize structure that will be passed to jlmt_proc_series
in_data = ensemble_init_data_struct;
in_data.vars = {'path', 'filename', 'path_no_ext', 'name_no_ext'};
ic = set_var_col_const(in_data.vars);
in_data.data{ic.path} = {};
in_data.data{ic.filename} = {};
in_data.data{ic.path_no_ext} = {};
in_data.data{ic.name_no_ext} = {};

%% iterate over datasets, collect stimuli to analyze
datasets = params.datasets;
nsets = length(datasets);
for iset = 1:nsets
	fprintf('\n\nDataset %d\n', iset);
	fprintf('Identifier: %s\n', datasets(iset).id);

	% Determine file path
	location_list = {};
	location_root = fullfile(params.paths.data_root, params.datasets(iset).id);
	
	if isfield(datasets(iset), 'subsets')
        nsub = length(datasets(iset).subsets);
		for isub = 1:nsub
			fprintf('\tSubset %d: %s\n', isub, datasets(iset).subsets(isub).id);
			location_list{isub} = fullfile(location_root, params.datasets(iset).subsets(isub).id);
		end
	else
		location_list = {location_root};
	end
	
	for iloc = 1:length(location_list)
		location = location_list{iloc};
		params.paths.stim_root = location;
		
		% Now create a data structure that is compatible with jlmt_proc_series, and
		% run jlmt_proc_series.
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

%% run jlmt_proc_series
jlmt_out = jlmt_proc_series(in_data, params.jlmt);

%% calculate metrics
% There are parameters within this script that you may want to experiment with.
CollinsEtAl_metric_compile(in_data, params.jlmt, jlmt_out); 
