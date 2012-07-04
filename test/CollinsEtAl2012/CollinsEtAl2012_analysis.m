% This script runs jlmt proc series on a selection of audio files from
% Collins, Tillmann, Delbe, Barrett, and Janata (2012). "title". journal.
% 
% Authors:
% October 2011, Tom Collins

% Load global analysis parameters
params = CollinsEtAl2012_globals;

list_datasets(params.datasets);  % list registered datasets

% Initialize structure that will be passed to ipem_proc_series
in_data = ensemble_init_data_struct;
in_data.vars = {'path', 'filename'};
ic = set_var_col_const(in_data.vars);
in_data.data{ic.path} = {};
in_data.data{ic.filename} = {};

% Loop over datasets and feed stuff into IPEM_proc_series
datasets = params.datasets;
nsets = length(datasets);

fprintf('Found %d datasets:\n', nsets);
for iset = 1:nsets
	fprintf('\n\nDataset %d\n', iset);
	fprintf('Identifier: %s\n', datasets(iset).id);

	% Determine file path
	location_list = {};
	location_root = fullfile(params.paths.project_root, params.datasets(iset).id);
	
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

% Run ipem_proc_series.
ipem_out = ipem_proc_series(in_data, params.ipem);

% calculate metrics
metric_readout_compile;

