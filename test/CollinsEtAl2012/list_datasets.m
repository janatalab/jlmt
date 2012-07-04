function list_datasets(datasets)

nsets = length(datasets);

fprintf('Found %d datasets:\n', nsets);
for iset = 1:nsets
	fprintf('\n\nDataset %d\n', iset);
	fprintf('Identifier: %s\n', datasets(iset).id);
	fprintf('Description: %s\n', datasets(iset).description);
	
    if isfield(datasets(iset), 'subsets')
        nsub = length(datasets(iset).subsets);
        if nsub
            for isub = 1:nsub
                fprintf('\tSubset %d: %s\n', isub, datasets(iset).subsets(isub).id);
            end
        end
    end
end