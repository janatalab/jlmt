function results = analyze_toract_pca_results(data, params)
% Summarizes results of PCA on torus activations across a series of input
% stimuli
%
% Inputs:
%   data: a cell array of strings with paths to toract matfiles
%   params: parameters guiding this function

% 17Jul2019 Petr Janata - started script

nstim = length(data);

for istim = 1:nstim
  fpath = data{istim};
  
  % Parse out the stimulus name
  tok = regexp(fpath,'/([\w\.-&]+)/toract','tokens');
  stimname = tok{1}{1};
  fprintf('Obtaining pca data for stimulus %d/%d: %s\n', istim, nstim, stimname);

  % Load the data
  load(fpath,'toract')
  
  % Convert to a table
  toract = ensemble_datastruct2table(toract);
  
  % Figure out how many time constants we have activations for
  numTimeConstants = length(toract.time_constants);
  
  % Get the pca data structure
  pcadata = toract.pca;
  
  % Determine which data sources we calculated PCA on
  sources = pcadata.vars;
  numSources = length(sources);
  
  % Loop over sources and time constants, initializing the output variables
  % if necessary
  for isrc = 1:numSources
    currSrc = sources{isrc};
    currSrcData = pcadata.data{isrc};
    
    for itc = 1:numTimeConstants
      currdata = currSrcData{itc};
      currTCLabel = toract.labels{itc};
      
      if isempty(currdata)
        continue
      end
      
      % Determine how many PCs were required to explain the criterion
      % amount of variance
      cumPctVarExplained = cumsum(currdata.explained);
      
      numPCs = find(cumPctVarExplained >= params.crit.pctVarExplained, 1, 'first');
      
      % Store the result in an output structure
      results.(currSrc).(currTCLabel)(istim,1) = numPCs;
      
    end % for itc = 1:numTimeConstants
  end % for isrc = 1:length(sources)
end % for istim=

% Generate summary statistics and plots
figure(1), clf
for isrc = 1:numSources
  currSrc = sources{isrc};
  
  for itc = 1:numTimeConstants
    currTCLabel = toract.labels{itc};
    if ~isfield(results.(currSrc),currTCLabel)
      continue
    end
    currdata = results.(currSrc).(currTCLabel);
    
    subplot(numSources, numTimeConstants,(isrc-1)*numTimeConstants + itc)
    bins = 1:18;
    hist(currdata,bins)
    set(gca,'xlim',[0 max(bins)])
    title(strrep(sprintf('%s: %s', currSrc,currTCLabel),'_','\_'),'fontsize',18)
    xlabel('#PCs')
    ylabel('#songs')
    
    xoffset = 0.05;
    align = 'left';
    text(xoffset, 0.95, sprintf('Mean: %.1f', mean(currdata)), ...
      'horizontalalign',align, ...
      'units','normalized');
    text(xoffset, 0.90, sprintf('Min: %.1f', min(currdata)), ...
      'horizontalalign',align, ...
      'units','normalized');  
    text(xoffset, 0.85, sprintf('Max: %.1f', max(currdata)), ...
      'horizontalalign',align, ...
      'units','normalized');
  end
end

% Save the figure if desired
if isfield(params,'output') && isfield(params.output, 'savefig') && params.output.savefig
  if ~isfield(params.output,'figstub')
    error('Name of file to which to save the figure was not provided!')
  end
  
  figname = fullfile(params.paths.figures, params.output.figstub);
  fprintf('Saving figure to %s\n', figname);
  print(figname, '-depsc')
end

end