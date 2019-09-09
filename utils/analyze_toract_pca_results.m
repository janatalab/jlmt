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
      results.(currSrc).(currTCLabel).numPCs(istim,1) = numPCs;
      
      % Calculate the mean score for each of the PCs
      results.(currSrc).(currTCLabel).mean_scores{istim} = mean(currdata.score(:,1:numPCs));
          
    end % for itc = 1:numTimeConstants
  end % for isrc = 1:length(sources)
end % for istim=

end