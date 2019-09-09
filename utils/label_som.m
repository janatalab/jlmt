function hitdata = label_som(trained_som, labeling_data, params)
% Labels the SOM defined in trained_som with the labels associated with
% input data provided in labeling_data.
%
% Inputs:
% trained_som: a structure returned by the somtoolbox containing
%              information about the SOM
% labeling_data: a structure array in which each element has the following
%                fields:
%   .data - matrix in which each row is an input
%   .label - the category label to with all of the rows in .data belong
%
% params:
%   .hit_type: one of 'kernel','crisp', or 'fuzzy'. default='kernel'

% 08Sep2019 Petr Janata - adapted from calc_som_loadings.m

% Determine how many categories we have
categories = {labeling_data.label};
numCategories = length(categories);

% Initialize the output variables
hits = nan(numCategories, size(trained_som.codebook,1));
prop_best = nan(numCategories,1);
best_unit = nan(numCategories,1);

for icat = 1:numCategories
  % Get the current category
  currCat = labeling_data(icat);
  
  % Get the number of observations
  nobs = size(currCat.data,1);
  
  % Convert the input data to an appropriate struct
  sD = som_data_struct(currCat.data);
  
  % Calculate the proportion of times each output unit was hit
  hits(icat,:) = som_hits(trained_som,sD,params.hit_type)'/nobs;
  
  % Get information about the proportion of time the best output unit was
  % activated, and its index
  [prop_best(icat),best_unit(icat)] = max(hits(icat,:));
end % for icat=

% Assign the new labels
new_labels = cell(size(trained_som.codebook,1),1);
new_labels(best_unit) = deal(categories);
hit_strength(best_unit) = prop_best;

% Find the best category for those units that haven't yet been assigned a
% category
switch params.hit_type
  case {'kernel','fuzzy'}
    unassigned_idx = find(~ismember(1:size(hits,2),best_unit));
    for iunit = 1:length(unassigned_idx)
      unit_idx = unassigned_idx(iunit);
      [hit_strength(unit_idx),best_cat] = max(hits(:,unit_idx));
      new_labels(unit_idx) = categories(best_cat);
    end
end

% Assign the labels to the output struct
trained_som.labels = new_labels;

% Add a couple of fields to the output struct
hitdata.categories = categories;
hitdata.hits = hits;
hitdata.hit_strength = hit_strength;
hitdata.labels = new_labels;

return