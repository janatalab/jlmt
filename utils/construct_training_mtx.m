function [training_mtx, Fs] = construct_training_mtx(datalist,params)
% Samples data vectors from the list of inputs in datalist and
% concatenates them into an output matrix. The items in datalist should
% typically be the leaky-integrated images obtained from periodicity pitch
% images.
%
% params is a structure that specifies which data should be extracted from
% the input data structure.
%
% params.representation - The representational space from which to extract
%                         the data. default='li' 
% params.time_constant - The leaky integration time constant to use.
%                        default=2
% params.prop_of_data - The proportion of data to sample from the file.
%                       default = 1.0
% params.shuffle - if false prop_of_data == 1 then this results in a training
%                  matrix with veridically ordered samples of the entire
%                  input. if false and prop_of_data < 1, then the correct
%                  number of samples is read, starting at the beginning.
%                  default = true  

% 07Sep2019 PJ - wrote initial script

training_mtx = [];
Fs = nan;

% Determine whether we are dealing with a list of files or structures
if ~iscell(datalist) && ~ischar(datalist{1})
  error('Only cell arrays of filenames are currently supported');
end

targSpace = params.representation;
targTimeConstant = params.time_constant;

% Determine how many input files we have
numFiles = length(datalist);

% Loop over the files and extract data
for ifile = 1:numFiles
  currFile = datalist{ifile};
  
  % Load the file
  currdata = load(currFile);
  
  % Make sure we have the target space
  if ~isfield(currdata, targSpace)
    error('Desired representational space (%s) missing', targSpace)
  end
  
  % Extract the data for the representational space
  repspace = currdata.(targSpace);
  cols = set_var_col_const(repspace.vars);
  
  if ~params.shuffle
    if isnan(Fs)
      Fs = repspace.data{cols.Fs};
    else
      if Fs ~= repspace.data{cols.Fs}
        error('Sampling rate (%.3f) does not match previous sampling rate (%.3f)',repspace.data{cols.Fs}, Fs)
      end
    end
  end
  
  % Make sure we have the desired time constant
  timeConstantStr = get_tc_names(targTimeConstant);
  tcIdx = find(ismember(timeConstantStr, repspace.data{cols.names}), 1);
  if isempty(tcIdx)
    error('Signal for desired time constant (%.2f) not found', targTimeConstant)
  end
  
  % Extract the signal matrix
  signal = repspace.data{cols.signals}{tcIdx};
  
  % Determine how many points we want to sample
  numOrigSamps = size(signal,2);
  nsamp = fix(numOrigSamps*params.prop_of_data);
  
  % Randomly sample this many points and concatenate with the training
  % matrix
  if params.shuffle
    idxs = randperm(numOrigSamps);
  else
    idxs = 1:numOrigSamps;
  end
  
  training_mtx = [training_mtx signal(:,idxs(1:nsamp));];
end % for ifile

return