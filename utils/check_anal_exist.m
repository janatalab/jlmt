function [foundFile,firstViolation,violationReason] = check_anal_exist(dirPath, params)
% Finds a saved calculation with matching parameters in a specified directory
%
% fname = check_anal_exist(dirpath,params);
%
% Checks the param structure in all of the MAT files in the directory specified
% by dirpath to see if they match the params structure. If a match is found,
% the name of the file is returned in fname.
%
% This function is used to check whether or not a particular calculation has
% already been performed on a particular stimulus.
%
% Accepts two params: params.varName, and params.matchParams, where
% params.varName is the name of the variable in the mat files to
% check, and params.matchParams is the parameter structure
% containing the parameter values we want to find.
%
% Copyright (c) 2007-2012 The Regents of the University of California
% All Rights Reserved.
%
% Author(s):
% 2007 Stefan Tomic
% 10/9/2008, Stefan Tomic - moved compare_structs to its own function,
%                           added functionality to check for a
%                           parameter subset
% 01/20/2011, Petr Janata - greatly reduced search times by having function
%                           return in the case that not all matches are
%                           requested and also searching existing files
%                           backward from most recent, rationale being that
%                           most recent calculations reflect desired state
%                           of the art.
% 01Nov2011, Petr Janata - added ability to request analyses matching subsets of values.  
%                          See compare_structs() for more info.

foundFile = '';
firstViolation = {};
violationReason = '';

varName = params.varName;

if exist(dirPath,'dir')
  % a search directory has been given
  fullpath = fullfile(dirPath,'*.mat');
elseif exist(fileparts(dirPath),'dir')
  % a search directory with a file stub has been given
  fullpath = [dirPath '*.mat'];
  dirPath = fileparts(dirPath);
else
  firstViolation = {'fileStruct'};
  violationReason = sprintf('input path %s not found',dirPath);
  return
end

fileStruct = dir(fullpath);
nFiles = length(fileStruct);

if ~nFiles
  firstViolation = {'fileStruct'};
  violationReason = sprintf('no files found for %s',fullpath');
  return
end

try
  ignore = params.ignore;
catch
  ignore = {};
end

try
  params.subset;
catch
  params.subset = 0;
end

try
  params.allMatches;
catch
  params.allMatches = false;
end

try
	params.subsetIsOK;
catch
	params.subsetIsOK = {};
end

paramFind = params.paramFind;
paramFind = removeIgnored(paramFind,ignore);
foundFile = {};

for iFile = nFiles:-1:1

  filename = fileStruct(iFile).name;
  filePath = fullfile(dirPath,filename);
  loadVars = load(filePath);

  %the variable to check params against
  if(isfield(loadVars,varName))
    checkVar = loadVars.(varName);
    
    if(isfield(checkVar,'data') && isfield(checkVar,'vars') &&  ismember('params',checkVar.vars))
      varConst = set_var_col_const(checkVar.vars);
      checkParams = checkVar.data{varConst.params};
      
      %deal with the storage of the param struct into a single cell
      %there are a few cases when this was not done properly
      %so this should be compatible with both (struct or cell struct)
      if(iscell(checkParams))
        checkParams = checkParams{1};
      end
      
      checkParams = removeIgnored(checkParams,ignore);
      
			
      [paramsAreEqual,firstViolation,violationReason] = compare_structs(paramFind,checkParams, ...
				'values',1, ...
				'substruct',params.subset, ...
				'types',1, ...
				'ignore_fieldnames', ignore, ...
				'subsetIsOK', params.subsetIsOK);
      
      
      
      if(paramsAreEqual)
        if(params.allMatches)
          foundFile = {foundFile{:} filename};
        else
          foundFile = {filename};
          return
        end
      end
    end
    
  else %the variable wasn't found in the file, so paramsAreEqual is false
    paramsAreEqual = 0;
  end
  
end


return

function params = removeIgnored(params,ignore)


for idx = 1:length(ignore)

  ignoreString = ignore{idx};
  
  %find out if we are ignoring a lower level param field (a field
  %of a field). Right now, only two field levels are supported.
  dotLoc = strfind(ignoreString,'.');
 
  
  if ~isempty(dotLoc) & isfield(params,ignoreString(1:dotLoc-1))
    if(isfield(params.(ignoreString(1:dotLoc-1)), ignoreString(dotLoc+1:end)))
      params.(ignoreString(1:dotLoc-1)) = rmfield(params.(ignoreString(1:dotLoc-1)), ...
		       ignoreString(dotLoc+1:end));
    end
  else
    if(isfield(params,ignoreString))
      params = rmfield(params,ignoreString);
    end
  end
end


return




