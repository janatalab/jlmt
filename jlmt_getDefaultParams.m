function params = jlmt_getDefaultParams(params,calc_names)

% return default params not otherwise specified for the given calc steps
% 
%   params = jlmt_getDefaultParams(params,calc_names)
% 
% REQUIRES
%   params (optional) - a structure containing parameter settings that you
%       want to assign specific (non-default) values to. These will be
%       passed as-is to the output
%   calc_names - this can be a string indicating one calc step to return
%       default parameters for, or it can be a cell array of strings
%       containing names of multiple calc steps to return default
%       parameters for. 
% 
% RETURNS
%   If calc_names contains a string or a cell array of string containing
%   only one calc name, the returned params will be the direct params for
%   that calc step. If params are requested for multiple calc steps at a
%   time, the return variable will contain one fieldname per requested calc
%   step, and the default params for that calc step will be contained
%   within.
% 
%       examples:
% 
%           jlmt_getDefaultParams('','ani') = this will return a params
%               struct containing the output of the params_ani() function
% 
%           jlmt_getDefaultParams('',{'ani','pp','li'}) = this will return
%               a structure with the fieldnames 'ani','pp', and 'li', each
%               containing default params for their respective calc step
% 
% FB 2012.07.03

def.glob.force_recalc = {};
def.glob.ignore = {};
def.glob.save_calc = {'ani','pp','li','toract'};
def.glob.process = {{'ani','pp','li','toract'},...
    {'ani','pp','pc','li','toract'}};
def.glob.outputType = 'data';

if ~iscell(calc_names)
  calc_names = {calc_names};
end % if ~iscell(calc_names

for k=1:length(calc_names)
  fh = parse_fh(['params_' calc_names{k}]);
  if isstruct(params) && isfield(params,calc_names{k})
    def.(calc_names{k}) = fh('getDefaultParams',params.(calc_names{k}));
  else
    def.(calc_names{k}) = fh('getDefaultParams');
  end
end % for k=1:length(calc_names

if (nargin)
  %only assign defaults to the processes that we are performing
  if(isfield(params,'glob') && isfield(params.glob,'process'))
    defTypes = [params.glob.process{:} {'glob'}];
  else
    defTypes = fieldnames(def);    
  end
    
  for iType = 1:length(defTypes)

    type = defTypes{iType};  

    if ~isfield(params,type)
      params.(type) = struct();
    end
    
    %see if there are any fields in params that are not in
    %defs. This should be reported since it may have been caused by
    %a misspelling. All possible fields should also have defaults
    %(even if the default is simply empty)
    
    % should we be checking sub-structs of the main analysis parameters? If
    % we do this, then we need to make sure that for every toract metric
    % function we create, we update the defaults for params.metrics in
    % params_toract.m to include the defaults for the given metric. If we
    % do not want to do this, we will have to either not descend into
    % parameters for each ipem_proc_series function, or we will have to
    % find another solution to checking the defaults. - FB 01/10/2011
%     [validParams,badFields,reason] = compare_structs(params.(type),def.(type),'values',0,'substruct',0,'types',0);
    [validParams,badFields,~] = compare_structs(params.(type),def.(type),'values',0,'substruct',1,'types',0);
    if(~validParams)
      if iscell(badFields), badFields = cell2str(badFields,', '); end
      warning(['Field(s) %s for ''%s'' calculation is not specified in '...
			 'default parameters. Check spelling or add'...
			 ' the field(s) to the list of defaults\n\n'],badFields,type);
    end
    
    %populate missing param fields with defaults
    params.(type) = struct_union(params.(type),def.(type));

  end
  
else
  params = def;
end
