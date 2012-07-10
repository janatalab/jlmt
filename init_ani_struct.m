function as = init_ani_struct(varargin)
% Initializes an ani data structure
%
% fieldnames:
%            params: []
%               sig: []
%                Fs: []
%       filterFreqs: []
%    start_time_sec: []
%           dur_sec: []
%
% Copyright (c) 2007-2012 The Regents of the University of California
% All Rights Reserved.
%
% Author(s):
% 04 Dec, 2006 - Stefan Tomic, Now accepts parameter,value tag pairs
% 12/09/06 Petr Janata - added params field for storing parameters that lead to the
%                given results
% March 2007 - Stefan Tomic, Changed data structure to comply with
%                     Ensemble data struct. No longer accepts arguments to initialize
%                     values of data (it doesn't seem like this was
%                     used anywhere).

as = ensemble_init_data_struct;
as.type = 'ani';
as.vars = {'params','sig','Fs','filterFreqs','start_time_sec','dur_sec'};
asCols = set_var_col_const(as.vars);

for iTag = 1:2:length(varargin)
  
  tag = varargin{iTag};
 
  if(~ismember(tag,as.vars))
    error(sprintf('Unknown tag: %s\n',tag));
  end
  
  as.data{asCols.(tag)} = varargin{iTag+1};
  
  
  
end
