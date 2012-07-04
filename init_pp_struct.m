function pp = init_pp_struct(varargin)
% pp = init_pp_struct(varargin)
%
% periodicity pitch data structure
%
%     params: []
%        sig: []
%         Fs: []
%    periods: []
%    filtANI: []
%
% Varargin can be a tag/value list of fieldnames and associated values
%
% Copyright (c) 2007 The Regents of the University of California
% All Rights Reserved.
%
%
%  12/04/06 PJ
%  12/09/06 PJ - added params field for storing parameters that lead to the
%                given results
%  4/25/07  ST - reformatted to comply with ensemble data struct

pp = ensemble_init_data_struct;
pp.type = 'pp';
pp.vars = {...
      'params', ...
      'sig', ...
      'Fs', ...
      'periods', ...
      'filtANI'};



return
