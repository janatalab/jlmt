function ci = init_ci_struct()
% contextuality index data structure, for use with calc_context2.m
% 
%   ci = init_ci_struct()
%
% RETURNS
%   params:     []
%   Fs: []
%   periods:    []
%   names:      {}
%   signals:    {}
%
% Copyright (c) 2007-2013 The Regents of the University of California
% All Rights Reserved.
%
%  12/04/06 PJ
%  12/09/06 PJ - added params field for storing parameters that lead to the
%                given results
%  4/25/07  ST - reformatted to comply with ensemble data struct
% 2010.02.16 FB - adapted from init_ci_struct, to handle multiple signals
% calculated, as well as names that can be associated with each signal

ci = ensemble_init_data_struct;
ci.type = 'ci';
ci.vars = {...
      'params', ...
      'Fs',...
      'periods',...
      'names',...
      'signals',...
      };

return
