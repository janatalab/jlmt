function ci = calc_ci(inData,par)
% JLMT wrapper for calc_li, which is a wrapper for IPEMLeakyIntegration()
% 
% this function passes its arguments directly to calc_li. See calc_li()
% for more information on this function. calc_ci is a wrapper to provide
% backwards compatibility with earlier versions of calc_li (calc_context)
%
% Copyright (c) 2012 The Regents of the University of California
% All Rights Reserved.
%
% FB 2012.07.12

if nargin < 2
  ci = calc_li(inData);
else
  ci = calc_li(inData,par);
end