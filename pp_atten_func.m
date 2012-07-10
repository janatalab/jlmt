function fh = pp_atten_func(atten_type)
% fh = pp_atten_func(atten_type);
%
% Returns a function handle to a particular type of attenuation function
%
% atten_type options:
%
% 'none' - returns empty; no attentuation will be applied
% 'ipem' - 1-(((1:x)/x - 0.5).^2)
%
% Also accepts any of the standard Matlab window names, e.g. 'kaiser',
% 'hanning', etc
%
%
% Copyright (c) 2007-2012 The Regents of the University of California
% All Rights Reserved.
%
% 12/10/06 Petr Janata

fh = [];
if nargin < 1
  return
end

switch lower(atten_type)
  case {'none', ''}
    fh = [];
    
  case 'ipem'
    fh = @(x) 1-(((1:x)/x - 0.5).^2);
    
  case 'ipem_squash_hf'
    fh = @ipem_squash_hf;
    
  otherwise
    try
      fh = str2func(atten_type);
    catch
      fprintf('Unknown attentuation or window function: %s\n', atten_type);
      fh = [];
    end
end

return

end

% ipem_squash_hf - forces the values in the bins corresponding to period=0 and
% period=1/Fs to zero as the results at these are not defined.
function win = ipem_squash_hf(x)
  win = 1-(((1:x)/x - 0.5).^2);
  win(1:2) = 0;
end
