function tmpdir = create_jlmt_tmpdir(tmpdir)
% Creates a temporary directory in which to put intermediate files created by the JLMT
% 
%   tmpdir = create_jlmt_tmpdir(tmpdir);
%
% If no directory name is passed in as input, a subdirectory is created in
% /tmp/jlmt with a unique identifier.
%
% It is absolutely necessary to use these temporary directories if using the
% JLMT in a multi-user environment. Otherwise there is the risk of
% different jobs clobbering each other.
%
% Copyright (c) 2006-2012 The Regents of the University of California
% All Rights Reserved.
%
% Author(s):
% 12/10/06 Petr Janata
% 2012.07.05 FB - updated for JLMT

% Make sure we have the basic root directory
dir_root = '/tmp/jlmt/';
check_dir(dir_root);

if nargin < 1
  tmpdir = fullfile(dir_root,sprintf('tmp%s', num2hex(datenum(now))));
end

if check_dir(tmpdir)
  error(sprintf('Could not create directory: %s', tmpdir))
end
