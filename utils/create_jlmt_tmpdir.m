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
% Copyright (c) 2006-2013 The Regents of the University of California
% All Rights Reserved.
%
% Author(s):
% 12/10/06 Petr Janata
% 2012.07.05 FB - updated for JLMT
% 2013.01.09 FB - added retry mechanism with random pauses between tries to accomodate
%   parallel computing environments. Previously, the same clock time was generated
%   by multiple jobs that were distributed, leading this script to try to create an
%   already existing tmpdir. The re-try mechanism places randomly spaced pauses. These
%   random pauses do not collide in parallel computing space. Rand() (without seeding)
%   will generate different random numbers for each worker.
% 2014.01.08 PJ - dealt with persistent issue of jlmt directory permissions
%   being set by another user by creating the jlmt directory within a user
%   directory within tmp

% Make sure we have the basic root directory

% Get the name of the user
[status,user] = system('whoami');
user = strtrim(user);

dir_root = fullfile('/tmp',user,'jlmt');
check_dir(dir_root,'parents',1);

tries = 10;
curtry = 1;
while curtry < tries
  if nargin < 1
    tmpdir = fullfile(dir_root,sprintf('tmp%s', num2hex(datenum(now))));
  end

  if check_dir(tmpdir)
    warning('Could not create directory: %s\nattempt %d/%d', tmpdir, curtry, tries)
    pause(rand);
  else
    return
  end
end % while curtry < tries

error('could not create directory: %s\n',tmpdir);
