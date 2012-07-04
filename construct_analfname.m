function analfname = construct_analfname(fname,anal_stub)
% Constructs an analysis filename that contains a unique hash code
%
% analfname = construct_analfname(fname,anal_stub);
%
% Takes a filepath and name in fname and generates the name of an output file
% (currently only a .mat file) that will be placed in a sibling directory that
% has the name in anal_stub. The sibling directory is created if it doesn't exist.
%
% Copyright (c) 2006 The Regents of the University of California
% All Rights Reserved.
%
% Author:
% 12/6/06 Petr Janata
% 08Nov2011 PJ - fixed handling in the case that fname is a directory path
% to the stimulus directory

if isdir(fname)
	fpath = fname;
	[~,fstub] = fileparts(fname);
	analdir = fullfile(fpath, anal_stub);
else
  [fpath,fstub] = fileparts(fname);
	analdir = fullfile(fileparts(fpath),anal_stub);
end

check_dir(analdir);
hashstr = sprintf('_%s',num2hex(datenum(now)));
analfname = fullfile(analdir,sprintf('%s%s.mat',fstub,hashstr));
end
