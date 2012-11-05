function flist = jlmt_prepare_dirs(fpath,params)
% Given a list of audio files, determines whether the proper directory exists for each file, 
% and if not, creates the directory and the requisite sub-directory (by audio
% file type), and moves the audio file to that directory.
%
% flist = jlmt_prepare_dirs(fpath, params);
%
% 

% 29Oct2012 Petr Janata

% Check to see if path is a cell array of strings (paths to files), or a
% directory
flist = {};
if isdir(fpath)
  
  % If the present directory is specified, get the full path to it so that
  % we can check for proper nesting
  if ismember(fpath,{'.',['.' filesep]})
    fpath = pwd;
  end
  
  fprintf('Searching for audio files in %s ...\n', fpath);
  
  % Look for all .wav, .mp3, and .aif or .aiff files
  fmts = {'wav','mp3','aif','aiff','mat'};
  for ifmt = 1:length(fmts);
    cfmt = fmts{ifmt};
    fprintf('\t%s ...\n', cfmt);
    list = dir(fullfile(fpath,sprintf('*.%s', cfmt)));
    flist = [flist; {list.name}'];
  end % for ifmt
  
  % Convert flist to be full path
  flist = strcat(fpath,filesep,flist);
elseif iscell(fpath)
  flist = fpath(:);
else
  flist = {fpath};
end % if isdir

% Loop over input files
nfiles = size(flist,1);
for ifile = 1:nfiles
  currFile = flist{ifile};
  [fpath,fstub,fext] = fileparts(currFile);
  fext(1) = []; % remove dot
  
  % Check to see if this file is already properly nested
  [parentPath,parentStub] = fileparts(fpath);
  if strcmp(fext,parentStub)
    [~,grandparentStub] = fileparts(parentPath);
    if strcmp(fstub, grandparentStub)
      continue
    end
  end
  
  % Check to see whether directory for this file exists
  targdir = fullfile(fpath,fstub);
  check_dir(targdir);
  
  % Check to see whether a sub-directory for the particular audio format
  % exists
  targdir = fullfile(targdir, fext);
  check_dir(targdir);
  
  % Move the file over to the sub-directory
  fprintf('Moving file %s to %s ...\n', currFile, targdir);
  movefile(currFile,targdir)
  
  % Replace name in flist with new name
  flist{ifile} = fullfile(targdir,sprintf('%s.%s', fstub, fext));
  
end

return