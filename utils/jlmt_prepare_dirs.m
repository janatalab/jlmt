function flist = jlmt_prepare_dirs(fpath,params)
% Given a list of audio files, determines whether the proper directory exists for each file, 
% and if not, creates the directory and the requisite sub-directory (by audio
% file type), and moves the audio file to that directory.
%
% flist = jlmt_prepare_dirs(fpath, params);
%
% 

% 29Oct2012 Petr Janata
% 2013.02.06 FB - previously only checked to see if fpath pointed to a file
%                 in its nested position. This script now also checks to see
%                 if fpath points to a file that had been nested, but where
%                 fpath indicated its pre-nested position. In this case, it
%                 returns the proper nested fpath path. Also, given a
%                 directory, this will check the first-level sub-directories
%                 underneath to search for previously nested stimulus
%                 analysis directories.

% Check to see if path is a cell array of strings (paths to files), or a
% directory
flist = {};
if iscell(fpath)
  flist = fpath(:);
elseif isdir(fpath)
  
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

  % look for previously prepared/nested directories
  dirlist = dir(fpath);
  dirlist(~[dirlist.isdir]) = [];
  fprintf(1,'checking directories...\n');
  for idir=1:length(dirlist)
    if ~isempty(strmatch(dirlist(idir).name,{'.','..'}))
      continue
    end
    
    for ifmt=1:length(fmts)
      list = dir(fullfile(fpath,dirlist(idir).name,fmts{ifmt},...
          sprintf('*.%s',fmts{ifmt})));
      if ~length(list), continue, end
      flist = [flist; cellfun(@fullfile,...
          repmat({fullfile(dirlist(idir).name,fmts{ifmt})},1,...
          length(list)),{list.name}','UniformOutput',false)];
    end % ifmt=1:length(fmts
  end % for idir=1:length(dirlist
  % Convert flist to be full path
  flist = strcat(fpath,filesep,flist);
else
  flist = {fpath};
end % if isdir

% Loop over input files
nfiles = size(flist,1);
for ifile = 1:nfiles
  currFile = flist{ifile};
  [fpath,fstub,fext] = fileparts(currFile);
  
  % Check to see if this file is already properly nested
  [parentPath,parentStub] = fileparts(fpath);
  if strcmp(fext(2:end),parentStub)
    [~,grandparentStub,grandparentExt] = fileparts(parentPath);
    if strcmp(fstub, grandparentStub) || ...
            strcmp(fstub,[grandparentStub grandparentExt])
      continue
    end
  end
  
  % Check to see whether directory for this file exists
  targdir = fullfile(fpath,fstub);
  check_dir(targdir);
  
  % Check to see whether a sub-directory for the particular audio format
  % exists
  targdir = fullfile(targdir, fext(2:end));
  check_dir(targdir);

  if ~exist(fullfile(targdir,[fstub fext]))
    % Move the file over to the sub-directory
    fprintf('Moving file %s to %s ...\n', currFile, targdir);
    movefile(currFile,targdir)
  end % if ~exist(fullfile(targdir,currFile
  
  % Replace name in flist with new name
  flist{ifile} = fullfile(targdir,[fstub fext]);
end

return
