% script to copy jlmt files from the private to the jlmt repository
% 
% this script will pull in external files to the private jlmt repository,
% then rsync between the private and public jlmt repositories
% 
% after this script runs, you need to inspect the changes that have been
% made and manually commit any changes that have been made to the private
% and public jlmt repositories
% 
% this file previously had a flag that automatically openned and inspected
% all .m files in the distribution. Going forward, this would require
% manual listing and maintenance of a list of all included .m files within
% this file. A much simpler approach (on mac and unix machines) is to run
% the following command in a terminal window from the top-level directory
% of jlmt:
% 
%       $> find . -regex '.*m$' -exec emacs -nw {} \;
% 
% if you prefer an editor other than emacs, you can replace 'emacs -nw'
% with the proper command for that editor. The above command will open each
% .m file under the current working directory, one at a time, and as you
% exit one file, the next file will open. If you change the editor from
% emacs to something else (vi maybe?) make sure to keep '{} \;' at the end
% of the command. {} is a place-holder for each filename found by 'find',
% and \; is required to specify the end of the argument to '-exec'

% FB 2012.05.21 - adapted from private/matlab/projects/rhythm_profiler/dist/rp_package_dist.m
% FB 2012.11.06 - adapted from private/matlab/jlmt/jlmt_package_dist.m
% FB 2013.03.18 - updated for use with private and public jlmt distros

%% initial variables
fprintf(1,'initializing variables\n');

GENERATE_DOCUMENTATION = 0;
RSYNC = 1; % 0=dry run, will list all files it would have copied, 1=rsync

[~,user] = unix('whoami');
user = regexprep(user,'\s','');

rsync_excludes = '--exclude=dist --exclude=deprecated';
m2html_ignore = {'dist/','deprecated/'};

%% paths
home_path = fullfile('/home',user);
priv_path = fullfile(home_path,'svn','private','matlab');
pub_path = fullfile(home_path,'svn','public','matlab');
jlmt_path = fullfile(home_path,'svn','jlmt');
priv_jlmt = fullfile(priv_path,'jlmt');

if isempty(which('btb'))
  path(path,genpath(fullfile(priv_path,'projects','rhythm_profiler')));
end
if isempty(which('listFilesOfType'))
  path(path,genpath(fullfile(priv_path,'projects','tonmodcomp')));
end

dstr = datestr(now(),30);
log_file = fullfile(home_path,'data',sprintf('jlmt_files_%s.log',dstr));
fid = fopen(log_file,'w');

%% these files will be pulled into the private JLMT repository
data_list = {...
    };
    
%     '/data2/tonmodcomp/bigand_etal_jep_2003/B_I_NTC/audio/B_I_NTC.wav',...
%     '/data2/tonmodcomp/bigand_etal_jep_2003/B_IV_NTC/audio/B_IV_NTC.wav',...
%     '/data2/tonmodcomp/marmel_and_tillmann_mp_2009/Riii0/audio/Riii0.mp3',...
%     '/data2/tonmodcomp/marmel_and_tillmann_mp_2009/Rvii0/audio/Rvii0.mp3',...
%     '/data2/tonmodcomp/marmel_etal_jep_2010/piano1MiF/audio/piano1MiF.wav',...
%     '/data2/tonmodcomp/marmel_etal_jep_2010/piano2BivF/audio/piano2BivF.wav',...
%     '/data2/tonmodcomp/marmel_etal_jep_2010/pure1AiF/audio/pure1AiF.wav',...
%     '/data2/tonmodcomp/marmel_etal_jep_2010/pure2AivF/audio/pure2AivF.wav',...
%     '/data2/tonmodcomp/marmel_etal_pp_2008/Gi/audio/Gi.wav',...
%     '/data2/tonmodcomp/marmel_etal_pp_2008/Giv/audio/Giv.wav',...
%     '/data2/tonmodcomp/tillmann_etal_cbr_2003/1RCinBB/audio/1RCinBB.wav',...
%     '/data2/tonmodcomp/tillmann_etal_cbr_2003/1RCinCB/audio/1RCinCB.wav',...
%     '/data2/tonmodcomp/tillmann_etal_cbr_2003/3UCinBB/audio/3UCinBB.wav',...
%     '/data2/tonmodcomp/tillmann_etal_cbr_2003/3UCinCB/audio/3UCinCB.wav',...
%     '/data2/tonmodcomp/tillmann_etal_jep_2003/C_I_baseline_cons/audio/C_I_baseline_cons.wav',...
%     '/data2/tonmodcomp/tillmann_etal_jep_2003/C_I_cons/audio/C_I_cons.wav',...
%     '/data2/tonmodcomp/tillmann_etal_jep_2003/C_IV_baseline_cons/audio/C_IV_baseline_cons.wav',...
%     '/data2/tonmodcomp/tillmann_etal_jep_2003/C_IV_cons/audio/C_IV_cons.wav',...
%     '/data2/tonmodcomp/tillmann_etal_jep_2008/CD/audio/CD.wav',...
%     '/data2/tonmodcomp/tillmann_etal_jep_2008/CDb/audio/CDb.wav',...
%     '/data2/tonmodcomp/tillmann_etal_jep_2008/CS/audio/CS.wav',...
%     '/data2/tonmodcomp/tillmann_etal_jep_2008/CSb/audio/CSb.wav',...
%     '/data2/tonmodcomp/tillmann_etal_jep_2008/CT/audio/CT.wav',...
%     '/data2/tonmodcomp/tillmann_etal_jep_2008/CTb/audio/CTb.wav',...
%     fullfile(pub_path,'colormaps','new_seismic'),...
%     '/usr/local/bin/mpg123',...
%     '/usr/local/bin/mp3info',...
%     '/data2/tonmodcomp/closure_work/closure_final/closure_struct_empr.mat',...
%     '/data2/tonmodcomp/closure_work/closure_final/closure_struct_hypo.mat',...
%     '/data2/tonmodcomp/closure_work/closure_final/closure_struct_smth.mat',...
%     '/data2/modulation/modpitchclass/pp2pitchclass_nnet_full_20120611T113
%     342.mat',...

%% these path replacements will be made on any files that are being pulled
% into the private repository
pathRep = {...
    '/data2/tonmodcomp','/data/CollinsEtAl';...
    };

%        fullfile(pub_path,'database'),'utils';...
% 	   fullfile(pub_path,'utils'),'utils';...
%        fullfile(pub_path,'database'),'utils';...
%        fullfile(priv_path,'analysis/tonality/metrics'),'utils';...
% 	   fullfile(priv_path,'data_types/event_objects'),'event_objects';...
%        fullfile(priv_path,'database_private'),'utils';...
% 	   fullfile(priv_path,'jlmt/rp_modules'),'proc/rp_modules';...
%        fullfile(priv_path,'jlmt/test'),'test';...
%        fullfile(priv_path,'jlmt/jlmt_startup.m'),'jlmt_startup.m';...
% 	   fullfile(priv_path,'jlmt/jlmt_proc_series.m'),'jlmt_proc_series.m';...
%        fullfile(priv_path,'jlmt/includes/mp3read_jlmt.m'),'utils/mp3read.m';...
% 	   fullfile(priv_path,'jlmt'),'proc';...
% 	   fullfile(priv_path,'projects/rhythm_profiler/params'),'proc/params';...
% 	   fullfile(priv_path,'projects/rhythm_profiler'),'.';...
%        fullfile(priv_path,'projects/tonmodcomp'),'utils';...
% 	   fullfile(priv_path,'utils/converters'),'utils';...
% 	   fullfile(priv_path,'utils/signal'),'utils';...
% 	   fullfile(priv_path,'utils/stats'),'utils';...
%      '/home/matrpcuser/svn/thirdparty/matlab/miditoolbox','midi';...
%        '/home/matrpcuser/svn/thirdparty/matlab/stats','utils';...
%        '/usr/local/bin','utils';...
%         ...%'audio','';...
%        '/data2/modulation/cont_popout','data/maps';...
%        '/data2/tonmodcomp/modorig/som','data/maps';...
%        '/data2/modulation/modpitchclass','data/maps';...
%        '/data2/tonmodcomp/closure_work/closure_final','data/CollinsEtAl';...
%        '/data2/tonmodcomp','data';...
%        '/data2/tontraj/orig_maps','data/maps';...
%        '/nfs/stimuli/audio/tension','data/test_jlmt';...
%        '/nfs/stimuli/audio/groove/polystream_move','data/test_jlmt';...

%% transfer files
nfiles = length(data_list);
if nfiles
  fprintf(fid,'copying %d files into %s\n',nfiles,priv_jlmt);
  fprintf(1,'copying %d files into %s\n',nfiles,priv_jlmt);
  
  dest_list = data_list;
  for iPath = 1:size(pathRep,1)
    dest_list = strrep(dest_list,pathRep{iPath,1},pathRep{iPath,2});  
  end
  dest_list = strcat(priv_jlmt,dest_list);

  for ifile = 1:nfiles
    check_dir(fileparts(dest_list{ifile}),0,1);
    copyfile(data_list{ifile},dest_list{ifile})
  end
end

if RSYNC
  fprintf(fid,'rsyncing: %s to %s\n',priv_jlmt,jlmt_path);
  fprintf(1,'rsyncing: %s to %s\n',priv_jlmt,jlmt_path);
  
  % rsync command
  rcmd = ['rsync -Cav ' rsync_excludes ' ' priv_jlmt ' ' fileparts(jlmt_path)];
else
  fprintf(fid,'rsync DRY RUN: %s to %s\n',priv_jlmt,jlmt_path);
  fprintf(1,'rsync DRY RUN: %s to %s\n',priv_jlmt,jlmt_path);
  
  % dry run command, which lists all files that would have been moved
  rcmd = ['rsync -Cavn ' rsync_excludes ' ' priv_jlmt ' ' fileparts(jlmt_path)];
end % if RSYNC
fprintf(fid,'executing: %s\n',rcmd);
[status,msg] = unix(rcmd);
if status
  fprintf(fid,['\n\n ---- RSYNC ERROR ---- \n\n']);
  fprintf(1,['\n\n ---- RSYNC ERROR ---- \n\n']);
end
fprintf(fid,['\n\n' msg '\n\n']);

%create documentation
if(GENERATE_DOCUMENTATION)
  fprintf(fid,'generating documentation\n');
  fprintf(1,'generating documentation\n');
  m2html('mfiles',jlmt_path,'ignore',m2html_ignore,...
      'htmldir',fullfile(jlmt_path,'documentation'),'recursive','on');
end

fprintf(fid,'DONE\n');
fprintf(1,'DONE\n');
fclose(fid);
