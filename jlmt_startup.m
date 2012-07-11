% startup file entries for the JLMT
% 
% If you do not currently have a 'startup.m' file in your
% /home/%username%/matlab directory (in Linux) or your
% C:\Documents\MATLAB\work directory (in Windows), or if that file exists
% and is empty, you can simply copy this file there and rename it
% 'startup.m'. It will then be run every time you open MATLAB, and it will
% set these paths and make the JLMT, IPEMToolbox, and mp3read available. If
% a 'startup.m' file exists and has paths in it already, you can simply
% copy the following lines into your existing 'startup.m' file and edit
% them accordingly. Either way, you must restart MATLAB for the changes to
% take effect. To test whether your paths have been set properly, execute:
%   >> which jlmt_proc_series.m
% This should return the path to your JLMT installation, with the 'proc'
% directory appended to it ('/path/to/jlmt/proc').
% 
% FB 2012.07.11

% % PATH TO JLMT : replace '/path/to/jlmt_dist' with the path to the JLMT
path(path,genpath('/path/to/jlmt_dist'));

% % PATH TO IPEM Toolbox : replace '/path/to/IPEMToolbox' with the path to
% % your IPEM Toolbox installation
path(path,genpath('/path/to/IPEMToolbox'));

% % PATH TO MP3READ : if you wish to process MP3 files, you need to
% % download 'mp3read' from the LabROSA web site:
% %     http://labrosa.ee.columbia.edu/matlab/mp3read.html
% % Then, uncomment the following line and change the path to add the
% % directory that holds your copy of mp3read:
% addpath('/path/to/directory/holding/mp3/read/');
