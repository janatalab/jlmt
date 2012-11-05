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
% and uncomment the following line
% path(path,genpath('/path/to/jlmt_dist'));
jlmtpath = fileparts(which('jlmt_proc_series'));
addpath(fullfile(jlmtpath,'maps'));
addpath(fullfile(jlmtpath,'params'));
addpath(fullfile(jlmtpath,'includes'));
addpath(fullfile(jlmtpath,'rp_modules'));

% % PATH TO IPEM Toolbox : replace '/path/to/IPEMToolbox' with the path to
% % your IPEM Toolbox installation and uncomment the following line
%path(path,genpath('/path/to/IPEMToolbox'));

