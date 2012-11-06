% startup file entries for the JLMT
% 
% After checking out the JLMT repository, add a line to your startup.m
% file adding the top-level directory of the repository to the path, and
% then add this script on the following line and restart Matlab. This
% will set the paths necessary to run JLMT. You may also use the MATLAB user
% interface to add JLMT paths by choosing:
%     File -> Set Path
% then clicking "Add with Subfolders" and selecting your JLMT directory.
% Make sure to include your JLMT paths AFTER you include the paths to your
% installation of the IPEM toolbox.
% 
% FB 2012.11.06

if isempty(which('IPEMCalcANI'))
  error(['YOU MUST ADD IPEM TOOLBOX TO YOUR PATHS!! '...
	 'jlmt will not run if you do not have IPEM Toolbox installed']);
end

fprintf(1,'Adding JLMT paths ...');

jlmtpath = fileparts(which('jlmt_proc_series'));

addpath(fullfile(jlmtpath,'event_objects'));
addpath(fullfile(jlmtpath,'includes'));
addpath(fullfile(jlmtpath,'maps'));
addpath(fullfile(jlmtpath,'midi'));
addpath(fullfile(jlmtpath,'params'));
addpath(fullfile(jlmtpath,'plots'));
addpath(fullfile(jlmtpath,'proc'));
addpath(fullfile(jlmtpath,'rp_modules'));
addpath(fullfile(jlmtpath,'test'));
addpath(fullfile(jlmtpath,'utils'));

fprintf(1,'\n');

