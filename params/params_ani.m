function params = params_ani(varargin)
% Returns a paramneter struct for use with calc_ani
%
% params = params_ani(varargin);
%
% This function accepts a variable number of tag/value pairs. If no
% argument is passed in, then a parameter structure with empty values
% is returned.
%
% Parameter structure for controlling ANI calculation using IPEM toolbox
%        start_time_sec: How far into the input signal to start
%                        processing ANI
%               dur_sec: The duration of the input signal to process
%         normalize_wav: Whether or not to normalize the input
%                        signal before processing.
%               aniPath: directory used for temporarily storing ANI
%                        calculations while calc_ani is processing.
%              PlotFlag: if non-zero, plots ANI (see IPEMCalcANI)
%      DownSampleFactor: integer factor by which ANI outcome is
%                        downsampled (see IPEMCalcANI)
%         NumOfChannels: number of channels (see IPEMCalcANI)
%              FirstCBU: frequency of first channel (in critical
%                        band units, see IPEMCalcANI)
%               CBUStep: frequency difference between channels (in CBU)
%
% Initialization values must be passed in as parameter/value pairs
%
% Copyright (c) 2006 The Regents of the University of California
% All Rights Reserved.
%
% Author(s):
% Stefan Tomic - first version
% 12/4/06 Petr Janata - made handling of varargin more dynamic

%this is the parameter list
%if you need to add more parameters, add them here
fields = {'start_time_sec','dur_sec','normalize_wav','normalize_maxVal',...
	  'aniPath','PlotFlag','DownSampleFactor',...
	  'NumOfChannels','FirstCBU','CBUStep'};

params = mkstruct(fields,varargin);
   
% If no path is specified in which to perform the ANI calculation, we have to
% add a unique path.  Otherwise, the calculations will occur in the default
% temp directory and run the risk of being clobbered if multiple jobs are
% running with default values
if isempty(params.aniPath)
  params.aniPath = create_ipem_tmpdir;
end

% set defaults if not otherwise specified
def.PlotFlag = 0;
def.DownSampleFactor = 4;
def.NumOfChannels = 40;
def.FirstCBU = 2.0;
def.CBUStep =0.5;
def.start_time_sec = [];
def.dur_sec = [];
def.aniPath = create_ipem_tmpdir;  % '/tmp/IPEM'
def.normalize_wav = 1;
def.normalize_maxVal = 0.95;

for ifld = 1:length(fields)
  if isempty(params.(fields{ifld}))
    params.(fields{ifld}) = def.(fields{ifld});
  end
end

end % params_ani
 
