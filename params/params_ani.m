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
%            inDataType: specify the type of input for calc_ani (optional)
%            prev_steps: specify the analysis steps previous to calc_li
%                        that should be encountered when using the given
%                        parameters (optional)
%
% Initialization values must be passed in as parameter/value pairs
%
% Copyright (c) 2006-2012 The Regents of the University of California
% All Rights Reserved.
%
% Author(s):
% Stefan Tomic - first version
% 12/4/06 Petr Janata - made handling of varargin more dynamic

%this is the parameter list
%if you need to add more parameters, add them here
fields = {'start_time_sec','dur_sec','normalize_wav','normalize_maxVal',...
	  'aniPath','PlotFlag','DownSampleFactor',...
	  'NumOfChannels','FirstCBU','CBUStep','prev_steps','inDataType'};

params = mkstruct(fields,varargin);
   
% set defaults if not otherwise specified
def.PlotFlag = 0;
def.DownSampleFactor = 4;
def.NumOfChannels = 40;
def.FirstCBU = 2.0;
def.CBUStep =0.5;
def.start_time_sec = [];
def.dur_sec = [];
def.aniPath = create_jlmt_tmpdir;  % '/tmp/jlmt'
def.normalize_wav = 1;
def.normalize_maxVal = 0.95;
def.prev_steps = [];
def.inDataType = [];

for ifld = 1:length(fields)
  if isempty(params.(fields{ifld}))
    params.(fields{ifld}) = def.(fields{ifld});
  end
end

end % params_ani
 
