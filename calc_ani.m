function outData = calc_ani(inData,ap)
% Calls the IPEM auditory nerve image calculation
%
% outData = calc_ani(inData,ap);
%
%
%
% INPUTS
% inData: an aud structure that complies with ensemble data
%         structs. The signal specified by the aud struct is passed
%         throught the IPEM auditory peripheral module
% ap:     a parameter structure for ani processing, see params_ani.m
%
% OUTPUTS
% outData: an ani data structure
%
%
% Copyright (c) 2007 The Regents of the University of California
% All Rights Reserved.
%
% Author(s):
% 12/04/06 Petr Janata - substantially rewrote code.
% 4/26/07  Stefan Tomic - edited to comply with ensemble data structs

outData = [];

%see if inData was a 'get_default_params' tag
if(ischar(inData) && strcmp(inData,'getDefaultParams'))
  if ~exist('ap','var'), ap = ''; end
  outData = getDefaultParams(ap);
  return
else
  aud = inData;
end

ani = init_ani_struct; 

aniCols = set_var_col_const(ani.vars);
audCols = set_var_col_const(aud.vars);

%populate missing params with defaults
ap = getDefaultParams(ap);

if ~isempty(aud.data{audCols.wvf})
  wvf = aud.data{audCols.wvf};
  
  %for compatibility with data cells
  if(iscell(wvf))
    wvf = wvf{1};
  end
  
  wavFs = aud.data{audCols.Fs};
  sigsize = size(wvf);
  
  %if start time and duration were set, then process only this
  %section of the waveform
  if(~isempty(ap.start_time_sec) & ~isempty(ap.dur_sec))
    startstop = [ap.start_time_sec, (ap.start_time_sec+ap.dur_sec)]*wavFs;
    check_startstop(startstop, sigsize);
    wvf = wvf(startstop,:);
  end
  
else
  warning('calc_ani: No waveform in aud structure');
  return;
end
   
%convert to mono if the waveform is multi-channel
if(size(wvf,2) > 1)
  wvf = mean(wvf,2);
end

%normalize the signal if this was set in params
if(ap.normalize_wav)
  wvf = wvf./max(abs(wvf))*ap.normalize_maxVal;
end
      
% Perform the actual calculation
try
  tmpDirStatus = check_dir(ap.aniPath,1);
  if(tmpDirStatus)
    error('Could not create temp directory for ANI calculation');
  end
  [ani.data{aniCols.sig},ani.data{aniCols.Fs},ani.data{aniCols.filterFreqs}] = IPEMCalcANI(...
      wvf, ...
      wavFs, ...
      ap.aniPath,...
      ap.PlotFlag, ...
      ap.DownSampleFactor,...
      ap.NumOfChannels, ...
      ap.FirstCBU, ...
      ap.CBUStep); 
catch
  fprintf('ANI calculation failed ...\n');
  if(exist('tmpDirStatus','var') & ~tmpDirStatus)
    fprintf('Cleaning up temp directory...');
    rmdir(ap.aniPath);
  end
  ani = [];
  return
end
    
if(isempty(ap.start_time_sec))
  ani.data{aniCols.start_time_sec} = 0;
else
  ani.data{aniCols.start_time_sec} = ap.start_time_sec;  
end

if(isempty(ap.dur_sec))
   end_time_sec = size(ani.data{aniCols.sig},2)/ani.data{aniCols.Fs};
   ani.data{aniCols.dur_sec} = end_time_sec - ...
       ani.data{aniCols.start_time_sec};
else
  ani.data{aniCols.dur_sec} = ap.dur_sec;
end

saveParams.ani = ap;
%save the params in the ani struct
ani.data{aniCols.params} = saveParams;

outData = ani;

return

  
function check_startstop(startstop,sigsize)
  
  if any(startstop) > sigsize(1);
    error('calc_ani: Specified ANI start time and/or duration out of bounds.');
    return;
  end

return
  


%makes sure that all necessary param fields are populated and
%non-empty. Sets empty or non-existent values to default values
function params = getDefaultParams(params)

def = params_ani;

fnames = fieldnames(def);

if(~exist('params','var'))
  params = def;
  return
else

  for ifld = 1:length(fnames)
    if (~isfield(params,fnames{ifld}) || isempty(params.(fnames{ifld})))
      params.(fnames{ifld}) = def.(fnames{ifld});
    end
  end

end

return
