function pp = jlmt_preproc_params(project_name)
%  pp = jlmt_preproc_params;
%
%  Initializes various pre-processing parameters for the JLMT
%
% Copyright (c) 2007-2013 The Regents of the University of California
% All Rights Reserved.
%
% Authors: Petr Janata, Stefan Tomic, Fred Barrett

if ~exist('project_name','var')
  project_name = '';
end

switch project_name
  case {'modulation','modulation_orig'}
    pp.notes_per_measure = 12;
    pp.s_per_note = 0.2;
    pp.wav_chunk_s = pp.notes_per_measure*pp.s_per_note;  
    pp.chunk_ramp_ms = 5;
    pp.nchan_ani = 40;
    pp.downsample_factor = 4;
    pp.normalize_wav = 1;
    pp.nwin_ani_chunk_overlap = 10;
    
    % periodicity pitch parameters
    pp.nwin_per_measure = 63;
    pp.frame_width = pp.wav_chunk_s/pp.nwin_per_measure;
    pp.periodpitch_prop_overlap = 0; % only guaranteed to work for zero
    pp.frame_stepsize = pp.frame_width * (1-pp.periodpitch_prop_overlap);

    % Parameters for calculating chord and tonal context images
    switch project_name
      case 'modulation_orig'
        pp.inHalfDecayChords = 0.2;   % 200 ms decay time for forming chord images (one note)
        pp.inHalfDecayToneCenters = 2;  % decay time (s) of tone center images
      otherwise
        pp.inHalfDecayTimes = [0.2 2 10];
    end
 case 'rhythm'
    pp.notes_per_measure = 12;
    pp.s_per_note = 0.2;
    pp.measures_per_chunk = 1;
    pp.wav_chunk_s = pp.measures_per_chunk*pp.notes_per_measure*pp.s_per_note;  
    pp.chunk_ramp_ms = 5;
    pp.nchan_ani = 40;
    pp.downsample_factor = 5;
    pp.normalize_wav = 1;
    pp.nwin_ani_chunk_overlap = 10;
    
     
  otherwise
    pp.notes_per_measure = 0;
    pp.s_per_note = 0;
    pp.wav_chunk_s = pp.notes_per_measure*pp.s_per_note;  
    pp.chunk_ramp_ms = 5;
    pp.nchan_ani = 40;
    pp.downsample_factor = 4;
    pp.nwin_ani_chunk_overlap = 10;
    
    % periodicity pitch parameters
    pp.nwin_per_measure = 1;
    pp.frame_width = pp.wav_chunk_s/pp.nwin_per_measure;
    pp.periodpitch_prop_overlap = 0; % only guaranteed to work for zero
    pp.frame_stepsize = pp.frame_width * (1-pp.periodpitch_prop_overlap);

    % Parameters for calculating chord and tonal context images
%     pp.inHalfDecayChords = 0.2;   % 200 ms decay time for forming chord images (one note)
%     pp.inHalfDecayToneCenters = 2;  % decay time (s) of tone center images
    pp.inHalfDecayTimes = [0.2 2 10];

end
