function ani = ani_paramGroups
% Returns ani parameter structure commonly used for BTB algorithm
%
% Copyright (c) 2007-2013 The Regents of the University of California
% All Rights Reserved
%
% Author(s):
% Stefan Tomic

% 13Nov2012 PJ - added initial call to params_ani to populate with global
% ani defaults

ani = params_ani;

ani.PlotFlag = 0;
ani.DownSampleFactor = 5;
ani.normalize_wav = 1;
ani.start_time_sec = [];
ani.dur_sec = [];