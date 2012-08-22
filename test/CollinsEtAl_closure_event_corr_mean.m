function ecm = CollinsEtAl_closure_event_corr_mean(corr_timecourse, Fs,...
    rp, params)

% Copyright (c) 2007-2012 The Regents of the University of California
% All Rights Reserved.

% This function takes a correlation time course, a rhythm profile, and
% some other parameters as input. The rhythm profile is accessed for
% onset estimates of the final two events in a stimulus. The correlation
% timecourse is then averaged in a region of these timepoints (window
% specified in params), and the two correlations are returned.

% INPUT
%  corr_timecourse is a vector of the correlation between local and global
%   activations over time.
%  Fs is the sampling rate of the correlation time course.
%  rp is a rhythm profile, containing estimated onsets among other results.
%  params is a parameter structure obtained from CollinsEtAl_globals and
%   modified in CollinsEtAl_calc_attributes.

% Tom Collins, 2012.06.26.

% Sampling rates for toract and rp variables (not the same).
Fs_t = Fs;
% Resonator band for accessing penultimate and final event time samples.
iband = params.closure.resonator_band;
% nstim = size(jlmt_out.data{1}, 1);
rc = set_var_col_const(rp.vars);
Fs_r = rp.data{rc.params}.rp.Fs;
% Onset information.
onsetInfo = rp.data{rc.onsetInfo};
oc = set_var_col_const(onsetInfo.vars);
onSampsBB = onsetInfo.data{oc.onsetTimesSampsByBand};
event_times = onSampsBB{iband}(end-1:end)/Fs_r;
% Retrieve the time windows.
win = params.closure.post_event_window/1000;
idx_s = ceil(Fs_t*(event_times + win(1))); % Start indices.
idx_e = floor(Fs_t*(event_times + win(2))); % End indices.
nevents = size(onSampsBB{iband}, 1);
nsamp = size(corr_timecourse, 1);
iback = 1;
% If the analysis window extends beyond the end of the audio, revert to
% the latest events times such that this is not the case.
while idx_e(2) > nsamp && iback < nevents + 1
    % idx_e(2) = nsamp;
    event_times = onSampsBB{iband}(end-(iback+1):end-iback)/Fs_r;
    idx_s = ceil(Fs_t*(event_times + win(1))); % Start indices.
    idx_e = floor(Fs_t*(event_times + win(2))); % End indices.
    iback = iback + 1;
end
if iback > 1
    fprintf('Time window exceeds length of stimulus.\n')
    fprintf('Reverted to earlier events.\n')
end
% In another version, the output was stored in an array, with rows and
% columns for time constants, and the third dimension for penultimate/
% final events.
event_corr_mtx = nan(2, 2, 2);
for kev = 1:2 % Indices for penultimate and final events.
    if iback > nevents
        fprintf('Could not define a sensible time window.\n')
        event_corr_mtx(1, 2, kev) = 0;
    else
        event_corr_mtx(1, 2, kev) =...
            mean(corr_timecourse(idx_s(kev):idx_e(kev)));
    end
end
ecm = struct;
ecm.type = 'ecm';
ecm.vars = {'params' 'event_corr_mtx'};
ecm.data{1} = params.closure;
ecm.data{2} = event_corr_mtx;
ecm.report = struct;
ecm.meta = struct;

end
