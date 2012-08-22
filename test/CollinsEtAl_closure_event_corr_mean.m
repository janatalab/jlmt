function ecm = CollinsEtAl_closure_event_corr_mean(corr_timecourse,...
    Fs, rp, params)

% Copyright (c) 2007-2012 The Regents of the University of California
% All Rights Reserved.

% This function takes a cell argument output by the function
% jlmt_proc_series. The cell contains the toract variables of interest
% and the rhythm profiles of each training stimulus. The rhythm profile is
% accessed to estimate the timepoints of the final two events in a
% stimulus. The toracts are then retrieved for these timepoints. For
% stored time-constant pairs, the toract correlations are then calculated
% at these timepoints, appended to the cell argument, and returned.

% INPUT
%  corr_timecourse is .
%  rp is .
%  Fs is .
%  params is a parameter structure obtained from CollinsEtAl_globals.

% Tom Collins, 2012.06.26.


% Column names.
% tc = set_var_col_const(jlmt_out.data{ic.toract}{1}.vars);
% Cell containing the time constants used in analysis.
% t_const = jlmt_out.data{ic.toract}{1}.data{tc.labels};
% Sampling rates for toract and rp variables (not the same).
Fs_t = Fs;
% Resonator band for accessing penultimate and final event time samples.
iband = params.closure.resonator_band;
% nstim = size(jlmt_out.data{1}, 1);

% Get a local variable of ordered fNames. NOT NECESSARY?
% dc = set_var_col_const(in_data.vars);
% fNames = cell(1, nstim);
% fPaths = cell(1, nstim);
% for istim = 1:nstim
%     [~, stim_name, ~] = fileparts(in_data.data{dc.filename}{istim});
%     fNames{istim} = stim_name;
%     [path_str, ~, ~] = fileparts(in_data.data{dc.path}{istim});
%     fPaths{istim} = path_str;
% end

% Access rhythm profile to estimate the timepoints of the final two
% events.
rc = set_var_col_const(rp.vars);
Fs_r = rp.data{rc.params}.rp.Fs;
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
% Store the penultimate and final correlations in an array, with rows
% and columns for time constants, and the third dimension for
% penultimate/final.
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
% ecm.data{2} = t_const;
% ecm.data{3} = tc_pair;
ecm.data{2} = event_corr_mtx;
ecm.report = struct;
ecm.meta = struct;

end
