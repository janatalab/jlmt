function ecm = CollinsEtAl_closure_event_corr_mean(ipem_out, in_data,...
    params, save_ecm)
% function ipem_out_all = closure_event_corr(ipem_out, in_data, params,...
%     save_ecm)

% This function takes a cell argument output by the function
% closure_run_ipem.m. The cell contains the toract variables of interest
% and the rhythm profiles of each training stimulus. The rhythm profile is
% accessed to estimate the timepoints of the final two events in a
% stimulus. The toracts are then retrieved for these timepoints. For
% stored time-constant pairs, the toract correlations are then calculated
% at these timepoints, appended to the cell argument, and returned.

% TC 2012.06.19

% Column names.
ic = set_var_col_const(ipem_out.vars);
tc = set_var_col_const(ipem_out.data{ic.toract}{1}.vars);
% Cell containing the time constants used in analysis.
t_const = ipem_out.data{ic.toract}{1}.data{tc.labels};
% Sampling rates for toract and rp variables (not the same).
Fs_t = ipem_out.data{ic.toract}{1}.data{tc.params}.toract.Fs;
% Resonator band for accessing penultimate and final event time samples.
% (Not sure how to pick a band from which to extract onsets.)
iband = params.closure_matrix.resonator_band;
nstim = size(ipem_out.data{1}, 1);

% Get a local variable of ordered fNames.
% if params.closure_matrix.annotated_onsets
dc = set_var_col_const(in_data.vars);
fNames = cell(1, nstim);
fPaths = cell(1, nstim);
for istim = 1:nstim
    [~, stim_name, ~] = fileparts(in_data.data{dc.filename}{istim});
    fNames{istim} = stim_name;
    [path_str, ~, ~] = fileparts(in_data.data{dc.path}{istim});
    fPaths{istim} = path_str;
end

% Iterate over stimuli.
ecm = cell(nstim, 1); % Variable to hold event correlation matrices.
for istim = 1:nstim
    % Access rhythm profile to estimate the timepoints of the final two
    % events.
    rp = ipem_out.data{ic.rp}{istim};
    rc = set_var_col_const(rp.vars);
    Fs_r = rp.data{rc.params}.rp.Fs;
    onsetInfo = rp.data{rc.onsetInfo};
    oc = set_var_col_const(onsetInfo.vars);
    onSampsBB = onsetInfo.data{oc.onsetTimesSampsByBand};
    event_times = onSampsBB{iband}(end-1:end)/Fs_r;
    % event_times = [3.720 4.340];
    % Retrieve the time windows.
    win = params.closure_matrix.win/1000;
    idx_s = ceil(Fs_t*(event_times + win(1))); % Start indices.
    idx_e = floor(Fs_t*(event_times + win(2))); % End indices.
    toract = ipem_out.data{ic.toract}{istim};
    toract_corr = toract.data{tc.toract_corr};
    tcc = set_var_col_const(toract_corr.vars);
    tc_pair = toract_corr.data{tcc.tc_pair};
    % This is a p x q matrix, where p is the number of time-constant pairs
    % and q is the number of samples:
    toract_corr_data = toract_corr.data{tcc.toract_corr};
    nsamp = size(toract_corr_data, 2);
    nevents = size(onSampsBB{iband}, 1);
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
        fprintf('Time window exceeds length of stim %s.\n', fNames{istim})
        fprintf('Reverted to earlier events.\n')
    end
    % Store the penultimate and final correlations in an array, with rows
    % and columns for time constants, and the third dimension for
    % penultimate/final.
    ntc = size(t_const, 2);
    event_corr_mtx = nan(ntc, ntc, 2);
    for kev = 1:2 % Indices for penultimate and final events.
        irow = 1; % Increment over the rows of toract_corr_data.
        for itc = 1:ntc-1 % Time constant 1.
            for jtc = itc+1:ntc % Time constant 2.
                if iback > nevents
                    fprintf('Could not define a sensible time window.\n')
                    event_corr_mtx(itc, jtc, kev) = 0;
                else
                    event_corr_mtx(itc, jtc, kev) =...
                        mean(toract_corr_data(irow,...
                        idx_s(kev):idx_e(kev)));
                end
                irow = irow + 1;
            end
        end
    end
    ecmi = struct;
    ecmi.type = 'ecm';
    ecmi.vars = {'params' 'time_constants' 'tc_pair'...
        'event_corr_mtx'};
    ecmi.data{1} = params.closure_matrix;
    ecmi.data{2} = t_const;
    ecmi.data{3} = tc_pair;
    ecmi.data{4} = event_corr_mtx;
    ecmi.report = struct;
    ecmi.meta = struct;
    % Save the ecm variable.
    if save_ecm
        save_dir = fullfile(fPaths{istim}, 'ecm');
        if exist(save_dir, 'dir') ~= 7
            mkdir(save_dir);
        end
        f_name = ['ecm_', datestr(now, 30)];
        ecmi.name = f_name;
        ecm{istim} = ecmi;
        save(fullfile(save_dir, [f_name '.mat']), 'ecmi');
    end
end

% Define the output variable.
ipem_out_all = ipem_out;
nvars = size(ipem_out_all.data, 2);
ipem_out_all.vars{nvars+1} = 'ecm';
ipem_out_all.data{nvars+1} = ecm;

end
