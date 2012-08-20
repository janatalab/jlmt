function p = CollinsEtAl_closure_probability(ecm, tc_pair,...
    distribution_str, params)

% Copyright (c) 2007-2012 The Regents of the University of California
% All Rights Reserved.

% This function retrieves the event correlations for a time-constant pair
% and calculates the closure probability for a specified distribution.

% Tom Collins, 2012.06.26

% Convert tc_pair to strings.
tc_pair = calc_li('calc_names', tc_pair);
ec = set_var_col_const(ecm.vars);
event_corr_mtx = ecm.data{ec.event_corr_mtx};
ecm_t_const = ecm.data{ec.time_constants};
ecm_tc_pair = ecm.data{ec.tc_pair};
rel_idx = find(strcmp(tc_pair(1), ecm_tc_pair(:,1)).*...
    strcmp(tc_pair(2), ecm_tc_pair(:,2)), 1);
if isempty(rel_idx)
    error('Could not find time-constant pair amongst existing analysis.')
end
itc = strcmp(tc_pair(1), ecm_t_const);
jtc = strcmp(tc_pair(2), ecm_t_const);
EC = zeros(1, 2);
for kev = 1:2
    EC(kev) = event_corr_mtx(itc, jtc, kev);
end

% Load the closure matrix and calculate probabilities of event
% correlations. Theoretical:
closure_struct = load(fullfile(params.paths.closure_struct,...
    ['closure_struct_' distribution_str '.mat']));
closure_struct = closure_struct.closure_struct;
closure_mtx = closure_struct.closure_mtx;
bin_size = closure_struct.bin_size;
vals = -1:bin_size:1;
r_idx = zeros(1, 2);
for kev = 1:2
    r_raw = EC(istim, kev);
    [~, idx] = min(abs(r_raw - vals));
    r_idx(kev) = idx;
    % May need to correct bin index:
    if idx > 1 && r_raw <= vals(idx)
        r_idx(kev) = idx - 1;
    end
end
p = closure_mtx(r_idx(1), r_idx(2));

end
