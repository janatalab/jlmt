function p = CollinsEtAl_closure_probability(ecm, closure_struct)

% Copyright (c) 2012 The Regents of the University of California
% All Rights Reserved.

% This function takes two values (correlations, stored in an event
% correlation matrix) and a two-dimensional probability density function
% (pdf, stored in a structure) as input. It locates the coordinate of the
% values in the pdf and returns the associated probability.

% INPUT
%  ecm is a 2x2x2 array, assumed to contain two correlation values.
%  closure_struct is a struct containing a probability distribution, as
%   well as parameters that were used to calculate this distribution.

% Tom Collins, 2012.06.26.

% Correlation values.
ec = set_var_col_const(ecm.vars);
event_corr_mtx = ecm.data{ec.event_corr_mtx};
EC = zeros(1, 2);
for kev = 1:2
    EC(kev) = event_corr_mtx(1, 2, kev);
end

% Locate and calculate probability.
closure_mtx = closure_struct.closure_mtx;
bin_size = closure_struct.bin_size;
vals = -1:bin_size:1;
r_idx = zeros(1, 2);
for kev = 1:2
    r_raw = EC(kev);
    [~, idx] = min(abs(r_raw - vals));
    r_idx(kev) = idx;
    % May need to correct bin index:
    if idx > 1 && r_raw <= vals(idx)
        r_idx(kev) = idx - 1;
    end
end
p = closure_mtx(r_idx(1), r_idx(2));

end
