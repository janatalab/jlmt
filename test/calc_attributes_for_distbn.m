
integration_constants = [.1 2 4]; 

representational_space = {'PP' 'CV' 'TS'};
calculation_type = {'MC' 'MV'};
window_comparison = {'abs' 'rel'};
post_target_window = {[0 200] [201 600]};

calculate_closure_tf = 1;
closure_post_event_window = [100 300];

making sure they are in
% ascending order.
integration_constants = sort(integration_constants, 'ascend');
nconstant = size(integration_constants, 2);
nspace = size(representational_space, 2);
ntype = size(calculation_type, 2);
ncomparison = size(window_comparison, 2);
nwindow = size(post_target_window, 2);

%% Iterate over stimuli.
for istim = 1:nstim
    %% Iterate over pairs of integration constants.
    for iconstant = 1:nconstant
        for jconstant = iconstant+1:nconstant
            %% Iterate over representational space.
            for ispace = 1:nspace
                % Get variable representing this space.
                
                %% Iterate over calculation type.
                for itype = 1:ntype
                    %% Iterate over window comparison.
                    for icomparison = 1:ncomparison
                        %% Iterate over post-target window.
                        for iwindow = 1:nwindow
                            % Do stuff!
                            
                        end % for iwindow
                    end % for icomparison
                end % for itype
                % Calculate closure attributes if desired.
            end % for ispace
        end % for jconstant
    end % for iconstant
end % for istim










