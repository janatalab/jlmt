% Analysis of example files from Collins et all 2012
% 
% 

audio_root = fullfile(fileparts(fileparts(which('jlmt_proc_series'))),...
    'data');

expmt_names = {'bigand_etal_jep_2003' 'bigand_etal_jep_2003' ...
    'marmel_and_tillmann_mp_2009' 'marmel_and_tillmann_mp_2009'};

audio_names = {'B_I_NTC' 'B_IV_NTC' 'Riii0' 'Rvii0'};

% Initialize structure that will be passed to jlmt_proc_series
in_data = ensemble_init_data_struct;
in_data.vars =...
    {'path', 'filename', 'path_no_ext', 'name_no_ext'};
ic = set_var_col_const(in_data.vars);
in_data.data{ic.path} = {};
in_data.data{ic.filename} = {};
in_data.data{ic.path_no_ext} = {};
in_data.data{ic.name_no_ext} = {};

% Loop over music stimuli.
nstim = size(audio_names, 2);
for istim = 1:nstim
    curr_path = fullfile(audio_root, expmt_names{istim},...
        audio_names{istim}, 'audio');
    file_list = listFilesOfType(curr_path, {'wav','mp3'});
    if length(file_list) > 1
        error('Found %d files in %s\n', length(file_list), curr_path);
    end
    in_data.data{ic.path}{istim} = curr_path;
    in_data.data{ic.filename}{istim} = file_list;
    [path_no_ext, ~, ~] = fileparts(curr_path);
    in_data.data{ic.path_no_ext}{istim} = path_no_ext;
    in_data.data{ic.name_no_ext}{istim} = audio_names{istim};
end

params = tonmodcomp_globals;
jlmt_out = jlmt_proc_series(in_data, params);



