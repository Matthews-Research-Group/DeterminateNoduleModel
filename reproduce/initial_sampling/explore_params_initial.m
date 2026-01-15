%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains supporting information of the following publication:
%
% Mechanistic modeling of the determinate nodule metabolism reveals enzymatic 
% influences on improving nitrogen fixation efficiency
%
% by Rourou Ji, Joshua A.M. Kaste and Megan L. Matthews
%
% This MATLAB script allows to explore the initial parameter space.
%
% Author: rourouj2@illinois.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function explore_params_initial(batch_index)
    % Validate input
    if nargin < 1
        error('Batch index must be provided as an input argument.');
    end

    main_dir = 'ODE_results';
    if ~exist(main_dir, 'dir')
        mkdir(main_dir); 
    end
    
    output_dir1 = fullfile(main_dir, 'LHS_Results_1e_6');
    if ~exist(output_dir1, 'dir')
        mkdir(output_dir1); 
    end
    
    output_dir2 = fullfile(main_dir, 'LHS_Results_4e_6');
    if ~exist(output_dir2, 'dir')
        mkdir(output_dir2);
    end

    output_dir5 = fullfile(main_dir, 'LHS_Results_2e_7');
    if ~exist(output_dir5, 'dir')
        mkdir(output_dir5); 
    end    

    % Directory to save log_paramvalues
    output_dir4 = 'log_param';
    if ~exist(output_dir4, 'dir')
        mkdir(output_dir4); 
    end

    % Add log file
    log_file = fullfile(output_dir4, sprintf('simulation_batch_%d_log.log', batch_index));
    fid = fopen(log_file, 'a'); 
    if fid == -1
        error('Cannot open log file for writing.');
    end

    fprintf(fid, 'Starting batch %d processing at %s\n', batch_index, datestr(now));

    %% Basic setting
    rng(42);

    addpath('../../kineticmodel/');
    [metabolite_names, reaction_names, x0, tspan, options] = Parameter_settings();
    Vmax=Vmaxvalues(); %get baseline values for Vmax
    param_names = fieldnames(Vmax);  % Get field names
    num_params = numel(param_names);  % Number of parameters in E


    %% Sample settings
    num_samples=6e5;%Sample size
    batch_size=6000; % Number of samples per batch
    num_batches = num_samples / batch_size; % Number of batches
    
    %Convert LHS matrix into an array of structs
    lhs_struct_array(num_samples) = Vmax; 

    %% LHS matrix
    % tic;
    [paramAll_min, paramAll_max] = param_ranges_initial(Vmax, param_names);
    [P_scaled]=lhsdesign_loguniform(num_samples, paramAll_min, paramAll_max);

    fprintf(fid, 'P_scaled is loaded at %s with size [%d %d]\n', datestr(now), size(P_scaled,1), size(P_scaled,2));

    save('./P_scaled_initial.mat','P_scaled');
     
    %% run LHS
    for i = 1:num_samples
        for j = 1:num_params
            % Assign scaled values to the corresponding fields in the struct
            lhs_struct_array(i).(param_names{j}) = P_scaled(i, j);
        end
    end

    %% parallel setting
    % Determine start and end index for this batch
    start_idx = (batch_index - 1) * batch_size + 1;
    end_idx = start_idx + batch_size - 1;
    fprintf(fid, 'Processing samples %d to %d\n', start_idx, end_idx);

    % Use the SLURM-provided TMPDIR
    tmp_dir = getenv('TMPDIR');
    cluster = parcluster('local');
    cluster.JobStorageLocation = tmp_dir;
    num_workers = str2double(getenv('SLURM_CPUS_PER_TASK'));
    if isempty(gcp('nocreate'))
        parpool(cluster, num_workers);
    end

    %% Run simulation for this batch
    parfor i = start_idx:end_idx
        
        nodLHSmatrix=lhs_struct_array(i);%Get the parameter set for this run
        
        % Solve ODE
        [ T, X, FLUXss ] = nodrunSimulation(nodLHSmatrix,tspan,x0,options);
        Xss = X(end,:);

        % Check if the last time point equals tmax
        if ~isempty(T) && T(end) == tspan(end)  
            slope=nodModel_ODE(T(end),X(end,:),nodLHSmatrix);
            
            if all(abs(slope) < 4e-6)
                % Choose output directory by slope threshold
                if all(abs(slope) < 2e-7)
                    save_dir = output_dir5;
                elseif all(abs(slope) < 1e-6)
                    save_dir = output_dir1;
                else
                    save_dir = output_dir2;
                end  

                % Generate output filename
                timestamp = datestr(now, 'yyyymmdd_HHMMSS');
                filename_base = fullfile(save_dir, sprintf('batch_%d_run_%d_%s_', batch_index, i, timestamp));
                   
                % Fill NaNs if needed
                num_rxn = length(reaction_names);
                Xss_full = [Xss, NaN(1, num_rxn - length(Xss))];
                slope_full = [slope', NaN(1, num_rxn - length(slope))];
                Vmax_LHS = cell2mat(struct2cell(nodLHSmatrix))';
    
                combined_table = table(Xss_full', slope_full', FLUXss', Vmax_LHS', ...
                    'VariableNames', {'Xss', 'Slope', 'FLUXss', 'Vmax_LHS'});
                writetable(combined_table, strcat(filename_base, 'combined_output.csv'));

            end
        end
    end

    %% all code finishes
    % Create a summary log file for the batch
    summary_log_file = fullfile(output_dir4, sprintf('batch_%d_summary.log', batch_index));
    fid_summary = fopen(summary_log_file, 'a');
    
    if fid_summary ~= -1
        fprintf(fid_summary, 'Batch %d completed at %s\n', batch_index, datestr(now));
        fclose(fid_summary);
    else
        warning('Could not write summary log for batch %d', batch_index);
    end


    rmpath('../../kineticmodel/');
end
    

