%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains supporting information of the following publication:
%
% Mechanistic modeling of the determinate nodule metabolism reveals enzymatic 
% influences on improving nitrogen fixation efficiency
%
% by Rourou Ji, Joshua A.M. Kaste and Megan L. Matthews
%
% This MATLAB script allows to determine parameter ranges for LHS 
%
% Author: rourouj2@illinois.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [paramAll_min, paramAll_max] = param_ranges(param_names)
%Vmax: a struct
%param_names: fieldnames of Vmax

    load('../../data/Vmax_range.mat');

    num_fields = length(param_names);
    mixVmaxE_num = num_fields-1;

    param_names_new = param_names;
    param_names_new(19,:)=[];

    % initialize min and max
    param_min = zeros(mixVmaxE_num, 1);
    param_max = zeros(mixVmaxE_num, 1);
    param_min(:,1)=1;
    param_max(:,1)=1000; 
   
    %make a cell
    range=cell(mixVmaxE_num,3);
    range(:,1)=param_names_new;
    range(:,2)=num2cell(param_min);
    range(:,3)=num2cell(param_max);

    idx_TKT1=find(strcmp(range(:,1), 'TKT1'));

    for i = 1:mixVmaxE_num
        fname = range{i,1};

        match_idx = find(strcmp(Vmax_range{:,1}, fname), 1);%if fname is in Vmax_range{:,1}

        % Special case for i=18
        if i == idx_TKT1
            range{i,2} = Vmax_range{match_idx, 7};   % TKT_E_min
            range{i,3} = Vmax_range{match_idx, 9};   % TKT_E_max
            continue;% skip
        end

        if ~isempty(match_idx)  
            range{i,2} = Vmax_range{match_idx, 8};  %min
            range{i,3} = Vmax_range{match_idx, 10};  %max
        end
    end
    
    %% output
    paramAll_min = cell2mat(range(:,2));
    paramAll_max = cell2mat(range(:,3));

end



