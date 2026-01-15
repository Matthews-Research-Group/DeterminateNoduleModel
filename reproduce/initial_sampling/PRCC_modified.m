%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains supporting information of the following publication:
%
% Mechanistic modeling of the determinate nodule metabolism reveals enzymatic
% influences on improving nitrogen fixation efficiency
%
% by Rourou Ji, Joshua A.M. Kaste and Megan L. Matthews
%
% This MATLAB script allows to calculate prcc values.
%
% Author: rourouj2@illinois.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [prcc, prcc_sign,sign_indices,sign_param,sign_prccvalues] = PRCC_modified(LHSmatrix, Y, PRCC_var, alpha)
    % Function to calculate Partial Rank Correlation Coefficients (PRCC) and significance values
    % Inputs:
    %   N: sampling times
    %   k: param numbers
    %   LHSmatrix: N x k matrix from Latin Hypercube Sampling
    %   Y: Output matrix (N*num_outputs)
    %   s: Row vector of time points to test
    %   PRCC_var: Cell array of parameter names
    %   alpha: Significance threshold
    % Outputs:
    %   prcc: PRCC values 
    %   sign: p-values (uncorrected)
    %   sign_label: Structure containing significant PRCC indices, labels, and values

    % Extract the subset of Y corresponding to the given time points
    % newY = Y';  % newY is N*1

    % Get matrix dimensions
    [N, k] = size(LHSmatrix);%sampling_times*param_nums
    num_outputs = size(Y, 2);  % Number of outputs

    % Initialize output matrices
    prcc = zeros(num_outputs, k);
    prcc_sign = zeros(num_outputs, k);

    % Loop through each parameter to compute PRCC
    for i = 1:k
        % Remove the current parameter from LHSmatrix to form the control set Z
        Z = LHSmatrix;
        Z(:, i) = []; % Remove the i-th parameter

        % Compute PRCC for each output (each time point)
        for j = 1:num_outputs
            A1=LHSmatrix(:, i);
            A2=Y(:, j);
            A=[A1,A2];
            [rho, p] = partialcorr(A, Z, 'type', 'Spearman');
            prcc(j, i) = rho(1, 2);  % Store PRCC value
            prcc_sign(j, i) = p(1, 2); % Store p-value

        end
    end

    % Identify significant PRCCs based on alpha threshold
    for r = 1:num_outputs
        sign_indices = find(prcc_sign(r, :) < alpha);
        sign_param=PRCC_var(sign_indices);%significant param names
        sign_prccvalues=prcc(r, sign_indices);
    end

end
