%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains supporting information of the following publication:
%
% Kinetic model of a determinate nodule reveals plant metabolic characteristics 
% for more efficient nitrogen fixation symbiosis
%
% by Rourou Ji, Joshua A.M. Kaste and Megan L. Matthews
%
% This MATLAB script allows to generate lhs matrix
%
% Author: rourouj2@illinois.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [X_scaled]=lhsdesign_loguniform(p_num, p_min,p_max)
%Inputs: 
%       p_num: number of radomly generated data points
%       min_ranges_p: [1xp] or [px1] vector that contains p values that correspond to the minimum value of each variable
%       max_ranges_p: [1xp] or [px1] vector that contains p values that correspond to the maximum value of each variable
%Outputs
%       X_scaled: [nxp] matrix of randomly generated variables within the
%       min/max range that the user specified
%       X_normalized: [nxp] matrix of randomly generated variables within the
%       0/1 range 

    dim=length(p_min);%num of variables

    %transpose if necessary
    [M1,N1]=size(p_min);
    if M1>N1
        p_min=p_min';
    end
        
    [M2,N2]=size(p_max);
    if M2>N2
        p_max=p_max';
    end

    lhs_raw = lhsdesign(p_num, dim);  
    log_min = log10(p_min);
    log_max = log10(p_max);
    
    log_samples = log_min + (log_max - log_min) .* lhs_raw;
    
    X_scaled = 10.^log_samples;  

end

