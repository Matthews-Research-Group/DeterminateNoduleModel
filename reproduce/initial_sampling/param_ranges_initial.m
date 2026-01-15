%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains supporting information of the following publication:
%
% Mechanistic modeling of the determinate nodule metabolism reveals enzymatic 
% influences on improving nitrogen fixation efficiency
%
% by Rourou Ji, Joshua A.M. Kaste and Megan L. Matthews
%
% This MATLAB script allows to determine parameter ranges for initial sampling 
%
% Author: rourouj2@illinois.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [paramAll_min, paramAll_max] = param_ranges_initial(Vmax, param_names)

    numAll = numel(fieldnames(Vmax));

    %% general Vmax range of all rxns
    paramAll_min=1*ones(numAll,1);%min Vmax mM/s
    paramAll_max=1000*ones(numAll,1);%max Vmax mM/s












 
