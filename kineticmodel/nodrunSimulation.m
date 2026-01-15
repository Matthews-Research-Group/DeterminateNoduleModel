%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains supporting information of the following publication:
%
% Mechanistic modeling of the determinate nodule metabolism reveals enzymatic
% influences on improving nitrogen fixation efficiency
%
% by Rourou Ji, Joshua A.M. Kaste and Megan L. Matthews
%
% This MATLAB script includes the function to run simulation.
%
% Author: rourouj2@illinois.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ T, X, FLUXss ] = nodrunSimulation(E,tspan,x0,options)

    funode = @(t,x) nodModel_ODE(t,x,E);
    [T, X] = ode15s(funode,tspan,x0,options);
    
    FLUXss= nodModel_Flux(T(end), X(end,:), E);

end
