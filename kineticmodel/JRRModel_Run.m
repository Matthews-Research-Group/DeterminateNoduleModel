%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains supporting information of the following publication:
%
% Mechanistic modeling of the determinate nodule metabolism reveals enzymatic 
% influences on improving nitrogen fixation efficiency
%
% by Rourou Ji, Joshua A.M. Kaste and Megan L. Matthews
%
% This MATLAB script allows to run ODE and get metabolite concentration and flux values
%
% Author: rourouj2@illinois.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;
rng(123);

% initial values and load parameter settings
[metabolite_names, reaction_names, x0, tspan, options] = Parameter_settings();
Vmax = Vmaxvalues();

Vmaxnew=Vmax;

% run ODE
[ T, X, FLUXss ] = JRRrunSimulation(Vmaxnew,tspan,x0,options);

%compute slope
dt=100;
Tindex_total=length(T);
tpoints_end=Tindex_total;
tpoints_last100=Tindex_total-dt;

slope_last100=JRRModel_ODE(T(tpoints_last100),X(tpoints_last100,:),Vmaxnew);
slopeend=JRRModel_ODE(T(end),X(end,:),Vmaxnew);
slopechange_end=(slopeend-slope_last100)/100;

slopedata=cell(size(X,2),2);
slopedata(:,1)=metabolite_names;
slopedata(:,2)=num2cell(slopeend);

% steady state concentration and flux
Xss = X(end,:)'; 
FLUXss = FLUXss';
concdata = cell(length(Xss),2);
concdata(:,1) = metabolite_names;
concdata(:,2) = num2cell(Xss);
fluxdata = cell(length(FLUXss),2);
fluxdata(:,1) = reaction_names;
fluxdata(:,2) = num2cell(FLUXss);

%% Plots
figure();
plot(T,X);
xlabel('Time');
ylabel('Metabolite Concentration(mM)');
title('ODE Simulation Results');

%%
gCgN = gcgn(FLUXss);
params = kineticparams(Vmaxnew);


%% Functions
function z=gcgn(v)
    %gCgN
    CO2fluxout = v(15,:) + v(12,:) + v(30,:) + v(31,:) + ...
                 v(24,:).*4 + v(33,:).*4 + v(35,:).*4 + v(61,:).*2 + v(67,:);
    CO2fluxin = v(54,:) + v(22,:);
    CO2efflux = CO2fluxout - CO2fluxin;
    NH4fluxout = v(68,:);
    z = 12 * CO2efflux ./ (14 * 4 * NH4fluxout);
end

