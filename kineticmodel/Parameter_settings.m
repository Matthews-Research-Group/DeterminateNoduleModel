%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains supporting information of the following publication:
%
% Mechanistic modeling of the determinate nodule metabolism reveals enzymatic 
% influences on improving nitrogen fixation efficiency
%
% by Rourou Ji, Joshua A.M. Kaste and Megan L. Matthews
%
% This MATLAB script includes the variable names.
%
% Author: rourouj2@illinois.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% This script contains all parameters to run ode
function [metabolite_names, reaction_names, x0, tspan, options] = Parameter_settings()
    % Variables Labels
    metabolite_names = {
        'GLC', 'G6P', 'F6P', 'FDP', 'GAP',...%5
        'DAP', 'BPG', 'PGA3','PGA2', 'PEP',...%10
        'PYR', 'GL6P','PGN','RB5P','R5P',... %15
        'X5P','E4P','S7P','OAA','MAL',... %20
        'ACCOA','CIT','ACO','ICIT','AKG',... %25
        'GLX','SUCCOA','SUC','FUM','B', ... %30
        'NH4','GLU','GLN','PRPP','PHP3', ...%35
        'SRN','GLY','PRA','GAR','FGAR' ...%40
        'ALN','THF','CH2THF','CHTHF','CHOTHF',...%45
        'FGAM','AIR','CAIR','ASP','SAICAR',...%50
        'AICAR','FAICAR','IMP','XMP','XAO',... %55
        'XAN','URATE','HIUH','OHCU','ALTN',%60
        };
        
    
    reaction_names = {
        'HKI','PGI','PFK', 'FBP', 'FBA',...%5
        'TPI','GDH','PGK','GPM','ENO',...%10
        'PYK','PDH','ZWF','PGL','GND',...%15
        'RPE','RPI','TKT1','TKT2','TAL',... %20
        'PRS','PEPC','MDH','MALOut','GLT',...%25
        'ACN1','ACN2','ICL','MALS','ICDH',...%30
        'LPD','SK','SUCOut','SDH','FUMOut',... %35
        'FUMA','ToNH4', 'GS','GOGAT','GD',...%40
        'PRAT','PGDH','PSTS','SHMT','GARS',...%45
        'GARTF','FGAMS','ALTSS','AGT','MTHFR',...%50
        'MTHFD','GLYDC','AIRS','CAIRS','SS',...%55
        'ATS','ADSL','AICARTF','IMPCH','IMPDH',...%60
        'GLYsink','NT','PNP','XOR','UOD',...%65
        'HIUHS','OHCUD','ALTNOut'%68
        };

    %% Initial conditions for the ode model
    num_met = length(metabolite_names);
    x0=28.89*ones(num_met,1);

    %% time span of the simulation
    tmax=20000; 
    t_step=0.01;
    tspan=[0:t_step:tmax]; 
    options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'Stats', 'on');

end
