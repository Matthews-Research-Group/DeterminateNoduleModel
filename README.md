# Determinate nodule kinetic model

This GitHub contains supporting information for the following publication:

Kinetic modeling of the determinate nodule metabolism reveals enzymatic influences on improving nitrogen fixation efficiency

by Rourou Ji, Joshua A.M. Kaste, and Megan L. Matthews
   
Author: rourouj2@illinois.edu

# Description  #
This archive contains all the scripts required to construct the model and perform the simulations described in the publication, and reproduce figures in the main context and supplemental information. 

- **💻 kinetic model**: Source code and parameters for the original model.
- **💾 data**: Contains all datasets required to replicate the results presented in the paper.
- **📊 reproduce**: This directory contains scripts for generating figures and performing statistical analyses.
  Parameter Space Exploration: To explore the parameter space of the model featured in Figure 1, run the command "sbatch run_explore_params.sh". The script will automatically create a directory named "ODE_results" to store all sampling and simulation data.


# Requirements #
MATLAB_R2023a, MATLAB toolbox, and Biocluster high-performance computing system (v3, Champaign, IL, USA) are required. If you are using a different version or operating system, some functions or toolbox features might be incompatible. Please adjust the source code accordingly to match your specific environment.


