# Determinate nodule kinetic model

Welcome! 👋  

This GitHub contains supporting information for the following publication:

### Kinetic model of a determinate nodule reveals plant metabolic characteristics for more efficient nitrogen fixation symbiosis

by **Rourou Ji, Joshua A.M. Kaste, and Megan L. Matthews**
   
Author: rourouj2@illinois.edu

---

## 📂 Description 
This archive contains all the scripts required to construct the model and perform the simulations described in the publication, and reproduce figures in the main context and supplemental information. 

**💻 kinetic model**: Source code and parameters for the original model.
  
**💾 data**: Contains all datasets required to replicate the results presented in the paper.
  
**📊 reproduce**: This directory contains scripts for generating figures and performing statistical analyses.
- **Parameter space exploration**: To explore the parameter space of the model featured in **Figure S1**, run the command `sbatch run_explore_params.sh` in each sampling repository. The script will automatically create a directory named **ODE_results** to store all sampling and simulation data.
- **Main_figures.m**: Generates all figures presented in the main text of the paper.
- **SI_figures.m**: Generates all figures included in the Supplementary Information.

---

## 🛠️ System requirements
MATLAB_R2023a, MATLAB toolbox, and HPC (Biocluster is used here) are required. If you are using a different version or operating system, some functions or toolbox features might be incompatible. Please adjust the source code accordingly to match your specific environment.


