## Kinetic Model Description

This folder contains all scripts required for numerical simulations of the kinetic model.

1. **Reaction Definition**  
   Reaction equations are generated in **rxns.m** based on reaction rules defined in **funs12.m**  

2. **Flux Construction**  
   Flux equations are assembled in **nodModel_Flux.m** based on the generated reactions

3. **ODE System Assembly**  
   The system of ordinary differential equations (ODEs) is defined in **nodModel_ODE.m**, incorporating all flux terms

4. **Parameters**  
   Kinetic parameters are loaded from **kineticparams.m**  
   Representative Vmax values are recorded in **Vmaxvalues.m**

5. **Model Simulation**  
   The main script **nodModel_Run.m** calls the ODE solver (*ode15s*) to simulate system dynamics and compute nitrogen fixation rate and efficiency

---

