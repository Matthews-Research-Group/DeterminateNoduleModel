function [ T, X, FLUXss ] = nodrunSimulation(E,tspan,x0,options)

    funode = @(t,x) nodModel_ODE(t,x,E);
    [T, X] = ode15s(funode,tspan,x0,options);
    
    FLUXss= nodModel_Flux(T(end), X(end,:), E);

end
