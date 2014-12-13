from __future__ import print_function, division
import numpy as n
from scipy.integrate import ode as scipy_ode


def odesolve(odes, boundaries, time, events=None):
    """A simple ode solver that uses scipy's 4th/5th order
    Runge-Kutta integrator as it's workhorse.
    """
    boundary_vars, boundary_vals = zip(*boundaries)
    
    #Set up the simulation object
    sim = scipy_ode(odes).set_integrator('dopri5')
    sim.set_initial_value(boundary_vals,time[0])

    ts = [sim.t]
    ys = [sim.y]

    if callable(events):
        while sim.t < time[1] and sim.successful() and not events(sim.y):
            sim.integrate(sim.t + time[2])
            ys.append(sim.y)
            ts.append(sim.t)
    
    else:
        while sim.t < time[1] and sim.successful():
            sim.integrate(sim.t + time[2])
            ys.append(sim.y)
            ts.append(sim.t)
    
    ys = n.array(ys).transpose()
    out = {'t':ts}
    for variable, values in zip(boundary_vars, ys):
        out[str(variable)] = values

    return out



