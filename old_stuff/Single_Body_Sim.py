from __future__ import print_function, division
from sympy import *
import numpy as n
from sympy.utilities import lambdify
from scipy.integrate import ode as scipy_ode
from matplotlib import pyplot as p
from IPython.display import display as disp
init_printing()


from tools import extract_variables

def gen_eom(r,m,f,vs):
    
    
    #Generate a list of  1st and 2nd derivatives
    dvs = [Symbol(r'\dot{' + var.name + r'}') for var in vs]
    ddvs = [Symbol(r'\ddot{' + var.name + r'}') for var in vs]
    

    #Will be needed for some manipulations
    t = Symbol('t')
    
    #Now find the velocity vector in terms of those variables
    v = Matrix([0]*len(r))
    for i in range(0,len(vs)):
        subbed = r.subs(vs[i],dvs[i]*t)
        diffed = Matrix([diff(component, t) for component in subbed])
        v += diffed.subs(dvs[i]*t,vs[i])
    
    #Now do the same for acceleration
    a = Matrix([0]*len(r))
    for i in range(0,len(vs)):
        subbed = v.subs(vs[i],dvs[i]*t)
        diffed = Matrix([diff(component, t) for component in subbed])
        a += diffed.subs(dvs[i]*t,vs[i])
        subbed = v.subs(dvs[i],ddvs[i]*t)
        diffed = Matrix([diff(component, t) for component in subbed])
        a += diffed.subs(ddvs[i]*t,dvs[i])
    
    #Now find the allowed directions of motion
    dirs = Matrix([[0]*len(r)]*len(vs))
    for i,var in enumerate(vs):
        direction = Matrix([diff(component, var) for component in r])
        mag = sqrt(simplify(direction.dot(direction)))
        unit = simplify(direction/mag)
        dirs[i,:] = unit.transpose()

    #Now we construct the LHS and RHS of the equation
    LHS = simplify(dirs*f)
    RHS = simplify(dirs*m*a)
    
    #Next we try to write it in the form LHS = mat*difs + vec
    #We extract the matrix (mat) multiplying the second derivatives (diff)
    diffs = Matrix(ddvs)
    mat = Matrix([[0]*len(vs)]*len(vs))
    for i, ddvar in enumerate(ddvs):
        mat[i,:] = Matrix([simplify(Poly(component,ddvar).coeff_monomial(ddvar))
                           for component in RHS]).transpose()

    #And the vector offset
    vec = simplify(RHS - mat * diffs)
    
    #Now we subtract vec from the LHS to be left with LHS = mat * diffs
    LHS = LHS - vec
    
    return mat, LHS, dvs


def make_odes(mat,LHS,vs,dvs):
    allvars = vs+dvs
    div = len(vs)
    #this is done to avoid problems with backslashes
    newvars = symbols(['x_'+str(i+1) for i in range(0,len(allvars))])
    nmat = lambdify(newvars, mat.subs(zip(allvars,newvars)),modules='numpy')
    nLHS = lambdify(newvars, LHS.subs(zip(allvars,newvars)),modules='numpy')

    def f(t,y):
        ddvs = n.linalg.inv(nmat(*y)) * nLHS(*y)
        return n.hstack((y[div:],n.ravel(ddvs)))
    
    return f

def odesolve(odes, boundaries, time, events=None):
    
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


class SingleBodySim(object):
    
    def __init__(self, m,r,f):
        self.m = m
        self.r = r
        self.f = f
        self.variables = list(extract_variables(r))
        self.compile_sim()


    def compile_sim(self):
        self.mat, self.LHS, self.dvariables \
            = gen_eom(self.r,self.m,self.f,self.variables)
    
    def sim(self, boundaries, time, events=None):
        odes = make_odes(self.mat, self.LHS, self.variables, self.dvariables)
        return odesolve(odes, boundaries, time, events=events)
        

#Define Parameters
# l = 0.5; m = 1; g = 9.81

# theta = Symbol(r'\theta',real=True)

# r = l*Matrix([sin(theta),-cos(theta)])
# f =  Matrix([0,-m*g])

# sim = SingleBodySim(m,r,f)
# y = sim.sim((('theta',0.8),('dtheta',0)), [0,10,0.01])

# p.plot(y['t'],y['theta'])
# p.show()


#Roller Coaster Problem

B=0; c=0;
h=1; L=1; m=1;
x = Symbol('x')

r = Matrix([x, h * (1 - B* x / L) * cos(3*n.pi*x / (2*L))**2])
f = Matrix([0, -m*9.81])
sim = SingleBodySim(m,r,f)
y = sim.sim((('x',0),('dx',0.01)), [0,10,0.01], lambda y: y[0] >= L)

p.plot(y['t'],y['x'])
p.show()
