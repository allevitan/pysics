from __future__ import print_function, division
import sympy as s
from sympy.utilities import lambdify
import numpy as n
from numbers import Number
from tools import *
from matplotlib.pyplot import *


class PointMass2D(object):
    
    def __init__(self, name, m, r):
        self.name = name
        self.m = m
        self.inits = None
        self.r = s.Matrix(r)
        self.basis = list(extract_variables(self.r))
        self.places = None


    def place(self, initial_condition):
        self.places = initial_condition
                        


class RigidBody2D(object):

    def __init__(self, name, m, I, r, ang):
        self.name = name
        self.m = m
        self.I = I
        self.r = s.Matrix(r)
        self.ang = ang
        self.places = None
    
    
    def place(self, initial_condition):
        self.places = initial_condition



class Force2D(object):
    
    def __init__(self, F=None, M=None):
        
        self.F = F
        self.M = M
    

    def __call__(self, bod):
        if self.F is not None:
            F = self.F(bod)
        else:
            F = [0,0]

        try:
            I = bod.I #Just to test if it has a defined angle
            if self.M is not None:
                M = self.M(bod)
            else:
                M = 0
            return s.Matrix(list(F) + [M])

        except AttributeError:
            pass
        
        return s.Matrix(F)
       
 

class Sim2D(object):
    
    def __init__(self):
        self.bods = []
        self.forces = []
    

    def PointMass(self, name, m, r=None):
        PM = PointMass2D(name,m,r)
        self.bods.append(PM)
        return PM
    

    def RigidBody(self, name, m, I, r=None, ang=None):
        RB = RigidBody2D(name,m,I,r)
        self.bods.append(RB)
        return RB


    def Force(self, F=None, M=None):
        F = Force2D(F,M)
        self.forces.append(F)
        return F


    def compile(self):
        
        self.r = []
        self.M = []

        for bod in self.bods:
            self.r.extend(list(bod.r))
            self.M.extend([bod.m]*len(bod.r))
            try:
                self.r.append(bod.ang)
                self.M.append(bod.I)
            except AttributeError:
                pass

        self.r = s.Matrix(self.r)
        M = s.Matrix([[0]*len(self.r)]*len(self.r))
        for i,m in enumerate(self.M):
            M[i,i] = m
        self.M = M

        
        self.basis = list(extract_variables(self.r))
        self.dbasis = [s.Symbol('d' + var.name) for var in self.basis]
        self.ddbasis = [s.Symbol('dd' + var.name) for var in self.basis]
        
        self.v, self.a = self._gen_derivs()
        self.dirs = self._gen_free_dirs()
        
        self.F = s.Matrix([0]*len(self.r))
        for force in self.forces:
            F = []
            for bod in self.bods:
                F.extend(force(bod))
            self.F = self.F + s.Matrix(F)
        
        #Now we write F=Ma in our new coordinate system
        LHS = s.simplify(self.dirs * self.F)
        RHS = s.simplify(self.dirs * self.M * self.a)
        
        #And put it into Vec = Mat*Diffs + Vec form
        self.dds = s.Matrix(self.ddbasis)
        self.RHMat = s.Matrix([[0]*len(self.basis)]*len(self.basis))
        for i, ddvar in enumerate(self.ddbasis):
            coeffs = [s.Poly(component,ddvar).coeff_monomial(ddvar)
                      for component in RHS]
            self.RHMat[i,:] = s.Matrix([s.simplify(coeff)
                                   for coeff in coeffs]).transpose()

        RHVec = s.simplify(RHS - self.RHMat * self.dds)
        
        self.LHS = LHS - RHVec
        
        self.odes = self._gen_odes()


    def run(self, time, events=None):
        self.compile()

        places = dict(sum([bod.places.items()
                           for bod in self.bods if bod.places], []))
        init_pos = []
        init_vel = []
        for bvar, dbvar in zip(self.basis,self.dbasis):
            init_pos.append((bvar.name, places[bvar][0]))
            init_vel.append((dbvar.name, places[bvar][1]))
        inits = init_pos + init_vel
        
        return odesolve(self.odes, inits, time, events)
        

    def _gen_derivs(self):
    
        r = self.r
        locs = self.basis
        vels = self.dbasis
        accs = self.ddbasis
        
        t = s.Symbol('t')
        
        v = s.Matrix([0]*len(r))
        for i in range(0,len(locs)):
            subbed = r.subs(locs[i],vels[i]*t)
            diffed = s.Matrix([s.diff(component, t) for component in subbed])
            v += diffed.subs(vels[i]*t,locs[i])

        a = s.Matrix([0]*len(r))
        for i in range(0,len(locs)):
            subbed = v.subs(locs[i],vels[i]*t)
            diffed = s.Matrix([s.diff(component, t) for component in subbed])
            a += diffed.subs(vels[i]*t,locs[i])
            subbed = v.subs(vels[i],accs[i]*t)
            diffed = s.Matrix([s.diff(component, t) for component in subbed])
            a += diffed.subs(accs[i]*t,vels[i])
        
        return (v,a)
    
    
    def _gen_free_dirs(self):
        
        r = self.r
        basis = self.basis
        
        dirs = s.Matrix([[0]*len(r)]*len(basis))
        for i,var in enumerate(basis):
            direction = s.Matrix([s.diff(component, var) for component in r])
            dirs[i,:] = direction.transpose()
            
            #Use if you want to have the directions normalized
            # mag = s.sqrt(simplify(direction.dot(direction)))
            # unit = s.simplify(direction/mag)
            # dirs[i,:] = unit.transpose()
        
        return dirs
    
    
    def _gen_odes(self):
        
        input_vars = self.basis + self.dbasis
        divider = len(self.basis) #between basis and dbasis
        
        #this is done to avoid problems with backslashes
        newvars = s.symbols(['x_'+str(i+1) for i in range(0,len(input_vars))])
        nRHMat = lambdify(newvars, self.RHMat.subs(zip(input_vars,newvars))
                        ,modules='numpy')
        nLHS = lambdify(newvars, self.LHS.subs(zip(input_vars,newvars))
                        ,modules='numpy')

        def odes(t,y):
            ddvs = n.linalg.inv(nRHMat(*y)) * nLHS(*y)
            return n.hstack((y[divider:],n.ravel(ddvs)))
        
        return odes

    

if __name__ == '__main__':

    sim = Sim2D()
    sim.Force(lambda bod: [0, -bod.m*9.81])
    
    th = s.Symbol('th', real=True)
    pend = sim.PointMass('pend', 2, 0.3*s.Matrix([s.sin(th),-s.cos(th)]))
    pend.place({th:(0.1,0)})
    y = sim.run([0,5,0.001])
    plot(y['t'],y['th'])
    show()
    
    #rpend = sim.RigidBody('rpend', 2, 3,
    #                      0.3*s.Matrix([s.sin(th),-s.cos(th)]), ang=th)
