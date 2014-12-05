from __future__ import print_function, division
import sympy as s
from sympy.utilities import lambdify
import numpy as n
from numbers import Number
from tools import *
from IPython.display import display as disp
s.init_printing()


t = s.Symbol('t')


#The Body Stuff

class PointMass2D(object):
    
    def __init__(self, name, m, r):
        self.name = name
        self.m = m

        self.r = s.Matrix(r)
        self.v = s.Matrix([s.diff(component, t) for component in self.r])
        self.a = s.Matrix([s.diff(component, t) for component in self.v])

        self.forces = []
        self.places = None


    def Force(self, F=None, M=None):
        F = Force2D(F,M)
        self.forces.append(F)
        return F


    def Gravity(self, g):
        G = Gravity2D(g)
        self.forces.append(G)
        return G
    

    def Drag(self, Cd, power):
        D = Drag2D(Cd, power)
        self.forces.append(D)
        return D
    

    def place(self, initial_condition):
        self.places = initial_condition
                        


class RigidBody2D(object):

    def __init__(self, name, m, I, r, ang):
        self.name = name
        self.m = m
        self.I = I

        self.r = s.Matrix(r)
        self.v = s.Matrix([s.diff(component, t) for component in self.r])
        self.a = s.Matrix([s.diff(component, t) for component in self.v])

        self.ang = ang
        self.omega = s.diff(self.ang, t)
        self.alpha = s.diff(self.omega, t)

        self.forces = []
        self.places = None
    

    def Force(self, F=None, M=None):
        F = Force2D(F,M)
        self.forces.append(F)
        return F


    def Gravity(self, g):
        G = Gravity2D(g)
        self.forces.append(G)
        return G
    

    def Drag(self, Cd, power):
        D = Drag2D(Cd, power)
        self.forces.append(D)
        return D
    

    def place(self, initial_condition):
        self.places = initial_condition



# The force stuff

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
       


class Gravity2D(Force2D):
    
    def __init__(self, g):
        super(Gravity2D, self).__init__(lambda bod: bod.m * s.Matrix(g))



class Drag2D(Force2D):
    
    def __init__(self, Cd, power):
        
        super(Drag2D, self).__init__(lambda bod: -Cd * bod.v * 
                                        bod.v.norm()**(power - 1))



#The Sim Stuff


class Sim2D(object):
    
    def __init__(self):
        self.bods = []
        self.forces = []
        self.basis = []
        self.dbasis = []
        self.ddbasis = []
        self.places = {}
    

    def DOF(self, *args):
        if len(args) == 1:
            dof = s.Function(args[0])(t)
            self.basis.append(dof)
            self.dbasis.append(s.diff(dof,t))
            self.ddbasis.append(s.diff(dof, t, 2))
            return dof
        else:
            dofs = [s.Function(name)(t) for name in args]
            self.basis.extend(dofs)
            self.dbasis.extend([s.diff(dof, t) for dof in dofs])
            self.ddbasis.extend([s.diff(dof, t, 2) for dof in dofs])
            return tuple(dofs)



    def PointMass(self, name, m, r=None):
        PM = PointMass2D(name,m,r)
        self.bods.append(PM)
        return PM
    

    def RigidBody(self, name, m, I, r=None, ang=None):
        RB = RigidBody2D(name,m,I,r,ang)
        self.bods.append(RB)
        return RB


    def Force(self, F=None, M=None):
        F = Force2D(F,M)
        self.forces.append(F)
        return F
        
    def Gravity(self, g):
        G = Gravity2D(g)
        self.forces.append(G)
        return G
    
    def Drag(self, Cd, power):
        D = Drag2D(Cd, power)
        self.forces.append(D)
        return D
    
    def place(self, initial_condition):
        self.places = initial_condition

    def compile(self):
        
        self.r, self.v, self.a  = [], [], []
        self.M = []

        for bod in self.bods:
            self.r.extend(list(bod.r))
            self.v.extend(list(bod.v))
            self.a.extend(list(bod.a))
            self.M.extend([bod.m]*len(bod.r))
            try:
                self.r.append(bod.ang)
                self.v.append(bod.omega)
                self.a.append(bod.alpha)
                self.M.append(bod.I)
            except AttributeError:
                pass

        self.r = s.Matrix(self.r)
        self.v = s.Matrix(self.v)
        self.a = s.Matrix(self.a)
        M = s.Matrix([[0]*len(self.r)]*len(self.r))
        for i,m in enumerate(self.M):
            M[i,i] = m
        self.M = M

        
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
                           for bod in self.bods + [self]
                           if bod.places], []))
            
        init_pos = []
        init_vel = []
        for bvar, dbvar in zip(self.basis,self.dbasis):
            init_pos.append((str(type(bvar)), places[bvar][0]))
            init_vel.append(('d' + str(type(dbvar)), places[bvar][1]))
        inits = init_pos + init_vel
        
        print('Simulation Compiled, now beginning to run')
        
        return odesolve(self.odes, inits, time, events)
        

    def _gen_free_dirs(self):
        
        r = self.r
        basis = self.basis
        
        dirs = s.Matrix([[0]*len(r)]*len(basis))
        for i,var in enumerate(basis):
            direction = s.Matrix([s.diff(component, var) for component in r])
            dirs[i,:] = direction.transpose()
            
        return dirs
    
    
    def _gen_odes(self):
        
        input_vars = self.basis + self.dbasis
        divider = len(self.basis) #between basis and dbasis
        
        #the replacement is done to avoid problems with backslashes
        newvars = s.symbols(['x_'+str(i+1) for i in range(0,len(input_vars))])
        #the flipping of the replacement vectors is to avoid problems with
        #variables being replaced before their derivatives
        nRHMat = lambdify(newvars, self.RHMat.subs(zip(input_vars[::-1],newvars[::-1]))
                        ,modules='numpy')
        nLHS = lambdify(newvars, self.LHS.subs(zip(input_vars[::-1],newvars[::-1]))
                        ,modules='numpy')


        def odes(t,y):
            ddvs = n.linalg.inv(nRHMat(*y)) * nLHS(*y)
            return n.hstack((y[divider:],n.ravel(ddvs)))
        
        return odes




print('Pysics imported! You now can simulate rigid body dynamics.')


if __name__ == '__main__':

    #Single Pendulum
    sim = Sim2D()
    sim.Force(lambda bod: [0, -bod.m*9.81])
    th = s.Symbol('th', real=True)
    pend = sim.PointMass('pend', 2, 0.3*s.Matrix([s.sin(th),-s.cos(th)]))
    pend.place({th:(0.1,0)})
    y = sim.run([0,5,0.001])
    from matplotlib.pyplot import *
    plot(y['t'],y['th'])
    show()
    
