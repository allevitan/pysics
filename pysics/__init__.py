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


    def Force(self, F=[0,0], M=0):
        F = Force2D(F,M)
        self.forces.append(F)


    def Gravity(self, g):
        G = Gravity2D(self, g)
        self.forces.append(G)
    

    def Drag(self, LCd,  power):
        D = Drag2D(self, LCd=LCd, power=power)
        self.forces.append(D)
    

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
    

    def Force(self, F=[0,0], M=0):
        F = Force2D(F,M)
        self.forces.append(F)


    def Gravity(self, g):
        G = Gravity2D(self, g)
        self.forces.append(G)
    

    def Drag(self, LCd=None, RCd=None, power=1):
        D = Drag2D(self, LCd, RCd, power)
        self.forces.append(D)
    

    def place(self, initial_condition):
        self.places = initial_condition



# The force stuff

def Force2D(F=[0,0], M=0):
    
    return s.Matrix(list(F) + [M])


def Gravity2D(bod, g):
    return Force2D(F=bod.m * s.Matrix(g))



def Drag2D(bod, LCd=None, RCd=None, power=1):
        
        #Don't let anything dependent on omega be created if bod
        #isn't a rigid body.
        try:
            I = bod.I
        except AttributeError:
            return Force2D(F= -LCd * bod.v * bod.v.norm()**(power - 1))

        #Now we know bod is a rigid body
        if RCd and LCd:
            return Force2D(F= - LCd * bod.v * bod.v.norm()**(power - 1),
                           M= -RCd * bod.omega**power)

        elif LCd:
            return Force2D(F= -LCd * bod.v * bod.v.norm()**(power - 1))

        elif RCd:
            return Force2D(M= -RCd * bod.omega**power)
        
        else:
            return Force2D()
            



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


    def Force(self, F=lambda bod: [0,0], M=lambda bod: 0):
        #IMPORTANT NOTE: sim.Force takes a lambda function,
        # Everything else just takes an expression
        self.forces.append((Force2D,F,M))

        
    def Gravity(self, g):
        self.forces.append((Gravity2D,g))

    
    def Drag(self, LCd=None, RCd=None, power=1):
        self.forces.append((Drag2D,LCd,RCd,power))

    
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
        
        
        for force in self.forces:
            #If it's from force2D, which for generality asks for a function
            # that generates an expression
            if force[0] == Force2D:
                for bod in self.bods:
                    bod.forces.append(force[0](*[arg(bod) for arg in force[1:]]))
            else:
                for bod in self.bods:
                    bod.forces.append(force[0](bod,*force[1:]))

        self.F = []
        for bod in self.bods:
            F = s.Matrix([0]*len(bod.r) + [0])
            for force in bod.forces:
                F += force
            try:
                I = bod.I
            except AttributeError:
                F = s.Matrix(F[0:2])
            self.F.extend(list(F))
        self.F = s.Matrix(self.F)
                
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
        
        input_vars = [t] + self.basis + self.dbasis
        divider = len(self.basis) #between basis and dbasis
        
        #the replacement is done to avoid problems with backslashes
        newvars = s.symbols(['x_'+str(i) for i in range(0,len(input_vars))])
        #the flipping of the replacement vectors is to avoid problems with
        #variables being replaced before their derivatives
        nRHMat = lambdify(newvars, self.RHMat.subs(zip(input_vars[::-1],newvars[::-1]))
                        ,modules='numpy')
        nLHS = lambdify(newvars, self.LHS.subs(zip(input_vars[::-1],newvars[::-1]))
                        ,modules='numpy')

        def odes(t,y):
            ddvs = n.linalg.inv(nRHMat(t,*y)) * nLHS(t,*y)
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
    
