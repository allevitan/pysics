from __future__ import print_function, division
import sympy as s
from sympy.utilities import lambdify
import numpy as n
from numbers import Number
from tools import *
from IPython.display import display as disp
s.init_printing()


# We definine a universal time symbol so the end user can
# define time-dependent constraints and generate time
# derivatives of the degrees of freedom.
t = s.Symbol('t', real=True)

#
# The following section defines the objects that represent bodies
# in the simulation. There are two types of bodies:
# 1. PointMass2D - a point mass, which has a mass and location
# 2. RigidBody2D - a rigid body, which has a mass, mass moment of
#         inertia, center of mass location, and angular displacement
#


class PointMass2D(object):
    """Represents a 2 dimensional point mass. Takes:
    name: a string to identify the mass
    m: the mass of the object
    r: the location of the object as a function of the 
    simulation's degrees of freedom.
    """

    
    def __init__(self, name, m, r):
        """Sets up the point mass object, and generates the velocity and
        acceleration vectors symbolically.
        """
        self.name = name
        self.m = m

        self.r = s.Matrix(r)
        self.v = s.Matrix([s.diff(component, t) for component in self.r])
        self.a = s.Matrix([s.diff(component, t) for component in self.v])

        self.forces = []
        self.places = None
    
    
    def Force(self, F=[0,0], M=0):
        """Adds an arbitrary force to the point mass. Syntax
        is the same as the syntax for creation of a Force2D.
        """
        F = Force2D(F,M)
        self.forces.append(F)


    def Gravity(self, g):
        """A convenience method to add gravity to the point mass.
        """
        G = Gravity2D(self, g)
        self.forces.append(G)
    

    def Drag(self, TCd,  power):
        """A convenience method to add drag to the point mass.
        """
        D = Drag2D(self, TCd=TCd, power=power)
        self.forces.append(D)
    



class RigidBody2D(object):
    """Represents a 2 dimensional rigid body. Takes:
    name: a string to identify the mass
    m: the mass of the object
    I: the mass moment of inertia of the object
    r: the location of the object as a function of the 
    simulation's degrees of freedom.
    ang: the angular position of the object as a function
    of the simulation's degrees of freedom
    """

    def __init__(self, name, m, I, r, ang):
        """Sets up the rigid object, and generates the translational
        and rotational velocity and acceleration symbolically.
        """
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
        """Adds an arbitrary force to the rigid body. Syntax
        is the same as the syntax for creation of a Force2D.
        """        
        F = Force2D(F,M)
        self.forces.append(F)


    def Gravity(self, g):
        """A convenience method to add gravity to the rigid body.
         """
        G = Gravity2D(self, g)
        self.forces.append(G)
    

    def Drag(self, TCd=None, RCd=None, power=1):
        """A convenience method to add drag to the point mass.
        """
        D = Drag2D(self, TCd, RCd, power)
        self.forces.append(D)
    



#
# The following section defines the constructors for forces
# in the simulation. These are mostly conveniences as forces
# are stored internally as vectors of length 3 - the first two
# spots being the force, and the third being the moment.
#


def Force2D(F=[0,0], M=0):
    """Returns a general force and/or moment. Takes:
    F: the force as a function of the simulation's degrees
    of freedom
    M: the moment as a function of the simulation's degrees
    of freedom"""
    return s.Matrix(list(F) + [M])



def Gravity2D(bod, g):
    """A convenience to generate gravitational forces. Takes:
    g: acceleration due to gravity"""
    return Force2D(F=bod.m * s.Matrix(g))



def Drag2D(bod, TCd=[0,0], RCd=0, power=1):
    """A convenience to generate drag forces. Takes:
    bod: The body experiencing the drag
    TCd: The translational drag coefficient
    RCd: The rotational drag coefficient
    power: The power dependence on velocity (and omega)
    and generates a force of the form
    F = TCd * (bod.v)**power, M = RCd * (bod.omega)**power.
    """    
    try:
        I = bod.I
    except AttributeError:
        # If the body is a point mass, don't make anything
        # dependent on bod.omega.
        return Force2D(F= -TCd * bod.v * bod.v.norm()**(power - 1))
        
    # If it's not omega, do what you want
    return Force2D(F= - TCd * bod.v * bod.v.norm()**(power - 1),
                   M= -RCd * bod.omega**power)




#
# The following section defines the object that represents an
# entire simulation, including the tools it uses to generate
# equations of motion.
#

class Sim2D(object):
    """The Sim2D object doesn't take anything in it's constructor,
    but for ease of use all simulation objects should be defined using
    its methods, including degrees of freedom and bodies.
    """
    
    def __init__(self):
        self.bods = []
        self.forces = []
        self.basis = []
        self.dbasis = []
        self.ddbasis = []
        self.places = {}
    

    def DOF(self, *args):
        """Generates a degree of freedom for the simulation. All
        position vectors and forces should be defined using degrees
        of freedom generated here.
        Takes an arbitrary number of names for degrees of freedom"""
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
        """Generates and returns a point mass in the simulation.
        Works just like the constructor for PointMass2d"""
        PM = PointMass2D(name,m,r)
        self.bods.append(PM)
        return PM
    

    def RigidBody(self, name, m, I, r=None, ang=None):
        """Generates and returns a point mass in the simulation.
        Works just like the constructor for PointMass2d"""
        RB = RigidBody2D(name,m,I,r,ang)
        self.bods.append(RB)
        return RB


    def Force(self, F=lambda bod: [0,0], M=lambda bod: 0):
        """Creates a force that acts on all bodies in the simulation.
        IMPORTANT NOTE: sim.Force lambda functions that take in the
        body to be acted on and return the appropriate force.
        For example:
        sim.Force(F=lambda bod: [0, -bod.m * 9,81])
        would add gravity to the simulation"""
        self.forces.append((Force2D,F,M))

        
    def Gravity(self, g):
        """Convenience function that adds gravity to the simulation.
        Takes:
        g: the acceleration due to gravity (i.e. [0, -9.81])"""
        self.forces.append((Gravity2D,g))

    
    def Drag(self, TCd=None, RCd=None, power=1):
        """Convenience function that adds drag to every body in the
        simulation."""
        
        self.forces.append((Drag2D,TCd,RCd,power))

    
    def place(self, initial_condition):
        """Sets the initial condition of the simulation. Takes:
        initial_condition: a dictionary whose keys are the system's
        degrees of freedom, and whose values are tuples corresponding
        to the initial value and time derivative of that degree of
        freedom.
        i.e., sim.place({s: (0.3,0)}) would place the degree of freedom
        s at 0.3, and set it's initial time derivative to 0."""
        self.places = initial_condition

    def compile(self):
        """Generates the equations of motion for the simulation and
        generates a numeric function that will serve as the ode for
        ode45"""
        
        
        # First, we generate a "system position", "system velocity",
        # and "system acceleration" vector, which contain all the
        # information on each object's location and rotation.
        # At the same time, we generate a "system mass" matrix which
        # is diagonal and serves to keep track of the corresponding
        # mass (or moment of inertia) for each term in the system
        # position vector.
        self.r, self.v, self.a  = [], [], []
        self.M = []; self.dM = []

        for bod in self.bods:
            self.r.extend(list(bod.r))
            self.v.extend(list(bod.v))
            self.a.extend(list(bod.a))
            self.M.extend([bod.m]*len(bod.r))
            self.dM.extend([s.diff(bod.m,t)]*len(bod.r))
            try:
                self.r.append(bod.ang)
                self.v.append(bod.omega)
                self.a.append(bod.alpha)
                self.M.append(bod.I)
                self.dM.append(s.diff(bod.I,t))
            except AttributeError:
                pass

        self.r = s.Matrix(self.r)
        self.v = s.Matrix(self.v)
        self.a = s.Matrix(self.a)
        M = s.Matrix([[0]*len(self.r)]*len(self.r))
        dM = s.Matrix([[0]*len(self.r)]*len(self.r))
        for i,m in enumerate(self.M):
            M[i,i] = m
            dM[i,i] = dM[i]
        self.M = M
        self.dM = dM

        
        # Next, we generate the free directions of motion. One free
        # direction of motion exists for each degree of freedom, and
        # is found by taking the partial derivative of self.r with
        # respect to that degree of freedom. 
        self.dirs = s.Matrix([[0]*len(self.r)]*len(self.basis))
        for i,var in enumerate(self.basis):
            direction = s.Matrix([s.diff(component, var)
                                  for component in self.r])
            self.dirs[i,:] = direction.transpose()
            
        
        # Now we compile all the forces and moments into a "system 
        # force vector". First, we generate body-specific forces from
        # each system force, and then we compile the body-specific
        # forces
        for force in self.forces:
            #If it's from force2D, which for generality asks for a function
            # that generates an expression
            if force[0] == Force2D:
                for bod in self.bods:
                    bod.forces.append(force[0](*[arg(bod) for arg in force[1:]]))
            else:
                for bod in self.bods:
                    bod.forces.append(force[0](bod,*force[1:]))
        self.forces = [] #clear forces in case the sim is run twice
        
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
        
        
        # Now we write F=d(mv)/dt in our new coordinate system LHS is the
        # F part, RHS is the d(mv) part.
        LHS = s.simplify(self.dirs * self.F)
        RHS = s.simplify(self.dirs * (self.M * self.a + self.dM * self.v)) 
        
        # And put it into Vec = Mat*Diffs + Vec form
        self.dds = s.Matrix(self.ddbasis)
        self.RHMat = s.Matrix([[0]*len(self.basis)]*len(self.basis))
        for i, ddvar in enumerate(self.ddbasis):
            collected = [component.expand().collect(ddvar)
                         for component in RHS]
            coeffs = [component.coeff(ddvar)
                      for component in collected]
            remainder = [coll - coeff*ddvar for coll, coeff 
                         in zip(collected, coeffs)]

            RHS = s.simplify(s.Matrix(remainder))
            self.RHMat[i,:] = s.Matrix([s.simplify(coeff)
                                   for coeff in coeffs]).transpose()

        
        # Now we simplify it to LHS = RHMat * dds. This is our
        # final equation of motion
        self.LHS = LHS - RHS
        
        
        #
        # Now we generate the numeric ode function        
        #
        
        input_vars = [t] + self.basis + self.dbasis
        divider = len(self.basis) #between basis and dbasis
        
        #the replacement is done to avoid problems with backslashes
        newvars = s.symbols(['x_'+str(i) for i in range(0,len(input_vars))])
        #the flipping of the replacement variables is to avoid problems with
        #variables being replaced before their derivatives
        
        
        # Generate a numeric version of RHMat and LHS
        nRHMat = lambdify(newvars,
                          self.RHMat.subs(zip(input_vars[::-1],newvars[::-1])),
                          modules='numpy')
        nLHS = lambdify(newvars,
                        self.LHS.subs(zip(input_vars[::-1],newvars[::-1])),
                        modules='numpy')
        # And at the same time we generate numeric r & v functions
        self.nr = lambdify(newvars,
                           self.r.subs(zip(input_vars[::-1],newvars[::-1])),
                           modules='numpy')
        self.nv = lambdify(newvars,
                           self.v.subs(zip(input_vars[::-1],newvars[::-1])),
                           modules='numpy')
        self.nM = lambdify(newvars,
                           self.M.subs(zip(input_vars[::-1],newvars[::-1])),
                           modules='numpy')
        
        def odes(t,y):
            ddvs = n.linalg.inv(nRHMat(t,*y)) * nLHS(t,*y)
            return n.hstack((y[divider:],n.ravel(ddvs)))
        
        self.odes = odes


    def run(self, time, events=None):
        """Compiles and runs the simulation. Takes:
        time: a list of [t_start, t_end, t_step]
        events: an optional function that returns true to signal the
        end of the simulation"""

        self.compile()
                
        # Now we put the initial conditions into a computer-friendly
        # (instead of user-friendly) format
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
        
        # And run the simulation
        self.y = odesolve(self.odes, inits, time, events)
        
        # Now we calculate the position and velocity vectors as a nicety
        values = [self.y['t']]
        for varname, init_val in inits:
            values.append(self.y[varname])

        r = zip(*[self.nr(*value).transpose().tolist()[0]
                  for value in zip(*values)])
        v = zip(*[self.nv(*value).transpose().tolist()[0]
                  for value in zip(*values)])
        M = [self.nM(*value).tolist() for value in zip(*values)]
        M = zip(*[[m[i][i] for i in range(0,len(m))] for m in M])
        
        i = 0
        self.y['KE'] = n.zeros(len(self.y['t']))
        for bod in self.bods:
            self.y['r_'+bod.name] = [n.array(r[i]),n.array(r[i+1])]
            self.y['v_'+bod.name] = [n.array(v[i]),n.array(v[i+1])]
            self.y['KE_'+bod.name] = 0.5 * n.array(M[i]) * \
                                     (n.array(v[i])**2 + n.array(v[i+1])**2)
            i = i+2
            
            try:
                I = bod.I
                self.y['ang_'+bod.name] = n.array(r[i])
                self.y['omega_'+bod.name] = n.array(v[i])
                self.y['KE_'+bod.name] += n.array(M[i]) * n.array(v[i])**2
                i = i + 1
            except AttributeError:
                pass
            
            self.y['KE'] += self.y['KE_'+bod.name]
                
        return self.y



print('Pysics imported! You now can simulate rigid body dynamics.')
