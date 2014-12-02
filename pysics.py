from __future__ import print_function, division
import sympy as s
import numpy as n
from tools import *


class PointMass2D(object):
    
    def __init__(self, name, m, r=None):
        self.name = name
        self.m = m

        if r is None:
            r = s.Matrix([s.Symbol('x_' + name),
                          s.Symbol('y_' + name)])
        
        self.constrain(r)

    
    def constrain(self, r):
        self.r = s.Matrix(r)


class RigidBody2D(object):
    def __init__(self, name, m, I, r=None, ang=None):
        self.name = name
        self.m = m
        self.I = I
        
        if r is None:
            r = s.Matrix([s.Symbol('x_' + name),
                          s.Symbol('y_' + name)])
        if ang is None:
            ang = s.Symbol('th_' + name)

        self.constrain(r, ang)

    def constrain(self, r=None, ang=None):
        if r is not None:
            self.r = s.Matrix(r)
        elif ang is not None:
            self.ang = ang

        
        # self.basis = list(extract_variables(self.r))
        # self.dbasis = [s.Symbol('d' + var.name) for var in self.basis]
        # self.ddbasis = [s.Symbol('dd' + var.name) for var in self.basis]
        
        # self.v, self.a = self._gen_derivs()
        # self.dirs = self._gen_free_dirs()
        

    # def _gen_derivs(self):
        
    #     r = self.r
    #     locs = self.basis
    #     vels = self.dbasis
    #     accs = self.ddbasis
        
    #     v = Matrix([0]*len(r))
    #     for i in range(0,len(locs)):
    #         subbed = r.subs(locs[i],vels[i]*t)
    #         diffed = Matrix([diff(component, t) for component in subbed])
    #         v += diffed.subs(vels[i]*t,locs[i])

    #     a = Matrix([0]*len(r))
    #     for i in range(0,len(locs)):
    #         subbed = v.subs(locs[i],vels[i]*t)
    #         diffed = Matrix([diff(component, t) for component in subbed])
    #         a += diffed.subs(vels[i]*t,locs[i])
    #         subbed = v.subs(vels[i],accs[i]*t)
    #         diffed = Matrix([diff(component, t) for component in subbed])
    #         a += diffed.subs(accs[i]*t,vels[i])
        
    #     return (v,a)
    
    
    # def _gen_free_dirs(self):
        
    #     r = self.r
    #     basis = self.basis
        
    #     dirs = Matrix([[0]*len(r)]*len(basis))
    #     for i,var in enumerate(basis):
    #         direction = Matrix([diff(component, var) for component in r])
    #         mag = sqrt(simplify(direction.dot(direction)))
    #         unit = simplify(direction/mag)
    #         dirs[i,:] = unit.transpose()
        
    #     return dirs


class Sim2D(object):
    
    def __init__(self):
        pass
    
    def PointMass(self, name, m, r=None):
        self.__dict__[name] = PointMass2D(name,m,r)
    
    def RigidBody(self, name, m, I, r=None, ang=None):
        self.__dict__[name] = RigidBody2D(name,m,r)
            

if __name__ == '__main__':

    sim = Sim2D()
    th = s.Symbol('th')
    #sim.PointMass('pend', 2)
    sim.RigidBody('rpend', 2, 3)
    #sim.pend.constrain(r=0.3*s.Matrix([s.sin(th),-s.cos(th)]))
    sim.rpend.constrain(r=0.3*s.Matrix([s.sin(th),-s.cos(th)]), ang=th)
