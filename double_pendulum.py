from pysics import *

#Create the simulation object and define gravity
sim = Sim2D()
sim.Gravity([0, -9.81])

#Define the constraints
th1, th2 = sim.DOF('th1','th2')
r1 = 0.4*s.Matrix([s.sin(th1),-s.cos(th1)])
r2 = r1 + 0.3*s.Matrix([s.sin(th2),-s.cos(th2)])

#And make the objects in the simulation
pend1 = sim.PointMass('pend1', 1, r1)
pend2= sim.PointMass('pend2', 0.7, r2)
                      
#Set the initial conditions
sim.place({th1: (n.pi/3,0), th2: (n.pi/4,0)})

#And simulate! Badabing badaboom
y = sim.run([0,5,0.01])

sim.analyze()

