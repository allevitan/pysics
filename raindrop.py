from pysics import *

#Create the simulation object
sim = Sim2D()
sim.Gravity([0,-9.81])

#Define the constraints
x,y = sim.DOF('x','y')
pend = sim.PointMass('pend', 0.01 + 0.01*t, [x,y])
                      
#Set the initial conditions
sim.place({x: (0, 0), y: (0,0)})

y = sim.run([0,7,0.01])

sim.analyze()
