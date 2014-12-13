from pysics import *

w = 3 #rad/s
R = 0.5 #m

#Create the simulation object
sim = Sim2D()

#Define the constraints
th = sim.DOF('th')
r = R*s.Matrix([s.sin(w*t),-s.cos(w*t)]) + \
    R*s.Matrix([s.sin(th + w*t),-s.cos(th + w*t)])

#And make the objects in the simulation
pend = sim.PointMass('pend', 1, r)
                      
#Set the initial conditions
sim.place({th: (1, 0.01)})

#And simulate! Badabing badaboom
y = sim.run([0,5,0.01])

sim.analyze()
