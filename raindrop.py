from pysics import *

sim = Sim2D()
sim.Gravity([0,-9.81])

x,y = sim.DOF('x','y')
m = 0.01 + 0.01*t
pend = sim.PointMass('pend', m, [x,y])
                      
sim.place({x: (0, 0), y: (0,0)})

y = sim.run([0,7,0.01])

sim.analyze()
