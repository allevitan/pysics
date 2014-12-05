from pysics import *
from matplotlib.pyplot import *

sim = Sim2D()
sim.Gravity([0,-9.81])
sim.Drag(1,1)
sim.Force(lambda bod: [s.cos(t),0])


th = sim.DOF('th')
pend = sim.PointMass('pend', 2, 0.3*s.Matrix([s.sin(th),-s.cos(th)]))
sim.place({th:(0.1,0)})

y = sim.run([0,20,0.01])


plot(y['t'],y['th'])
show()

