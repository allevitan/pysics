from pysics import *
from matplotlib.pyplot import *

sim = Sim2D()
sim.Force(lambda bod: [0, -bod.m*9.81])

th = sim.DOF('th')
pend = sim.PointMass('pend', 2, 0.3*s.Matrix([s.sin(th),-s.cos(th)]))
pend.place({th:(0.1,0)})

y = sim.run([0,5,0.01])


plot(y['t'],y['th'])
show()

