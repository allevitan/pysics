from pysics import *
from matplotlib.pyplot import *


sim = Sim2D()
sim.Force(lambda bod: [0, -bod.m*9.81])

th1 = DOF('th1')
th2 = DOF('th2')
pend1 = sim.PointMass('pend1', 1, 0.4*s.Matrix([s.sin(th1),-s.cos(th1)]))
pend2= sim.PointMass('pend2', 0.7, pend1.r + \
                      0.3*s.Matrix([s.sin(th2),-s.cos(th2)]))
pend1.place({th1:(n.pi/3,0)})
pend2.place({th2:(n.pi/4,0)})

y = sim.run([0,5,0.01])


plot(y['t'],y['th1'])
plot(y['t'],y['th2'])
show()
