from pysics import *
from matplotlib.pyplot import *


sim = Sim2D()
sim.Force(lambda bod: [0, -bod.m*9.81])

th1 = DOF('th1')
th2 = DOF('th2')
th3 = DOF('th3')
pend1 = sim.PointMass('pend1', 1, s.Matrix([s.sin(th1),-s.cos(th1)]))
pend2= sim.PointMass('pend2', 1, pend1.r + \
                      s.Matrix([s.sin(th2),-s.cos(th2)]))
pend3= sim.PointMass('pend3', 1, pend2.r + \
                      s.Matrix([s.sin(th3),-s.cos(th3)]))
pend1.place({th1:(0.1,0)})
pend2.place({th2:(0.1,0)})
pend3.place({th3:(0.1,0)})

y = sim.run([0,5,0.01])


plot(y['t'],y['th1'])
plot(y['t'],y['th2'])
plot(y['t'],y['th3'])
show()
