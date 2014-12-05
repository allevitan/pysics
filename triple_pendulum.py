from pysics import *
from matplotlib.pyplot import *


sim = Sim2D()
sim.Gravity([0,-9.81])

th1,th2,th3 = sim.DOF('th1','th2','th3')
pend1 = sim.PointMass('pend1', 1, s.Matrix([s.sin(th1),-s.cos(th1)]))
pend2= sim.PointMass('pend2', 1, pend1.r + \
                    s.Matrix([s.sin(th2),-s.cos(th2)]))
pend3= sim.PointMass('pend3', 1, pend2.r + \
                      s.Matrix([s.sin(th3),-s.cos(th3)]))
sim.place({th1:(0.8,0), th2:(0.7,0), th3:(0.9,0)})

y = sim.run([0,5,0.01])


plot(y['t'],y['th1'])
plot(y['t'],y['th2'])
plot(y['t'],y['th3'])
show()
