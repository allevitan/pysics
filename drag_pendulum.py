from pysics import *
from matplotlib.pyplot import *

sim = Sim2D()
sim.Gravity([0,-9.81])

th = sim.DOF('th')
pend = sim.PointMass('pend', 2, 0.3*s.Matrix([s.sin(th),-s.cos(th)]))
pend.Drag(LCd=1, power=1) #viscous hinge drag
#pend.Drag(LCd=10,power=2) #aerodynamic body drag

sim.place({th:(0.1,0)})

y = sim.run([0,20,0.01])


plot(y['t'],y['th'])
show()

