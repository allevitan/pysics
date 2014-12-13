from pysics import *
from matplotlib.pyplot import *

sim = Sim2D()
sim.Gravity([0,-9.81])

th = sim.DOF('th')
pend = sim.PointMass('pend', 2, 0.3*s.Matrix([1+s.sin(th),1-s.cos(th)]))
sim.place({th:(0.1,0)})

y = sim.run([0,5,0.01])


sim.analyze()

