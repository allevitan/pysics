from pysics import *
from matplotlib.pyplot import *

m = 1
B = 0
h = 1
L = 1

sim = Sim2D()
sim.Force(lambda bod: [0, -bod.m*9.81])

x = DOF('x')
pend = sim.PointMass('pend', m, [x,h*(1-B* x/L) * s.cos(3*n.pi*x/(2*L))**2])
pend.place({x:(0,0.1)})

y = sim.run([0,5,0.01], events=lambda y: y[0] >= L)


plot(y['t'],y['x'])
show()

