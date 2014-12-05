from pysics import *
from matplotlib.pyplot import *

m = 1
B = 0
h = 1
L = 1

sim = Sim2D()
sim.Gravity([0, -9.81])

x = sim.DOF('x')
cart = sim.PointMass('cart', m, [x,h*(1-B* x/L) * s.cos(3*n.pi*x/(2*L))**2])
sim.place({x:(0,0.1)})

y = sim.run([0,5,0.01], events=lambda y: y[0] >= L)


plot(y['t'],y['x'])
show()

