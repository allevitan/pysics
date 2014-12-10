from __future__ import division, print_function
from pysics import *
from matplotlib.pyplot import *

h = 1
B = 0.3
L = 1

sim = Sim2D()
sim.Gravity([0,-9.81])

x = sim.DOF('x')
cart = sim.PointMass('cart', 1, [x, h * (1- B*x / L) * s.cos(3*n.pi*x/(2*L))**2])
cart.Drag(TCd=0.05, power=2)
sim.place({x: (0.01,0)})

y = sim.run([0,8,0.01], lambda y: y[0] > L)


#plot(y['t'],y['x'])
plot(*y['r_cart'])
show()

