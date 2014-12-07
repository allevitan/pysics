from pysics import *
from matplotlib.pyplot import *

mc = 1
mp = 1
r = 1

#Create the simulation object
sim = Sim2D()
sim.Gravity([0,-9.81])

#
x, th1, th2 = sim.DOF('x','th1','th2')

cart = sim.PointMass('cart', mc, [x,0])
p1 = sim.PointMass('p1', mp, cart.r + r*s.Matrix([s.sin(th1),-s.cos(th1)]))
p2 = sim.PointMass('p2', mp, p1.r + r*s.Matrix([s.sin(th2),-s.cos(th2)]))

sim.place({x: (0,0), th1: (1,0), th2: (0.5,0)})

y = sim.run([0, 5, 0.01])

plot(y['t'],y['x'])
plot(y['t'],y['th1'])
plot(y['t'],y['th2'])
show()
