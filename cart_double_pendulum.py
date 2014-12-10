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
#sim.Drag(TCd=0.02,power=2)

sim.place({x: (0,0), th1: (1,0), th2: (0.5,0)})

y = sim.run([0, 1.45, 0.01])


plot(y['r_cart'][0], y['r_cart'][1], 'b--')
plot(y['r_cart'][0][-1], y['r_cart'][1][-1], 'bo')
plot(y['r_p1'][0], y['r_p1'][1], 'g--')
plot(y['r_p1'][0][-1], y['r_p1'][1][-1], 'go')
plot((y['r_cart'][0][-1],y['r_p1'][0][-1]),
     (y['r_cart'][1][-1],y['r_p1'][1][-1]), 'k-')
plot(y['r_p2'][0], y['r_p2'][1], 'r--')
plot(y['r_p2'][0][-1], y['r_p2'][1][-1], 'ro')
plot((y['r_p1'][0][-1],y['r_p2'][0][-1]),
     (y['r_p1'][1][-1],y['r_p2'][1][-1]), 'k-')
title('Double Pendulum on a Track with Drag')
ylim([-2,0.2])
xlim([-0.4,1.8])
figure()


#plot(y['t'],y['x'])
#plot(y['t'],y['th1'])
#plot(y['t'],y['th2'])
show()
