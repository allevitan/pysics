from pysics import *

mc = 1
mp = 1
r = 1

#Create the simulation object
sim = Sim2D()
sim.Gravity([0,-9.81])

x, th1, th2 = sim.DOF('x','th1','th2')

cart = sim.PointMass('cart', mc, [x,0])
p1 = sim.PointMass('p1', mp, cart.r + r*s.Matrix([s.sin(th1),-s.cos(th1)]))
p2 = sim.PointMass('p2', mp, p1.r + r*s.Matrix([s.sin(th2),-s.cos(th2)]))
#sim.Drag(TCd=0.02,power=2)

sim.place({x: (0,0), th1: (1,0), th2: (0.5,0)})

y = sim.run([0, 1.45, 0.01])

sim.analyze()


# This plots a slightly prettier version of the trajectory
# p.plot(y['r_cart'][0], y['r_cart'][1], 'b--')
# p.plot(y['r_cart'][0][-1], y['r_cart'][1][-1], 'bo')
# p.plot(y['r_p1'][0], y['r_p1'][1], 'g--')
# p.plot(y['r_p1'][0][-1], y['r_p1'][1][-1], 'go')
# p.plot((y['r_cart'][0][-1],y['r_p1'][0][-1]),
#      (y['r_cart'][1][-1],y['r_p1'][1][-1]), 'k-')
# p.plot(y['r_p2'][0], y['r_p2'][1], 'r--')
# p.plot(y['r_p2'][0][-1], y['r_p2'][1][-1], 'ro')
# p.plot((y['r_p1'][0][-1],y['r_p2'][0][-1]),
#      (y['r_p1'][1][-1],y['r_p2'][1][-1]), 'k-')
# p.title('Double Pendulum on a Track with Drag')
# p.ylim([-2,0.2])
# p.xlim([-0.4,1.8])
# p.show()
