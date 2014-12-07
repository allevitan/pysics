from pysics import *
from matplotlib.pyplot import *

sim = Sim2D()
sim.Gravity([0,-9.81])

th1, th2 = sim.DOF('th1', 'th2')
pend1 = sim.RigidBody('pend1', 1, 0.2, s.Matrix([0.4*s.sin(th1),-0.4*s.cos(th1)]), th1)
pend2 = sim.RigidBody('pend2', 0.7, 0.1, pend1.r + 
                      s.Matrix([0.3*s.sin(th2),-0.3*s.cos(th2)]), th2)
sim.place({th1: (n.pi/3,0), th2: (n.pi/4,0)})


y = sim.run([0,5,0.01])


#plot(y['t'],y['th1'])
#plot(y['t'],y['th2'])
plot(*y['r_pend2'])
show()
