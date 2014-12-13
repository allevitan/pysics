from pysics import *

sim = Sim2D()
sim.Gravity([0,-9.81])

th = sim.DOF('th')
pend = sim.RigidBody('pend', 2, 0.5, [0.3*s.sin(th),-0.3*s.cos(th)], th)
sim.place({th:(0.1,0)})

y = sim.run([0,5,0.01])

sim.analyze()
