from pysics import *
from matplotlib.pyplot import *


sim = Sim2D()
sim.Gravity([0,-9.81])

x,y = sim.DOF('x','y')
ball = sim.PointMass('ball', 1, [x,y])
ball.Drag(TCd=1, power=2)
sim.place({x: (0,10), y:(0,20)})

y = sim.run([0,5,0.01],events=lambda y: y[1] < 0)

sim.analyze()
#plot(y['x'],y['y'])
#show()
