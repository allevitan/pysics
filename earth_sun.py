from pysics import *

#Create the simulation object
sim = Sim2D()

xs,ys,xe,ye = sim.DOF('x_s','y_s','x_e','y_e')
earth = sim.PointMass('earth', 5.972e24, [xe,ye])
sun = sim.PointMass('sun', 1.989e30, [xs,ys])
sim.Potential(-6.67e-11 * sun.m * earth.m / 
                s.sqrt((earth.r - sun.r).dot(earth.r - sun.r)))
sim.place({xs: (0,0), ys: (0,0),
           xe: (147098070000, 0),
           ye: (0,30300)})
y = sim.run([0, 31536000, 3600*24])

sim.analyze()
