## Synopsis

Pysics is a simple python module capable of simulating arbitrary 2D multibody mechanical systems obeying holonomic constraints, including systems with time-dependent constraints and forcing. It was written as a final project for the sophomore-level Dynamics class at Olin College of Engineering. Included in this repository are approximately 10 example simulations demonstrating the range of systems that pysics is capable of simulating.

## Dependencies

Pysics is dependent on several scientific computing packages for python:

* **Numpy**, for efficient array math
* **Sympy**, for symbolic manipulations
* **Scipy**, for numeric integration
* **Matplotlib**, for creationg of figures

Pysics should work both in Python 2.7 and 3.x, provided the most recent versions of the above packages are installed.

## Code Example

Pysics takes great pains to only request the minimum possible amount of information about the system. For example, the snippet of code below defines a double solid pendulum with aerodynamic drag.

```python
from pysics import *

#Create the simulation object and define gravity
sim = Sim2D()
sim.Gravity([0, -9.81])

#Define the constraints
th1, th2 = sim.DOF('th1','th2')
r1 = 0.4*s.Matrix([s.sin(th1),-s.cos(th1)])
r2 = r1 + 0.3*s.Matrix([s.sin(th2),-s.cos(th2)])

#And define the bodies and forces
pend1 = sim.RigidBody('pend1', 1, 0.3, r1, th1)
pend2= sim.RigidBody('pend2', 0.7, 0.15, r2, th2)
sim.Drag(TCd=0.1, power=2)
                      
#Set the initial conditions
sim.place({th1: (n.pi/3,0), th2: (n.pi/4,0)})

#And simulate! Badabing badaboom
y = sim.run([0,5,0.01])

#Now we generate some figures to look at
sim.analyze()
```

## Installation

Prior to installation of pysics, Python must be installed with numpy, sympy, scipy, and matplotlib. Pysics may be installed system-wide by moving the pysics folder into your system's pythonpath, or simply used on a project-by-project basis by placing a copy of the folder into each project's directory.
