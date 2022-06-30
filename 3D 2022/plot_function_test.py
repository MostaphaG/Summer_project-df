import Vector_field_object as vf_obj

import numpy as np
import matplotlib.pyplot as plt
from mayavi import mlab

grid = np.linspace(-5, 5,10)


xg, yg, zg = np.meshgrid(grid, grid, grid)




Fx = xg
Fy = yg
Fz = zg

fx = yg/np.sqrt(xg**2+yg**2)
fy = -xg/np.sqrt(xg**2+yg**2)
fz = 0*zg


F3 = vf_obj.vector_field3(xg, yg, zg, Fx, Fy, Fz)
f3 = vf_obj.vector_field3(xg, yg, zg, fx, fy, fz)


#No need to define axis!

F3.plot()
f3.plot()


