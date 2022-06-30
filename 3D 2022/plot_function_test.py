from regex import F
import Vector_field_object as vf_obj

import numpy as np
import matplotlib.pyplot as plt
from mayavi import mlab

grid = np.linspace(-5,5,11)


xg, yg, zg = np.meshgrid(grid, grid, grid)


Fx = yg/np.sqrt(xg**2+yg**2)
Fy = -xg/np.sqrt(xg**2+yg**2)
Fz = 0*zg


fx = xg/np.sqrt(xg**2+yg**2+zg**2)
fy = yg/np.sqrt(xg**2+yg**2+zg**2)
fz = zg/np.sqrt(xg**2+yg**2+zg**2)

Ffx = xg
Ffy = yg
Ffz = zg


f3 = vf_obj.vector_field3(xg, yg, zg, fx, fy, fz)
F3 = vf_obj.vector_field3(xg, yg, zg, Fx, Fy, Fz)
Ff3 = vf_obj.vector_field3(xg, yg, zg, Ffx, Ffy, Ffz)

f3.plot()
F3.plot()
Ff3.plot()





