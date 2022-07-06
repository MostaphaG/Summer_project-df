from regex import F
import dformpy3D as vf_obj

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
f3.give_eqn('x/sqrt(x**2+y**2+z**2)','(-y)/sqrt(x**2+y**2+z**2)','z/sqrt(x**2+y**2+z**2)')
C3 = f3.curl()


f3.plot(add_curl='yeah')
C3.plot()


'''
Ff3 = vf_obj.vector_field3(xg, yg, zg, Ffx, Ffy, Ffz)

F3 = vf_obj.vector_field3(xg, yg, zg, Fx, Fy, Fz)


f3.plot()

Ff3.plot()


'''







