import Vector_field_object as vf_obj

import numpy as np
import matplotlib.pyplot as plt
from mayavi import mlab

grid = np.linspace(-5,5,5)


xg, yg, zg = np.meshgrid(grid, grid, grid)

potential = np.sqrt(xg**2 + yg**2 + zg**2)
expr = 'sqrt(x**2 + y**2 + z**2)'

f0 = vf_obj.form_0_3d(xg, yg, zg, potential, expr)
f0.set_density(100)
f0.log_scaling()
f0.plot(cross_sec_plane = 'y')
f0.plot()


