
from regex import F
import dformpy3D as df3

import numpy as np
import matplotlib.pyplot as plt
from mayavi import mlab

#grid = np.linspace(-5,5,11)
grid = np.linspace(-10,10,17)
#grid_0 = np.linspace(-150,150,25)

xg, yg, zg = np.meshgrid(grid, grid, grid)
#x0, y0, z0 = np.meshgrid(grid, grid, grid)


#fx = -yg/np.sqrt(xg**2+yg**2-zg**2)
#fy = xg/np.sqrt(xg**2+yg**2-zg**2)
#fz = zg/np.sqrt(xg**2+yg**2-zg**2)

#fx = -yg
#fy = xg
#fz = np.ones((np.size(xg)))

#potential = np.sqrt(x0**2 + y0**2 - z0**2)

#f3 = df3.vector_field3(xg, yg, zg, fx, fy, fz)
#f3.give_eqn('-y/sqrt(x**2+y**2-z**2)','x/sqrt(x**2+y**2-z**2)','z/sqrt(x**2+y**2-z**2)')
#c3 = f3.curl()
#c3.log_scaling()
#f3.plot(add_curl='y', scaling=0.3, arrow_colour='c', arrow_cmap='viridis', opacity=1.0, curl_opacity=0.2)
#c3.plot(scaling=1.0)

#f3.div(0, -2, 3)

#f0 = df3.form_0_3d(xg, yg, zg, potential)
#f0.give_eqn('sqrt(x**2 + y**2 - z**2)')
#f0.set_density(100)
#f0.levels(15)
#f0.log_scaling()
#f0.plot(cross_sec_plane='n')
#f0.plot(cross_sec_plane='y')

#form_1 = df3.form_1_3d(xg, yg, zg, fx, fy, fz)
#form_1.log_scaling()
#form_1.plot()

#fz = zg
#fx = xg
#fy = yg

#f2 = df3.form_2_3d(xg, yg, zg, Fz=fx, Fx=fy, Fy=fz)
#f2.plot()


pot = xg/np.sqrt(xg**2 + yg**2 + zg**2)
f3 = df3.form_3_3d(xg, yg, zg, pot)
f3.plot()

'''

Fx = yg/np.sqrt(xg**2+yg**2)
Fy = -xg/np.sqrt(xg**2+yg**2)
Fz = 0*zg

F3 = vf_obj.vector_field3(xg, yg, zg, Fx, Fy, Fz)
F3.give_eqn('y/sqrt(x**2+y**2)','-x/sqrt(x**2+y**2)','0')
C3 = F3.curl()


F3.plot(add_curl='y',arrow_colour='c',scaling=2.0)
C3.plot(arrow_cmap='jet', arrow_colour='c')




Ffx = xg
Ffy = yg
Ffz = zg



Ff3 = vf_obj.vector_field3(xg, yg, zg, Ffx, Ffy, Ffz)
F3 = vf_obj.vector_field3(xg, yg, zg, Fx, Fy, Fz)
f3.plot()
Ff3.plot()
'''