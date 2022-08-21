
from regex import F
import dformpy3D as df3

import numpy as np
import matplotlib.pyplot as plt
from mayavi import mlab

#grid = np.linspace(-5,5,11)
grid = np.linspace(-2,2,8)
#grid_0 = np.linspace(-150,150,25)

xg, yg, zg = np.meshgrid(grid, grid, grid)
#x0, y0, z0 = np.meshgrid(grid, grid, grid)


fx = xg/np.sqrt(xg**2+yg**2-zg**2)
fy = yg/np.sqrt(xg**2+yg**2-zg**2)
fz = zg/np.sqrt(xg**2+yg**2-zg**2)


potential = np.sqrt(xg**2 + yg**2 - zg**2)

vf = df3.vector_field3(xg, yg, zg, fx, fy, fz)
#vf.autoscale()
vf.give_eqn('x/sqrt(x**2 + y**2 - z**2)','y/sqrt(x**2 + y**2 - z**2)','z/sqrt(x**2 + y**2 - z**2)')
#cov_f1 = vf.covariant()
#vff = vf.zoom(mag = 5, target=[0,0,-8], dpd = 6)
#cvf = vf.curl()
#cvf.log_scaling()
#vf.plot(add_curl='n', scaling=0.1, arrow_colour='c', arrow_cmap='viridis', opacity=1.0, curl_opacity=0.2)
#cov_f1.plot()
#cvf.plot(scaling=1.0)
#deriv_field = vf.deriv(target=[2,2,0], mag=2, dpd=5)
#deriv_field.plot()

#f3.div(0, -2, 3)

f0 = df3.form_0_3d(xg, yg, zg, potential)
f0.give_eqn('sqrt(x**2 + y**2 - z**2)')
#f00 = f0.zoom(mag=2, target=[-5,0,0], dpd=100)
#f0.set_density(10)
#f0.levels(15)
#f00.log_scaling()
#f0.plot(cross_sec_plane='n')
#f0_ext_d = f0.ext_d()
#f0_ext_d_num = f0.num_ext_d()
#f0_ext_d.set_density(8)
#f0_ext_d.plot()
#f0_ext_d_num.plot()
#f0.plot(cross_sec_plane='y')
#f00.plot(cross_sec_plane='y')

form_1 = df3.form_1_3d(xg, yg, zg, fx, fy, fz)
form_1.give_eqn('x/sqrt(x**2+y**2-z**2)','y/sqrt(x**2+y**2-z**2)','z/sqrt(x**2+y**2-z**2)')
#form_11 = form_1.zoom(mag = 2.8, target=[0,0,10], dpd=4)
#form_1.log_scaling()
#contravariant_field = form_1.contravariant()
#contravariant_field.plot(scaling=0.1)
#f1_extd = form_1.ext_d()
#f1_extd_num = form_1.num_ext_d()
#f1_extd.plot()
#f1_extd_num.plot()
#hodged_f1 = form_1.hodge()
#hodged_f1.plot()
#hodged_f1_num = form_1.num_hodge()
#hodged_f1_num.plot()
#f1_intd = form_1.interior_d(vf)
#f1_intd.plot()
#f1_intd_num = form_1.num_interior_d(vf)
#f1_intd_num.plot()


f2 = df3.form_2_3d(xg, yg, zg, Fz=fz, Fx=fx, Fy=fy)
f2.give_eqn('x/sqrt(x**2+y**2-z**2)','y/sqrt(x**2+y**2-z**2)','z/sqrt(x**2+y**2-z**2)')
#f2.log_scaling()
#f22 = f2.zoom(mag = 5, target=[2,2,2], dpd=8)
#f2.plot()


#pot = (xg**2 + yg**2 + zg**2)
f3 = df3.form_3_3d(xg, yg, zg, potential)
f3.give_eqn('sqrt(x**2+y**2-z**2)')
#f33 = f3.zoom(mag = 5, target=[2,2,2], dpd=8)
#f3.log_scaling()
#f33.plot()

#f0_wedge_f0 = f0.wedge(f0)
#f0_wedge_f0.plot()
#f0_wedge_f1 = f0.wedge(form_1)
#f0_wedge_f1.plot()
#f0_wedge_f2 = f0.wedge(f2)
#f0_wedge_f2.plot()
#f0_wedge_f3 = f0.wedge(f3)
#f0_wedge_f3.plot()

#f0_wedge_f0 = f0.num_wedge(f0)
#f0_wedge_f0.plot()
#f0_wedge_f1 = f0.num_wedge(form_1)
#f0_wedge_f1.plot()
#f0_wedge_f2 = f0.num_wedge(f2)
#f0_wedge_f2.plot()
#f0_wedge_f3 = f0.num_wedge(f3)
#f0_wedge_f3.plot()



#f1_wedge_f0 = form_1.wedge(f0)
#f1_wedge_f0.plot()
#f1_wedge_f1 = form_1.wedge(form_1)
#f1_wedge_f1.plot()
#f1_wedge_f2 = form_1.wedge(f2)
#f1_wedge_f2.plot()
#f1_wedge_f3 = form_1.wedge(f3)

#f1_num_wedge_f0 = form_1.num_wedge(f0)
#f1_num_wedge_f0.plot()
#f1_num_wedge_f1 = form_1.num_wedge(form_1)
#f1_num_wedge_f1.plot()
#f1_num_wedge_f2 = form_1.num_wedge(f2)
#f1_num_wedge_f2.plot()
#f1_num_wedge_f3 = form_1.num_wedge(f3)


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