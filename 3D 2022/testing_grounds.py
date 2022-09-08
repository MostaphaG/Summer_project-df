
from regex import F
import dformpy3D as df3
import dformpy as df

import numpy as np
import matplotlib.pyplot as plt
from mayavi import mlab

#grid = np.linspace(-5,5,11)
grid = np.linspace(-1,1,3)
#grid_0 = np.linspace(-150,150,25)

xg, yg, zg = np.meshgrid(grid, grid, grid)
xg2, yg2= np.meshgrid(grid, grid)
#x0, y0, z0 = np.meshgrid(grid, grid, grid)


fx = -yg/np.sqrt(xg**2+yg**2)
fy = xg/np.sqrt(xg**2+yg**2)
fz = 0*zg


fx2 = -yg2/np.sqrt(xg2**2+yg2**2)
fy2 = xg2/np.sqrt(xg2**2+yg2**2)


potential = np.sqrt(xg**2 + yg**2 + zg**2)
potential2 = np.sqrt(xg2**2 + yg2**2)


#vf = df3.vector_field3(xg, yg, zg, fx, fy, fz)
#vf.autoscale()
#vf.give_eqn('x/sqrt(x**2 + y**2 - z**2)','y/sqrt(x**2 + y**2 - z**2)','z/sqrt(x**2 + y**2 - z**2)')
#cov_f1 = vf.covariant()
#vff = vf.zoom(mag = 5, target=[0,0,-8], dpd = 6)
#cvf = vf.curl()
#cvf.log_scaling()
#vf.plot(add_curl='n', scaling=0.1, arrow_colour='c', arrow_cmap='viridis', opacity=1.0, curl_opacity=0.2)
#cov_f1.plot()
#cvf.plot(scaling=1.0)
#deriv_field = vf.deriv(target=[2,2,0], mag=2, dpd=5)
#deriv_field.plot()
#vf.div(0, -2, 3)







#3f0 = df3.form_0_3d(xg, yg, zg, potential)
#f0.give_eqn('sqrt(x**2 + y**2 + z**2)')
#f00 = f0.zoom(mag=2, target=[-5,0,0], dpd=100)
#f0.set_density(100)
#f0.levels(15)
#f00.log_scaling()
#f0.plot(cross_sec_plane='y')
#f0_ext_d = f0.ext_d()
#f0_ext_d_num = f0.num_ext_d()
#f0_ext_d.set_density(8)
#f0_ext_d.plot()
#f0_ext_d_num.plot()
#f0.plot(cross_sec_plane='y')
#f00.plot(cross_sec_plane='y')

#----compare difference between analytic/num 2d & 3d-----------------------

#f0_2d = df.form_0(xg2, yg2, potential2)
#f0_2d.give_eqn('sqrt(x**2 + y**2)')
#f0_2d_ext_d = f0_2d.ext_d()
#f0_2d_ext_d_num = f0_2d.num_ext_d()

#fig = plt.figure()
#ax = fig.add_subplot(121)
#ax1 = fig.add_subplot(122)
#ax.set_aspect('equal')
#ax1.set_aspect('equal')


#f0_2d_ext_d.plot(axis=ax)
#f0_2d_ext_d_num.plot(axis=ax1)
#plt.show()

#--------------------------------------------------------------------------







#form_1 = df3.form_1_3d(xg, yg, zg, fx, fy, fz)
#form_1.give_eqn('-y/sqrt(x**2+y**2)','x/sqrt(x**2+y**2)','0')
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
#form_1.plot()

#----compare difference between analytic/num 2d & 3d-----------------------

#f1_2d = df.form_1(xg2, yg2, fx2, fy2)
#f1_2d.give_eqn('-y/sqrt(x**2 + y**2)','x/sqrt(x**2 + y**2)')
#f1_2d_ext_d = f1_2d.ext_d()
#f1_2d_ext_d_num = f1_2d.num_ext_d()

#fig = plt.figure()
#ax = fig.add_subplot(121)
#ax1 = fig.add_subplot(122)
#ax.set_aspect('equal')
#ax1.set_aspect('equal')

#f1_2d_ext_d.plot(axis=ax)
#f1_2d_ext_d_num.plot(axis=ax1)
#plt.show()

#--------------------------------------------------------------------------














#f2 = df3.form_2_3d(xg, yg, zg, Fz=fz, Fx=fx, Fy=fy)
#f2.give_eqn('x/sqrt(x**2+y**2-z**2)','y/sqrt(x**2+y**2-z**2)','z/sqrt(x**2+y**2-z**2)')
#f2.log_scaling()
#f22 = f2.zoom(mag = 5, target=[2,2,2], dpd=8)
#f2.plot()
#f2_ext_d = f2.ext_d()
#f2_ext_d.plot()
#f2_num_ext_d = f2.num_ext_d()
#f2_num_ext_d.plot()
#hodged_f2 = f2.hodge()
#hodged_f2.plot()
#hodged_f2_num = f2.num_hodge()
#hodged_f2_num.plot()
#f2_intd = f2.interior_d(vf)
#f2_intd.plot()
#f2_intd_num = f2.num_interior_d(vf)
#f2_intd_num.plot()


pot = (xg**2 + yg**2 + zg**2)
f3 = df3.form_3_3d(xg, yg, zg, potential)
f3.give_eqn('sqrt(x**2+y**2-z**2)')
#f33 = f3.zoom(mag = 5, target=[2,2,2], dpd=8)
#f3.log_scaling()
#f33.plot()
#f3_intd = f3.interior_d(vf)
#f3_intd.plot()
#f3_intd_num = f3.num_interior_d(vf)
#f3_intd_num.plot()
f3_f0 = f3.hodge()
f3_f0.plot()


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



#f2_wedge_f0 = f2.wedge(f0)
#f2_wedge_f0.plot()
#f2_wedge_f1 = f2.wedge(form_1)
#f2_wedge_f1.plot()
#f2_wedge_f2 = f2.wedge(f2)
#f2_wedge_f3 = f2.wedge(f3)

#f2_wedge_f0_num = f2.num_wedge(f0)
#f2_wedge_f0_num.plot()
#f2_wedge_f1_num = f2.num_wedge(form_1)
#f2_wedge_f1_num.plot()
#f2_wedge_f2_num = f2.num_wedge(f2)
#f2_wedge_f3_num = f2.num_wedge(f3)



#f3_wedge_f0 = f3.wedge(f0)
#f3_wedge_f0.plot()
#f3_wedge_f1 = f3.wedge(form_1)
#f3_wedge_f2 = f3.wedge(f2)
#f3_wedge_f3 = f3.wedge(f3)

#f3_wedge_f0_num = f3.num_wedge(f0)
#f3_wedge_f0_num.plot()
#f3_wedge_f1_num = f3.num_wedge(form_1)
#f3_wedge_f2_num = f3.num_wedge(f2)
#f3_wedge_f3_num = f3.num_wedge(f3)



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