from turtle import color
from matplotlib import projections, scale
import numpy as np
from Vector_Fields_Library import vector_field as vf 
import matplotlib.pyplot as plt



crd_cart = vf.coordinate_list_cartesian(-4, 4, 2)
crd_sph = vf.coordinate_list_spherical(15,6,6,4)
crd_cyl = vf.coordinate_list_cylindrical(15, 4, 8, -5, 5 ,3)

F1 = vf.vfield('1','1','0', crd_cart)
C1 = vf.curl('1','1','0', crd_cart)

F = vf.vfield('y/sqrt(x**2+y**2)','-x/sqrt(x**2+y**2)','0', crd_cart)
C = vf.curl('y/sqrt(x**2+y**2)','-x/sqrt(x**2+y**2)','0', crd_cart)




fig = plt.figure()
ax = fig.add_subplot(121, projection='3d')
ax.quiver(crd_cart[0], crd_cart[1], crd_cart[2], F[0], F[1], F[2],
             color = 'r')
ax.quiver(crd_cart[0], crd_cart[1], crd_cart[2], C[0], C[1], C[2],
             color = 'blue')
ax1 = fig.add_subplot(122, projection='3d')
ax1.quiver(crd_cart[0], crd_cart[1], crd_cart[2], F1[0], F1[1], F1[2],
             color = 'r')
ax1.quiver(crd_cart[0], crd_cart[1], crd_cart[2], C1[0], C1[1], C1[2],
             color = 'blue')
plt.show()









