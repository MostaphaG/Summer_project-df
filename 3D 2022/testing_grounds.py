from turtle import color
from matplotlib import projections, scale
import numpy as np
from Vector_Fields_Library import vector_field as vf 
import matplotlib.pyplot as plt
from mayavi.mlab import *

crd_cart = vf.coordinate_list_cartesian(-10, 10, 3)
crd_sph = vf.coordinate_list_spherical(15,6,6,3)
crd_cyl = vf.coordinate_list_cylindrical(15, 4, 8, -15, 15 ,3)

F1 = vf.vfield('(x)/sqrt(x**2+y**2+z**2)',
               '(y)/sqrt(x**2+y**2+z**2)',
               '(z)/sqrt(x**2+y**2+z**2)', crd_sph)
C1 = vf.curl('(x)/sqrt(x**2+y**2+z**2)',
             '(y)/sqrt(x**2+y**2+z**2)',
             '(z)/sqrt(x**2+y**2+z**2)', crd_sph)

F = vf.vfield('y/sqrt(x**2+y**2)','-x/sqrt(x**2+y**2)','0', crd_cyl)
C = vf.curl('y/sqrt(x**2+y**2)','-x/sqrt(x**2+y**2)','0', crd_cyl)



vf.plot(F, C)









