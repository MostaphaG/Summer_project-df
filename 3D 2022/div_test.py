from turtle import color
from matplotlib import projections, scale
import numpy as np
from Vector_Fields_Library import vector_field as vf 
import matplotlib.pyplot as plt



crd_cart = vf.coordinate_list_cartesian(-10, 10, 5)
crd_sph = vf.coordinate_list_spherical(15,6,6,4)
crd_cyl = vf.coordinate_list_cylindrical(15, 4, 8, -5, 5 ,3)


D1 = vf.div('(x)/sqrt(x**2+y**2+z**2)',
            '(y)/sqrt(x**2+y**2+z**2)',
            '(z)/sqrt(x**2+y**2+z**2)', crd_cart, 2, 3, 1)


