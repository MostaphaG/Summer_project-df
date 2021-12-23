# -*- coding: utf-8 -*-
"""
Created on Thu Dec 23 10:47:41 2021

@author: single user
"""

import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(8,8))
ax = fig.gca()
r = np.linspace(-1,5,21)
xg, yg = np.meshgrid(r,r)
u = xg*np.cos(yg)
v = -yg*np.sin(xg)
ax.set_aspect('equal')

plt.quiver(xg,yg,u,v)

x_m, y_m = (1,1)
x_range = xg[0,-1] - xg[0,0]
y_range = yg[-1,0] - yg[0,0]

xi = (x_m - xg[0,0])/x_range
yi = (y_m - yg[0,0])/y_range

inax = ax.inset_axes([xi, yi, 0.3, 0.3])