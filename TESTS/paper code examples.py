import dformpy as fp
import numpy as np
import matplotlib.pyplot as plt

# Example code for Vector field zooming

# Set up the Vector Field
r = np.linspace(-5, 5, 21)
xg, yg = np.meshgrid(r, r)  # axis grids
u = yg*np.sin(xg)
v = -xg*np.cos(yg)   # components
VF = fp.vector_field(xg, yg, u, v)  # initialise instance
VF.give_eqn('y*sin(x)','-x*cos(y)')  # supply equations to it

# set up figure
fig = plt.figure(figsize=(10, 10))
ax = fig.gca()
ax.set_aspect('equal')
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')

# plot field
VF.plot(ax)

# zoom onto it
zoomed_ax, zoomed_VF = VF.zoom(target=[2, 3], mag=1.25, dpd=7,
                               inset=True, axis=ax, insize=0.25)

# customise insset plot
zoomed_VF.colour('r')
zoomed_ax.set_yticks(np.linspace(zoomed_VF.yg[0, 0], zoomed_VF.yg[-1, 0], 5))
zoomed_ax.set_xticks(np.linspace(zoomed_VF.xg[0, 0], zoomed_VF.xg[0, -1], 5))

# replot it
zoomed_VF.plot(zoomed_ax)


# %%

import dformpy as fp
import numpy as np
import matplotlib.pyplot as plt

# Example code for 1-form zooming

# Set up the 1-form
r = np.linspace(-5, 5, 21)
xg, yg = np.meshgrid(r, r)  # axis grids
u = yg*np.sin(xg)
v = -xg*np.cos(yg)   # components
f1 = fp.form_1(xg, yg, u, v)  # initialise instance
f1.give_eqn('y*sin(x)', '-x*cos(y)')  # supply equations to it

# set up figure
fig = plt.figure(figsize=(6, 6))
ax = fig.gca()
ax.set_aspect('equal')
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')

# plot the 1-form
f1.plot(ax)

# zoom onto it
ax_zoom, f_zoom = f1.zoom(target=[2, 3], mag=1, dpd=7,
                             inset=True, axis=ax, insize=0.3)

# customise insset plot
f_zoom.colour('r')
ax_zoom.set_yticks(np.linspace(f_zoom.yg[0, 0], f_zoom.yg[-1, 0], 5))
ax_zoom.set_xticks(np.linspace(f_zoom.xg[0, 0], f_zoom.xg[0, -1], 5))

# replot it
f_zoom.plot(ax_zoom)

# %%

import dformpy as fp
import numpy as np
import matplotlib.pyplot as plt

# Example code for exterior derivative of Yukawa 0-form

# Set up the axis grids
r = np.linspace(-1, 1, 27)
xg, yg = np.meshgrid(r, r)
rg = np.sqrt(xg**2 + yg**2)

# set up 0-form instance
scalar = 1/rg * np.exp(-rg)
f0 = fp.form_0(xg, yg, scalar)

# supply equation
f0.give_eqn('1/sqrt(x**2 + y**2) * e**(-sqrt(x**2 + y**2))')

# get 1-form via exterior derivative and customise
dphi = f0.ext_d()
dphi.log_scaling()
dphi.sheet_size(0.04)

# customise 0-from to plot
f0.levels(50)  # number of level lines
f0.density_increase(3)  # more definition

# set up figure and axis
fig = plt.figure(figsize=(10, 10))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax1.set_aspect('equal')
ax2.set_aspect('equal')

# plot
f0.plot(ax1)
dphi.plot(ax2)

# %%

import dformpy as fp
import numpy as np
import matplotlib.pyplot as plt

# 2D Black hole

# set up grids
x = np.linspace(-4, 4, 26)
t = np.linspace(0, 4, 19)
xg, tg = np.meshgrid(x, t)

# set up 1-form
u = np.ones(np.shape(xg))
v = np.tanh(xg)*(np.cosh(xg))**(2/3)
f1 = fp.form_1(xg, tg, u, v)

# set up figure and axis
fig = plt.figure(figsize=(6, 12))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
ax1.set_ylabel(r'$t$')
ax2.set_ylabel(r'$t$')
ax2.set_xlabel(r'$x$')

# Hodge the result
f1_star = f1.num_hodge()

# plot each
f1.plot(ax1)
f1_star.plot(ax2)


