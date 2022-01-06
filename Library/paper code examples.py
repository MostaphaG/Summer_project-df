import formpy as fp
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

import formpy as fp
import numpy as np
import matplotlib.pyplot as plt

# Example code for exterior derivative of Yukawa 0-form

# Set up the 0-form
r = np.linspace(-1, 1, 27)
xg, yg = np.meshgrid(r, r)  # axis grids
rg = np.sqrt(xg**2 + yg**2)
scalar = 1/rg * np.exp(-rg)  # scalar funtion
f0 = fp.form_0(xg, yg, scalar)  # initialise instance
f0.give_eqn('1/sqrt(x**2 + y**2) * e**(-sqrt(x**2 + y**2))')  # supply equation

# get 1-form via exterior derivative and customise
dphi = f0.ext_d()
dphi.log_scaling()
dphi.sheet_size(0.04)

# customise 0-from to plot
f0.lines_number(50)  # number of level lines
f0.density_increase(3)

# set up figure and axis
fig = plt.figure(figsize=(10, 10))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax1.set_aspect('equal')
ax2.set_aspect('equal')

# plot
f0.plot(ax1)
dphi.plot(ax2)















