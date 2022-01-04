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
