'''

Copies of code that produces some figures not used in paper

'''

import formpy as fp
import numpy as np
import matplotlib.pyplot as plt

# %%

# Vector field with insets: (zoom, DF, div and curl), as subplots


# Set up the Vector Field
r = np.linspace(-5, 5, 31)
xg, yg = np.meshgrid(r, r)

u = yg*np.sin(xg)
v = -xg*np.cos(yg)

VF = fp.vector_field(xg, yg, u, v)
VF.give_eqn('y*sin(x)','-x*cos(y)')

# set up subplots for different zooming in values
fig = plt.figure(figsize=(10, 10))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

ax1.set_aspect('equal')
ax1.set_xlabel(r'$x$')
ax1.set_ylabel(r'$y$')
ax1.set_title(r'$VF \ with \ Zoom \ Inset$')
ax2.set_aspect('equal')
ax2.set_xlabel(r'$x$')
ax2.set_ylabel(r'$y$')
ax2.set_title(r'$VF \ with \ derivative \ Inset$')
ax3.set_aspect('equal')
ax3.set_xlabel(r'$x$')
ax3.set_ylabel(r'$y$')
ax3.set_title(r'$VF \ with \ div \ Inset$')
ax4.set_aspect('equal')
ax4.set_xlabel(r'$x$')
ax4.set_ylabel(r'$y$')
ax4.set_title(r'$VF \ with \ curl \ Inset$')

# plot VF on all
VF.plot(ax1)
VF.plot(ax2)
VF.plot(ax3)
VF.plot(ax4)


# zoom on each and plot

# ax1
zoomed_ax1, zoomed_VF = VF.zoom(target=[2, 3], mag=2, dpd=7, inset=True, axis=ax1)
zoomed_VF.colour('r')
zoomed_VF.plot(zoomed_ax1)

# ax2
df_ax1, df_VF = VF.deriv(target=[2, 3], mag=2, dpd=7, inset=True, axis=ax2)
df_VF.colour('r')
df_VF.plot(df_ax1)

# ax3
div_ax1, div_VF = VF.div(target=[2, 3], mag=2, dpd=7, inset=True, axis=ax3)
div_VF.colour('r')
div_VF.plot(div_ax1)

# ax4
curl_ax1, curl_VF = VF.curl(target=[2, 3], mag=2, dpd=7, inset=True, axis=ax4)
curl_VF.colour('r')
curl_VF.plot(curl_ax1)