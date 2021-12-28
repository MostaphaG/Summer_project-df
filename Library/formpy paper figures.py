import formpy as fp
import numpy as np
import matplotlib.pyplot as plt

# All use the examples:
# u = yg*np.sin(xg)
# v = -xg*np.cos(yg)
# scalar values = u*v

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
ax2.set_title(r'$VF \ with \ DF \ Inset$')
ax3.set_aspect('equal')
ax3.set_xlabel(r'$x$')
ax3.set_ylabel(r'$y$')
ax3.set_title(r'$VF \ with \ Div \ Inset$')
ax4.set_aspect('equal')
ax4.set_xlabel(r'$x$')
ax4.set_ylabel(r'$y$')
ax4.set_title(r'$VF \ with \ Curl \ Inset$')

# plot VF on all
VF.plot(ax1)
VF.plot(ax2)
VF.plot(ax3)
VF.plot(ax4)


# zoom on each and plot

# ax1
zoomed_ax1, zoomed_VF = VF.zoom(target=[2, 3], zoom=2, dpd=7, inset=True, axis=ax1)
zoomed_VF.colour('r')
zoomed_VF.plot(zoomed_ax1)

# ax2
df_ax1, df_VF = VF.DF(target=[2, 3], zoom=2, dpd=7, inset=True, axis=ax2)
df_VF.colour('r')
df_VF.plot(df_ax1)

# ax3
div_ax1, div_VF = VF.Div(target=[2, 3], zoom=2, dpd=7, inset=True, axis=ax3)
div_VF.colour('r')
div_VF.plot(div_ax1)

# ax4
curl_ax1, curl_VF = VF.Curl(target=[2, 3], zoom=2, dpd=7, inset=True, axis=ax4)
curl_VF.colour('r')
curl_VF.plot(curl_ax1)

# %%

# 1-form plot and zooming on it

r = np.linspace(-5, 5, 31)
xg, yg = np.meshgrid(r, r)

u = yg*np.sin(xg)
v = -xg*np.cos(yg)

form1 = fp.form_1(xg, yg, u, v)
form1.sheet_size(0.04)
form1.give_eqn('y*sin(x)','-x*cos(y)')

# set up siubplots for different zooming in values
fig = plt.figure(figsize=(13, 7))
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)
ax1.set_aspect('equal')
ax1.set_xlabel(r'$x$')
ax1.set_ylabel(r'$y$')
ax2.set_aspect('equal')
ax2.set_xlabel(r'$x$')
ax2.set_ylabel(r'$y$')
ax3.set_aspect('equal')
ax3.set_xlabel(r'$x$')
ax3.set_ylabel(r'$y$')

# plot 1-form on all
form1.plot(ax1)
form1.plot(ax2)
form1.plot(ax3)

# zoom on each and plot

# ax1
zoomed_ax, zoomed_form = form1.zoom(target=[2, 3.5], zoom=1, dpd=7, inset=True, axis=ax1)
zoomed_form.colour('r')
zoomed_form.plot(zoomed_ax)

# ax2
zoomed_ax, zoomed_form = form1.zoom(target=[2, 3.5], zoom=10, dpd=7, inset=True, axis=ax2)
zoomed_form.colour('r')
zoomed_form.plot(zoomed_ax)

# ax3
zoomed_ax, zoomed_form = form1.zoom(target=[2, 3.5], zoom=100, dpd=7, inset=True, axis=ax3)
zoomed_form.colour('r')
zoomed_form.plot(zoomed_ax)

# %%

# 2-form plot with zoom inset


# %%

# Example of exterior derivative



# %%

# Example of interior derivative
# Perhaps Lorentz force if it ends up working


# %%

# Hodge Example


# %%

# Wedge example


# %%

# Subplot example
# Could do:
# Combo of ext. alg. operations
# Comparison between numerical and analaytical operations for some method
    # (Not great figure though, showing two same plots)





# %%

# Example of metric:
# polar or hyperbolic metric, VF, 1-form and the interior deriv of the 2


# %%

'''
Showcasing, Perhaps most for an Appendix?
'''

# %%

# Example showing proof of Stokes theorem


# %%

# electric/magnetic fields examples


# %%

# 2D BH example (?)


# %%

# d^2 (any 0-form) = 0




