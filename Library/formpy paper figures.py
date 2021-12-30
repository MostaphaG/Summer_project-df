import formpy as fp
import numpy as np
import matplotlib.pyplot as plt

# All use the examples:
# u = yg*np.sin(xg)
# v = -xg*np.cos(yg)
# scalar values = u*v

# %%

# Vector field with insets


# Set up the Vector Field
r = np.linspace(-5, 5, 21)
xg, yg = np.meshgrid(r, r)

u = yg*np.sin(xg)
v = -xg*np.cos(yg)

VF = fp.vector_field(xg, yg, u, v)
VF.give_eqn('y*sin(x)','-x*cos(y)')

# set up subplots for different zooming in values
fig = plt.figure(figsize=(10, 10))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

ax1.set_aspect('equal')
ax1.set_xlabel(r'$x$', fontsize=18)
ax1.set_ylabel(r'$y$', fontsize=18)
ax1.set_title(r'$VF \ with \ Zoom \ Inset$', fontsize=18)
ax2.set_aspect('equal')
ax2.set_xlabel(r'$x$', fontsize=18)
ax2.set_title(r'$VF \ with \ derivative \ Inset$', fontsize=18)

# plot VF on all
VF.plot(ax1)
VF.plot(ax2)


# zoom on each and plot

# ax1
zoomed_ax1, zoomed_VF = VF.zoom(target=[2, 3], zoom=2, dpd=7, inset=True, axis=ax1)
zoomed_VF.colour('r')
zoomed_ax1.set_yticks(np.linspace(zoomed_VF.yg[0, 0], zoomed_VF.yg[-1, 0], 5))
zoomed_ax1.set_xticks(np.linspace(zoomed_VF.xg[0, 0], zoomed_VF.xg[0, -1], 5))
zoomed_VF.plot(zoomed_ax1)

# ax2
df_ax1, df_VF = VF.deriv(target=[2, 3], zoom=2, dpd=7, inset=True, axis=ax2)
df_VF.colour('r')
df_ax1.set_yticks(np.linspace(df_VF.yg[0, 0], df_VF.yg[-1, 0], 5))
df_ax1.set_xticks(np.linspace(df_VF.xg[0, 0], df_VF.xg[0, -1], 5))
df_VF.plot(df_ax1)

# %%
'''

Code to present example in paper

'''


# Set up the Vector Field
r = np.linspace(-5, 5, 21)
xg, yg = np.meshgrid(r, r)
u = yg*np.sin(xg)  # components
v = -xg*np.cos(yg)
VF = fp.vector_field(xg, yg, u, v)  # initialise instance
VF.give_eqn('y*sin(x)','-x*cos(y)')  # supply equations

# set figure
fig = plt.figure(figsize=(10, 10))
ax = fig.gca()
ax.set_aspect('equal')
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')

# plot field
VF.plot(ax)

# zoom onto it
zoomed_ax, zoomed_VF = VF.zoom(target=[2, 3], zoom=1.25, dpd=7,
                               inset=True, axis=ax, insize=0.25)

# customise
zoomed_VF.colour('r')
zoomed_ax.set_yticks(np.linspace(zoomed_VF.yg[0, 0], zoomed_VF.yg[-1, 0], 5))
zoomed_ax.set_xticks(np.linspace(zoomed_VF.xg[0, 0], zoomed_VF.xg[0, -1], 5))

# replot inset
zoomed_VF.plot(zoomed_ax)



# %%

# 1-form plot and zooming on it

r = np.linspace(-5, 5, 23)
xg, yg = np.meshgrid(r, r)

u = yg*np.sin(xg)
v = -xg*np.cos(yg)

form1 = fp.form_1(xg, yg, u, v)
form1.sheet_size(0.04)
form1.give_eqn('y*sin(x)','-x*cos(y)')

# set up siubplots for different zooming in values
fig = plt.figure(figsize=(13, 7))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax1.set_aspect('equal')
ax1.set_xlabel(r'$x$', fontsize=18)
ax1.set_ylabel(r'$y$', fontsize=18)
ax1.set_title(r'$1-form \ with \ zoom=1 \ inset$', fontsize=18)
ax1.tick_params(labelsize=16)
ax2.set_aspect('equal')
ax2.set_xlabel(r'$x$', fontsize=18)
ax2.set_ylabel(r'$y$', fontsize=18)
ax2.set_title(r'$1-form \ with \ zoom=2 \ inset$', fontsize=18)
ax2.tick_params(labelsize=16)

# plot 1-form on all
form1.plot(ax1)
form1.plot(ax2)

# zoom on each and plot

# ax1
zoomed_ax, zoomed_form = form1.zoom(target=[2, 3], zoom=1, dpd=7, inset=True, axis=ax1)
zoomed_form.colour('r')
zoomed_form.plot(zoomed_ax)

# ax2
zoomed_ax, zoomed_form = form1.zoom(target=[2, 3], zoom=2, dpd=7, inset=True, axis=ax2)
zoomed_form.colour('r')
zoomed_form.plot(zoomed_ax)

# %%

# 2-form plot with zoom inset

v = np.linspace(-6, 6, 31)
xg, yg = np.meshgrid(v, v)
zooming = 1.2
u = yg*np.sin(xg)
v = -xg*np.cos(yg)
form = u*v
form2 = fp.form_2(xg, yg, form)
form2.give_eqn('-y*x*sin(x)*cos(y)')

# set up a figure and axis to plot it and its inset
fig = plt.figure()
ax = fig.gca()
ax.set_aspect('equal')
ax.tick_params(labelsize=14)
ax.set_xlabel(r'$x$', fontsize=16)
ax.set_ylabel(r'$y$', fontsize=16)
ax.set_title(r'$2-form \ with \ zoom={:.1f} \ inset$'.format(zooming), fontsize=16)

form2.plot(ax)

zoomed_ax, form2_zoomed = form2.zoom(target=[2.5, 3], zoom=zooming, dpd=9, inset=True, axis=ax)

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




