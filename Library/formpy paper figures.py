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
fig = plt.figure(figsize=(6, 12))
fig.tight_layout()
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

ax1.set_aspect('equal')
ax1.set_ylabel(r'$y$', fontsize=16)
ax1.set_title(r'$VF \ with \ Zoom \ Inset$', fontsize=16)
ax2.set_aspect('equal')
ax2.set_xlabel(r'$x$', fontsize=16)
ax2.set_ylabel(r'$y$', fontsize=16)
ax2.set_title(r'$VF \ with \ derivative \ Inset$', fontsize=16)

plt.subplots_adjust(left=0.06, bottom=0.06, right=0.94, top=0.94, wspace=0.2)

# plot VF on all
VF.plot(ax1)
VF.plot(ax2)


# zoom on each and plot

# ax1
zoomed_ax1, zoomed_VF = VF.zoom(target=[2, 3], mag=2, dpd=7, inset=True, axis=ax1, insize=0.25)
zoomed_VF.colour('r')
zoomed_ax1.set_yticks(np.linspace(zoomed_VF.yg[0, 0], zoomed_VF.yg[-1, 0], 4))
zoomed_ax1.set_xticks(np.linspace(zoomed_VF.xg[0, 0], zoomed_VF.xg[0, -1], 4))
#zoomed_ax1.tick_params(axis='x', rotation=45)
plt.setp(zoomed_ax1.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
zoomed_VF.plot(zoomed_ax1)

# ax2
df_ax1, df_VF = VF.deriv(target=[2, 3], mag=2, dpd=7, inset=True, axis=ax2, insize=0.25)
df_VF.colour('r')
df_ax1.set_yticks(np.linspace(df_VF.yg[0, 0], df_VF.yg[-1, 0], 4))
df_ax1.set_xticks(np.linspace(df_VF.xg[0, 0], df_VF.xg[0, -1], 4))
#df_ax1.tick_params(axis='x', rotation=45)
plt.setp(df_ax1.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
df_VF.plot(df_ax1)

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
ax1.set_title(r'$1-form \ with \ mag=1 \ inset$', fontsize=18)
ax1.tick_params(labelsize=16)
ax2.set_aspect('equal')
ax2.set_xlabel(r'$x$', fontsize=18)
ax2.set_ylabel(r'$y$', fontsize=18)
ax2.set_title(r'$1-form \ with \ mag=2 \ inset$', fontsize=18)
ax2.tick_params(labelsize=16)

# plot 1-form on all
form1.plot(ax1)
form1.plot(ax2)

# zoom on each and plot

# ax1
zoomed_ax, zoomed_form = form1.zoom(target=[2, 3], mag=1, dpd=7, inset=True, axis=ax1)
zoomed_form.colour('r')
zoomed_form.plot(zoomed_ax)

# ax2
zoomed_ax, zoomed_form = form1.zoom(target=[2, 3], mag=2, dpd=7, inset=True, axis=ax2)
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
ax.set_title(r'$2-form \ with \ mag={:.1f} \ inset$'.format(zooming), fontsize=16)

form2.plot(ax)

zoomed_ax, form2_zoomed = form2.zoom(target=[2.5, 3], mag=zooming, dpd=9, inset=True, axis=ax)

# %%

# Example of exterior derivative



# %%

# Example of interior derivative
# Perhaps Lorentz force if it ends up working

# NB working with muI/2pi = 1

# set up the 2-form
v = np.linspace(-5, 5, 23)
xg, yg = np.meshgrid(v, v)

form = -xg/xg**2 + yg**2
form2 = fp.form_2(xg, yg, form)
form2.give_eqn('-x/(x**2 + y**2)')

# define the VF definining the velocity
u = np.zeros(np.shape(xg))
v = np.ones(np.shape(xg))
VF = fp.vector_field(xg, yg, u, v)
VF.give_eqn('0', '1')

# set up figure and axis
fig = plt.figure()
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)
for i in [1, 2, 3, 4]:
    exec('ax' + str(i) + '.set_aspect(\'equal\')')
    exec('ax' + str(i) + '.set_xlabel(r\'$x$\')')
    exec('ax' + str(i) + '.set_ylabel(r\'$y$\')')

# plot form and VF
form2.plot(ax1)
VF.plot(ax2)

# find numerical and analytical interior derivative and plot
num_int = form2.interior_d(VF, numerical_only=True)

# plot these
num_int.plot(ax3)

# use cross product:
ax4.quiver(xg, yg, -xg/(xg**2 + yg**2), 1/(xg**2 + yg**2))

# %%

# Trying Lorentz with treating d\rho as normal 1-form, not radial for cont. phi

rho = np.linspace(0, 5, 23)
z = np.linspace(0, 5, 23)
rhog, zg = np.meshgrid(rho, z)

form = 1/rhog
form2 = fp.form_2(rhog, zg, form)
form2.give_eqn('1/x')

# define the VF definining the velocity
u = np.zeros(np.shape(zg))
v = np.ones(np.shape(zg))
VF = fp.vector_field(rhog, zg, u, v)
VF.give_eqn('0', '1')

fig = plt.figure(figsize=(12, 12))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

for i in [1, 2, 3, 4]:
    exec('ax' + str(i) + '.set_aspect(\'equal\')')
    if i == 3 or i == 4:
        exec('ax' + str(i) + '.set_xlabel(r\'$x$\')')
    if i == 1 or i == 3:
        exec('ax' + str(i) + '.set_ylabel(r\'$y$\')')

ax1.set_title(r'$Magnetic \ 2-form \ from \ wire$')
ax2.set_title(r'$veclocity \ vector \ field$')
ax3.set_title(r'$1-form \ Lorentz \ Force$')
ax4.set_title(r'$Vector \ field \ force$')

form2.plot(ax1)

VF.plot(ax2)

# find numerical and analytical interior derivative and plot
num_int = form2.interior_d(VF, numerical_only=True)

# plot these
num_int.plot(ax3)

# use cross product:
VF_c = fp.vector_field(rhog, zg, -1/rhog, np.zeros(np.shape(zg)))
VF.give_eqn('-1/x', '0')

VF_c.log_scaling()
VF_c.plot(ax4)

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




