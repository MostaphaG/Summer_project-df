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
fig = plt.figure(figsize=(12, 6))
fig.tight_layout()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

ax1.set_aspect('equal')
ax1.set_xlabel(r'$x$', fontsize=16)
ax1.set_ylabel(r'$y$', fontsize=16)
ax1.set_title(r'$VF \ with \ Zoom \ Inset$', fontsize=16)
ax2.set_aspect('equal')
ax2.set_xlabel(r'$x$', fontsize=16)
ax2.set_title(r'$VF \ with \ derivative \ Inset$', fontsize=16)

# plot VF on all
VF.plot(ax1)
VF.plot(ax2)

# zoom on each and plot

# ax1
zoomed_ax1, zoomed_VF = VF.zoom(target=[2, 3], mag=1.25, dpd=7, inset=True, axis=ax1, insize=0.25)
zoomed_VF.colour('r')
zoomed_ax1.set_yticks(np.linspace(zoomed_VF.yg[0, 0], zoomed_VF.yg[-1, 0], 5))
zoomed_ax1.set_xticks(np.linspace(zoomed_VF.xg[0, 0], zoomed_VF.xg[0, -1], 5))
#zoomed_ax1.tick_params(axis='x', rotation=45)
plt.setp(zoomed_ax1.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
zoomed_VF.plot(zoomed_ax1)

# ax2
df_ax1, df_VF = VF.deriv(target=[2, 3], mag=1.25, dpd=7, inset=True, axis=ax2, insize=0.25)
df_VF.colour('r')
df_ax1.set_yticks(np.linspace(df_VF.yg[0, 0], df_VF.yg[-1, 0], 5))
df_ax1.set_xticks(np.linspace(df_VF.xg[0, 0], df_VF.xg[0, -1], 5))
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
fig = plt.figure(figsize=(6, 12))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
ax1.set_aspect('equal')
ax1.set_ylabel(r'$y$', fontsize=16)
ax1.set_title(r'$1-form \ with \ mag=1 \ inset$', fontsize=16)
ax1.tick_params(labelsize=14)
ax2.set_aspect('equal')
ax2.set_xlabel(r'$x$', fontsize=16)
ax2.set_ylabel(r'$y$', fontsize=16)
ax2.set_title(r'$1-form \ with \ mag=2 \ inset$', fontsize=16)
ax2.tick_params(labelsize=14)

plt.subplots_adjust(left=0.06, bottom=0.06, right=0.94, top=0.94, wspace=0.2)

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

# Exterior derivative 

# Set up figures and axes for plotting
fig1 = plt.figure(figsize=(8,8))
fig2 = plt.figure(figsize=(8,8))
fig3 = plt.figure(figsize=(8,8))
fig4 = plt.figure(figsize=(8,8))

ax1 = fig1.gca()
ax2 = fig2.gca()
ax3 = fig3.gca()
ax4 = fig4.gca()

# 1F setup
r = np.linspace(-5, 5, 21)
xg, yg = np.meshgrid(r, r)
u = yg*np.sin(xg)
v = -xg*np.cos(yg)

# Create object and provide component expressions
alpha = fp.form_1(xg, yg, u, v)
alpha.give_eqn('y*sin(x)','-x*cos(y)')

# Exterior derivative (analytical)
d_alpha1 = alpha.ext_d()

# Exterior derivative (numerical)
d_alpha2 = alpha.num_ext_d()

# Plot
alpha.plot(ax1)
d_alpha1.plot(ax2)
d_alpha2.plot(ax3)

# Continue - use contravariant with flat metric to get the vector field equivalent 
# and show the regions of +ve -ve curl agree with the exterior derivative
F = alpha.contravariant()

F.plot(ax4)
 
z1 = F.curl(target=(2,0), mag=10, dpd=7, inset=True, axis=ax4, insize=0.3)
z2 = F.curl(target=(-2,3), mag=10, dpd=7, inset=True, axis=ax4, insize=0.3)


# %%

# Alternative ext deriv example

 # Set up figures and axes for plotting
fig = plt.figure(figsize=(6, 12))
plt.subplots_adjust(left=0.06, bottom=0.06, right=0.94, top=0.94, wspace=0.2)
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

ax1.set_yticks(np.linspace(-1, 1, 5))
ax2.set_yticks(np.linspace(-1, 1, 5))

ax1.tick_params(labelsize=14)
ax1.set_aspect('equal')
ax1.set_ylabel(r'$y$', fontsize=16)
ax1.set_title(r'$0-form \ \Phi(x, y) \ Yukawa \ Potential$', fontsize=16)
ax2.tick_params(labelsize=14)
ax2.set_aspect('equal')
ax2.set_xlabel(r'$x$', fontsize=16)
ax2.set_ylabel(r'$y$', fontsize=16)
ax2.set_title(r'$1-form \ d \Phi$', fontsize=16)

# 0F setup
r = np.linspace(-1, 1, 21)
xg, yg = np.meshgrid(r, r)
phi = 1/np.sqrt(xg**2 + yg**2) * np.exp(-1*np.sqrt(xg**2 + yg**2))

# Create object and provide component expressions
form0 = fp.form_0(xg, yg, phi)
form0.give_eqn('1/sqrt(x**2 + y**2) * e**(-sqrt(x**2 + y**2))')

# Exterior derivative (numerical)
#d_f0_n = form0.num_ext_d()
#figure = plt.figure()
#axis = figure.gca()
#d_f0_a.plot(axis)

# Exterior derivative (analytical)
d_f0_a = form0.ext_d()
d_f0_a.log_scaling()
d_f0_a.sheet_size(0.04)

# Plot
form0.lines_number(40)
form0.density_increase(3)
form0.plot(ax1)
d_f0_a.plot(ax2)


# %%

# Trying Lorentz with treating d\rho as normal 1-form, not radial for cont. phi

rho = np.linspace(0, 1, 23)
z = np.linspace(0, 1, 23)
rhog, zg = np.meshgrid(rho, z)

form = 1/rhog
form2 = fp.form_2(rhog, zg, form)
form2.give_eqn('1/x')
form2.max_sheets(8)
form2.log_scaling()

# define the VF definining the velocity
u = np.zeros(np.shape(zg))
v = np.ones(np.shape(zg))
VF = fp.vector_field(rhog, zg, u, v)
VF.give_eqn('0', '1')
VF.log_scaling()

fig = plt.figure(figsize=(12, 12))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

for i in [1, 2, 3, 4]:
    exec('ax' + str(i) + '.set_aspect(\'equal\')')
    if i == 3 or i == 4:
        exec('ax' + str(i) + '.set_xlabel(r\'$ \\rho $\')')
    if i == 1 or i == 3:
        exec('ax' + str(i) + '.set_ylabel(r\'$z$\')')

ax1.set_title(r'$Magnetic \ 2-form \ from \ wire$')
ax2.set_title(r'$veclocity \ \vec{v}$')
ax3.set_title(r'$1-form \ Lorentz \ Force$')
ax4.set_title(r'$\vec{F}$')

form2.plot(ax1)

VF.plot(ax2)

# find numerical and analytical interior derivative and plot
num_int = form2.interior_d(VF)

# plot these
num_int.log_scaling()
num_int.plot(ax3)

# use cross product:
VF_c = fp.vector_field(rhog, zg, -1/rhog, np.zeros(np.shape(zg)))
VF.give_eqn('-1/x', '0')

VF_c.log_scaling()
VF_c.plot(ax4)

# %%

# Wedge example

# set up grids
r = np.linspace(-2, 2, 23)
xg, yg = np.meshgrid(r, r)

# set up figure and subplots
fig = plt.figure(figsize=(18, 6))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

#fig = plt.figure(figsize=(6, 12))
#plt.subplots_adjust(left=0.06, bottom=0.06, right=0.94, top=0.94, wspace=0.3)
#ax1 = fig.add_subplot(211)
#ax2 = fig.add_subplot(212)

ax1.set_aspect('equal')
ax1.set_title(r'$1-form \ \alpha$', fontsize=16)
ax1.set_xlabel(r'$x$', fontsize=16)
ax1.set_ylabel(r'$y$', fontsize=16)
ax1.tick_params(labelsize=14)

ax2.set_aspect('equal')
ax2.set_title(r'$Wedge \ \alpha \wedge \star \alpha$', fontsize=16)
ax2.set_xlabel(r'$x$', fontsize=16)
ax2.tick_params(labelsize=14)


# set up first 1-form
u1 = np.ones(np.shape(xg))
v1 = np.tanh(xg)*(np.cosh(xg))**(2/3)
form11 = fp.form_1(xg, yg, u1, v1)
form11.give_eqn('1', '-tanh(x)**2 * cosh(x)**(4/3)')

# and second via the Hodge
form12 = form11.hodge()

# wegde them:
form2 = form11.wedge(form_1_second=form12)

# plot
form11.plot(ax1)
form2.log_scaling()
form2.plot(ax2)


# %%

# Example of metric:
# polar or hyperbolic metric, VF, 1-form and the interior deriv of the 2

from matplotlib import cm

# set up grids
r = np.linspace(-2, 2, 21)
xg, yg = np.meshgrid(r, r)

# set up figure and subplots

fig = plt.figure(figsize=(12, 6))
# plt.subplots_adjust(left=0.06, bottom=0.06, right=0.94, top=0.94, wspace=0.3)
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

#fig = plt.figure(figsize=(6, 18))
#plt.subplots_adjust(left=0.06, bottom=0.06, right=0.94, top=0.94, wspace=0.3)
#ax1 = fig.add_subplot(311)
#ax2 = fig.add_subplot(312)
#ax3 = fig.add_subplot(313)

ax1.set_aspect('equal')
ax1.set_title(r'$Starting \ VF \ v^{i}$', fontsize=16)
ax1.set_xlabel(r'$x$', fontsize=16)
ax1.set_ylabel(r'$y$', fontsize=16)
ax1.set_yticks(np.linspace(-2, 2, 5))
ax1.set_xticks(np.linspace(-2, 2, 5))
ax1.tick_params(labelsize=14)

ax2.set_aspect('equal')
ax2.set_title(r'$v_{i} = g_{ij} v^{j}$', fontsize=16)
ax2.set_xlabel(r'$x$', fontsize=16)
ax2.set_yticks(np.linspace(-2, 2, 5))
ax2.set_xticks(np.linspace(-2, 2, 5))
ax2.tick_params(labelsize=14)

# set up the VF
u = np.ones(np.shape(xg))
v = np.ones(np.shape(xg))
VF = fp.vector_field(xg, yg, u, v)
VF.give_eqn('1', '1')

# plot it
VF.plot(ax1)

# set up the metric in strings
metric = [['1', '0'],
          ['0', '-tanh(x)**2 * cosh(x)**(4/3)']]

# find the 1-form
form1 = VF.covariant(g=metric)

# plot it
form1.plot(ax2)


# %%

# 2D BH example
p = np.linspace(-4, 4, 31)
q = np.linspace(0, 4, 19)
x, y = np.meshgrid(p, q)

u = np.ones(np.shape(x))
v = np.tanh(x)*(np.cosh(x))**(2/3)
f1 = fp.form_1(x, y, u, v)
f1.sheet_size(0.04)

fig = plt.figure(figsize=(6, 12))
ax1 = fig.add_subplot(211, adjustable='box')
ax2 = fig.add_subplot(212, adjustable='box')

ax1.tick_params(labelsize=14)
#ax1.set_aspect('equal')
ax1.set_ylabel(r'$t$', rotation=0, labelpad=10, fontsize=16)
ax1.set_title(r'$2D \ Black \ Hole \ frame \ field \ \omega$', fontsize=16)

ax2.tick_params(labelsize=14)
#ax2.set_aspect('equal')
ax2.set_ylabel(r'$t$', rotation=0, labelpad=10, fontsize=16)
ax2.set_xlabel(r'$x$', fontsize=16)
ax2.set_title(r'$\star \omega$', fontsize=16)

# Hodge the result
f1_star = f1.num_hodge(keep_object=False)
f1_star.sheet_size(0.04)

# plot
f1.plot(ax1)
f1_star.plot(ax2)

# %%


# Hodge Example previous option

# 1F setup
r = np.linspace(-5, 5, 21)
x, y = np.meshgrid(r, r)

# Black hole field

u = -2*x*((x**2+y**2)**(-1.5))*(1-(2/np.sqrt(x**2+y**2)))**(-2)
v = -2*y*((x**2+y**2)**(-1.5))*(1-(2/np.sqrt(x**2+y**2)))**(-2)

f1 = fp.form_1(x, y, u, v)
f1.log_scaling()
f1.max_sheets(10)

fig1 = plt.figure(figsize=(10,5))
ax1 = fig1.add_subplot(121)
ax2 = fig1.add_subplot(122)

f1.plot(ax1)

f1hodge = f1.num_hodge(keep_object=False)

f1hodge.plot(ax2)

# %%

# d^2 (any 0-form) = 0


# %%

# electric/magnetic fields examples


# %%

# Example showing proof of Stokes theorem


