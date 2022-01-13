'''

Copies of code that produces some figures not used in paper / Wrong

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