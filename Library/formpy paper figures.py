import formpy1 as fp
import numpy as np
import matplotlib.pyplot as plt


<<<<<<< HEAD
# 1-form plot and zooming on it


# %%

# 2-form plot with zoom inset


=======
# All use the examples:
# u = yg*np.sin(xg)
# v = -xg*np.cos(yg)
# scalar values = u*v

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

# ax1
zoomed_ax, zoomed_form = form1.zoom(target=[2, 3.5], zoom=10, dpd=7, inset=True, axis=ax2)
zoomed_form.colour('r')
zoomed_form.plot(zoomed_ax)

# ax1
zoomed_ax, zoomed_form = form1.zoom(target=[2, 3.5], zoom=100, dpd=7, inset=True, axis=ax3)
zoomed_form.colour('r')
zoomed_form.plot(zoomed_ax)

# %%

# 2-form plot with zoom inset


>>>>>>> Attempt-at-changing-all-figures-and-axis-supplies
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
Showcasing
'''

# %%

# Example showing proof of Stokes theorem
<<<<<<< HEAD
=======


# %%

# electric/magnetic fields examples


# %%

# 2D BH example (?)
>>>>>>> Attempt-at-changing-all-figures-and-axis-supplies


# %%

<<<<<<< HEAD
# electric/magnetic fields examples


# %%

# 2D BH example (?)


# %%

=======
>>>>>>> Attempt-at-changing-all-figures-and-axis-supplies
# d^2 (any 0-form) = 0




