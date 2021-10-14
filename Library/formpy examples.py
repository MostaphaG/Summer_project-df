# formpy examples
import formpy as fp
import numpy as np
import matplotlib.pyplot as plt

# %%

# plotting a 1-form

# set up needed parameters
v = np.linspace(-6, 6, 31)
xg, yg = np.meshgrid(v, v)
F_x = yg*np.sin(xg)
F_y = xg*np.cos(yg)

# PLOT, note, it will create a figure for user
# we probably don't want that, otherwise we would have to make this
# a method to the matplotlib object, which might mean we need to play with
# their library, which I suppose we can't.
form_obj = fp.form_1(xg, yg, F_x, F_y)

form_obj.axis.set_xlabel('x')
form_obj.axis.set_ylabel('y')
form_obj.arrows()
form_obj.colour('blue')
form_obj.fig_size(8, 9)
form_obj.head_width(0.3)
form_obj.orient('tail')
form_obj.max_sheets(8)
form_obj.sheet_size(0.04)
form_obj.surround_space(6)

form_obj.plot()


# %%

# 1 -form density change example

# set up params.
v = np.linspace(-6, 6, 31)
xg, yg = np.meshgrid(v, v)
F_x = yg*np.sin(xg)
F_y = xg*np.cos(yg)

# create object
form_obj = fp.form_1(xg, yg, F_x, F_y)

# plot
form_obj.plot()

# pause
plt.pause(3)

# supply equation
form_obj.give_eqn('y*sin(x)', 'x*cos(y)')

# change density
form_obj.same_range_density(11)

# replot
form_obj.plot(False)

# %%

# plotting a 2-form

v = np.linspace(-6, 6, 21)
xg, yg = np.meshgrid(v, v)

form_2 = xg*yg


form_obj = fp.form_2(xg, yg, form_2)
form_obj.plot()

plt.pause(3)

# redo this with some 'customisations' (bad ones for an example)
form_obj.give_eqn('x*y')
form_obj.same_range_density(18)
# form_obj.sheet_size(0.466)
form_obj.plot(False)


# %%

# plotting a 0-form

v = np.linspace(-4.5, 4.5, 11)
xg, yg = np.meshgrid(v, v)


form_0 = np.cos(xg*yg)

form_obj = fp.form_0(xg, yg, form_0)
form_obj.lines_number(4)
# form_obj.density_increase(20)  # demonstation of an error
form_obj.plot()

plt.pause(5)

# customising grids with an equation:

form_obj.give_eqn('cos(x*y)')

# change the density:
form_obj.density_increase(25)
form_obj.lines_number(15)
form_obj.plot(keep=False)

# %%

# testing 0-form exterior derivative

# set up grids
v = np.linspace(-4.5, 4.5, 11)
xg, yg = np.meshgrid(v, v)

# set up the 0 form object and plot it
form_0 = np.cos(xg*yg)
form_0_obj = fp.form_0(xg, yg, form_0)
form_0_obj.plot()

# try exterior derivative without having first given an equation in
# throws an error
form_1_obj = form_0_obj.ext_d()

# supply equation and complete ext. deriv.
form_0_obj.give_eqn('cos(x*y)')
form_1_obj = form_0_obj.ext_d()  # this supplies the 1-form with equations too

# plot that 1-form object
form_1_obj.plot()

plt.pause(2)

# change its density and replot
form_1_obj.same_range_density(26)
form_1_obj.sheet_size(0.04)
form_1_obj.plot(False)

print(form_1_obj.return_string())

