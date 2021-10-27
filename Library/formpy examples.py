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

plt.pause(2)

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
# form_1_obj = form_0_obj.ext_d()

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

# %%

# test of 0-form ext deriv with constant resulting 1-form


# set up grids
v = np.linspace(-4.5, 4.5, 11)
xg, yg = np.meshgrid(v, v)

# set up the 0 form object and plot it
form_0 = xg**2 + 3*yg
form_0_obj = fp.form_0(xg, yg, form_0)
form_0_obj.plot()

# supply equation and complete ext. deriv.
form_0_obj.give_eqn('x**2 + 3*y')
form_1_obj = form_0_obj.ext_d()  # this supplies the 1-form with equations too

# plot that 1-form object
form_1_obj.plot()

# %%

# Testing ext deriv of 1-form

# set up grids
v = np.linspace(-4.5, 4.5, 21)
xg, yg = np.meshgrid(v, v)

# set up the 1 form object and plot it
form_1_x = xg*np.cos(yg)
form_1_y = yg*np.sin(xg)

form_1_obj = fp.form_1(xg, yg, form_1_x, form_1_y)
form_1_obj.plot()

# supply equation and complete ext. deriv.
form_1_obj.give_eqn('x*cos(y)', 'y*sin(x)')

# compute the exterior derivative
form_2_obj = form_1_obj.ext_d()  # this supplies the 2-form with equations too

# plot that 2-form object
form_2_obj.plot()

plt.pause(2)

# change its density and replot
form_2_obj.same_range_density(26)
form_2_obj.max_sheets(10)
form_2_obj.plot(False)

print(form_2_obj.return_string())


# %%

# Testing the Hodge of a 1-form

# set up grids
v = np.linspace(-4.5, 4.5, 21)
xg, yg = np.meshgrid(v, v)

# set up the 1 form object and plot it
form_1_x = yg
form_1_y = -xg

form_1_obj = fp.form_1(xg, yg, form_1_x, form_1_y)
form_1_obj.plot()

# wait
plt.pause(2)

# supply equation and complete ext. deriv.
form_1_obj.give_eqn('y', '-x')

# compute the Hodge and replot
form_1_obj.Hodge(numerical_only=False, keep_object=True)

# replot and wait
form_1_obj.plot(keep=False)
plt.pause(2)

# compute the Hodge again, now set to a new form
form_1_obj_2H = form_1_obj.Hodge(numerical_only=False, keep_object=False)

form_1_obj_2H.plot(False)

# %%

# example of using subplots

# set up grids
v = np.linspace(-4.5, 4.5, 21)
xg, yg = np.meshgrid(v, v)

# set up the 1 form object and plot it
form_1_x = yg
form_1_y = -xg

# set up a figure
fig = plt.figure(figsize=(14, 7))

# create a form object using these
form_1_obj = fp.form_1(xg, yg, form_1_x, form_1_y, fig=fig, subplots=True)

form_1_obj.add_subplot(121)
form_1_obj.add_subplot(122)

form_1_obj.plot(keep=False, subplot_index=0)

form_1_obj.Hodge(numerical_only=True, keep_object=True)
form_1_obj.plot(keep=True, subplot_index=1)

# %%

# try it with 2-forms


# set up grids
v = np.linspace(-4.5, 4.5, 21)
xg, yg = np.meshgrid(v, v)

fig = plt.figure(figsize=(14, 7))

# set up the 1 form object and plot it
form_1_x = yg
form_1_y = -xg
form_1_obj = fp.form_1(xg, yg, form_1_x, form_1_y, fig=fig, subplots=True)

form_2 = xg*yg
form_2_obj = fp.form_2(xg, yg, form_2, fig=fig, subplots=True)

form_1_obj.add_subplot(121)
form_2_obj.add_subplot(121)
form_2_obj.add_subplot(122)

form_1_obj.plot(keep=False, subplot_index=0)
form_2_obj.plot(keep=True, subplot_index=1)


# %%


# Do the same but with exterior derivative to show that it can pass form
# function to function

# set up grids
v = np.linspace(-4.5, 4.5, 21)
xg, yg = np.meshgrid(v, v)

fig = plt.figure(figsize=(14, 7))

# set up the 1 form object and plot it
form_1_x = yg*np.cos(xg)
form_1_y = -xg
form_1_obj = fp.form_1(xg, yg, form_1_x, form_1_y, fig=fig, subplots=True)

# supply eqautions into the 1-form
form_1_obj.give_eqn('y*cos(x)', '-x')

form_2_obj = form_1_obj.ext_d(pass_on_figure=True)

form_1_obj.add_subplot(121)
form_2_obj.add_subplot(121)  # form_2 did not have the first one given, so indexes will be wrong if we don't supply it a dummy one
form_2_obj.add_subplot(122)

form_1_obj.plot(keep=False, subplot_index=0)
form_2_obj.plot(keep=True, subplot_index=1)

# %%

# testing proividing subplot axis straight in:


# set up grids
v = np.linspace(-4.5, 4.5, 11)
xg, yg = np.meshgrid(v, v)

fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

# set up the 0 form object and plot it
form_0 = xg**2 + 3*yg
form_0_obj = fp.form_0(xg, yg, form_0, fig=fig, subplots=True, sub_axis_list=[ax1, ax2])

form_0_obj.plot(keep=False, subplot_index=0)

# supply equation and change its density
form_0_obj.give_eqn('x**2 + 3*y')
form_0_obj.same_range_density(31)
form_0_obj.lines_number(20)

# plot that changed 0-form object
form_0_obj.plot(keep=True, subplot_index=1)

# %%

# example use of Wedge of two 1-forms, analytically

# set up grids
v = np.linspace(-4.5, 4.5, 26)
xg, yg = np.meshgrid(v, v)

# set up the 1 form object
form_1_x = yg*np.cos(xg)
form_1_y = -xg
form_1_obj = fp.form_1(xg, yg, form_1_x, form_1_y)

# change the stack sizes
form_1_obj.sheet_size(0.04)

# supply equations
form_1_obj.give_eqn('y*cos(x)', '-x')

# plot it
form_1_obj.plot()

# now find wedge of it with a different form
form_wedged_2 = form_1_obj.wedge_analytical('x*y', '2*y')

# plot it
form_wedged_2.plot()

# return string
print(form_wedged_2.return_string())

# %%

# Example of 2-form Hodge

v = np.linspace(-6, 6, 21)
xg, yg = np.meshgrid(v, v)

form_2 = xg*yg

fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

form_obj = fp.form_2(xg, yg, form_2, fig=fig, subplots=True, sub_axis_list=[ax1, ax2])

form_obj.plot(subplot_index=0)

# supply equations
form_obj.give_eqn('x*y')

# complete Hodge
form_0 = form_obj.Hodge(numerical_only=False, pass_on_figure=True)

form_0.plot(subplot_index=1)

#%%

r = np.linspace(-5, 5, 15)
xg, yg = np.meshgrid(r, r)

u = xg*cos(yg)
v = yg

form1 = fp.form_1(xg, yg, u, v)
form1.give_eqn('x*cos(y)', 'y')
form1.max_sheets(5)
form1.plot()
form1.return_string()

#%%
form1_zoom = form1.zoom((5,5), 10, 9)
form1_zoom.stacks()
form1_zoom.arrows()
form1_zoom.colour('red')
form1_zoom.plot()

