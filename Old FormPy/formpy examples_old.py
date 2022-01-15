# formpy examples
import formpy as fp
import numpy as np
import matplotlib.pyplot as plt
import timeit

# %%

tstart = timeit.default_timer()

# plotting a 1-form

# set up needed parameters
v = np.linspace(-6, 6, 41)
xg, yg = np.meshgrid(v, v)

F_x = 10*xg*yg
F_y = 1/(np.sin(yg))

# PLOT, note, it will create a figure for user
# we probably don't want that, otherwise we would have to make this
# a method to the matplotlib object, which might mean we need to play with
# their library, which I suppose we can't.
form_obj = fp.form_1(xg, yg, F_x, F_y)

form_obj.axis.set_xlabel('x')
form_obj.axis.set_ylabel('y')
form_obj.colour('blue')
form_obj.fig_size(8, 9)
form_obj.head_width(0.3)
form_obj.orient('tail')
form_obj.max_sheets(8)
form_obj.sheet_size(0.04)
form_obj.surround_space(6)

form_obj.plot()

tstop = timeit.default_timer()
print(tstop-tstart)

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

# trying subplots with 2-forms


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
form_1_x = xg**2
form_1_y = yg**2
form_1_obj = fp.form_1(xg, yg, form_1_x, form_1_y)

# change the stack sizes
form_1_obj.sheet_size(0.04)

# supply equations
form_1_obj.give_eqn('x**2', 'y**2')

# plot it
form_1_obj.plot()

# now find wedge of it with a different form
form_wedged_2 = form_1_obj.wedge_analytical('1/y**2', 'x**2')

# plot it
form_wedged_2.plot()

# return string
print(form_wedged_2.return_string())


# %%

# Example of Wedge done numerically

# set up grids
v = np.linspace(-4.5, 4.5, 26)
xg, yg = np.meshgrid(v, v)

# set up a figure
fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

# set up the 1 form object
form_1_x = 1/xg**3
form_1_y = yg**2
form_1_obj = fp.form_1(xg, yg, form_1_x, form_1_y, fig=fig, subplots=True, sub_axis_list=[ax1, ax2])

# plot it
form_1_obj.plot(keep=True, subplot_index=0)

# now find wedge of it with a different form
form_wedged_2 = form_1_obj.wedge_num(form_1_second=(1/yg**2, xg**2), pass_on_figure=True)

# plot it
form_wedged_2.plot(keep=True, subplot_index=1)

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

u = xg
v = yg

form1 = fp.form_1(xg, yg, u, v, 'x', 'y')
# form1.give_eqn('x*cos(y)', 'y')
form1.max_sheets(5)
form1.plot()
form1.return_string()

form1_zoom = form1.zoom((10,10), 10000, 9)
form1_zoom.colour('red')
form1_zoom.plot()

# %%

# Testing a vector field object, customisations and plotting


# set up needed parameters
v = np.linspace(-6, 6, 31)
xg, yg = np.meshgrid(v, v)

# set up the field
F_x = xg/(xg**2 + yg**2)**1.5
F_y = yg/(xg**2 + yg**2)**1.5

# set up a figure, with subplots
fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

# set up the object
field_obj = fp.vector_field(xg, yg, F_x, F_y, fig=fig, subplots=True, sub_axis_list=[ax1, ax2])

# on axis 1, plot default
field_obj.plot(keep=True, subplot_index=0)

# change some properties and plot the second subplot
field_obj.axis[1].set_xlabel('x')
field_obj.axis[1].set_ylabel('y')
field_obj.colour('blue')
field_obj.orient('tail')
field_obj.surround_space(6)

field_obj.plot(keep=True, subplot_index=1)

# %%

# Testing numerical ext deriv of 0-forms, using numpy gradient

# set up grids
v = np.linspace(-4.5, 4.5, 21)
xg, yg = np.meshgrid(v, v)

# set up a figure with subplots
fig = plt.figure()
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)
ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax3.set_aspect('equal')

# set up the 0 form object and plot it
form_0 = (xg+yg)/(xg**2 +yg**2)**1.5

# create an object with these axis in it
form_0_obj = fp.form_0(xg, yg, form_0, fig=fig, subplots=True, sub_axis_list=[ax1, ax2, ax3])

# plot it on first subplot
form_0_obj.plot(keep=True, subplot_index=0)

# compute the numerical ext deriv and plot it on second subplot
form_1_num_e = form_0_obj.num_ext_d(edge_order=1, pass_on_figure=True)  # pass figure to pass on subplot axis

# plot it on second axis set
form_1_num_e.plot(keep=True, subplot_index=1)

# supply equation and complete ext. deriv., then plot that on 3rd axis set
form_0_obj.give_eqn('(x+y)/(x**2 + y**2)**1.5')
form_1_obj_a = form_0_obj.ext_d(pass_on_figure=True)  # this supplies the 1-form with equations too

# plot that 1-form object
form_1_obj_a.plot(keep=True, subplot_index=2)

# plot the difference between trhe two (analytical and numerical)
#form_test = fp.form_1(xg, yg, form_1_num_e.F_x - form_1_obj_a.F_x, form_1_num_e.F_y - form_1_obj_a.F_y)
#form_test.plot()

# %%

# Testing Hodge of a 0-form

# set up grids
v = np.linspace(-4.5, 4.5, 11)
xg, yg = np.meshgrid(v, v)

# set up a figure with subplots
fig = plt.figure()
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

# set up the 0 form object and plot it
form_0 = xg**2 + 3*yg
form_0_obj = fp.form_0(xg, yg, form_0, fig=fig, subplots=True, sub_axis_list=[ax1, ax2, ax3, ax4])
form_0_obj.plot(keep=True, subplot_index=0)
form_0_obj.plot(keep=True, subplot_index=2)

# now compute the Hodge numerically only and plot the resulting 2-form
form_2_num = form_0_obj.Hodge(numerical_only=True, pass_on_figure=True)
form_2_num.plot(keep=True, subplot_index=1)

# supply equations and do it analytically
form_0_obj.give_eqn('x**2 + 3*y')
form_2_an = form_0_obj.Hodge(numerical_only=False, pass_on_figure=True)

form_2_an.plot(keep=True, subplot_index=3)

#%%

# Vector field set up and zooming example


r = np.linspace(-4.5, 4.5, 21)
xg, yg = np.meshgrid(r, r)

u = xg*np.cos(yg)
v = -yg*np.sin(xg)

vf1 = fp.vector_field(xg, yg, u, v)
vf1.give_eqn('x*cos(y)', '-y*sin(x)')
vf1.plot()

vfz = vf1.zoom((2, 2), zoom=3, dpd=9)
vfz.plot()

# %%

# Testing numerical 1-form exterior derivative


# set up grids
v = np.linspace(-4.5, 4.5, 19)
xg, yg = np.meshgrid(v, v)

# set up the 1 form object and plot it
form_1_x = -xg/(xg**2 + yg**2)**1.5
form_1_y = -yg/(xg**2 + yg**2)**1.5

# set up a figure with sublots
fig = plt.figure()
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)

ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax3.set_aspect('equal')

# set up a 1-form object with these:
form_1_obj = fp.form_1(xg, yg, form_1_x, form_1_y, fig=fig, subplots=True, sub_axis_list=[ax1, ax2, ax3])

# plot it
form_1_obj.plot(keep=True, subplot_index=0)

# complete the numerical exterior derivative
form_2_num = form_1_obj.num_ext_d(pass_on_figure=True)

# plot it
form_2_num.plot(keep=True, subplot_index=1)

# complete the analytical exterior derivative
form_1_obj.give_eqn('-x/(x**2 + y**2)**1.5', '-y/(x**2 + y**2)**1.5')
form_2_an = form_1_obj.ext_d(pass_on_figure=True)

# plot it
form_2_an.plot(keep=True, subplot_index=2)


# %%

# Testing zooming on 2-forms


v = np.linspace(-6, 6, 21)
xg, yg = np.meshgrid(v, v)

form_2 = xg*yg

fig = plt.figure()

fig.tight_layout()

form_obj = fp.form_2(xg, yg, form_2, fig=fig)

form_obj.axis.set_aspect('equal')

form_obj.plot()

# for zooming in, supply the equation
form_obj.give_eqn('x*y')

form_zoomed = form_obj.zooming(target=[2, 2], zoom=5, dpd=21)

# plot it
form_zoomed.plot()

# %%

# TESTING CONSTANT COMPONENTs

# 1-form

# set up needed parameters
v = np.linspace(-6, 6, 21)
xg, yg = np.meshgrid(v, v)

# set up form components
F_x = 1*np.ones(np.shape(xg))  # upto user to make sure their arrays are of correct size
F_y = 3*np.ones(np.shape(xg))

form_obj1 = fp.form_1(xg, yg, F_x, F_y)

# beofre plotting, give equations, to see if if constant forms are dealt with
# correctly there at least
form_obj1.give_eqn('1', '3')
form_obj1.same_range_density(16)

form_obj1.figure.tight_layout()
form_obj1.axis.set_aspect('equal')

form_obj1.plot()


# %%

# Testing analytical interior derivative of 1-form

# set up needed parameters
r = np.linspace(-6, 6, 21)
xg, yg = np.meshgrid(r, r)

# set up form components and the object
F_x = xg**2
F_y = yg*xg
form_obj1 = fp.form_1(xg, yg, F_x, F_y)

# give equations to it
form_obj1.give_eqn('x**2', 'y*x')

# plot it
form_obj1.plot()

# set up a vector field object:
u = np.cos(yg)
v = np.sin(xg)

vf = fp.vector_field(xg, yg, u, v)
vf.give_eqn('cos(y)', 'sin(x)')

# plot it
vf.plot()

# complete interior derivative, with object, with equations and without any input:
form_0_i_1 = form_obj1.interior_d(vector_field=vf, pass_on_figure=False, numerical_only=False)
form_0_i_2 = form_obj1.interior_d(pass_on_figure=False, numerical_only=False)
form_0_i_3 = form_obj1.interior_d(('cos(y)', 'sin(x)'), pass_on_figure=False, numerical_only=False)

# plot each
form_0_i_1.plot()
form_0_i_2.plot()
form_0_i_3.plot()

# %%

# Testing numerical interior derivative of a 1-form:

# set up needed parameters
r = np.linspace(-6, 6, 21)
xg, yg = np.meshgrid(r, r)

# set up form components and the object
F_x = xg**2
F_y = yg*xg
form_obj1 = fp.form_1(xg, yg, F_x, F_y)

# plot it
form_obj1.plot()

# set up a vector field object:
u = np.cos(yg)
v = np.sin(xg)

vf = fp.vector_field(xg, yg, u, v)

# plot it:
vf.plot()

# complete int deriv. with object, arrays and with no input:
form_0_i_1 = form_obj1.interior_d(vector_field=vf, pass_on_figure=False, numerical_only=True)
form_0_i_2 = form_obj1.interior_d(pass_on_figure=False, numerical_only=True)
form_0_i_3 = form_obj1.interior_d((u, v), pass_on_figure=False, numerical_only=True)

# plot each
form_0_i_1.plot()
form_0_i_2.plot()
form_0_i_3.plot()


# %%

# Testing numerical int deriv with constant fields

# set up needed parameters
r = np.linspace(-6, 6, 21)
xg, yg = np.meshgrid(r, r)

# set up form components and the object
F_x = 3* np.ones(np.shape(xg))
F_y = 1 * np.ones(np.shape(xg))
form_obj1 = fp.form_1(xg, yg, F_x, F_y)

# plot it
form_obj1.plot()

# set up a vector field object:
u = 1 * np.ones(np.shape(xg))
v = np.sin(xg)

vf = fp.vector_field(xg, yg, u, v)

# plot it:
vf.plot()

# complete int deriv. with object, arrays and with no input:
form_0_i_1 = form_obj1.interior_d(vector_field=vf, pass_on_figure=False, numerical_only=True)
form_0_i_2 = form_obj1.interior_d(pass_on_figure=False, numerical_only=True)
form_0_i_3 = form_obj1.interior_d((u, v), pass_on_figure=False, numerical_only=True)

# plot each
form_0_i_1.plot()
form_0_i_2.plot()
form_0_i_3.plot()

# NB no contours on plot 4 (form_0_i_2) as 0-form = 4. constant, so no contorus

# %%

# Testing numerical interior derivative of a 2-form

# set up needed parameters
r = np.linspace(-6, 6, 21)
xg, yg = np.meshgrid(r, r)

# set up form
form2_comp = xg*yg**2
form_obj2 = fp.form_2(xg, yg, form2_comp)

# plot it
form_obj2.plot()

# set up a vector field object:
u = np.cos(yg)
v = np.sin(xg)

vf = fp.vector_field(xg, yg, u, v)

# plot it:
vf.plot()

# complete int deriv. with object, arrays and with no input:
form_1_i_1 = form_obj2.interior_d(vector_field=vf, pass_on_figure=False, numerical_only=True)
form_1_i_2 = form_obj2.interior_d(pass_on_figure=False, numerical_only=True)
form_1_i_3 = form_obj2.interior_d((u, v), pass_on_figure=False, numerical_only=True)

# plot each
form_1_i_1.plot()
form_1_i_2.plot()
form_1_i_3.plot()


# %%


# Testing analytical interior derivative of a 2-form

# set up needed parameters
r = np.linspace(-6, 6, 21)
xg, yg = np.meshgrid(r, r)

# set up form
form2_comp = xg*yg**2
form_obj2 = fp.form_2(xg, yg, form2_comp)

# give it the equation for analytical calculations
form_obj2.give_eqn('x*y**2')

# plot it
form_obj2.plot()

# set up a vector field object:
u = np.cos(yg)
v = np.sin(xg)

vf = fp.vector_field(xg, yg, u, v)

# give it the equation
vf.give_eqn('cos(y)', 'sin(x)')

# plot it:
vf.plot()

# complete int deriv. with object, arrays and with no input:
form_1_i_1 = form_obj2.interior_d(vector_field=vf, pass_on_figure=False, numerical_only=False)
form_1_i_2 = form_obj2.interior_d(pass_on_figure=False, numerical_only=False)
form_1_i_3 = form_obj2.interior_d(('cos(y)', 'sin(x)'), pass_on_figure=False, numerical_only=False)

# plot each
form_1_i_1.plot()
form_1_i_2.plot()
form_1_i_3.plot()

# %%

# Testing interior derivative of a 2-form with constant 2-form and/or vf

# set up needed parameters
r = np.linspace(-6, 6, 21)
xg, yg = np.meshgrid(r, r)

# set up form
form2_comp = xg
form_obj2 = fp.form_2(xg, yg, form2_comp)

# give it the equation for analytical calculations
form_obj2.give_eqn('x')

# plot it
form_obj2.plot()

# set up a vector field object:
u = np.cos(yg)
v = np.sin(xg)

vf = fp.vector_field(xg, yg, u, v)

# give it the equation
vf.give_eqn('1', '3')

# plot it:
vf.plot()

# complete int deriv. with object, arrays and with no input:
form_1_i_1 = form_obj2.interior_d(vector_field=vf, pass_on_figure=False, numerical_only=False)
form_1_i_2 = form_obj2.interior_d(pass_on_figure=False, numerical_only=False)
form_1_i_3 = form_obj2.interior_d(('10', 'pi'), pass_on_figure=False, numerical_only=False)
form_1_i_4 = form_obj2.interior_d(vector_field=vf, pass_on_figure=False, numerical_only=True)

# plot each
form_1_i_1.plot()
form_1_i_2.plot()
form_1_i_3.plot()
form_1_i_4.plot()


# %%

# Testing interior derivative of a 2-form with singularities

# set up needed parameters
r = np.linspace(-6, 6, 21)
xg, yg = np.meshgrid(r, r)

# set up form
form2_comp = 1/np.sqrt(xg**2 + yg**2)
form_obj2 = fp.form_2(xg, yg, form2_comp)

# give it the equation for analytical calculations
form_obj2.give_eqn('1/sqrt(x**2 + y**2)')

# plot it
form_obj2.plot()

# set up a vector field object:
u = np.cos(yg)
v = np.sin(xg)

vf = fp.vector_field(xg, yg, u, v)

# give it the equation
vf.give_eqn('1', 'x')  # Note, this overrites the u and v in the vf object too

# plot it:
vf.plot()

# complete int deriv. with object, arrays and with no input:
form_1_i_1 = form_obj2.interior_d(vector_field=vf, pass_on_figure=False, numerical_only=False)
form_1_i_2 = form_obj2.interior_d(pass_on_figure=False, numerical_only=False)
form_1_i_3 = form_obj2.interior_d(('10', 'pi'), pass_on_figure=False, numerical_only=False)
form_1_i_4 = form_obj2.interior_d(vector_field=vf, pass_on_figure=False, numerical_only=True)

# plot each
form_1_i_1.plot()
form_1_i_2.plot()
form_1_i_3.plot()
form_1_i_4.plot()

# %%

# Testing conversion form vector field to 1-form with metric as 1, 1 for now:

r = np.linspace(-4.5, 4.5, 21)
xg, yg = np.meshgrid(r, r)

u = xg*np.cos(yg)
v = -yg*np.sin(xg)

vf1 = fp.vector_field(xg, yg, u, v)
vf1.plot()

# change it to a 1-form
form_1_obj = vf1.formalise()

# plot that
form_1_obj.plot()

# %%

# Testing conversion form vector field to 1-form with string based metric

r = np.linspace(-4.5, 4.5, 21)
xg, yg = np.meshgrid(r, r)

# set up the field and vector field objecct
u = xg
v = -yg
vf1 = fp.vector_field(xg, yg, u, v)
vf1.plot()

# give equations
vf1.give_eqn('x', '-y')

# set up a metric
metric = [['1', '0'],
          ['0', 'y**2']]

# via this, set up a 1-form
form_1_obj = vf1.formalise(g=metric)

# plot it
form_1_obj.plot()

# create a comparison 1-form with correct components already given
form_1_correct = fp.form_1(xg, yg, u, v*yg**2)

# plot it
form_1_correct.plot()

# %%

# Testing conversion form vector field to 1-form with array based metric

r = np.linspace(-4.5, 4.5, 21)
xg, yg = np.meshgrid(r, r)

# set up the field and vector field objecct
u = xg
v = -yg
vf1 = fp.vector_field(xg, yg, u, v)
vf1.plot()

# set up a metric
metric = [[np.ones(np.shape(xg)), np.zeros(np.shape(xg))],
          [np.zeros(np.shape(xg)), yg**2]]

# via this, set up a 1-form
form_1_obj = vf1.formalise(g=metric)

# plot it
form_1_obj.plot()

# create a comparison 1-form with correct components already given
form_1_correct = fp.form_1(xg, yg, u, v*yg**2)

# plot it
form_1_correct.plot()

#%% Testing derivative of vector field

#Set up grids
r = np.linspace(-5, 5, 21)
xg, yg = np.meshgrid(r, r)

u = yg
v = -xg

# Set up subplots
fig1 = plt.figure(figsize=(10,5))
ax1 = fig1.add_subplot(121)
ax2 = fig1.add_subplot(122)


# Create vector field
vf1 = fp.vector_field(xg, yg, u, v, fig = fig1, subplots=True, sub_axis_list=[ax1, ax2])
vf1.give_eqn('y','-x')
vf1.plot(keep=False, subplot_index=0)

# Problem: As the method creates a new vector field
D_vf1 = vf1.DF((3,3), 5, 9, pass_on_figure=True)
D_vf1.plot(keep=False, subplot_index=1)

ax1.set_xlabel('Vector Field')
ax2.set_xlabel('Derivative Vector Field')

#%% Test zooming of vector field with subplot

#Set up grids
r = np.linspace(-5, 5, 21)
xg, yg = np.meshgrid(r, r)

u = yg*np.cos(xg)
v = -xg

# Set up subplots
fig1 = plt.figure(figsize=(10,5))
ax1 = fig1.add_subplot(121)
ax2 = fig1.add_subplot(122)

# Create vector field
vf1 = fp.vector_field(xg, yg, u, v, fig = fig1, subplots=True, sub_axis_list=[ax1, ax2])
vf1.give_eqn('y*cos(x)','-x')
vf1.plot(keep=False, subplot_index=0)

zoom_field = vf1.zoom((2,-2), 20, 13, pass_on_figure=True)
zoom_field.plot(keep=False, subplot_index=1)

ax1.set_xlabel('Vector Field')
ax2.set_xlabel('Zoomed Vector Field')

#%% Test divergence of vector field with subplot

#Set up grids
r = np.linspace(-5, 5, 21)
xg, yg = np.meshgrid(r, r)

u = yg*np.cos(xg)
v = -xg

# Set up subplots
fig1 = plt.figure(figsize=(10,5))
ax1 = fig1.add_subplot(121)
ax2 = fig1.add_subplot(122)

# Create vector field
vf1 = fp.vector_field(xg, yg, u, v, fig = fig1, subplots=True, sub_axis_list=[ax1, ax2])
vf1.give_eqn('y*cos(x)','-x')
vf1.plot(keep=False, subplot_index=0)

div_field = vf1.Div((5,0), 100, 13, pass_on_figure=True)
div_field.autoscale()
div_field.plot(keep=False, subplot_index=1)

ax1.set_xlabel('Vector Field')
ax2.set_xlabel('Div Vector Field')

#%% Test curl of vector field with subplot

#Set up grids
r = np.linspace(-5, 5, 21)
xg, yg = np.meshgrid(r, r)

u = yg*np.cos(xg)
v = -xg

# Set up subplots
fig1 = plt.figure(figsize=(10,5))
ax1 = fig1.add_subplot(121)
ax2 = fig1.add_subplot(122)

# Create vector field
vf1 = fp.vector_field(xg, yg, u, v, fig = fig1, subplots=True, sub_axis_list=[ax1, ax2])
vf1.give_eqn('y*cos(x)','-x')
vf1.plot(keep=False, subplot_index=0)

curl_field = vf1.Curl((0,0), 100, 13, pass_on_figure=True)
# curl_field.autoscale()
curl_field.plot(keep=False, subplot_index=1)

ax1.set_xlabel('Vector Field')
ax2.set_xlabel('Curl Vector Field')

# %%

# Testing the metric with the inetrior derivative (numerically)

r = np.linspace(-4.5, 4.5, 21)
xg, yg = np.meshgrid(r, r)

# set up the field and vector field objecct
u = np.ones(np.shape(xg))
v = np.ones(np.shape(xg))
vf1 = fp.vector_field(xg, yg, u, v)
vf1.plot()

# set up a metric
metric = [[np.ones(np.shape(xg)), np.zeros(np.shape(xg))],
          [np.zeros(np.shape(xg)), xg**2 + yg**2]]  # polar transformation

# via this, set up a 1-form
form_1_obj = vf1.formalise(g=metric)

# plot it
form_1_obj.plot()

# get the interioir deriv of that
zero_form = form_1_obj.interior_d(vector_field=vf1, pass_on_figure=False, numerical_only=True)

# plot it
zero_form.plot()


# %%

# Testing 1-form to VF via the metric


# ####################################################
# with string based metric:
# ####################################################


r = np.linspace(-4.5, 4.5, 21)
xg, yg = np.meshgrid(r, r)

# set up the 1-form
u = xg
v = -yg
form1 = fp.form_1(xg, yg, u, v)
form1.plot()

# give equations
form1.give_eqn('x', '-y')

# set up a metric
metric = [['1', '0'],
          ['0', 'y**2']]

# via this, get the VF
vf = form1.vectorise(g=metric)

# plot it
vf.plot()


# ####################################################
# with array based metric:
# ####################################################


# set up a metric
metric = [[np.ones(np.shape(xg)), np.zeros(np.shape(xg))],
          [np.zeros(np.shape(xg)), yg**2]]

# via this, get the VF
vf2 = form1.vectorise(g=metric)

# plot it
vf2.plot()


# ##################################################
# Compare to know result
# ##################################################


# create a comparison 1-form with correct components already given
vf_correct = fp.vector_field(xg, yg, u, v*yg**2)

# plot it
vf_correct.plot()


# %%

# Testing the interior derivative of 1-form wrt to VF from it via the metric

r = np.linspace(-4.5, 4.5, 21)
xg, yg = np.meshgrid(r, r)

# set up the 1-form
u = np.ones(np.shape(xg))
v = np.ones(np.shape(xg))
form1 = fp.form_1(xg, yg, u, v)
form1.plot()

# give equations
form1.give_eqn('1', '1')

# set up the metric
metric = [['1', '0'],
          ['0', 'x**2 + y**2']]

# find the VF
vf = form1.vectorise(g=metric)

# plot it
vf.plot()

# find the interior derivative of 1-form wrt that VF
zero_form = form1.interior_d(vector_field=vf, pass_on_figure=False, numerical_only=False)

# plot it:
zero_form.plot()

# %%

# Test EM potential 0-form ext deriv twice to 2-form (that should be 0)

# set up grids
v = np.linspace(-4.5, 4.5, 21)
xg, yg = np.meshgrid(v, v)

# set up a figure with subplots
fig = plt.figure()
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)
ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax3.set_aspect('equal')

# set up the 0 form object and plot it
form_0 = xg + yg

# create an object with these axis in it
form_0_obj = fp.form_0(xg, yg, form_0, fig=fig, subplots=True, sub_axis_list=[ax1, ax2, ax3])

# compute the numerical ext deriv
form_1_num = form_0_obj.num_ext_d(edge_order=1, pass_on_figure=True)  # pass figure to pass on subplot axis

# supply equation and complete ext. deriv. analytically
form_0_obj.give_eqn('x + y')
form_1_ana = form_0_obj.ext_d(pass_on_figure=True)  # this supplies the 1-form with equations too

# find the numerical exterior derivaitve of numerical 1-form and plot
form_2_num = form_1_num.num_ext_d(pass_on_figure=True)
form_2_num.plot(keep=True, subplot_index=0)

# find numerical exterior derivative of analytical 1-form and plot
form_2_numana = form_1_ana.num_ext_d(pass_on_figure=True)
form_2_numana.plot(keep=True, subplot_index=0)

# find analytical ext deriv. of analytical 1-form and plot
form_2_ana = form_1_ana.ext_d(pass_on_figure=True)
form_2_ana.plot(keep=True, subplot_index=0)


# %%

# Test 0-form for EM potential

# set up grids
v = np.linspace(-4.5, 4.5, 21)
xg, yg = np.meshgrid(v, v)

# set up a figure with subplots
fig = plt.figure()
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)
ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax3.set_aspect('equal')

# set up the 0 form object and plot it
form_0 = 1/(xg**2 + yg**2)

# create an object with these axis in it
form_0_obj = fp.form_0(xg, yg, form_0, fig=fig, subplots=True, sub_axis_list=[ax1, ax2, ax3])

# plot it on first subplot
form_0_obj.plot(keep=True, subplot_index=0)

# compute the numerical ext deriv and plot it on second subplot
form_1_num_e = form_0_obj.num_ext_d(edge_order=1, pass_on_figure=True)  # pass figure to pass on subplot axis

# plot it on second axis set
form_1_num_e.plot(keep=True, subplot_index=1)

# supply equation and complete ext. deriv., then plot that on 3rd axis set
form_0_obj.give_eqn('1/(x**2 + y**2)')
form_1_obj_a = form_0_obj.ext_d(pass_on_figure=True)  # this supplies the 1-form with equations too

# plot that 1-form object
form_1_obj_a.plot(keep=True, subplot_index=2)

# %%

# Further tests with numerical exterior derivative of 1-form

np.set_printoptions(True)

# set up grids
v = np.linspace(-0.1, 0.1, 21)
xg, yg = np.meshgrid(v, v)

# set up the 1 form object and plot it
#form_1_x = -xg/((xg**2 + yg**2)**1.5)
#form_1_y = -yg/((xg**2 + yg**2)**1.5)

#form_1_x = xg/((xg**2 + yg**2)**1.5)
#form_1_y = -yg/((xg**2 + yg**2)**1.5)

#form_1_x = 1/(xg**2 + yg**2)
#form_1_y = 1

#form_1_x = -xg/(xg**2 + yg**2)
#form_1_y = -yg/(xg**2 + yg**2)

#form_1_x = np.ones(np.shape(xg))
#form_1_y = yg

#form_1_x = -yg
#form_1_y = np.cos(xg)

form_1_x = xg
form_1_y = yg

# set up a 1-form object with these:
form_1_obj = fp.form_1(xg, yg, form_1_x, form_1_y)
plt.close(form_1_obj.figure)  # not plotting the 1-form anyway

# complete the numerical exterior derivative
form_2_num = form_1_obj.num_ext_d()

# plot it
form_2_num.plot()

# complete the analytical exterior derivative
#form_1_obj.give_eqn('-x/((x**2 + y**2)**1.5)', '-y/((x**2 + y**2)**1.5)')

#form_1_obj.give_eqn('x/((x**2 + y**2)**1.5)', '-y/((x**2 + y**2)**1.5)')

#form_1_obj.give_eqn('1/(x**2 + y**2)', '1')

#form_1_obj.give_eqn('-x/(x**2 + y**2)', '-y/(x**2 + y**2)')

#form_1_obj.give_eqn('1', 'y')

#form_1_obj.give_eqn('-y', 'cos(x)')

form_1_obj.give_eqn('x', 'y')

form_2_an = form_1_obj.ext_d()

# plot it
form_2_an.plot()


# %%

# 1-form ext deriv testing for mag field

# set up grids
v = np.linspace(-2, 2, 22)
xg, yg = np.meshgrid(v, v)

# set up the 1 form object and plot it
form_1_x = -xg/np.sqrt(xg**2 + yg**2)**3
form_1_y = -yg/np.sqrt(xg**2 + yg**2)**3

#form_1_x = 2*xg**np.exp(xg*yg) + xg**2*yg*np.exp(xg*yg)  # IMPORTANT EXAMPLE
#form_1_y = xg**3 * np.exp(xg*yg) + 2*yg

# set up a 1-form object with these:
form_1_obj = fp.form_1(xg, yg, form_1_x, form_1_y)
#form_1_obj.plot()
plt.close(form_1_obj.figure)  # not plotting the 1-form anyway

# complete the numerical exterior derivative
form_2_num = form_1_obj.num_ext_d()

# plot it
form_2_num.plot()

# complete the analytical exterior derivative
#form_1_obj.give_eqn('-x/((x**2 + y**2)**1.5)', '-y/((x**2 + y**2)**1.5)')

form_1_obj.give_eqn('2*x*e**(x*y) + x**2*y*e**(x*y)', 'x**3 * e**(x*y) + 2*y')

form_2_an = form_1_obj.ext_d()

# plot it
form_2_an.plot()

#%% Testing numpy exterior derivative

import formpy as fp
import numpy as np
import matplotlib.pyplot as plt

r = np.linspace(-5,5,15)
xg,yg = np.meshgrid(r,r)

# u = -xg/(xg**2 + yg**2)**1.5
# v = -yg/(xg**2 + yg**2)**1.5

fig1 = plt.figure(figsize=(15,5))
ax1 = fig1.add_subplot(131)
ax2 = fig1.add_subplot(132)
ax3 = fig1.add_subplot(133)
ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax3.set_aspect('equal')
ax1.set_xlabel('1-Form')
ax2.set_xlabel('Analytical d')
ax3.set_xlabel('Numerical d')

# Field is the exterior derivative of function x*cos(y), so 2-form result should be 0
u = np.cos(yg)
v = -xg*np.sin(yg)

#u = xg*np.cos(yg)
#v = yg*np.sin(xg)

field = fp.form_1(xg, yg, u, v, fig=fig1, subplots=True, sub_axis_list=[ax1, ax2, ax3])
field.give_eqn('cos(y)','-x*sin(y)')
field.plot(keep=True, subplot_index = 0)

ext_derivative = field.ext_d(pass_on_figure=True)
ext_derivative_np = field.num_ext_d(pass_on_figure=True)

ext_derivative.plot(keep=True, subplot_index = 1)
ext_derivative_np.plot(keep=True, subplot_index = 2)

#%%
ana_2f = ext_derivative.form_2()
num_2f = ext_derivative_np.form_2()

#%% Zooming inset axis

r = np.linspace(-5,5,15)
xg,yg = np.meshgrid(r,r)

fig1 = plt.figure(figsize=(7,7))
ax1 = fig1.gca()
ax1.set_aspect('equal')

u = xg*np.cos(yg)
v = yg*np.sin(xg)

field = fp.vector_field(xg, yg, u, v, fig=fig1)
field.give_eqn('x*cos(y)','y*sin(x)')

field.plot()

# The method zoom_inset creates the zoomed field, creates the inset axis and plots the field on the inset axis.
# Target, zoom and dpd carry over from original zoom.
field.zoom_inset((3, -3), 10, 9)

# %%

r = np.linspace(-2, 2, 30)
xg, yg = np.meshgrid(r, r)

u = yg**2
v = xg

dy_u = np.zeros(shape=(len(r), len(r)))
dx_v = np.zeros(shape=(len(r), len(r)))

# # Need to take the gradients of the rows of v and the cols of u

# u0 = u[:,0]
# v0 = v[0,:]

# dy_u0 = np.gradient(u0)
# dx_v0 = np.gradient(v0)

# # Store in new arrays

# dy_u[:,0] = dy_u0
# dx_v[0,:] = dx_v0

for i in range(len(r)):
    dy_u[:, i] = np.gradient(u[:, i], edge_order=2)
    dx_v[i, :] = np.gradient(v[i, :], edge_order=2)

ed = dx_v - dy_u

# Result appears to be half of what the numerical ext derivative function does.

xg2 = xg[1:-1, 1:-1]
yg2 = yg[1:-1, 1:-1]
ed2 = ed[1:-1, 1:-1]

f2 = fp.form_2(xg2, yg2, ed2)

f2.plot()

f3 = fp.form_2(xg,yg,ed)

f3.plot()

#%%

r = np.linspace(-2, 2, 30)
xg, yg = np.meshgrid(r, r)

form0 = 1/(xg**2 + yg**2)**(0.5)

f0 = fp.form_0(xg, yg, form0)

f0.plot()

f1 = f0.num_ext_d()

f1.plot()

f2 = f1.num_ext_d()

f2.plot()

# %%

r = np.linspace(-2, 2, 30)
xg, yg = np.meshgrid(r, r)

form_1_x = -0.5/(xg**2 + yg**2)
form_1_y = -0.5/(xg**2 + yg**2)

f1 = fp.form_1(xg, yg, form_1_x, form_1_y)
f1.plot()

f0 = f1.interior_d((1, 1), numerical_only=True)
f0.plot()

f1new = f0.num_ext_d()
f1new.plot()

f2 = f1new.num_ext_d()

f2.plot()



# %%

tstart = timeit.default_timer()

# plotting a 1-form

# set up needed parameters
v = np.linspace(-6, 6, 31)
xg, yg = np.meshgrid(v, v)

F_x = yg*np.sin(xg)
F_y = -xg*np.cos(yg)

# PLOT, note, it will create a figure for user
# we probably don't want that, otherwise we would have to make this
# a method to the matplotlib object, which might mean we need to play with
# their library, which I suppose we can't.
form1 = fp.form_1(xg, yg, F_x, F_y)
form1.sheet_size(0.03)
form1.plot()

form0 = fp.form_0(xg, yg, F_x*F_y)
form0.plot()

form2 = fp.form_2(xg, yg, F_x*F_y)
form2.plot()

tstop = timeit.default_timer()
print(tstop-tstart)

# %%


r = np.linspace(-5, 5, 21)
xg, yg = np.meshgrid(r, r)

fig = plt.figure(figsize=(7, 7))
ax = fig.gca()
ax.set_aspect('equal')

u = yg*np.sin(xg)
v = -xg*np.cos(yg)

field = fp.vector_field(xg, yg, u, v, fig=fig, ax=ax)
field.give_eqn('y*sin(x)','x*cos(y)')

field.plot()

# zoom
zoomed = field.zoom(target=[2, 2], zoom=1, dpd=9, fig=fig, ax=ax, inset=True)
zoomed.plot()


# %%

v = np.linspace(-6, 6, 31)
xg, yg = np.meshgrid(v, v)

form = -yg*np.sin(xg) * xg*np.cos(yg)

form2 = fp.form_2(xg, yg, form)

form2.give_eqn('-y*sin(x)*x*cos(y)')

form2_zoomed = form2.zooming(target=[-2, 2], zoom=2, dpd=9)

form2_zoomed.plot()
