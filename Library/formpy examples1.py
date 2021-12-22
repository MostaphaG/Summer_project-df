# formpy examples
import formpy1 as fp
import numpy as np
import matplotlib.pyplot as plt
import timeit

# %%

'''

PLOTTING FORMS AND VF

'''

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

form_obj.colour('blue')
form_obj.head_width(0.3)
form_obj.orient('tail')
form_obj.max_sheets(6)
form_obj.sheet_size(0.03)
form_obj.surround_space(6)

# set up a plot to put these on:
fig = plt.figure(figsize=(6, 6))
ax = fig.gca()
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_aspect('equal')

form_obj.plot(ax)

tstop = timeit.default_timer()
print(tstop-tstart)

# %%

# plotting a 2-form

# set up 2-form
v = np.linspace(-6, 6, 21)
xg, yg = np.meshgrid(v, v)
form_2 = xg*yg
form_obj = fp.form_2(xg, yg, form_2)

# Create a figure and axis to plot it on
fig = plt.figure()
ax = fig.gca()
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_aspect('equal')

form_obj.plot(ax)

# wait, then change some properties and redraw
plt.pause(3)
form_obj.give_eqn('x*y')
form_obj.same_range_density(18)
ax.clear()
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_aspect('equal')
form_obj.plot(ax)

# %%

# plotting a 0-form

# set up
v = np.linspace(-4.5, 4.5, 11)
xg, yg = np.meshgrid(v, v)
form_0 = np.cos(xg*yg)
form_obj = fp.form_0(xg, yg, form_0)
form_obj.lines_number(4)
# form_obj.density_increase(20)  # demonstation of an error

# Create a figure and axis to plot it on
fig = plt.figure()
ax = fig.gca()
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_aspect('equal')

form_obj.plot(ax)


# wait, recustomise and redraw
plt.pause(2)
form_obj.give_eqn('cos(x*y)')
form_obj.density_increase(25)
form_obj.lines_number(15)
ax.clear()
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_aspect('equal')
form_obj.plot(ax)

# %%

# Plotting a vector field

# set up needed parameters
v = np.linspace(-6, 6, 31)
xg, yg = np.meshgrid(v, v)

# set up the field grids
F_x = xg/(xg**2 + yg**2)**1.5
F_y = yg/(xg**2 + yg**2)**1.5

# set up a figure, with subplots
fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax1.set_aspect('equal')
ax2.set_aspect('equal')

# set up the object
field_obj = fp.vector_field(xg, yg, F_x, F_y)

# on axis 1, plot default
field_obj.plot(ax1)

# change some properties and plot the second subplot
field_obj.colour('blue')
field_obj.orient('tail')
field_obj.autoscale()
field_obj.surround_space(6)

field_obj.plot(ax2)

# %%

# More complex example of using subplots

# set up grids
v = np.linspace(-4.5, 4.5, 21)
xg, yg = np.meshgrid(v, v)

# set up the 1 form object and plot it
form_1_x = yg
form_1_y = -xg

# set up a figure, with subplots
fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax1.set_aspect('equal')
ax2.set_aspect('equal')

# create a form object using these
form_1_obj = fp.form_1(xg, yg, form_1_x, form_1_y)

form_1_obj.plot(ax1)

star_1_form = form_1_obj.Hodge(numerical_only=True, keep_object=False)
star_1_form.plot(ax2)


# %%

# subplots with 2-forms

# set up grids
v = np.linspace(-4.5, 4.5, 21)
xg, yg = np.meshgrid(v, v)


# set up the 1 form object
form_1_x = yg
form_1_y = xg
form_1_obj = fp.form_1(xg, yg, form_1_x, form_1_y)

# set up the 2 form
form_2 = xg*yg
form_2_obj = fp.form_2(xg, yg, form_2)

# set up a figure, with subplots
fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax1.set_aspect('equal')
ax2.set_aspect('equal')

form_1_obj.plot(ax1)
form_2_obj.plot(ax2)


# %%

# Further subplot example

# set up grids
v = np.linspace(-4.5, 4.5, 11)
xg, yg = np.meshgrid(v, v)

# set up the 0 form object
form_0 = xg**2 + 3*yg
form_0_obj = fp.form_0(xg, yg, form_0)

# set up axis
fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax1.set_aspect('equal')
ax2.set_aspect('equal')

form_0_obj.plot(ax1)

# supply equation and change its density
form_0_obj.give_eqn('x**2 + 3*y')
form_0_obj.same_range_density(31)
form_0_obj.lines_number(20)

# plot that changed 0-form object on second axis set
form_0_obj.plot(ax2)


# %%

'''

Zooming

'''


#%%

# Vector field set up and zooming example

# set up VF
r = np.linspace(-4.5, 4.5, 21)
xg, yg = np.meshgrid(r, r)
u = xg*np.cos(yg)
v = -yg*np.sin(xg)
vf1 = fp.vector_field(xg, yg, u, v)
vf1.give_eqn('x*cos(y)', '-y*sin(x)')

# set up figure and axis to plot it on
fig = plt.figure()
ax = fig.gca()
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_aspect('equal')

vf1.plot(ax)

zomm_axis, zoom_vf = vf1.zoom((2, 2), zoom=3, dpd=9, inset=True, axis=ax)

#%%

# Zooming with inset is false

# set up VF
r = np.linspace(-5,5,15)
xg, yg = np.meshgrid(r, r)
u = xg*np.cos(yg)
v = yg*np.sin(xg)
field = fp.vector_field(xg, yg, u, v)
field.give_eqn('x*cos(y)','y*sin(x)')

# set up figure and axis
fig = plt.figure(figsize=(7, 7))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax1.set_aspect('equal')
ax2.set_aspect('equal')

# plot
field.plot(ax1)

# zoom in, to get field, but not inset, plot on subplot ax2
zoomed_VF = field.zoom((3, -3), 10, 9)
zoomed_VF.surround_space(4)
zoomed_VF.colour('r')

zoomed_VF.plot(ax2)

# %%

# VF zooming more

r = np.linspace(-5, 5, 21)
xg, yg = np.meshgrid(r, r)

u = yg*np.sin(xg)
v = -xg*np.cos(yg)

field = fp.vector_field(xg, yg, u, v)
field.give_eqn('y*sin(x)','x*cos(y)')

fig = plt.figure(figsize=(7, 7))
ax = fig.gca()
ax.set_aspect('equal')

field.plot(ax)

# zoom, but put down as inset, with customisations.
zoomed_axis, zoomed_field = field.zoom(target=[2, 2], zoom=10, dpd=9, inset=True, axis=ax)
zoomed_field.colour('r')
zoomed_field.plot(zoomed_axis)

# %%

# Testing 1-form zooming

r = np.linspace(-5, 5, 21)
xg, yg = np.meshgrid(r, r)

u = yg*np.sin(xg)
v = -xg*np.cos(yg)

form1 = fp.form_1(xg, yg, u, v)
form1.give_eqn('y*sin(x)','x*cos(y)')


# set up axis and plot it
fig = plt.figure(figsize=(7, 7))
ax = fig.gca()
ax.set_aspect('equal')

form1.plot(ax)

# zoom
zoomed_ax, zoomed_form = form1.zoom(target=[2, 2], zoom=10, dpd=7, inset=True, axis=ax)


# %%

# Testing 2-form zooming

# set up the 2-form
v = np.linspace(-6, 6, 31)
xg, yg = np.meshgrid(v, v)
#form = -yg*np.sin(xg) * xg*np.cos(yg)
#form2 = fp.form_2(xg, yg, form)
#form2.give_eqn('-y*sin(x)*x*cos(y)')
form = xg*yg**2
form2 = fp.form_2(xg, yg, form)
form2.give_eqn('x**2*y')

# set up a figure and axis to plot it and its inset
fig = plt.figure()
ax = fig.gca()
ax.set_aspect('equal')

form2.plot(ax)

zoomed_ax, form2_zoomed = form2.zooming(target=[2, -3], zoom=10, dpd=9, inset=True, axis=ax)

# %%

'''

VF insets with calculus

'''

# %%

# Testing derivative of vector field

# Set up field
r = np.linspace(-5, 5, 21)
xg, yg = np.meshgrid(r, r)
u = yg
v = -xg
vf1 = fp.vector_field(xg, yg, u, v)
vf1.give_eqn('y','-x')

# Set up axis
fig = plt.figure()
ax = fig.gca()
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_aspect('equal')

# plot field
vf1.plot(ax)

# find local deriv and plot as inset
D_inset_ax, D_vf1 = vf1.DF((3, 3), 5, 9, inset=True, axis=ax)

#%%

# Test divergence of vector field with subplot

# Set up fields
r = np.linspace(-5, 5, 21)
xg, yg = np.meshgrid(r, r)
u = yg*np.cos(xg)
v = -xg
vf1 = fp.vector_field(xg, yg, u, v)
vf1.give_eqn('y*cos(x)','-x')

# Set up axis
fig = plt.figure()
ax = fig.gca()
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_aspect('equal')

# plot VF
vf1.plot(ax)

# plot inset with its div, with new changed properties.
dif_axis, div_field = vf1.Div((0,2), 100, 13, inset=True, axis=ax)
#div_field.autoscale()
div_field.plot(dif_axis)

# %%

# Test curl of vector field with inset

# Set up field
r = np.linspace(-5, 5, 21)
xg, yg = np.meshgrid(r, r)
u = yg*np.cos(xg)
v = -xg

# Create vector field
vf1 = fp.vector_field(xg, yg, u, v)
vf1.give_eqn('y*cos(x)','-x')

# Set up axis
fig = plt.figure()
ax = fig.gca()
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_aspect('equal')

vf1.plot(ax)

curl_inset_ax, curl_field = vf1.Curl((0,2), 4, 13, inset=True, axis=ax)
# curl_field.autoscale()
curl_field.plot(curl_inset_ax)

# %%

'''

Exterior Alg.

'''


# %%

# testing 0-form exterior derivative

# set up grids
v = np.linspace(-4.5, 4.5, 11)
xg, yg = np.meshgrid(v, v)

# set up the 0 form object
form_0 = np.cos(xg*yg)
form_0_obj = fp.form_0(xg, yg, form_0)
form_0_obj.give_eqn('cos(x*y)')
form_0_obj.density_increase(25)
form_0_obj.lines_number(15)

# set up figure to plot on
fig = plt.figure(figsize=(7, 7))
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)
ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax3.set_aspect('equal')

# plot the 0-form
form_0_obj.plot(ax1)

# try exterior derivative without having first given an equation in
# throws an error
# form_1_obj = form_0_obj.ext_d()


# complete analytical ext. deriv.
form_1_obj = form_0_obj.ext_d()  # this supplies the 1-form with equations too

# plot that 1-form object
form_1_obj.plot(ax2)

# change its density and plot on last axis pair
form_1_obj.same_range_density(26)
form_1_obj.sheet_size(0.04)
form_1_obj.plot(ax3)

# return its equation
print(form_1_obj.return_string())

# %%

# Testing ext deriv of 1-form

# set up grids and 1 form object
v = np.linspace(-4.5, 4.5, 21)
xg, yg = np.meshgrid(v, v)
form_1_x = xg*np.cos(yg)
form_1_y = yg*np.sin(xg)
form_1_obj = fp.form_1(xg, yg, form_1_x, form_1_y)

# set up figure to plot on
fig = plt.figure(figsize=(7, 7))
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)
ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax3.set_aspect('equal')


form_1_obj.plot(ax1)

# supply equation and complete ext. deriv.
form_1_obj.give_eqn('x*cos(y)', 'y*sin(x)')

# compute the exterior derivative
form_2_obj = form_1_obj.ext_d()  # this supplies the 2-form with equations too

# plot that 2-form object
form_2_obj.plot(ax2)

# change its density and plot on last axis set
form_2_obj.same_range_density(26)
form_2_obj.max_sheets(10)
form_2_obj.plot(ax3)

# print 2-form equation
print(form_2_obj.return_string())


# %%

# exterior derivaitve of 1-form on subplots

# set up grids and 1-form
v = np.linspace(-4.5, 4.5, 21)
xg, yg = np.meshgrid(v, v)
form_1_x = yg*np.cos(xg)
form_1_y = -xg
form_1_obj = fp.form_1(xg, yg, form_1_x, form_1_y)
form_1_obj.give_eqn('y*cos(x)', '-x')

# set up figure to plot on
fig = plt.figure(figsize=(7, 7))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
ax1.set_aspect('equal')
ax2.set_aspect('equal')

# complete ext deriv.
form_2_obj = form_1_obj.ext_d()

# plot these:
form_1_obj.plot(ax1)
form_2_obj.plot(ax2)

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
fig = plt.figure()
ax = fig.add_subplot(121)
axw = fig.add_subplot(122)
ax.set_xlabel(r'$x$')
ax.set_xlabel(r'$x$')
ax.set_aspect('equal')

form_1_obj.plot(ax)

# now find wedge of it with a different form
form_wedged_2 = form_1_obj.wedge_analytical(('1/y**2', 'x**2'))

# plot it on separate axis
axw.set_xlabel(r'$x$')
axw.set_ylabel(r'$y$')
axw.set_aspect('equal')
form_wedged_2.plot(axw)

# return string
print(form_wedged_2.return_string())


# %%

# Example of Wedge done numerically

# set up grids and 1-form
v = np.linspace(-4.5, 4.5, 26)
xg, yg = np.meshgrid(v, v)
form_1_x = 2/xg**3
form_1_y = yg**2
form_1_obj = fp.form_1(xg, yg, form_1_x, form_1_y)

# set up a figure and axis
fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

# plot it
form_1_obj.plot(ax1)

# now find wedge of it with a different form
form_wedged_2 = form_1_obj.wedge_num(form_1_second=(1/yg**2, xg**3))

# plot it
form_wedged_2.plot(ax2)

# %%

# Testing the Hodge of a 1-form

# set up grids and 1-form
v = np.linspace(-4.5, 4.5, 21)
xg, yg = np.meshgrid(v, v)
form_1_x = yg
form_1_y = -xg
form_1_obj = fp.form_1(xg, yg, form_1_x, form_1_y)
form_1_obj.give_eqn('y', '-x')

# set up a figure and axis
fig = plt.figure()
ax = fig.gca()
ax.set_xlabel(r'$x$')
ax.set_xlabel(r'$x$')
ax.set_aspect('equal')

# plot the 1-form
form_1_obj.plot(ax)

# compute the Hodge and replot
form_1_obj.Hodge(numerical_only=False, keep_object=True)

plt.pause(1)
ax.clear()
ax.set_xlabel(r'$x$')
ax.set_xlabel(r'$x$')
ax.set_aspect('equal')

form_1_obj.plot(ax)

# compute the Hodge again, now set to a new form and noe numerically
form_1_obj_2H = form_1_obj.Hodge(numerical_only=True, keep_object=False)

# replot
plt.pause(1)
ax.clear()
ax.set_xlabel(r'$x$')
ax.set_xlabel(r'$x$')
ax.set_aspect('equal')
form_1_obj_2H.plot(ax)

# %%

# Example of 2-form Hodge

# set up 2-form
v = np.linspace(-6, 6, 21)
xg, yg = np.meshgrid(v, v)
form_2 = xg*yg
form_obj = fp.form_2(xg, yg, form_2)

# set up figures and axis
fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

form_obj.plot(ax1)

# supply equations
form_obj.give_eqn('x*y')

# complete analytical Hodge and plot
form_0 = form_obj.Hodge(numerical_only=False)
form_0.plot(ax2)

# complete numerical Hodge, wait and replot
form_0 = form_obj.Hodge(numerical_only=True)
plt.pause(2)
ax2.clear()
form_0.plot(ax2)  # should observe no difference (and do)

# %%

# Testing Hodge of a 0-form

# set up grids and 0-form
v = np.linspace(-4.5, 4.5, 11)
xg, yg = np.meshgrid(v, v)
form_0 = xg**2 + 3*yg
form_0_obj = fp.form_0(xg, yg, form_0)
form_0_obj.give_eqn('x**2 + 3*y')

# set up a figure with subplots
fig = plt.figure()
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

# plot it on two of the axis
form_0_obj.plot(ax1)
form_0_obj.plot(ax3)

# compute the Hodge analytically and plot the resulting 2-form
form_2_ana = form_0_obj.Hodge(numerical_only=True)
form_2_ana.plot(ax2)

# supply equations and do it analytically
form_2_num = form_0_obj.Hodge(numerical_only=False)
form_2_num.plot(ax4)


# %%

# Testing analytical interior derivative of 1-form

# set up needed parameters
r = np.linspace(-6, 6, 21)
xg, yg = np.meshgrid(r, r)

# set up form components and the object
F_x = xg**2
F_y = yg*xg
form_obj1 = fp.form_1(xg, yg, F_x, F_y)
form_obj1.give_eqn('x**2', 'y*x')

# set up components and VF object
u = np.cos(yg)
v = np.sin(xg)
vf = fp.vector_field(xg, yg, u, v)
vf.give_eqn('cos(y)', 'sin(x)')


# set up figure and axis
fig = plt.figure()
ax1 = fig.add_subplot(231)
ax2 = fig.add_subplot(232)
ax3 = fig.add_subplot(233)
ax4 = fig.add_subplot(234)
ax5 = fig.add_subplot(235)

# plot the 1-form
form_obj1.plot(ax1)

# plot it
vf.plot(ax2)

# complete interior derivative, with object, with equations and without any input:
form_0_i_1 = form_obj1.interior_d(vector_field=vf, numerical_only=False)
form_0_i_2 = form_obj1.interior_d(numerical_only=False)
form_0_i_3 = form_obj1.interior_d(('cos(y)', 'sin(x)'), numerical_only=False)

# plot each
form_0_i_1.plot(ax3)
form_0_i_2.plot(ax4)
form_0_i_3.plot(ax5)

# %%

# Testing numerical interior derivative of a 1-form:

# set up needed parameters
r = np.linspace(-6, 6, 21)
xg, yg = np.meshgrid(r, r)

# set up form object
F_x = xg**2
F_y = yg*xg
form_obj1 = fp.form_1(xg, yg, F_x, F_y)

# set up a vector field object:
u = np.cos(yg)
v = np.sin(xg)
vf = fp.vector_field(xg, yg, u, v)

# set up figure and axis
fig = plt.figure()
ax1 = fig.add_subplot(231)
ax2 = fig.add_subplot(232)
ax3 = fig.add_subplot(233)
ax4 = fig.add_subplot(234)
ax5 = fig.add_subplot(235)

# plot form
form_obj1.plot(ax1)

# plot VF:
vf.plot(ax2)

# complete int deriv. with object, arrays and with no input:
form_0_i_1 = form_obj1.interior_d(vector_field=vf, numerical_only=True)
form_0_i_2 = form_obj1.interior_d(numerical_only=True)
form_0_i_3 = form_obj1.interior_d((u, v), numerical_only=True)

# plot each
form_0_i_1.plot(ax3)
form_0_i_2.plot(ax4)
form_0_i_3.plot(ax5)

# %%

# Testing numerical interior derivative of a 2-form

# set up needed parameters
r = np.linspace(-6, 6, 21)
xg, yg = np.meshgrid(r, r)

# set up form
form2_comp = xg*yg**2
form_obj2 = fp.form_2(xg, yg, form2_comp)

# create figure with subplots to put these all on:
fig = plt.figure()
ax1 = fig.add_subplot(151)
ax2 = fig.add_subplot(152)
ax3 = fig.add_subplot(153)
ax4 = fig.add_subplot(154)
ax5 = fig.add_subplot(155)
ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax3.set_aspect('equal')
ax4.set_aspect('equal')
ax5.set_aspect('equal')

# plot it
form_obj2.plot(ax1)

# set up a vector field object:
u = np.cos(yg)
v = np.sin(xg)

vf = fp.vector_field(xg, yg, u, v)

# plot it:
vf.plot(ax2)

# complete int deriv. with object, arrays and with no input:
form_1_i_1 = form_obj2.interior_d(vector_field=vf, numerical_only=True)
form_1_i_2 = form_obj2.interior_d(numerical_only=True)
form_1_i_3 = form_obj2.interior_d((u, v), numerical_only=True)

# plot each
form_1_i_1.plot(ax3)
form_1_i_2.plot(ax4)
form_1_i_3.plot(ax5)


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


# create figure with subplots to put these all on:
fig = plt.figure()
ax1 = fig.add_subplot(151)
ax2 = fig.add_subplot(152)
ax3 = fig.add_subplot(153)
ax4 = fig.add_subplot(154)
ax5 = fig.add_subplot(155)
ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax3.set_aspect('equal')
ax4.set_aspect('equal')
ax5.set_aspect('equal')


# plot it
form_obj2.plot(ax1)

# set up a vector field object:
u = np.cos(yg)
v = np.sin(xg)

vf = fp.vector_field(xg, yg, u, v)

# give it the equation
vf.give_eqn('cos(y)', 'sin(x)')

# plot it:
vf.plot(ax2)

# complete int deriv. with object, arrays and with no input:
form_1_i_1 = form_obj2.interior_d(vector_field=vf, numerical_only=False)
form_1_i_2 = form_obj2.interior_d(numerical_only=False)
form_1_i_3 = form_obj2.interior_d(('cos(y)', 'sin(x)'), numerical_only=False)

# plot each
form_1_i_1.plot(ax3)
form_1_i_2.plot(ax4)
form_1_i_3.plot(ax5)

# %%

# Test Double ext deriv giving zero as 2-form
# Analytically

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

'''

Curved Spacetime

'''


# %%

# Testing conversion form vector field to 1-form with string based metric

r = np.linspace(-4.5, 4.5, 21)
xg, yg = np.meshgrid(r, r)

# set up the field and vector field objecct
u = np.ones(np.shape(xg))
v = yg
vf1 = fp.vector_field(xg, yg, u, v)
vf1.give_eqn('1', 'y')

# set up figure and axis to plot on
fig = plt.figure()
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)
ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax3.set_aspect('equal')

# plot the VF
vf1.plot(ax1)

# set up a metric with strings
metric = [['1', '0'],
          ['0', '(x**2 + y**2)']]

# via this, set up a 1-form
form_1_obj = vf1.formalise(g=metric)

# plot it
form_1_obj.plot(ax2)

# create a comparison 1-form with correct components already given
form_1_correct = fp.form_1(xg, yg, u, v*(xg**2+yg**2))

# plot it
form_1_correct.plot(ax3)

# %%

# Testing conversion form vector field to 1-form with array based metric

r = np.linspace(-4.5, 4.5, 21)
xg, yg = np.meshgrid(r, r)

# set up the field and vector field objecct
u = np.ones(np.shape(xg))
v = yg
vf1 = fp.vector_field(xg, yg, u, v)

# set up figure and axis to plot on
fig = plt.figure()
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)
ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax3.set_aspect('equal')

# plot VF
vf1.plot(ax1)

# set up a metric with arrays
metric = [[np.ones(np.shape(xg)), np.zeros(np.shape(xg))],
          [np.zeros(np.shape(xg)), yg**2 + xg**2]]

# via this, set up a 1-form
form_1_obj = vf1.formalise(g=metric)

# plot it
form_1_obj.plot(ax2)

# create a comparison 1-form with correct components already given
form_1_correct = fp.form_1(xg, yg, u, v*(xg**2 + yg**2))

# plot it
form_1_correct.plot(ax3)

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


'''

Tricky examples - const comps. singularities etc

'''

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

