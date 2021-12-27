'''

testing Exterior derivative numerical problems

'''


import formpy as fp
import numpy as np
import matplotlib.pyplot as plt
import timeit

# %%

# Test Double ext deriv giving zero as 2-form
# Numerically

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


# %%

# Testing numpy exterior derivative

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


