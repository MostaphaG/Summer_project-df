'''

File for testing, separate to any official codes.

'''

import numpy as np
import formpy as fp
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import patches as patch
from math import isnan
import timeit
import time

# input many numpy functions to deal with user input
from numpy import sin, cos, tan, sqrt, log, arctan, arcsin, arccos, tanh
from numpy import sinh, cosh, arcsinh, arccosh, arctanh, exp, pi, e

# %%

'''
Testing insets
'''

fig = plt.figure(figsize=(8,8))
ax = fig.gca()
r = np.linspace(-5,5,21)
xg, yg = np.meshgrid(r,r)
u = xg*np.cos(yg)
v = -yg*np.sin(xg)
ax.set_aspect('equal')

insize = 0.3

plt.quiver(xg,yg,u,v)

x_m, y_m = (4,0)
x_range = xg[0,-1] - xg[0,0]
y_range = yg[-1,0] - yg[0,0]

xi = (x_m - xg[0,0])/x_range
yi = (y_m - yg[0,0])/y_range

q = 0.92
inax = ax.inset_axes([(xi - q*0.5*insize), (yi - q*0.5*insize), insize, insize])

# %%

'''
Central force check
'''

# set up needed parameters
v = np.linspace(-4, 4, 21)
xg, yg = np.meshgrid(v, v)

F_x = xg/((xg**2 + yg**2)**(1.5))
F_y = yg/((xg**2 + yg**2)**(1.5))
form_obj = fp.form_1(xg, yg, F_x, F_y)

# set up a plot to put these on:
fig = plt.figure(figsize=(6, 6))
ax = fig.gca()
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_aspect('equal')

# plot
form_obj.plot(ax)


# set up a similar VF to test


VF = fp.vector_field(xg, yg, F_x, F_y)

# set up a plot to put these on:
fig1 = plt.figure(figsize=(6, 6))
ax1 = fig1.gca()
ax1.set_xlabel(r'$x$')
ax1.set_ylabel(r'$y$')
ax1.set_aspect('equal')

# plot
VF.plot(ax1)


# Similarly with 2-form

form_obj2 = fp.form_2(xg, yg, F_x*F_y)

# set up a plot to put these on:
fig2 = plt.figure(figsize=(6, 6))
ax2 = fig2.gca()
ax2.set_xlabel(r'$x$')
ax2.set_ylabel(r'$y$')
ax2.set_aspect('equal')

# plot
form_obj2.plot(ax2)

# Similarly with 0-form

form_obj0 = fp.form_0(xg, yg, F_x*F_y)
form_obj0.density_increase(5)
form_obj0.give_eqn('x*y/(x**2 + y**2)**(1.5)')

# set up a plot to put these on:
fig3 = plt.figure(figsize=(6, 6))
ax3 = fig3.gca()
ax3.set_xlabel(r'$x$')
ax3.set_ylabel(r'$y$')
ax3.set_aspect('equal')

# plot
form_obj0.plot(ax3)


'''
Messy with all the figures, I know, but its just for quick tests for singularities
'''

# %%

'''

Attempts at optimising stackplot

'''

# ############################################################################
# Define stackplot function, equivalently to method in formpy
# ############################################################################


def G(s, n, c):
    if c == 0:
        return ((2*s + 1)/(2*(n-1)))
    else:
        return (s/(n-1))


def stackplot(axis, xg, yg, F_x, F_y, s_max=6, s_min=2, color='green', w_head=1/8, h_head=1/4, fract=0.05, delta_factor=10, logarithmic_scale_bool=0, arrowheads=True):
    global indexes
    # get the lengths of x and y from their grids
    x_len = len(xg[:, 0])
    y_len = len(yg[0, :])
    
    # Extract L from the x and y grids. Assumes they are square.
    L = 0.5*(xg[0, -1] - xg[0, 0])
    x0 = xg[0,0] + L
    y0 = yg[0,0] + L
    
    # scale to size accordinly to deltafactor
    ax_L = L + L/delta_factor
    axis.set_xlim(-ax_L + x0, ax_L + x0)
    axis.set_ylim(-ax_L + y0, ax_L + y0)
    
    # find the distance between neightbouring points on the grid
    dist_points = xg[0, 1] - xg[0, 0]
    
    # define an empty array of magnitudes, to then fill with integer rel. mags
    # these will dtermine, number of stacks to plot for each mag
    R_int = np.zeros(shape=((x_len), (y_len)))

    # #########################################################################
    # get variables needed for the initial, simplified stack plot
    # #########################################################################
    
    # set all insignificant values to zero:
    F_x[np.abs(F_x) < 1e-15] = 0
    F_x[np.abs(F_x) < 1e-15] = 0
    
    # find the arrow length corresponding to each point and store in mag array
    mag = np.sqrt(F_x**2 + F_y**2)
    
    # find direction of each arrow
    angles = np.arctan2(F_y, F_x)   # theta defined from positive x axis ccw
    
    
    isnan_arr = np.isnan(mag)
    for i in range(x_len):
        for j in range(y_len):
            # set to zero points that are not defined or inf
            if isnan_arr[i, j]:
                #colour this region as a shaded square
                rect = patch.Rectangle((xg[i, j] - dist_points/2, yg[i, j]  - dist_points/2), dist_points, dist_points, color='#B5B5B5')
                axis.add_patch(rect)
                mag[i, j] = 0
            if abs(mag[i, j]) == np.inf  or abs(mag[i, j]) > 1e15:
                # colour this point as a big red dot
                circ = patch.Circle((xg[i, j], yg[i, j]), L*fract/3, color='red')
                axis.add_patch(circ)
                mag[i, j] = 0

#    isnan_arr = np.isnan(mag)
#    indexes = np.vstack(np.nonzero(isnan_arr))
#    if indexes.size != 0:
#        for i in indexes[:, 0]:
#            for j in indexes[:, 1]:
#                rect = patch.Rectangle((xg[i, j] - dist_points/2, yg[i, j]  - dist_points/2), dist_points, dist_points, color='#B5B5B5')
#                axis.add_patch(rect)
#                mag[i, j] = 0
#    for i in range(x_len):
#        for j in range(y_len):
#            if abs(mag[i, j]) == np.inf  or abs(mag[i, j]) > 1e15:
#                # colour this point as a big red dot
#                circ = patch.Circle((xg[i, j], yg[i, j]), L*fract/3, color='red')
#                axis.add_patch(circ)
#                mag[i, j] = 0
    
    # #########################################################################
    # use the the direction of arrows to define stack properties
    # #########################################################################
    
    # s_L and sheet_L are defined the same? 
    
    # define length of sheet as a fraction of total graph scale
    # sheet_L = L * fract
    # set up the max, total height of stack (along arrow)
    
    s_L = fract * L
    
    # #########################################################################
    # define the stacks based on geometrical arguments
    # to be perp. to arrow. shifted parallel to it, their density porp to mag
    # of the arrow and with an arrowhead on top.
    # #########################################################################
    
    # find the maximum magnitude for scaling
    max_size = np.max(mag)   # careful with singularities, else ---> nan
    
    # Define scaling factor
    #ScaleFactor = max_size/(0.9*(2*L/pt_den))
    
    # find the relative magnitude of vectors to maximum, as an array
    R = mag/max_size
    
    # logarithmic attempt
    if logarithmic_scale_bool == 1:
        log_a = 1000000
        R = np.where(R<=1/log_a, 1/log_a, R)  # Remove the values less than critical for log
        R = log(log_a*R)/log(log_a)
    else:
        pass

    # define tigonometirc shifts
    I_sin = np.sin(angles)
    I_cos = np.cos(angles)
    
    # define the points that set out a line of the stack sheet (middle line)
    A_x = xg + (s_L/2)*I_sin
    A_y = yg - (s_L/2)*I_cos
    B_x = xg - (s_L/2)*I_sin
    B_y = yg + (s_L/2)*I_cos
    
    # define points of stack arrowheads as arrays for all stacks
    p_sh1x = xg + (s_L/2)*I_cos + (s_L*w_head)*I_sin
    p_sh1y = yg + (s_L/2)*I_sin - (s_L*w_head)*I_cos
    p_sh2x = xg + (s_L/2)*I_cos - (s_L*w_head)*I_sin
    p_sh2y = yg + (s_L/2)*I_sin + (s_L*w_head)*I_cos
    p_sh3x = xg + (s_L*0.5 + s_L*h_head)*I_cos
    p_sh3y = yg + (s_L*0.5 + s_L*h_head)*I_sin
    
    # define these for when there is only 1 line in the stack plot:
    P_sh1x = xg + (s_L*w_head)*I_sin
    P_sh1y = yg - (s_L*w_head)*I_cos
    P_sh2x = xg - (s_L*w_head)*I_sin
    P_sh2y = yg + (s_L*w_head)*I_cos
    P_sh3x = xg + (s_L*h_head)*I_cos
    P_sh3y = yg + (s_L*h_head)*I_sin
    
    # loop over each arrow coordinate in x and y
    for i in range(x_len):
        for j in range(y_len):
            
            # Label each element with the number of stacks required: linear scaling
            for t in range(1, s_max+1):
                if (t-1)/s_max <= R[i, j] <= t/s_max:
                    R_int[i, j] = t
            
            # set a varible for current considered magnitude as it is reused
            # avoids extracting from R many times.
            n = R_int[i, j]
            
            #if axis_check == 1 and click_opt_int > 1 and i == i_m and j == j_m:
            if mag[i,j] == 0:
                continue
            
            # deal with even number of sheets from magnitudes:
            if n % 2 == 0:
                # define a parameter to loop over in the recursion equation
                s = 0
                
                # Define the points for sheets required for the given magnitude
                # from these define all the needed lines and plot them
                while s <= 0.5*(n-2):  # maximum set by equations (documentation)
                    # define all the points for the 2 currently looped +- sheets in while loop
                    Ax1 = A_x[i, j] + G(s, n, 0)*s_L*I_cos[i, j]
                    Ay1 = A_y[i, j] + G(s, n, 0)*s_L*I_sin[i, j]
                    Bx1 = B_x[i, j] + G(s, n, 0)*s_L*I_cos[i, j]
                    By1 = B_y[i, j] + G(s, n, 0)*s_L*I_sin[i, j]
                    Ax2 = A_x[i, j] - G(s, n, 0)*s_L*I_cos[i, j]
                    Ay2 = A_y[i, j] - G(s, n, 0)*s_L*I_sin[i, j]
                    Bx2 = B_x[i, j] - G(s, n, 0)*s_L*I_cos[i, j]
                    By2 = B_y[i, j] - G(s, n, 0)*s_L*I_sin[i, j]
                    
                    # from these, define the 2 lines, for this run
                    axis.add_line(Line2D((Ax1, Bx1), (Ay1, By1), linewidth=1, color=color))
                    axis.add_line(Line2D((Ax2, Bx2), (Ay2, By2), linewidth=1, color=color))
                    
                    # update parameter to reapet and draw all needed arrows
                    s += 1
            # deal with the odd number of stacks:
            else:
                # Add the centre line for odd numbers of stacks
                axis.add_line(Line2D((A_x[i, j], B_x[i, j]), (A_y[i, j], B_y[i, j]), linewidth=1, color=color))
                
                # then loop over the remaining lines as per the recursion formula:
                s = 1  # change the looping parametr to exclude already completed 0 (corr. to middle sheet here)
                
                # define all remaining sheets for the magnitude:
                while s <= 0.5*(n-1):  # maximum set by equations (documentation)
                    # define all the points for the current +- displacement in while loop
                    Ax1 = A_x[i, j] + G(s, n, 1)*s_L*I_cos[i, j]
                    Ay1 = A_y[i, j] + G(s, n, 1)*s_L*I_sin[i, j]
                    Bx1 = B_x[i, j] + G(s, n, 1)*s_L*I_cos[i, j]
                    By1 = B_y[i, j] + G(s, n, 1)*s_L*I_sin[i, j]
                    Ax2 = A_x[i, j] - G(s, n, 1)*s_L*I_cos[i, j]
                    Ay2 = A_y[i, j] - G(s, n, 1)*s_L*I_sin[i, j]
                    Bx2 = B_x[i, j] - G(s, n, 1)*s_L*I_cos[i, j]
                    By2 = B_y[i, j] - G(s, n, 1)*s_L*I_sin[i, j]
                    
                    # from these, define the 2 displaced lines
                    axis.add_line(Line2D((Ax1,Bx1),(Ay1,By1), linewidth=1, color=color))
                    axis.add_line(Line2D((Ax2,Bx2),(Ay2,By2), linewidth=1, color=color))
                    
                    # change the parameter to loop over all changes in displacement for current magnitude
                    s += 1
            if arrowheads == True:
                # plot lines of arrowheads from central sheet for n = 1 or on top sheet for n>1 
                if n > 1:   # for all lines ubt the single sheet one
                    axis.add_line(Line2D((p_sh1x[i, j],p_sh3x[i, j]),(p_sh1y[i, j],p_sh3y[i, j]), linewidth=1, color = color))
                    axis.add_line(Line2D((p_sh2x[i, j],p_sh3x[i, j]),((p_sh2y[i, j],p_sh3y[i, j])), linewidth=1, color = color))
                # then define it for the stacks with only 1 sheet:
                else:
                    axis.add_line(Line2D((P_sh1x[i, j], P_sh3x[i, j]), (P_sh1y[i, j], P_sh3y[i, j]), linewidth=1, color = color))
                    axis.add_line(Line2D((P_sh2x[i, j], P_sh3x[i, j]), ((P_sh2y[i, j], P_sh3y[i, j])), linewidth=1, color = color))
            else:
                pass
# %%
# ############################################################################
# Run it and precisely time
# ############################################################################

# use a precision counter
start = time.perf_counter()

# call the fucntion for some parameters:
r = np.linspace(-5, 5, 41)
xg, yg = np.meshgrid(r, r)

F_x = 1/np.sin(xg+yg)
F_y = np.exp(xg) + 1/yg

fig = plt.figure()
ax = fig.gca()
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.tick_params(labelsize=20)

stackplot(ax, xg, yg, F_x, F_y, fract=0.03)

# stop and return precise time to run
stop = time.perf_counter()
print('time to run: {:f}'.format(stop-start))


# %%

# testing the red dots on 1-form

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
ana_int = form2.interior_d(VF, numerical_only=False)

# plot these
ana_int.plot(ax3)

# use cross product:
VF_c = fp.vector_field(rhog, zg, -1/rhog, np.zeros(np.shape(zg)))
VF.give_eqn('-1/x', '0')

VF_c.log_scaling()
VF_c.plot(ax4)


# %%


# Test s_min

v = np.linspace(-5, 5, 21)
xg, yg = np.meshgrid(v, v)

F_x = xg
F_y = yg

form_obj = fp.form_1(xg, yg, F_x, F_y)
form_obj.s_min = 4
form_obj.max_sheets(12)
form_obj.colour('#AFAFAF')
form_obj.give_eqn('x', 'y')

fig = plt.figure()
ax = fig.gca()
ax.set_aspect('equal')

form_obj.plot(ax)

# %%

# tets girds off origin and not square for 1-form


# set up needed parameters
x = np.linspace(-2, 6, 21)
y = np.linspace(-3, 2, 31)
xg, yg = np.meshgrid(x, y)

F_x = 10*xg*yg
F_y = 1/(np.sin(yg))

# PLOT, note, it will create a figure for user
# we probably don't want that, otherwise we would have to make this
# a method to the matplotlib object, which might mean we need to play with
# their library, which I suppose we can't.
form_obj = fp.form_1(xg, yg, F_x, F_y)

form_obj.colour('blue')
form_obj.head_width(0.3)
form_obj.s_min = 2
form_obj.max_sheets(7)
form_obj.sheet_size(0.03)
form_obj.surround_space(6)

# set up a plot to put these on:
fig = plt.figure(figsize=(6, 6))
ax = fig.gca()
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_aspect('equal')

form_obj.plot(ax)


# %%

# Changing colours for 2-form

v = np.linspace(-6, 6, 21)
xg, yg = np.meshgrid(v, v)
form_2 = xg*yg
form_obj = fp.form_2(xg, yg, form_2)
form_obj.colours(['orange', 'blue', 'black'])

# Create a figure and axis to plot it on
fig = plt.figure()
ax = fig.gca()
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_aspect('equal')

form_obj.plot(ax)

# %%

p = np.linspace(-4, 4, 23)
q = np.linspace(0, 4, 23)
x, y = np.meshgrid(p, q)


u = np.ones(np.shape(x))
v = np.tanh(x)*(np.cosh(x))**(2/3)
fp1 = fp.form_1(x,y,u,v)

plt.figure(1,figsize=(6,6))
ax = plt.axes()
fp1.plot(ax)

# %%

# Testing default colours in 2-forms

v = np.linspace(-6, 6, 21)
xg, yg = np.meshgrid(v, v)
form_2 = np.ones(np.shape(xg))
form_obj = fp.form_2(xg, yg, form_2)

# Create a figure and axis to plot it on
fig = plt.figure()
ax = fig.gca()
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_aspect('equal')

form_obj.plot(ax)

# Should be positive in dx/\dy therefore counterclockwise, this is red.
# And it is.

# %%

# plotting a 0-form with custom levels

# set up
v = np.linspace(-4.5, 4.5, 11)
xg, yg = np.meshgrid(v, v)
form_0 = 1/np.sqrt(xg**2 + yg**2) * np.exp(-1*np.sqrt(xg**2 + yg**2))
form_obj = fp.form_0(xg, yg, form_0)
form_obj.levels([0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2])
form_obj.give_eqn('1/sqrt(x**2 + y**2) * e**(-sqrt(x**2 + y**2))')
form_obj.density_increase(10)
form_obj.labels()

# Create a figure and axis to plot it on
fig = plt.figure()
ax = fig.gca()
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.set_aspect('equal')

form_obj.plot(ax)

# %%

# New 1-form Wedge

# set up 1-form
v = np.linspace(-2, 2, 23)
xg, yg = np.meshgrid(v, v)
u = xg*yg
v = np.sin(xg*yg)
f1 = fp.form_1(xg, yg, u, v, 'x*y', 'sin(x*y)')

# another one
f1a = fp.form_1(xg, yg, u*v, v, 'x*y*sin(x*y)', 'sin(x*y)')

# set up 0-form
f0 = fp.form_0(xg, yg, u, 'x*y')

# set up 2-form
f2 = fp.form_2(xg, yg, v, 'sin(x*y)')

# try many wedges

'''
1-form/\1-form
'''
f2a = f1.wedge(f1a)
f2b = f1.wedge(('x*y*sin(x*y)', 'sin(x*y)'))
f2c = f1.wedge()
#f2d = f1.wedge(('13'))

# plot the working ones
fig = plt.figure()
ax1 = fig.add_subplot(131)
ax2 = fig.add_subplot(132)
ax3 = fig.add_subplot(133)
ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax3.set_aspect('equal')

f2a.plot(ax1)
f2b.plot(ax2)
f2c.plot(ax3)

'''
1-form/\0-form
'''

f1a = f1.wedge(f0)
f1b = f1.wedge('x*y*sin(x*y)', degree=0)
#f1d = f1.wedge(('13', 'x'))
f1c = f1.wedge()

# plot the working ones
fig1 = plt.figure()
ax1 = fig1.add_subplot(131)
ax2 = fig1.add_subplot(132)
ax3 = fig1.add_subplot(133)
ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax3.set_aspect('equal')

f1a.plot(ax1)
f1b.plot(ax2)
f1c.plot(ax3)

'''
1-form/\2-form
'''

f3a = f1.wedge(f2)
f3b = f1.wedge('x*y*cos(x*y)', degree=2)


# %%

# Same numerically on 1-form

v = np.linspace(-2, 2, 23)
xg, yg = np.meshgrid(v, v)
u = xg*yg
v = np.sin(xg*yg)
f1 = fp.form_1(xg, yg, u, v)

f1a = fp.form_1(xg, yg, u*v, v)

f0 = fp.form_0(xg, yg, u)

f2 = fp.form_2(xg, yg, v)

'''
1-form/\1-form
'''

f2a = f1.num_wedge(f1a)
f2b = f1.num_wedge(('x*y*sin(x*y)', 'sin(x*y)'))
f2c = f1.num_wedge((xg*yg*np.sin(xg*yg), np.sin(xg*yg)))
f2d = f1.num_wedge()
#f2d = f1.wedge(('13'))

# plot the working ones
fig = plt.figure()
ax1 = fig.add_subplot(141)
ax2 = fig.add_subplot(142)
ax3 = fig.add_subplot(143)
ax4 = fig.add_subplot(144)
ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax3.set_aspect('equal')
ax4.set_aspect('equal')

f2a.plot(ax1)
f2b.plot(ax2)
f2c.plot(ax3)
f2d.plot(ax4)

'''
1-form/\0-form
'''

f1a = f1.num_wedge(f0)
f1b = f1.num_wedge('x*y*sin(x*y)', degree=0)
f1c = f1.num_wedge(xg*yg*np.sin(xg*yg), degree=0)
#f1d = f1.wedge(('13', 'x'))
f1d = f1.num_wedge()

# plot the working ones
fig1 = plt.figure()
ax1 = fig1.add_subplot(141)
ax2 = fig1.add_subplot(142)
ax3 = fig1.add_subplot(143)
ax4 = fig1.add_subplot(144)

ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax3.set_aspect('equal')
ax4.set_aspect('equal')

f1a.plot(ax1)
f1b.plot(ax2)
f1c.plot(ax3)
f1d.plot(ax4)

'''
1-form/\2-form
'''

f3a = f1.num_wedge(f2)
f3b = f1.num_wedge('x*y*cos(x*y)', degree=2)


# %%

# TESTING 2-form wegdges analytically

v = np.linspace(-2, 2, 23)
xg, yg = np.meshgrid(v, v)

f2 = fp.form_2(xg, yg, np.sin(xg*yg), 'sin(x*y)')

f0 = fp.form_0(xg, yg, xg*yg)
f0.give_eqn('x*y')

f2a = f2.wedge(f0)
f2b = f2.wedge('x**2')
f2c = f2.wedge()

fig1 = plt.figure()
ax1 = fig1.add_subplot(131)
ax2 = fig1.add_subplot(132)
ax3 = fig1.add_subplot(133)
ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax3.set_aspect('equal')

f2a.plot(ax1)
f2b.plot(ax2)
f2c.plot(ax3)

# %%

# TESTING 2-form wegdges numerically

v = np.linspace(-2, 2, 23)
xg, yg = np.meshgrid(v, v)

f2 = fp.form_2(xg, yg, np.sin(xg*yg))

f0 = fp.form_0(xg, yg, xg*yg)

f2a = f2.num_wedge(f0)
f2b = f2.num_wedge('x**2')
f2c = f2.num_wedge(yg**2)
f2d = f2.num_wedge()

fig1 = plt.figure()
ax1 = fig1.add_subplot(141)
ax2 = fig1.add_subplot(142)
ax3 = fig1.add_subplot(143)
ax4 = fig1.add_subplot(144)

ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax3.set_aspect('equal')
ax4.set_aspect('equal')

f2a.plot(ax1)
f2b.plot(ax2)
f2c.plot(ax3)
f2d.plot(ax4)

# %%

# TESTING 0-form wegdges analytically

v = np.linspace(-2, 2, 23)
xg, yg = np.meshgrid(v, v)

u = xg*yg
v = np.sin(xg*yg)

f0 = fp.form_0(xg, yg, v, 'sin(x*y)')

# another 0-form
f0a = fp.form_0(xg, yg, u, 'x*y')

f1 = fp.form_1(xg, yg, u, v, 'x*y', 'sin(x*y)')

# another one
f1a = fp.form_1(xg, yg, u*v, v, 'x*y*sin(x*y)', 'sin(x*y)')

# set up 2-form
f2 = fp.form_2(xg, yg, v, 'sin(x*y)')

f0a = f0.wedge(f0)
f0b = f0.wedge(f0a)
f1a = f0.wedge(f1)
f1b = f0.wedge(f1a)
f2a = f0.wedge(f2)

f0c = f0.wedge('x + y')
f0d = f0.wedge()
f1c = f0.wedge(('x**2', 'y**2'))
f1d = f0.wedge()
f2b = f0.wedge('x + 2*y')
f2c = f0.wedge()

# %%

# TESTING 0-form wegdges numerically

v = np.linspace(-2, 2, 23)
xg, yg = np.meshgrid(v, v)

u = xg*yg
v = np.sin(xg*yg)

f0 = fp.form_0(xg, yg, v, 'sin(x*y)')

# another 0-form
f0a = fp.form_0(xg, yg, u, 'x*y')

f1 = fp.form_1(xg, yg, u, v, 'x*y', 'sin(x*y)')

# another one
f1a = fp.form_1(xg, yg, u*v, v, 'x*y*sin(x*y)', 'sin(x*y)')

# set up 2-form
f2 = fp.form_2(xg, yg, v, 'sin(x*y)')

f0a = f0.wedge(f0)
f0b = f0.wedge(f0a)
f1a = f0.wedge(f1)
f1b = f0.wedge(f1a)
f2a = f0.wedge(f2)

f0c = f0.wedge('x + y')
f0d = f0.wedge()
f1c = f0.wedge(('x**2', 'y**2'))
f1d = f0.wedge()
f2b = f0.wedge('x + 2*y')
f2c = f0.wedge()




