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

Attempts at optimising stackplot

'''

# ############################################################################
# Define stackplot function, equivalently to method in formpy
# ############################################################################


def parity(x):
    if x % 2 == 1:
        result = False
    elif x % 2 == 0:
        result = True
    return result


def G(s, n, c):
    if c == 0:
        return ((2*s + 1)/(2*(n-1)))
    else:
        return (s/(n-1))


def stackplot(axis, xg, yg, F_x, F_y, s_max=6, s_min=2, color='green', w_head=1/8, h_head=1/4, fract=0.05, delta_factor=10, logarithmic_scale_bool=0, arrowheads=True):
    # get the lengths of x and y from their grids
    x_len = len(xg[:, 0])
    y_len = len(yg[0, :])
    
    # Extract L from the x and y grids. Assumes they are square.
    L = 0.5*(xg[0, -1] - xg[0, 0])
    x0 = xg[0,0] + L
    y0 = yg[0,0] + L
    
    ax_L = L + L/delta_factor
    axis.set_xlim(-ax_L + x0, ax_L + x0)
    axis.set_ylim(-ax_L + y0, ax_L + y0)
    
    # find the distance between neightbouring points on the grid
    dist_points = xg[0, 1] - xg[0, 0]
    
    # define an empty array of magnitudes, to then fill with integer rel. mags
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
    
    # find regions ON GRID that are nan or inf as a bool array
    #bool_array = undef_region(mag)
    
    # deal with infs and nans in mag
    for i in range(x_len):
        for j in range(y_len):
            # set to zero points that are not defined or inf
            if isnan(mag[i, j]) is True:
                #colour this region as a shaded square
                rect = patch.Rectangle((xg[i, j] - dist_points/2, yg[i, j]  - dist_points/2), dist_points, dist_points, color='#B5B5B5')
                axis.add_patch(rect)
                mag[i, j] = 0
            if abs(mag[i, j]) == np.inf  or abs(mag[i, j]) > 1e15:
                # colour this point as a big red dot
                circ = patch.Circle((xg[i, j], yg[i, j]), L*fract/3, color='red')
                axis.add_patch(circ)
                mag[i, j] = 0
    
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
            if parity(n) is True:
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
            elif parity(n) is False:
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

# ############################################################################
# Run it and precisely time
# ############################################################################

# use a precision counter
start = time.perf_counter()

# call the fucntion for some parameters:
r = np.linspace(-5, 5, 41)
xg, yg = np.meshgrid(r, r)

F_x = np.sin(xg+yg)**2
F_y = np.exp(xg) + yg

fig = plt.figure()
ax = fig.gca()
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$y$')
ax.tick_params(labelsize=20)

stackplot(ax, xg, yg, F_x, F_y, fract=0.03)

# stop and return precise time to run
stop = time.perf_counter()
print('time to run: {:f}'.format(stop-start))



