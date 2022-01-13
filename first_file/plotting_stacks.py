# Lats attempt to plot stack fileds - optimised
# import nec. modules
import timeit
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# %%

# start the timer
start = timeit.default_timer()


# deifne a fucntion to check number parity
def parity(x):
    '''
    return True when number provided is even and False when it is odd. Checks
    integers only, and not arrays.
    '''
    if x % 2 == 1:
        result = False
    elif x % 2 == 0:
        result = True
    return result


# define scale of the graph
L = 5
pt_den = 20   # number of points on each axis

# define x and y values
x = np.linspace(-L, L, pt_den)
y = np.linspace(-L, L, pt_den)

# create a grid on x-y plane
xg, yg = np.meshgrid(x, y)

# define an example vector field
a = 0.05
u = a*xg*(1-yg)  # x component
v = a*yg*(xg-1)  # y component

'''
To define it in polar:

hash out the above u and v, then:

# define polar axis
theta = np.arange(0, 361, 10) * np.pi/180
radius = np.arange(0.2, 1, 0.1)

# set up grid in polar coordinates
thetag, rg = np.meshgrid(theta, radius)

# set up wanted constants
a = 1

# define the field in polar coordiantes
Fr = thetag
Ftheta = - a/rg

# convert to cartesian
u = Fr*np.cos(thetag) - Ftheta*np.sin(thetag)  # x component
v = Fr*np.sin(thetag) + Ftheta*np.cos(thetag)  # y component

CONTINUE FROM THERE
'''

# for no dependance, use : np.zeros(np.shape(grid))

''' plot the initial quiver plot to work from '''

# create a figure
fig = plt.figure(figsize=(8, 8))

# set up axis
ax = fig.gca()

# set up quiver factors
orientation = 'mid'  # how the arrow rotates about its assigned grid point - options: tail, mid and tip as string
scale = 1  # the scale reduction factor, if None (as None-type, not str), automatically computed by average, if 1 = mag

# plot the quiver plot on grid points
ax.quiver(xg, yg, u, v, pivot=orientation, scale=scale, scale_units='xy')

# set up visuals - axis labels
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')

# Scaling of axes and setting equal proportions circles look like circles
ax.set_aspect('equal')
ax_L = L + L/10
ax.set_xlim(-ax_L, ax_L)
ax.set_ylim(-ax_L, ax_L)

''' get variables needed for the initial, simplified stack plot '''

# find the arrow length corresponding to each point and store in mag array
mag = np.sqrt(u**2 + v**2)

# find direction of each arrow
theta = np.arctan2(v, u)   # theta defined from positive x axis ccw

''' use the the direction of arrows to define stack properties '''

# define length of sheet as a fraction of total graph scale
fract = 0.05     # fraction of sheet length to graph length
sheet_L = L*fract

# set up the max, total height of stack (along arrow)
s_L = fract*L

''' define stacks at all positions depedning on rounded magnitude'''

# find the maximum magnitude for scaling
max_size = np.max(mag)   # careful with fields that inc. singularities --> nan

# find the relative magnitude of vectors to maximum, as an array
R = mag/max_size

# define the maximum number of stack to plot, dep. on magnitude
s_max = 20

# define tigonometirc shifts
I_sin = np.sin(theta)
I_cos = np.cos(theta)

# define the points that set out a line of the stack sheet (middle line)
A_x = xg + (sheet_L/2)*I_sin
A_y = yg - (sheet_L/2)*I_cos
B_x = xg - (sheet_L/2)*I_sin
B_y = yg + (sheet_L/2)*I_cos

# Arrow points
p_sh1x = xg + (s_L/2)*I_cos + (sheet_L/8)*I_sin
p_sh1y = yg + (s_L/2)*I_sin - (sheet_L/8)*I_cos
p_sh2x = xg + (s_L/2)*I_cos - (sheet_L/8)*I_sin
p_sh2y = yg + (s_L/2)*I_sin + (sheet_L/8)*I_cos
p_sh3x = xg + (3*s_L/4)*I_cos
p_sh3y = yg + (3*s_L/4)*I_sin


# Arrow points for single sheet scenario
P_sh1x = xg + (sheet_L/8)*I_sin
P_sh1y = yg - (sheet_L/8)*I_cos
P_sh2x = xg - (sheet_L/8)*I_sin
P_sh2y = yg + (sheet_L/8)*I_cos
P_sh3x = xg + (s_L/4)*I_cos
P_sh3y = yg + (s_L/4)*I_sin


# define function that sets the recursion constant for the loop to plot stacks
# pre-define the displacements from mid point needed
# c is the Truth value from parity (odd or even number n)
def G(s, n, c):
    if c == 0:
        return ((2*s + 1)/(2*(n-1)))
    else:
        return (s/(n-1))


# define an empty array of magnitudes, to then fill with integer rel. mags
R_int = np.zeros(shape=(len(x), len(y)))

# loop over each arrow coordinate in x and y
for i in range(len(x)):
    for j in range(len(y)):
        # Label each element with the number of stacks required: linear scaling
        for t in range(1, s_max+1):
            if (t-1)/s_max <= R[i, j] <= t/s_max:
                R_int[i, j] = t
        
        # set a varible for current considered magnitude as it is reused
        # avoids extracting from R many times.
        n = R_int[i, j]
        
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
                ax.add_line(Line2D((Ax1, Bx1), (Ay1, By1), linewidth=1, color='green'))
                ax.add_line(Line2D((Ax2, Bx2), (Ay2, By2), linewidth=1, color='green'))
                
                # update parameter to reapet and draw all needed arrows
                s += 1
        # deal with the odd number of stacks:
        elif parity(n) is False:
            # Add the centre line for odd numbers of stacks
            ax.add_line(Line2D((A_x[i, j], B_x[i, j]), (A_y[i, j], B_y[i, j]), linewidth=1, color='green'))
            
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
                ax.add_line(Line2D((Ax1,Bx1),(Ay1,By1), linewidth=1, color='green'))
                ax.add_line(Line2D((Ax2,Bx2),(Ay2,By2), linewidth=1, color='green'))
                
                # change the parameter to loop over all changes in displacement for current magnitude
                s += 1
        # plot lines of arrowheads from central sheet for n = 1 or on top sheet for n>1 
        if n > 1:   # for all lines ubt the single sheet one
            ax.add_line(Line2D((p_sh1x[i, j],p_sh3x[i, j]),(p_sh1y[i, j],p_sh3y[i, j]), linewidth=1, color='green'))
            ax.add_line(Line2D((p_sh2x[i, j],p_sh3x[i, j]),((p_sh2y[i, j],p_sh3y[i, j])), linewidth=1, color='green'))
        # then define it for the stacks with only 1 sheet:
        else:
            ax.add_line(Line2D((P_sh1x[i, j], P_sh3x[i, j]), (P_sh1y[i, j], P_sh3y[i, j]), linewidth=1, color='green'))
            ax.add_line(Line2D((P_sh2x[i, j], P_sh3x[i, j]), ((P_sh2y[i, j], P_sh3y[i, j])), linewidth=1, color='green'))

# return time to run
stop = timeit.default_timer()
print('Time: ', stop - start)

# %%
