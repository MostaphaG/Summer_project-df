# Visualising 2-forms on R^2

# import modules
import timeit
import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, diff
from sympy.parsing.sympy_parser import parse_expr
from matplotlib.lines import Line2D
from matplotlib import cm
from matplotlib import patches as patch

# %%

# start the timer
start = timeit.default_timer()

'''

Find the mathematical expression of the 2-form, for correct magnitudes

'''



# define a function to replace input string to be 'python understood'
def format_eq(string):
    # replace all the x and y with xg and yg:
    string = string.replace('x', 'xg')
    string = string.replace('y', 'yg')
    string = string.replace('z', 'zg')
    # where there are special functions, replace them with library directions
    string = string.replace('pi', 'np.pi')
    string = string.replace('sqrt', 'np.sqrt')
    string = string.replace('sin', 'np.sin')
    string = string.replace('cos', 'np.cos')
    string = string.replace('tan', 'np.tan')
    string = string.replace('^', '**')
    string = string.replace('ln', 'np.log')        
    return string


# define a function that takes input string that is python understood and turn into vector components:
def eq_to_comps(string_x, string_y, xg, yg, u, v):
    global equation_x, equation_y
    # use this fucntion to replace given string to python understood equations:
    equation_x = format_eq(string_x)
    equation_y = format_eq(string_y)
    # use these to define the field:
    # also: check if equation equals zero, to then replace it with an array and not just 0:
    if equation_x == '0':
        u = np.zeros(np.shape(xg))
        v = eval(equation_y)
    elif equation_y == '0':
        u = eval(equation_x)
        v = np.zeros(np.shape(yg))
    elif equation_x and equation_y == '0':
        u = np.zeros(np.shape(xg))
        v = np.zeros(np.shape(yg))
    else:
        u = eval(equation_x)
        v = eval(equation_y)
    # return these
    return u, v


# define scale of the graph
L = 5
pt_den = 26   # number of points on each axis
a = 0.05  # linear scaling factor

# define x and y values
x = np.linspace(-L, L, pt_den)
y = np.linspace(-L, L, pt_den)

# create the grids
xg, yg = np.meshgrid(x, y)

# define an example vector field, now - from string, even initially
string_x = 'x*sin(y)'  # x component
string_y = 'y*cos(x)'  # y component

# to define a 2 from, need to perform the exterior derrivative on
# the given 1 form (vector field).
# To do this, seed to define partial derrivatives.
# to make it general, need to define it through derrivatives w.r.t. not
# already present coordinates.

# set the dimensionality
m = 2

# take the input strings and turn them into sympy expressions to be able to
# use sympy's partial differentiation
sympy_expr_x = parse_expr(string_x, evaluate=False)
sympy_expr_y = parse_expr(string_y, evaluate=False)
# for m > 2, need more components, and need these in 'expressions' too!

# define a sympy expression for string 0
sympy_expr_zero = parse_expr('0*x', evaluate=False)

# combine the 2 into a list:
expressions = np.array([sympy_expr_x, sympy_expr_y])

# set up an array to store derrivatives.
ext_ds = np.empty((m, m), dtype='object')

# use sympy partial derrivatives on these, as to get a 2-form on R2:
# need to differentiate each component w.r.t the coordinates that it's
# elementary 1 form does not contain.

# set up an array of coordinates that need to be used (in standard order)
coords = ['x', 'y']

# set up an array to store the results
# in 2D only dx^dy, in 3D (in order): dx^dy, dx^dz, dy^dz
result = np.empty((int((m-1)*m/2), 1), dtype='object')
for i in range(int((m-1)*m/2)):
    result[i] = str(result[i])

# loop over differentiating each, when differentiating w.r.t its coord, set to 0
for coord_index in range(len(coords)):
    # loop over differentiating each component:
    for comp_index in range(len(expressions)):
        # when equal set to 0, when not-differentiate:
        if comp_index == coord_index:
            ext_ds[comp_index, coord_index] = str(sympy_expr_zero)
        elif comp_index != coord_index:
            ext_ds[comp_index, coord_index] = str(diff(expressions[comp_index], coords[coord_index]))
        # change the signs for wedges in wrong order
        if comp_index < coord_index:
            ext_ds[comp_index, coord_index] = ' - (' + str(ext_ds[comp_index, coord_index]) + ')'
        elif comp_index > coord_index:
            ext_ds[comp_index, coord_index] = ' + ' + str(ext_ds[comp_index, coord_index])

'''
merge the results into a 2 form (for 2-form on R^2, the result is a single component (dx^xy))
do so by adding opposite elements along the diagonal ( / ) components of ext_ds
this  includes taking elemets with switched i and j
'''

# set up a variable to count pairs (pairs because we are forming 2 forms):
pair = 0

# loop over opposing elements (matching elementary 2-forms)
for i in range(1, m):
    for j in range(i):
        # initially clear the element from its Nonetype placeholder
        result[pair, 0] = ''
        # extract opposing elements
        temp = ext_ds[i, j]
        temp1 = ext_ds[j, i]
        # check these against zero entries:
        if (temp == '0') or (temp == '-(0)') or (temp == '0*x'):
            pass
        else:
            result[pair, 0] += temp
        if (temp1 == '0') or (temp1 == '-(0)') or (temp1 == '0*x'):
            pass
        else:
            result[pair, 0] += temp1
        # update the result row counter
        pair += 1


# return them for user to inspect, note- objects thereofre cant do so in var. explorer.
print('\n')
print(' the resulting 2 form is :')
print(result)
print('\n')
print('the matrix of derrivatives was:')
print(ext_ds)
print('\n')

# format string in each result row
for d in range(pair):
    # format the result to be 'python understood' to be able to use the eval()
    result[d, 0] = format_eq(result[d, 0])

# set up a vector to store the 2 form numerically, from xg and yg
form_2 = np.empty((len(result[:, 0]), pt_den, pt_den))    # Note - need pt_den m times.

# evaluate the expressions again:
for d in range(0, len(result[:, 0])):
    # numerical evaluation of the 2 form on R^2, at all defined grid points in R^2
    form_2[d, :, :] = eval(result[d, 0])  # Note, need : m times


'''

Define a stack vector field and wedge them

'''

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


# define function that sets the recursion constant for the loop to plot stacks
# pre-define the displacements from mid point needed
# c is the Truth value from parity (odd or even number n)
def G(s, n, c):
    if c == 0:
        return ((2*s + 1)/(2*(n-1)))
    else:
        return (s/(n-1))


# define a function that will complete all stack plotting:
def stack_plot(xg, yg, ax, u, v, s_max, L, pt_den, fract, arrows='True', orientation='mid', scale=1):
    global s_L
    # get axis lengths:
    x_len = len(xg[:, 0])
    y_len = len(yg[:, 0])
    
    # Scaling of axes and setting equal proportions circles look like circles
    ax.set_aspect('equal')
    ax_L = L + L/delta_factor
    ax.set_xlim(-ax_L, ax_L)
    ax.set_ylim(-ax_L, ax_L)
    
    # define an empty array of magnitudes, to then fill with integer rel. mags
    R_int = np.zeros(shape=((x_len), (y_len)))
    
    # #########################################################################
    # plot the initial quiver plot to work from
    # #########################################################################
    
    # plot the quiver plot on grid points if chosen in original function
    if arrows is True:
        ax.quiver(xg, yg, u, v, pivot=orientation, scale=scale, scale_units='xy')
    else:
        pass
    
    # #########################################################################
    # get variables needed for the initial, simplified stack plot
    # #########################################################################
    # find the arrow length corresponding to each point and store in mag array
    mag = np.sqrt(u**2 + v**2)
    # find direction of each arrow
    theta = np.arctan2(v, u)   # theta defined from positive x axis ccw
    
    # #########################################################################
    # use the the direction of arrows to define stack properties
    # #########################################################################
    
    # define length of sheet as a fraction of total graph scale
    sheet_L = L*fract
    # set up the max, total height of stack (along arrow)
    s_L = fract*L
    
    # #########################################################################
    # define the stacks based on geometrical arguments
    # to be perp. to arrow. shifted parallel to it, their density porp to mag
    # of the arrow and with an arrowhead on top.
    # #########################################################################
    # find the maximum magnitude for scaling
    max_size = np.max(mag)   # careful with singularities, else ---> nan
    
    # find the relative magnitude of vectors to maximum, as an array
    R = mag/max_size
    
    # define tigonometirc shifts
    I_sin = np.sin(theta)
    I_cos = np.cos(theta)
    
    # define the points that set out a line of the stack sheet (middle line)
    A_x = xg + (sheet_L/2)*np.sin(theta)
    A_y = yg - (sheet_L/2)*np.cos(theta)
    B_x = xg - (sheet_L/2)*np.sin(theta)
    B_y = yg + (sheet_L/2)*np.cos(theta)
    
    # loop over each arrow coordinate in x and y
    for i in range(x_len):
        for j in range(y_len):
            # define it for all magnitudes. Separately for odd and even corr. number of sheets:
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
                    ax.add_line(Line2D((Ax1, Bx1), (Ay1, By1), linewidth=0.5, color='green'))
                    ax.add_line(Line2D((Ax2, Bx2), (Ay2, By2), linewidth=0.7, color='green'))
                    
                    # update parameter to reapet and draw all needed arrows
                    s += 1
            # deal with the odd number of stacks:
            elif parity(n) is False:
                # Add the centre line for odd numbers of stacks
                ax.add_line(Line2D((A_x[i, j], B_x[i, j]), (A_y[i, j], B_y[i, j]), linewidth=0.7, color='green'))
                
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
                    ax.add_line(Line2D((Ax1,Bx1),(Ay1,By1), linewidth=0.7, color='green'))
                    ax.add_line(Line2D((Ax2,Bx2),(Ay2,By2), linewidth=0.7, color='green'))
                    
                    # change the parameter to loop over all changes in displacement for current magnitude
                    s += 1


# evaluate the u and v given previously, with formating for them to be
# python understood
u = a* eval(format_eq(str(expressions[0])))
v = a* eval(format_eq(str(expressions[1])))

# set up a zero vector filed to plot x and y components as 2 separate fields:
zero_field = np.zeros(np.shape(xg))

# set up the delta_factor of additional axis space L/delta_factor gives extra space on axis
delta_factor = 10

# fraction of sheet length to graph length
fract = 0.05

# define the maximum number of stack to plot, dep. on magnitude
s_max = 10
    
# create a figure
fig = plt.figure(figsize=(8, 8))

# set up axis
ax = fig.gca()

# set up visuals - axis labels
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')

# Scaling of axes and setting equal proportions circles look like circles
ax.set_aspect('equal')
ax_L = L + L/delta_factor
ax.set_xlim(-ax_L, ax_L)
ax.set_ylim(-ax_L, ax_L)

# plot the fields with desired parameters as specidfied above
# arrow params not needed as arrows arent plotted
stack_plot(xg, yg, ax, u, zero_field, s_max, L, pt_den, fract, False, 'mid', 1)
stack_plot(xg, yg, ax, zero_field, v, s_max, L, pt_den, fract, False, 'mid', 1)
# these fields are the original components

'''

using the magnitude of the 2-form component (dx^dy), define
direction by finding the total 2-form over that area and checkiong signs

'''

# integration is not needed, as grid is evenly spaced
# as form_2 is defined in terms of xg and yg that are the same as the ones used
# when defining stacks:
# determine the sign of the 2 form at each grid point
# based on put a coloured region inside the stack ; blue (-ve) or red (+ive)

# change the 2-form array to 0s 1s and (-1)s:
form_2_sgn = np.sign(form_2[0])  # as it is a 2_form on R^2, it only has one component, only select that one

# loop over all grid points
for i in range(len(xg[:, 0])):
    for j in range(len(yg[0, :])):
        # depending on sign, colour the area of the stack
        if form_2_sgn[i, j] == -1:
            rect1 = patch.Rectangle((xg[i, j] - s_L/4, yg[i, j] - s_L/4), s_L/2, s_L/2, color='blue')
            ax.add_patch(rect1)
        elif form_2_sgn[i, j] == 1:
            rect2 = patch.Rectangle((xg[i, j] - s_L/4, yg[i, j] - s_L/4), s_L/2, s_L/2, color='red')
            ax.add_patch(rect2)
        else:  # for when it is exacly equal to zero (just in case)
            print('ZERO!')




# can add a background colourmap to see between gridpoints by the following:
'''
cp = ax.contourf(xg, yg, form_2[0], cmap='bwr')
fig.colorbar(cp)
'''

# return time to run
stop = timeit.default_timer()
print('Time to run: ', stop - start)
