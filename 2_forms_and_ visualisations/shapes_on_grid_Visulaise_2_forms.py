# Visualising 2-forms through shapes at grid points

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
L = 8
pt_den = 31   # number of points on each axis
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

# evaluate the expressions for each elemental 2-form:
for d in range(0, len(result[:, 0])):
    # numerical evaluation of the 2 form on R^2, at all defined grid points in R^2
    form_2[d, :, :] = eval(result[d, 0])  # Note, need : m times

'''

Use these to plot a square, upto some maximum size at each grid point

'''

# find the maximum of the 2-form, over the grid
max_val2 = np.max(form_2[0, :, :])

# get an array of relative magnitudes:
Mag = form_2[0, :, :]/max_val2

# make Mag unitary:
Mag = Mag/np.max(Mag)

# set a maximum side size as a fraction of the graph size
fract_s = 0.1
max_s = fract_s*L

# from the maximum size, set the Mag array to store shape sizes
Mag = Mag*max_s

# extract signs from that (+1, -1, and 0)
sign_mag = np.sign(Mag)

# clear the negative Magnitudes to be able to use these as shape size:
Mag = abs(Mag)

# create a figure
fig = plt.figure(figsize=(8, 8))

# set up axis
ax = fig.gca()

# set up visuals - axis labels
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')

# Scaling of axes and setting equal proportions circles look like circles
ax.set_aspect('equal')
delta_factor = 10   # sets denominator which established how much extra space to add to axis, past where field is defined
ax_L = L + L/delta_factor
ax.set_xlim(-ax_L, ax_L)
ax.set_ylim(-ax_L, ax_L)


# define integer values of magnitude, linearly scaled
# save these as the Mag array.
for i in range(len(Mag[:, 0])):
    for j in range(len(Mag[0, :])):
        # Now, at every grid point, plot a square of size given by magnitude
        # use colours to indicate orientation - clockwise (negative ==> blue)
        # counterclockwise (positive ==> red)
        if sign_mag[i, j] == -1:
            rect = patch.Rectangle((xg[i, j] - Mag[i, j]/4, yg[i, j] - Mag[i, j]/4), Mag[i, j]/2, Mag[i, j]/2, color='blue')
            ax.add_patch(rect)
        elif sign_mag[i, j] == 1:
            rect = patch.Rectangle((xg[i, j] - Mag[i, j]/4, yg[i, j] - Mag[i, j]/4), Mag[i, j]/2, Mag[i, j]/2, color='red')
            ax.add_patch(rect)
        else:
            rect = patch.Rectangle((xg[i, j] - Mag[i, j]/4, yg[i, j] - Mag[i, j]/4), Mag[i, j]/2, Mag[i, j]/2, color='grey')
            ax.add_patch(rect)  # not needed as zero sign will have zero magnitude therefore will not be seen anyway.
            print('can\'t establish orientation for a zero magnitude 2 form')

