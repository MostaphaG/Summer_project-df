# Plotting 2 form on R3, from given 1 form on R3.
# Displays chosen plane, and allows for sliding across the viewing plane

# import modules
import timeit
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from sympy.parsing.sympy_parser import parse_expr
from sympy import diff
from matplotlib.lines import Line2D
from matplotlib.backend_bases import key_press_handler
from sympy import simplify
from matplotlib import patches as patch

# %%

# start the timer
start = timeit.default_timer()

# start a tkinter window
root = tk.Tk()

# set its title
root.title('2 forms on R3')

# set window size
root.geometry('1100x900')

# define a frame to put the plot into:
plot_frame = tk.LabelFrame(root, text='plot Frame', padx=5, pady=5)
plot_frame.grid(row=0, column=0)


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
    string = string.replace('log', 'np.log')
    return string


# function to unformat equation to display back to user in same form
def unformat(string):
    # replace all the x and y with xg and yg:
    string = string.replace('xg', 'x')
    string = string.replace('yg', 'y')
    string = string.replace('zg', 'z')
    # where there are special functions, replace them with library directions
    string = string.replace('np.pi', 'pi')
    string = string.replace('np.sqrt', 'sqrt')
    string = string.replace('np.sin', 'sin')
    string = string.replace('np.cos', 'cos')
    string = string.replace('np.tan', 'tan')
    string = string.replace('np.arctan2', 'ARCTAN')
    string = string.replace('ARCSIN', 'np.arcsin')
    string = string.replace('np.arccos', 'ARCCOS')
    string = string.replace('np.tanh', 'TANH')
    string = string.replace('np.sinh', 'SINH')
    string = string.replace('np.cosh', 'COSH')
    string = string.replace('**', '^')
    string = string.replace('np.log()', 'ln')
    string = string.replace('np.exp', 'e**')
    return string


# define a function that takes input string that is python understood and turn into vector components:
def eq_to_comps(string_x, string_y, xg, yg):
    global equation_x, equation_y
    # use the format_eq fucntion to make given string readable by python
    equation_x = format_eq(string_x)
    equation_y = format_eq(string_y)
    # use these to define the field:
    # also: checking if equation equals zero, to then replace it with an array and not just 0:
    u = eval(equation_x)
    v = eval(equation_y)
    if equation_x.find('x') & equation_x.find('y') == -1:
        u = eval(equation_x)*np.ones(np.shape(xg))
    if equation_y.find('x') & equation_y.find('y') == -1:
        v = eval(equation_y)*np.ones(np.shape(yg))
    # deal properly with zero and constant fields too:
    # check for when the derivative is zero, do not plot it as nothing
    # if the other component is not zero.
    if equation_x == '0' and equation_y != '0':
        u = np.ones(np.shape(xg))
    if equation_y == '0' and equation_x != '0':
        v = np.ones(np.shape(yg))
    # return these
    return u, v

# define a function that will find the 2 form from given expressions
# in a given number of dimensions and in terms of given coordinate symbols
def find_2_form(expressions, coords, pt_den, m=3):
    global ext_ds, result
    # define a sympy expression for string 0
    sympy_expr_zero = parse_expr('0*x', evaluate=False)
    
    # set up an array to store derrivatives.
    ext_ds = np.empty((m, m), dtype='object')
    
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
                ext_ds[comp_index, coord_index] = ' +' + str(ext_ds[comp_index, coord_index])
    
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
            # clear the used ext_ds
            ext_ds[i, j] = ''
            ext_ds[j, i] = ''
            # check these against zero entries:
            if (temp == '0') or (temp == ' - (0)') or (temp == '0*x'):
                ext_ds[i, j] = '0'
            else:
                result[pair, 0] += temp
                ext_ds[i, j] = temp
            if (temp1 == '0') or (temp1 == ' - (0)') or (temp1 == '0*x'):
                ext_ds[j, i] = '0'
            else:
                result[pair, 0] += temp1
                ext_ds[j, i] = temp1
            # update the result row counter
            pair += 1
    # format string in each result row
    # making sure to format it correctly even if it contains constants or '0'
    # this is done if result is to be directly used later for any reason
    # instead of form_2 numerically.
    for d in range(pair):
        if result[d, 0].find('x') & result[d, 0].find('y') & result[d, 0].find('z') == -1:
            if result[d, 0] == '':
                result[d, 0] = '0*x'  # no need for unit field, as it is * by 0
                result[d, 0] = format_eq(result[d, 0])
            else:
                result[d, 0] = '(' + result[d, 0] + ')* field_unit'
                result[d, 0] = format_eq(result[d, 0])
        else:
            # format the result to be 'python understood' to be able to use the eval()
            result[d, 0] = format_eq(result[d, 0])
    
    # set up a vector to store the 2 form numerically, from xg and yg
    form_2 = np.empty((len(result[:, 0]), pt_den, pt_den, pt_den))    # Note - need pt_den m times.
    
    # evaluate the expressions again:
    for d in range(0, len(result[:, 0])):
        # numerical evaluation of the 2 form on R^2, at all defined grid points in R^2
        form_2[d, :, :, :] = eval(result[d, 0])  # Note, need : m times
    return form_2


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


# define a function that will plot stack components, coloured
# as per the orientation of the 2 form at that grid point
def form_2_components_plot_3(grid_x, grid_y, h_index, axis_view, u, v, s_max, L, pt_den, fract, colour_str):
    global s_L
    
    # depending on axis_view and h_index, get the planar u and v from the given
    # 3D ones and the grids sorted
    if axis_view == 'z':
        grid_x = grid_x[:, :, h_index]
        grid_y = grid_y[:, :, h_index]
        u = u[:, :, h_index]
        v = v[:, :, h_index]
    elif axis_view == 'y':
        grid_x = grid_x[h_index, :, :]
        grid_y = grid_y[h_index, :, :]
        u = u[h_index, :, :]
        v = v[h_index, :, :]
    elif axis_view == 'x':
        grid_x = grid_x[:, h_index, :]
        grid_y = grid_y[:, h_index, :]
        u = u[:, h_index, :]
        v = v[:, h_index, :]
    else:
        print('Error can\'t find this axis')

    # get axis lengths:
    x_len = len(grid_x[:, 0])  # no need to change with axis_view
    y_len = len(grid_y[0, :])  # if grids are all same size
    
    # Scaling of axes and setting equal proportions circles look like circles
    ax.set_aspect('equal')
    ax_L = L + L/delta_factor
    ax.set_xlim(-ax_L, ax_L)
    ax.set_ylim(-ax_L, ax_L)
    
    # define an empty array of magnitudes, to then fill with integer rel. mags
    R_int = np.zeros(shape=((x_len), (y_len)))
    
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
    A_x = grid_x + (sheet_L/2)*np.sin(theta)
    A_y = grid_y - (sheet_L/2)*np.cos(theta)
    B_x = grid_x - (sheet_L/2)*np.sin(theta)
    B_y = grid_y + (sheet_L/2)*np.cos(theta)
    
    # get the 2 form signs and
    # change the 2_form signs to just be in terms of the selected plane
    if axis_view == 'x':
        form_2_sgn = np.sign(form_2[2])
        form_2_sgn_planar = form_2_sgn[h_index, :, :]
    elif axis_view == 'y':
        form_2_sgn = np.sign(form_2[1])
        form_2_sgn_planar = form_2_sgn[:, h_index, :]
    elif axis_view == 'z':
        form_2_sgn = np.sign(form_2[0])
        form_2_sgn_planar = form_2_sgn[:, :, h_index]
    
    # loop over each arrow coordinate in x and y
    for i in range(x_len):
        for j in range(y_len):
            # define it for all magnitudes. Separately for odd and even corr. number of sheets:
            # Label each element with the number of stacks required: linear scaling
            
            if form_2_sgn_planar[i, j] == +1:
                color_index = 0
            elif form_2_sgn_planar[i, j] == -1:
                color_index = 1
            else:
                color_index = 2  # in case it is zero exactly
            
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
                    ax.add_line(Line2D((Ax1, Bx1), (Ay1, By1), linewidth=0.5, color=colour_str[color_index]))
                    ax.add_line(Line2D((Ax2, Bx2), (Ay2, By2), linewidth=0.7, color=colour_str[color_index]))
                    
                    # update parameter to reapet and draw all needed arrows
                    s += 1
            # deal with the odd number of stacks:
            elif parity(n) is False:
                # Add the centre line for odd numbers of stacks
                ax.add_line(Line2D((A_x[i, j], B_x[i, j]), (A_y[i, j], B_y[i, j]), linewidth=0.7, color=colour_str[color_index]))
                
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
                 
                    ax.add_line(Line2D((Ax1,Bx1),(Ay1,By1), linewidth=0.7, color=colour_str[color_index]))
                    ax.add_line(Line2D((Ax2,Bx2),(Ay2,By2), linewidth=0.7, color=colour_str[color_index]))
                    
                    # change the parameter to loop over all changes in displacement for current magnitude
                    s += 1
    plt.close()


# define a function to plot the simplified 2 forms, with coloured squares
def plot_form(form_2, grid_x, grid_y, fract_s):
    global Mag
    # redefine the axis limits
    ax.set_xlim(-ax_L, ax_L)
    ax.set_ylim(-ax_L, ax_L)
    # crop xg, yg and zg to size:
    if axis_view == 'z':
        grid_x = grid_x[:, :, h_index]
        grid_y = grid_y[:, :, h_index]
    elif axis_view == 'y':
        grid_x = grid_x[h_index, :, :]
        grid_y = grid_y[h_index, :, :]
    elif axis_view == 'x':
        grid_x = grid_x[:, h_index, :]
        grid_y = grid_y[:, h_index, :]
    # find the maximum of the 2-form, over the grid
    max_val2 = np.max(form_2)
    # get an array of relative magnitudes:
    Mag = form_2/abs(max_val2)
    # make Mag unitary:
    Mag = Mag/abs(np.max(Mag))
    # set a maximum side size as a fraction of the graph size
    max_s = fract_s*L
    # from the maximum size, set the Mag array to store shape sizes
    Mag = Mag*max_s
    # extract signs from that (+1, -1, and 0)
    sign_mag = np.sign(Mag)
    # clear the negative Magnitudes to be able to use these as shape size:
    Mag = abs(Mag)
    # define integer values of magnitude, linearly scaled
    # save these as the Mag array.
    for i in range(len(Mag[:, 0])):
        for j in range(len(Mag[0, :])):
            # Now, at every grid point, plot a square of size given by magnitude
            # use colours to indicate orientation - clockwise (negative ==> blue)
            # counterclockwise (positive ==> red)
            if sign_mag[i, j] == -1:
                rect = patch.Rectangle((grid_x[i, j] - Mag[i, j]/4, grid_y[i, j] - Mag[i, j]/4), Mag[i, j]/2, Mag[i, j]/2, color='blue')
                ax.add_patch(rect)
            elif sign_mag[i, j] == 1:
                rect = patch.Rectangle((grid_x[i, j] - Mag[i, j]/4, grid_y[i, j] - Mag[i, j]/4), Mag[i, j]/2, Mag[i, j]/2, color='red')
                ax.add_patch(rect)
            else:
                rect = patch.Rectangle((grid_x[i, j] - Mag[i, j]/4, grid_y[i, j] - Mag[i, j]/4), Mag[i, j]/2, Mag[i, j]/2, color='grey')
                ax.add_patch(rect)  # not needed as zero sign will have zero magnitude therefore will not be seen anyway.


'''

Set up all needed parameters, plots etc

'''

# define parameters of the axis
L = 5
pt_den = 21  # number of points on each axis

# define x, y and z values
x = np.linspace(-L, L, pt_den)
y = np.linspace(-L, L, pt_den)
z = np.linspace(-L, L, pt_den)

# create the grids
xg, yg, zg = np.meshgrid(x, y, z)

# to start with, set as viewing aling z axis onto x-y plane
axis_view = 'z'

# set the initial index of height of the viewed plane along viewing axis
h_index = 11

# set the dimensionality
m = 3

# set up a zero vector filed to plot x and y components as 2 separate fields:
zero_field = np.zeros(np.shape(xg))

# define an array of ones of the correct size in case result in needed with correct magnitude
field_unit = np.ones(np.shape(zg))

# set up the delta_factor of additional axis space L/delta_factor gives extra space on axis
delta_factor = 10

# fraction of sheet length to graph length
fract = 0.05

# same for block side size
fract_s = 0.15

# define the maximum number of stack to plot, dep. on magnitude
s_max = 5

# create a figure
fig = plt.figure(figsize=(8, 6))

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

# set up a colour string for the orientations
# first string defines colour for positive (ccw), second for negative (cw)
# last one is an in case, for when the magnitude is exactly zero.
colour_str = ['red', 'blue', 'grey']

# predefine strings for splitting in x and in y for 2 form stacks
xy_x_split_str = '1/2'
xy_y_split_str = '1/2'
xz_x_split_str = '1/2'
xz_y_split_str = '1/2'
yz_x_split_str = '1/2'
yz_y_split_str = '1/2'

# define the wanted 1 form on R3 in terms of each component:
string_x = 'x*y*z'  # x component
string_y = 'y*z*x'  # y component
string_z = 'z*x*y'  # z component

# INITIAL CALCULATIONS

# take the input strings and turn them into sympy expressions to be able to
# use sympy's partial differentiation
sympy_expr_x = parse_expr(string_x, evaluate=False)
sympy_expr_y = parse_expr(string_y, evaluate=False)
sympy_expr_z = parse_expr(string_z, evaluate=False)
# for m > 2, need more components, and need these in 'expressions' too!

# combine the 2 into a list:
expressions = np.array([sympy_expr_x, sympy_expr_y, sympy_expr_z])

# use sympy partial derrivatives on these, as to get a 2-form on R2:
# need to differentiate each component w.r.t the coordinates that it's
# elementary 1 form does not contain.

# set up an array of coordinates that need to be used (in standard order)
coords = ['x', 'y', 'z']

# from these, use the find_2_form function to get the 2 form
form_2 = find_2_form(expressions, coords, pt_den, m)

# get it's string
form_2_str_dxdy = str(simplify(str(unformat(result[0][0]))))
form_2_str_dxdz = str(simplify(str(unformat(result[1][0]))))
form_2_str_dydz = str(simplify(str(unformat(result[2][0]))))

# because on R3, the x, y and z components are made up of
# dy^dz, dx^dz and dx^dy
# need to superpose stacks from the x y and z fields (depending on viewing)
# the filed will point in the extra dimension, or by polar sense
# will be defined in plane based on orientation of stacks. Indicated by colour

# note that depending on the viewing axis, any result stacks
# cannot exist in that direction, as that will close the area off in the viewed
# plane
# the plane are done as before, by splitting and so aren't unique
# blocks are!
# therefore F_x. F_y and F_z need to be defined depending on viewed plane
# find all separately, to make switching axis and sliding faster

# define components in the x-y plane:
eq_1_xy = str(simplify('(' + form_2_str_dxdy + ')' + '*(' + xy_x_split_str + ')'))
eq_2_xy = str(simplify('(' + form_2_str_dxdy + ')' + '*(' + xy_y_split_str + ')'))
F_xy_x, F_xy_y = eq_to_comps(eq_1_xy, eq_2_xy, xg, yg)

# define them in x-z plane
eq_1_xz = str(simplify('(' + form_2_str_dxdz + ')' + '*(' + xz_x_split_str + ')'))
eq_2_xz = str(simplify('(' + form_2_str_dxdz + ')' + '*(' + xz_y_split_str + ')'))
F_xz_x, F_xz_z = eq_to_comps(eq_1_xz, eq_2_xz, xg, yg)


# define them in y-z plane
eq_1_yz = str(simplify('(' + form_2_str_dydz + ')' + '*(' + yz_x_split_str + ')'))
eq_2_yz = str(simplify('(' + form_2_str_dydz + ')' + '*(' + yz_y_split_str + ')'))
F_yz_y, F_yz_z = eq_to_comps(eq_1_yz, eq_2_yz, xg, yg)


# plot the starting field with desired parameters as specidfied above
# arrow params not needed as arrows arent plotted
# starting with viewing axis ='z' therefore xy plane
form_2_components_plot_3(xg, yg, h_index, axis_view, F_xy_x, zero_field, s_max, L, pt_den, fract, colour_str)
form_2_components_plot_3(xg, yg, h_index, axis_view, zero_field, F_xy_y, s_max, L, pt_den, fract, colour_str)

# reduce white space from the figure in the plot frame
fig.tight_layout()

# set up the space for the plot to be put into AKA the plot frame
canvas = FigureCanvasTkAgg(fig, master=plot_frame)
canvas.draw()
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

# put the default matplotlib toolbar, below the figure
toolbar = NavigationToolbar2Tk(canvas, plot_frame)
toolbar.update()  # allow the plot to update based on the toolbar
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)


# track the mouse presses for the toolbar to respond to
def on_key_press(event):
    print("you pressed {}".format(event.key))
    key_press_handler(event, canvas, toolbar)


# define a function to respond to new 1 forms being input
def form_1_to_2onR3():
    global string_x, string_y, string_z, form_2, F_xy_x, F_xy_y, F_xz_x, F_xz_z, F_yz_y, F_yz_z
    global eq_1_xy, eq_2_xy, eq_1_xz, eq_2_xz, eq_1_yz, eq_2_yz
    # take the inputs from user into strings
    string_x = str(simplify(form_1_x_entry.get()))
    string_y = str(simplify(form_1_y_entry.get()))
    string_z = str(simplify(form_1_z_entry.get()))
    # turn to sumpy expressions
    sympy_expr_x = parse_expr(string_x, evaluate=False)
    sympy_expr_y = parse_expr(string_y, evaluate=False)
    sympy_expr_z = parse_expr(string_z, evaluate=False)
    # combine the 2 into an array:
    expressions = np.array([sympy_expr_x, sympy_expr_y, sympy_expr_z])
    # ALL AS BEFORE:
    form_2 = find_2_form(expressions, coords, pt_den, m)
    form_2_str_dxdy = str(simplify(str(unformat(result[0][0]))))
    form_2_str_dxdz = str(simplify(str(unformat(result[1][0]))))
    form_2_str_dydz = str(simplify(str(unformat(result[2][0]))))
    eq_1_xy = str(simplify('(' + form_2_str_dxdy + ')' + '*(' + xy_x_split_str + ')'))
    eq_2_xy = str(simplify('(' + form_2_str_dxdy + ')' + '*(' + xy_y_split_str + ')'))
    F_xy_x, F_xy_y = eq_to_comps(eq_1_xy, eq_2_xy, xg, yg)
    eq_1_xz = str(simplify('(' + form_2_str_dxdz + ')' + '*(' + xz_x_split_str + ')'))
    eq_2_xz = str(simplify('(' + form_2_str_dxdz + ')' + '*(' + xz_y_split_str + ')'))
    F_xz_x, F_xz_z = eq_to_comps(eq_1_xz, eq_2_xz, xg, yg)
    eq_1_yz = str(simplify('(' + form_2_str_dydz + ')' + '*(' + yz_x_split_str + ')'))
    eq_2_yz = str(simplify('(' + form_2_str_dydz + ')' + '*(' + yz_y_split_str + ')'))
    F_yz_y, F_yz_z = eq_to_comps(eq_1_yz, eq_2_yz, xg, yg)
    # call the slide function to easily plot them depending on h and on axis_view
    slide()
    # display this new 2 form in its frame as labels
    Label_2_form_xy.configure(text=str(simplify(unformat(result[0][0]))) + '  dx^dy')
    Label_2_form_xz.configure(text=str(simplify(unformat(result[1][0]))) + '  dx^dz')
    Label_2_form_yz.configure(text=str(simplify(unformat(result[2][0]))) + '  dy^dz')


# define a function to update z Index and redraw the plot based on the slider
def slide():
    # remove the currently displayed plot
    ax.clear()
    # replot the graph with that new h_index
    # and change the label under the slider to be the value of the chosen axis at that h_index
    if stack_block_int == 1:
            # 1 was chosen therefore complete stacks
        if axis_view == 'z':
            form_2_components_plot_3(xg, yg, h_index, axis_view, F_xy_x, zero_field, s_max, L, pt_den, fract, colour_str)
            form_2_components_plot_3(xg, yg, h_index, axis_view, zero_field, F_xy_y, s_max, L, pt_den, fract, colour_str)
            ax.set_xlabel('$x$')
            ax.set_ylabel('$y$')
            # update the label based on that
            axis_height_txt.configure(text=str(z[h_index]))
        elif axis_view == 'y':
            form_2_components_plot_3(xg, zg, h_index, axis_view, F_xz_x, zero_field, s_max, L, pt_den, fract, colour_str)
            form_2_components_plot_3(xg, zg, h_index, axis_view, zero_field, F_xz_z, s_max, L, pt_den, fract, colour_str)
            ax.set_xlabel('$x$')
            ax.set_ylabel('$z$')
            # update the label based on that
            axis_height_txt.configure(text=str(y[h_index]))
        elif axis_view == 'x':
            form_2_components_plot_3(yg, zg, h_index, axis_view, F_yz_y, zero_field, s_max, L, pt_den, fract, colour_str)
            form_2_components_plot_3(yg, zg, h_index, axis_view, zero_field, F_yz_z, s_max, L, pt_den, fract, colour_str)
            ax.set_xlabel('$y$')
            ax.set_ylabel('$z$')
            # update the label based on that
            axis_height_txt.configure(text=str(x[h_index]))
    elif stack_block_int == 0:
        if axis_view == 'z':
            plot_form(form_2[0][:, :, h_index], xg, yg, fract_s)
        elif axis_view == 'y':
            plot_form(form_2[1][:, h_index, :], xg, zg, fract_s)
        elif axis_view == 'x':
            plot_form(form_2[2][:, :, h_index], yg, zg, fract_s)
    # draw that onto the screen
    canvas.draw()


# deifne a function that will update the label
def label_update(var):
    global h_index
    # update current height
    h_index += var
    # update the label based on that
    if axis_view == 'z':
        # update the label based on that
        axis_height_txt.configure(text=str(z[h_index]))
    elif axis_view == 'y':
        # update the label based on that
        axis_height_txt.configure(text=str(y[h_index]))
    elif axis_view == 'x':
        # update the label based on that
        axis_height_txt.configure(text=str(x[h_index]))


# define a function that will repond to changing axis view with radiobuttons
def view_response(view_var):
    # get and make global the axis view variable
    global axis_view
    axis_view = view_tk.get()
    # call slide function to plot based on the wanted axis view and height.
    slide()
    # draw that onto the screen
    canvas.draw()


# define a function to decide if stacks or blocks are being plotted
def plot_type_responseR3(type_stack_block):
    global stack_block_int
    stack_block_int = int(type_stack_block)


# define a label and arrows to change the axis value, and index
# instead of the previously tried slider and entry boxes:

# NOTE NOT GREAT I KNOW BUT TEMPORARY:
# define a frame for the buttons moving the view_axis value
height_frame = tk.LabelFrame(root, text='Fields frame', padx=32, pady=5)
height_frame.grid(row=1, column=0)

# Label to show current axis value
axis_height_txt = tk.Label(height_frame, text=str(z[11]))
axis_height_txt.grid(row=1, column=0)

# on the left, make a 'move down' button
down_height = tk.Button(height_frame, text=' \/ ', command=lambda: label_update(-1))
down_height.grid(row=2, column=0)

# on the right, make a 'move up' button
up_height = tk.Button(height_frame, text=' /\ ', command=lambda: label_update(1))
up_height.grid(row=0, column=0)

# define a button to submit the currently chosen value:
Submit_h_btn = tk.Button(height_frame, text='SUBMIT', padx=10, pady=50, command=slide)
Submit_h_btn.grid(row=0, column=1, rowspan=3, padx=20)


# define rediobuttons to chose from which axis the user is looking:
view_tk = tk.StringVar()
view_tk.set('z')
view_z_btn = tk.Radiobutton(height_frame, text='z', variable=view_tk, value='z', command=lambda: view_response(view_tk.get())).grid(row=0, column=2)
view_z_btn = tk.Radiobutton(height_frame, text='y', variable=view_tk, value='y', command=lambda: view_response(view_tk.get())).grid(row=1, column=2)
view_z_btn = tk.Radiobutton(height_frame, text='x', variable=view_tk, value='x', command=lambda: view_response(view_tk.get())).grid(row=2, column=2)


# NOTE  NOT GREAT I KNOW BUT TEMPORARY:
# define a new frame for the fields to be input
field_input_frame = tk.LabelFrame(root, text='Fields frame', padx=32, pady=5)
field_input_frame.grid(row=0, column=1)

# define radiobuttons to pick stacks or blocks
stack_block = tk.IntVar()  # Tkinter variable for Radiobuttons
stack_block.set(1)
stack_block_int = stack_block.get()
blocks_rb = tk.Radiobutton(field_input_frame, text='blocks', variable=stack_block, value=0, command=lambda: plot_type_responseR3(stack_block.get()))
blocks_rb.grid(row=0, column=1)
stacks_rb = tk.Radiobutton(field_input_frame, text='stacks', variable=stack_block, value=1, command=lambda: plot_type_responseR3(stack_block.get()))
stacks_rb.grid(row=0, column=2)

# define entry boxes for the three 1 forms that are being plotted
# define entries for a 1 form
tk.Label(field_input_frame, text='x component of 1 form').grid(row=1, column=1)
form_1_x_entry = tk.Entry(field_input_frame, width=20, borderwidth=2)
form_1_x_entry.grid(row=2, column=1, columnspan=2)
form_1_x_entry.insert(0, string_x)

tk.Label(field_input_frame, text='y component of 1 form').grid(row=3, column=1)
form_1_y_entry = tk.Entry(field_input_frame, width=20, borderwidth=2)
form_1_y_entry.grid(row=4, column=1, columnspan=2)
form_1_y_entry.insert(0, string_y)

tk.Label(field_input_frame, text='z component of 1 form').grid(row=5, column=1)
form_1_z_entry = tk.Entry(field_input_frame, width=20, borderwidth=2)
form_1_z_entry.grid(row=6, column=1, columnspan=2)
form_1_z_entry.insert(0, string_z)

# deifne a button to plot from these:
form_2_from_1_R3_btn = tk.Button(field_input_frame, text='2 form from 1 forms R3', padx=3, pady=5, command=form_1_to_2onR3)
form_2_from_1_R3_btn.grid(row=7, column=1, columnspan=2)

# define a frame where the resulting 2 form can display
form_2_frame = tk.LabelFrame(root, text='2 form frame', padx=32, pady=5)
form_2_frame.grid(row=1, column=1)

# in it, FOR NOW:
# put the labels of the 2 form components
# for the 2 form that is being plotted.
Label_2_form_xy = tk.Label(form_2_frame, text=str(simplify(unformat(result[0][0]))) + '  dx^dy')
Label_2_form_xy.grid(row=0, column=0)

Label_2_form_xz = tk.Label(form_2_frame, text=str(simplify(unformat(result[1][0]))) + '  dx^dz')
Label_2_form_xz.grid(row=1, column=0)

Label_2_form_yz = tk.Label(form_2_frame, text=str(simplify(unformat(result[2][0]))) + '  dy^dz')
Label_2_form_yz.grid(row=2, column=0)


# return time to run
stop = timeit.default_timer()
print('Time to run: ', stop - start)

# keep the program running until otherwise specified by user
tk.mainloop()
