# Drawing together all calculus on geom forms that I have into a single GUI
# All of which can later be added onto the main GUI
# through clearing the 'right_frame' when correct option is selected
# and inserting this.

# import all needed modules:
import timeit
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import tkinter as tk
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
from sympy import diff
from sympy.parsing.sympy_parser import parse_expr
from sympy import simplify, integrate, Symbol
from matplotlib import patches as patch

# %%

# start the timer
start = timeit.default_timer()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set up basic layout of the window
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# define an object tracker - for GUI set up
root = tk.Tk()

# set its title
root.title('Differential Forms Calculus')

# set a window size for it all to initially appear in
# do so by extracting the size of user's screen
# the method was found on blog.pythonlibrary
width = root.winfo_screenwidth()
height = root.winfo_screenheight()

root.geometry(str(width) + 'x' + str(height))

# set up frames for each of:
# bottom side (field, scalings etc) and the right side (with detailed options)
# and top left for plot

# right frame:
right_frame = tk.LabelFrame(root, text='Options Frame', padx=32, pady=5)
right_frame.grid(row=1, column=1)

# bot frame:
bot_frame = tk.LabelFrame(root, text='Inputs Frame', padx=30, pady=25)
bot_frame.grid(row=2, column=0)

# plot frame:
plot_frame = tk.LabelFrame(root, text='Graph Frame', padx=5, pady=5)
plot_frame.grid(row=1, column=0)

# plot characteristics frame and plot button
small_frame = tk.LabelFrame(root, text='Plot buttons Frame', padx=35, pady=5)
small_frame.grid(row=2, column=1)


'''

Define initially needed functions

'''


# function make variables python understood
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
    string = string.replace('ARCTAN', 'np.arctan')
    string = string.replace('ARCSIN', 'np.arcsin')
    string = string.replace('ARCCOS', 'np.arccos')
    string = string.replace('TANH', 'np.tanh')
    string = string.replace('SINH', 'np.sinh')
    string = string.replace('COSH', 'np.cosh')
    string = string.replace('^', '**')
    string = string.replace('ln', 'np.log')
    string = string.replace('e**', 'np.exp')
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


# define a function that will take care of constant and zero 2 forms
def form_2_constant_correction(form_2_eq):
    # want to check if it contains x or y and if not, make the shape correct
    # without using eq_to_comps becuase we do not change it to any components
    # here, although the process is very similar
    if form_2_eq.find('xg') & form_2_eq.find('yg') == -1:
        form_2_eq = str(form_2_eq) + '* np.ones(np.shape(xg))'
    else:
        pass
    return form_2_eq


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


# define a function to plot the simplified 2 forms, with coloured squares
def plot_form(form_2, fract_s):
    # celar the currently present plot
    ax.clear()
    # redefine the axis limits
    ax.set_xlim(-ax_L, ax_L)
    ax.set_ylim(-ax_L, ax_L)
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
                rect = patch.Rectangle((xg[i, j] - Mag[i, j]/4, yg[i, j] - Mag[i, j]/4), Mag[i, j]/2, Mag[i, j]/2, color='blue')
                ax.add_patch(rect)
            elif sign_mag[i, j] == 1:
                rect = patch.Rectangle((xg[i, j] - Mag[i, j]/4, yg[i, j] - Mag[i, j]/4), Mag[i, j]/2, Mag[i, j]/2, color='red')
                ax.add_patch(rect)
            else:
                rect = patch.Rectangle((xg[i, j] - Mag[i, j]/4, yg[i, j] - Mag[i, j]/4), Mag[i, j]/2, Mag[i, j]/2, color='grey')
                ax.add_patch(rect)  # not needed as zero sign will have zero magnitude therefore will not be seen anyway.


# #############################################################################
# supply canvas to the plot frame with a figure and connect it to mouse clicks
# #############################################################################

fig = plt.figure(figsize=(8, 6))
ax = fig.gca()
ax.set_xlabel('$x$')
ax.set_ylabel('$y$')
ax.set_aspect('equal')
fig.tight_layout()
plt.close()
canvas = FigureCanvasTkAgg(fig, master=plot_frame)
canvas.draw()
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
toolbar = NavigationToolbar2Tk(canvas, plot_frame)
toolbar.update()  # allow the plot to update based on the toolbar
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)


# bind it to the mouse events
def on_key_press(event):
    print("you pressed {}".format(event.key))
    key_press_handler(event, canvas, toolbar)


'''

Define all needed intial variables to start with

'''

# define scale of the graph
L = 5
pt_den = 31   # number of points on each axis

# with L, correct the axis and redraw
delta_factor = 10
ax_L = L + L/delta_factor
ax.set_xlim(-ax_L, ax_L)
ax.set_ylim(-ax_L, ax_L)
canvas.draw()

# define x and y values
x = np.linspace(-L, L, pt_den)
y = np.linspace(-L, L, pt_den)
# create a grid on x-y plane
xg, yg = np.meshgrid(x, y)

# set the dimensionality
m = 2

'''
start the program with a 2 form being supplied and plotted by blocks only
define intiial variables for these:
'''

# define the initial 2 form, in terms of a string, equation and numerically
form_2_str = 'x*y**2'  # dx^dy component
form_2_eq = format_eq(form_2_str)
form_2 = eval(form_2_eq)

# get the signs of the 2 form
form_2_sgn = np.sign(form_2)

# set a maximum side size for blocks as a fraction of the graph size
fract_s = 0.1

# put the initial plot onto the canvas
plot_form(form_2, fract_s)
canvas.draw()

'''
define initial variables for all pother operations that will be completed
 using the GUI by the user
'''

# ##### GENERAL ONES FOR STACKS #######

# define colours to use for 2 form plots with components
# first string defines colour for positive (ccw), second for negative (cw)
# last one is an in case, for when the magnitude is exactly zero.
colour_str = ['red', 'blue', 'grey']

# define fract of graph to set as size of stack
fract = 0.06

# define max number of sheets to start with
s_max = 5

# set up a zero vector filed to plot x and y components as 2 separate fields:
zero_field = np.zeros(np.shape(xg))

# ###### ONES NEEDED FOR 2 FORM FROM 1 FORMS #######

# define an example vector field, now - from string, initially
string_x = 'sin(x+y) - y'  # x component
string_y = 'sin(x+y)'  # y component

# take the input strings and turn them into sympy expressions to be able to
# use sympy's partial differentiation
sympy_expr_x = parse_expr(string_x, evaluate=False)
sympy_expr_y = parse_expr(string_y, evaluate=False)
# for m > 2, need more components, and need these in 'expressions' too!

# combine the 2 into a list:
expressions = np.array([sympy_expr_x, sympy_expr_y])

# use sympy partial derrivatives on these, as to get a 2-form on R2:
# need to differentiate each component w.r.t the coordinates that it's
# elementary 1 form does not contain.

# set up an array of coordinates that need to be used (in standard order)
coords = ['x', 'y']

# define a unit field so that code deals with constant fields
field_unit = np.ones(np.shape(xg))  # need it here cuz its used in eq_to_comps and in find_2_form

# set up initial strings for 2 forms window to display, for it to save properly after
to_wedge_x_1_str = ''
to_wedge_y_1_str = ''
to_wedge_x_2_str = ''
to_wedge_y_2_str = ''

# predefine vectors to store for the interior derivative
vector_ex_str = ''
vector_ey_str = ''

# predefine strings for splitting in x and in y for 2 form stacks
x_split_str = '1/2'
y_split_str = '1/2'


'''

Define other functions needed to be used in the responses to GUI interactions

'''


# define a function that will plot stack components, coloured
# as per the orientation of the 2 form at that grid point
def form_2_components_plot(xg, yg, u, v, form_2_sgn, s_max, L, fract, colour_str, arrowheads=False, w_head=1/8, h_head=1/4):
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
    
    # define points of stack arrowheads as arrays for all stacks
    p_sh1x = xg + (s_L/2)*I_cos + (sheet_L*w_head)*I_sin
    p_sh1y = yg + (s_L/2)*I_sin - (sheet_L*w_head)*I_cos
    p_sh2x = xg + (s_L/2)*I_cos - (sheet_L*w_head)*I_sin
    p_sh2y = yg + (s_L/2)*I_sin + (sheet_L*w_head)*I_cos
    p_sh3x = xg + (s_L*0.5 + s_L*h_head)*I_cos
    p_sh3y = yg + (s_L*0.5 + s_L*h_head)*I_sin
    
    # define these for when there is only 1 line in the stack plot:
    P_sh1x = xg + (sheet_L*w_head)*I_sin
    P_sh1y = yg - (sheet_L*w_head)*I_cos
    P_sh2x = xg - (sheet_L*w_head)*I_sin
    P_sh2y = yg + (sheet_L*w_head)*I_cos
    P_sh3x = xg + (s_L*h_head)*I_cos
    P_sh3y = yg + (s_L*h_head)*I_sin
    
    # loop over each arrow coordinate in x and y
    for i in range(x_len):
        for j in range(y_len):
            # define it for all magnitudes. Separately for odd and even corr. number of sheets:
            # Label each element with the number of stacks required: linear scaling
            
            if form_2_sgn[i, j] == +1:
                color_index = 0
            elif form_2_sgn[i, j] == -1:
                color_index = 1
            else:
                color_index = 2  # just in case
            
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
                    
                    ax.add_line(Line2D((Ax1, Bx1), (Ay1, By1), linewidth=0.7, color=colour_str[color_index]))
                    ax.add_line(Line2D((Ax2, Bx2), (Ay2, By2), linewidth=0.7, color=colour_str[color_index]))
                    
                    # change the parameter to loop over all changes in displacement for current magnitude
                    s += 1
            if arrowheads is True:
                # plot lines of arrowheads from central sheet for n = 1 or on top sheet for n>1 
                if n > 1:   # for all lines ubt the single sheet one
                    ax.add_line(Line2D((p_sh1x[i, j], p_sh3x[i, j]), (p_sh1y[i, j], p_sh3y[i, j]), linewidth=1, color=colour_str[color_index]))
                    ax.add_line(Line2D((p_sh2x[i, j], p_sh3x[i, j]), ((p_sh2y[i, j], p_sh3y[i, j])), linewidth=1, color=colour_str[color_index]))
                # then define it for the stacks with only 1 sheet:
                else:
                    ax.add_line(Line2D((P_sh1x[i, j], P_sh3x[i, j]), (P_sh1y[i, j], P_sh3y[i, j]), linewidth=1, color=colour_str[color_index]))
                    ax.add_line(Line2D((P_sh2x[i, j], P_sh3x[i, j]), ((P_sh2y[i, j], P_sh3y[i, j])), linewidth=1, color=colour_str[color_index]))
            else:
                pass


# define a function that will complete all stack plotting:
def stack_plot(xg, yg, axis, u, v, s_max, L, fract, arrows=False, orientation='mid', scale=1, w_head=1/8, h_head=1/4, axis_check=0, arrowheads=True, colour='green'):
    global s_L
    # get the lengths of x and y from their grids
    
    x_len = len(xg[:, 0])
    y_len = len(yg[0, :])
    
    # Scaling of axes and setting equal proportions circles look like circles
    # axis.set_aspect('equal')
    # ax_L = L + L/delta_factor
    # axis.set_xlim(-ax_L, ax_L)
    # axis.set_ylim(-ax_L, ax_L)
    ax_L = L + L/delta_factor
    axis.set_xlim(-ax_L, ax_L)
    axis.set_ylim(-ax_L, ax_L)
    
    # define an empty array of magnitudes, to then fill with integer rel. mags
    R_int = np.zeros(shape=((x_len), (y_len)))
    
    # #########################################################################
    # plot the initial quiver plot to work from
    # #########################################################################
    
    # plot the quiver plot on grid points if chosen in original function
    if arrows is True:
        axis.quiver(xg, yg, u, v, pivot=orientation, scale=scale, scale_units='xy')
    else:
        pass
    
    # #########################################################################
    # get variables needed for the initial, simplified stack plot
    # #########################################################################
    # find the arrow length corresponding to each point and store in mag array
    mag = np.sqrt(u**2 + v**2)
    # find direction of each arrow
    angles = np.arctan2(v, u)   # theta defined from positive x axis ccw
    
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
    I_sin = np.sin(angles)
    I_cos = np.cos(angles)
    
    # define the points that set out a line of the stack sheet (middle line)
    A_x = xg + (sheet_L/2)*I_sin
    A_y = yg - (sheet_L/2)*I_cos
    B_x = xg - (sheet_L/2)*I_sin
    B_y = yg + (sheet_L/2)*I_cos
    
    # define points of stack arrowheads as arrays for all stacks
    p_sh1x = xg + (s_L/2)*I_cos + (sheet_L*w_head)*I_sin
    p_sh1y = yg + (s_L/2)*I_sin - (sheet_L*w_head)*I_cos
    p_sh2x = xg + (s_L/2)*I_cos - (sheet_L*w_head)*I_sin
    p_sh2y = yg + (s_L/2)*I_sin + (sheet_L*w_head)*I_cos
    p_sh3x = xg + (s_L*0.5 + s_L*h_head)*I_cos
    p_sh3y = yg + (s_L*0.5 + s_L*h_head)*I_sin
    
    # define these for when there is only 1 line in the stack plot:
    P_sh1x = xg + (sheet_L*w_head)*I_sin
    P_sh1y = yg - (sheet_L*w_head)*I_cos
    P_sh2x = xg - (sheet_L*w_head)*I_sin
    P_sh2y = yg + (sheet_L*w_head)*I_cos
    P_sh3x = xg + (s_L*h_head)*I_cos
    P_sh3y = yg + (s_L*h_head)*I_sin
    
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
                    axis.add_line(Line2D((Ax1, Bx1), (Ay1, By1), linewidth=1, color=colour))
                    axis.add_line(Line2D((Ax2, Bx2), (Ay2, By2), linewidth=1, color=colour))
                    
                    # update parameter to reapet and draw all needed arrows
                    s += 1
            # deal with the odd number of stacks:
            elif parity(n) is False:
                # Add the centre line for odd numbers of stacks
                axis.add_line(Line2D((A_x[i, j], B_x[i, j]), (A_y[i, j], B_y[i, j]), linewidth=1, color=colour))
                
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
                    axis.add_line(Line2D((Ax1,Bx1),(Ay1,By1), linewidth=1, color=colour))
                    axis.add_line(Line2D((Ax2,Bx2),(Ay2,By2), linewidth=1, color=colour))
                    
                    # change the parameter to loop over all changes in displacement for current magnitude
                    s += 1
            if arrowheads == True:
                # plot lines of arrowheads from central sheet for n = 1 or on top sheet for n>1 
                if n > 1:   # for all lines ubt the single sheet one
                    axis.add_line(Line2D((p_sh1x[i, j],p_sh3x[i, j]),(p_sh1y[i, j],p_sh3y[i, j]), linewidth=1, color='green'))
                    axis.add_line(Line2D((p_sh2x[i, j],p_sh3x[i, j]),((p_sh2y[i, j],p_sh3y[i, j])), linewidth=1, color='green'))
                # then define it for the stacks with only 1 sheet:
                else:
                    axis.add_line(Line2D((P_sh1x[i, j], P_sh3x[i, j]), (P_sh1y[i, j], P_sh3y[i, j]), linewidth=1, color='green'))
                    axis.add_line(Line2D((P_sh2x[i, j], P_sh3x[i, j]), ((P_sh2y[i, j], P_sh3y[i, j])), linewidth=1, color='green'))
            else:
                pass
    plt.close()


# define a function that will find the 2 form from given expressions
# in a given number of dimensions and in terms of given coordinate symbols
def find_2_form(expressions, coords, m=2):
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
    
    return form_2


'''

Define response functions to GUI interactions

'''


# gets 2 form from entry box and plots it as coloured blocks only
def form_2_response():
    global form_2_str, form_2_eq, form_2_sgn, form_2, comp_x, comp_y, u, v
    # get the input 2 form
    form_2_str = str(simplify(form_2_entry.get()))
    # format it to be python understood
    form_2_eq = format_eq(form_2_str)
    # check against constant and zero 2 forms being supplied
    form_2_eq = form_2_constant_correction(form_2_eq)
    # get the numerical evaluation of it
    form_2 = eval(form_2_eq)
    # get the signs of thsi new 2 form
    form_2_sgn = np.sign(form_2)
    # depending on chosen setting, complete the calcualtions and plots
    if stack_block_int == 0:
        # plot the new form using the previously define funtion
        plot_form(form_2, fract_s)
        canvas.draw()
    elif stack_block_int == 1:
        # split this equation into two ARBITRARY parts
        # supposing that the 2 form came from those via the exterior derivative
        # for the example, split it into two EQUAL components
        # use the inut splittings:
        eq_1 = str(simplify('(' + form_2_str + ')' + '*(' + x_split_str + ')'))
        eq_2 = str(simplify('(' + form_2_str + ')' + '*(' + y_split_str + ')'))
        # the above is needs the input as a string, and symbol, that needs to be changed
        # to a string to be able to add minus to it from dy^dx ---> dx^dy
        # and that needs to be simplified
        # after which, it needs to be changed to a string again
        # NOTE:
        # It can be kept as a sympy type expression, but its not as easily manipulated
        # and won't show in variable explorer.
        # turn these equatiosn to components to use in plotting the 2 form from stacks
        u, v = eq_to_comps(eq_1, eq_2, xg, yg)
        # clear the current plot
        ax.clear()
        # use plotting stacks to display these
        form_2_components_plot(xg, yg, u, zero_field, form_2_sgn, s_max, L, fract, colour_str)
        form_2_components_plot(xg, yg, zero_field, v, form_2_sgn, s_max, L, fract, colour_str)
        # display the new plot
        canvas.draw()
        # define as string message to show in the message box
        uniqueness_message = '''
        THIS RESULT IN NOT UNIQUE!!! \n
        Depeding on the way that the 2 form
        is split into dx^dy and dy^dx, the result will be
        different. \n
        All of them are quivalent in giving the same 2 form
        in the end, but local details are not consistant. \n
        This result is one possible representation, out of all which
        belong to a rotational group. The superposition of all of of its
        elements gives the total representation - blocks, not local.
        For a consistant expression, need to use the blocks method,
        lossing local detail
        '''
        # display it
        tk.messagebox.showinfo('Uniqueness of this result', uniqueness_message)
    # display a background green on the 2 form entry to show that
    # this entry is being displayed now.
    form_2_entry.configure(bg='#C0F6BB')
    # undo it for 1 forms
    form_1_x_entry.configure(bg='#FFFFFF')
    form_1_y_entry.configure(bg='#FFFFFF')


# plots the vetor field with stacks only
def form_1_stacks_response():
    global u, v, stack_block_int
    # clear the current axis
    ax.clear()
    # get the supplied 1 forms from entry boxes:
    x_comp_str = str(simplify(form_1_x_entry.get()))
    y_comp_str = str(simplify(form_1_y_entry.get()))
    # take all these values, and the input from field component bnoxes to set up the field:
    u, v = eq_to_comps(x_comp_str, y_comp_str, xg, yg)
    # plot the new field
    stack_plot(xg, yg, ax, u, v, s_max, L, fract)
    # put it onto the screen
    canvas.draw()
    # display the fact that this is a 1 form therefore always shows stacks
    stack_block.set(1)
    # and update its parameter also
    stack_block_int = 1
    # display a background green on 1 form components and get rid of the 2 form
    # colour, to show which one is being plotted
    form_1_x_entry.configure(bg='#C0F6BB')
    form_1_y_entry.configure(bg='#C0F6BB')
    # get rid of the 2 form colour:
    form_2_entry.configure(bg='#FFFFFF')


# performs the interior derivative on supplied 2 form and plots it as stacks
# for 2 form only, not including the 1 form.
# if combined, use different fucntion.
def Int_deriv_2_form():
    global u, v, u_str, v_str, vector_ex_str, vector_ey_str, vector_ex_eq, vector_ey_eq, vector_ex, vector_ey
    # take the supplied componments and save them globally
    vector_ex_str = str(simplify(int_vect_ex_entry.get()))
    vector_ey_str = str(simplify(int_vect_ey_entry.get()))
    # turn these into equations
    vector_ex_eq = format_eq(vector_ex_str)
    vector_ey_eq = format_eq(vector_ey_str)
    # check against zeros
    vector_ex_eq = form_2_constant_correction(vector_ex_eq)
    vector_ey_eq = form_2_constant_correction(vector_ey_eq)
    # find numerical evaluation of it
    vector_ex = eval(vector_ex_eq)
    vector_ey = eval(vector_ey_eq)
    # get the input 2 form
    form_2_str = str(simplify(form_2_entry.get()))
    # format it to be python understood
    form_2_eq = format_eq(form_2_str)
    # check against constant and zero 2 forms being supplied
    form_2_eq = form_2_constant_correction(form_2_eq)
    # get the numerical evaluation of it
    form_2 = eval(form_2_eq)
    # get the signs of thsi new 2 form
    form_2_sgn = np.sign(form_2)
    # using interior product, get the u and v (dx and dy) components
    # of the resulting 1 form
    u = -form_2 * vector_ey
    v = form_2 * vector_ex
    # to be usable in ext_deriv, define strings of these variables
    # later to put them into their entry boxes.
    u_str = str(simplify('-(' + form_2_str + ')*(' + vector_ey_str + ')' ))
    v_str = str(simplify( '(' + form_2_str + ')*(' + vector_ex_str + ')' ))
    # use the stacks plotter to present this
    form_2_components_plot(xg, yg, u, v, form_2_sgn, s_max, L, fract, colour_str, arrowheads=True)


# define a function that will find the interior derivative of a given 1 form
def Int_deriv_1_form():
    global zero_form_str, zero_form, vector_ex_str, vector_ey_str, vector_ex_eq, vector_ey_eq, vector_ex, vector_ey, x_comp, y_comp, x_comp_eq, y_comp_eq
    # take the supplied componments and save them globally
    vector_ex_str = str(simplify(int_vect_ex_entry.get()))
    vector_ey_str = str(simplify(int_vect_ey_entry.get()))
    # turn these into equations
    vector_ex_eq = format_eq(vector_ex_str)
    vector_ey_eq = format_eq(vector_ey_str)
    # check against zeros
    vector_ex_eq = form_2_constant_correction(vector_ex_eq)
    vector_ey_eq = form_2_constant_correction(vector_ey_eq)
    # find numerical evaluation of it
    vector_ex = eval(vector_ex_eq)
    vector_ey = eval(vector_ey_eq)
    # get the input 1 forms
    x_comp_str = str(simplify(form_1_x_entry.get()))
    y_comp_str = str(simplify(form_1_y_entry.get()))
    # format them to be python understood
    x_comp_eq = format_eq(x_comp_str)
    y_comp_eq = format_eq(y_comp_str)
    # check against constant and zeros being supplied
    x_comp_eq = form_2_constant_correction(x_comp_eq)
    y_comp_eq = form_2_constant_correction(y_comp_eq)
    # get the numerical evaluations
    x_comp = eval(x_comp_eq)
    y_comp = eval(y_comp_eq)
    # as per the ineterior product, get the scalar function from this
    # first, numerically
    zero_form = x_comp*vector_ex + y_comp*vector_ey
    # then as a string to display in an extra window
    zero_form_str = str(simplify('(' + x_comp_str + ')*(' + vector_ex_str + ')' + ' + (' + y_comp_str + ')*(' + vector_ey_str + ')'))
    # plot the zero_form as contours with labeled levels
    contour_levels = np.linspace(np.min(zero_form), np.max(zero_form), 10)
    ax.contour(xg, yg, zero_form, contour_levels)


# define a function that will respond to plotting the 2 form only
def Int_deriv_22_form():
    global stack_block_int
    # clear the axis
    ax.clear()
    # call the asked fucntion
    Int_deriv_2_form()
    # draw its result
    canvas.draw()
    # display the fact that this is a 1 form therefore always shows stacks
    stack_block.set(1)
    # and update its parameter also
    stack_block_int = 1
    # change the entry box 1 form to the calculated ones
    form_1_x_entry.delete(0, 'end')
    form_1_y_entry.delete(0, 'end')
    form_1_x_entry.insert(0, u_str)
    form_1_y_entry.insert(0, v_str)
    # display that the 1 form is now being plotted, therefore get rid of
    # 2 form colour and show 1 form components in green:
    form_1_x_entry.configure(bg='#C0F6BB')
    form_1_y_entry.configure(bg='#C0F6BB')
    form_2_entry.configure(bg='#FFFFFF')


# define a function that will find the interior derivative of both the 2 form
# and a 1 form, merged.
def Int_deriv_21_form():
    global stack_block_int
    # clear the axis
    ax.clear()
    # first call the function to plot a 2 form and plot it
    Int_deriv_2_form()
    # then call the function that will do this for the 1 form
    # plot these together
    Int_deriv_1_form()
    # draw its result
    canvas.draw()
    # make a label for the found 0 form
    Label_zero_form.configure(text=zero_form_str)
    # display the fact that this is a 1 form therefore always shows stacks
    stack_block.set(1)
    # and update its parameter also
    stack_block_int = 1
    # change the entry box 1 form to the calculated ones
    form_1_x_entry.delete(0, 'end')
    form_1_y_entry.delete(0, 'end')
    form_1_x_entry.insert(0, u_str)
    form_1_y_entry.insert(0, v_str)
    # display that the 1 form is now being plotted, therefore get rid of
    # 2 form colour and show 1 form components in green:
    form_1_x_entry.configure(bg='#C0F6BB')
    form_1_y_entry.configure(bg='#C0F6BB')
    form_2_entry.configure(bg='#FFFFFF')
    return


# Interior derivative response for 2 form only
# asks to give vectors w.r.t to which perform iota.
def Int_deriv_2_form_response(type_form):
    global int_vector_window, int_vect_ex_entry, int_vect_ey_entry
    # open a titled new window
    int_vector_window = tk.Toplevel()
    int_vector_window.title('input a vector for the interior derivative')
    # define entry boxes for e^x and e^y components of desired vector
    tk.Label(int_vector_window, text='vector component e^x').grid(row=0, column=0)
    int_vect_ex_entry = tk.Entry(int_vector_window, width=30, borderwidth=1)
    int_vect_ex_entry.insert(0, vector_ex_str)
    int_vect_ex_entry.grid(row=1, column=0)
    tk.Label(int_vector_window, text='vector component e^y').grid(row=2, column=0)
    int_vect_ey_entry = tk.Entry(int_vector_window, width=30, borderwidth=1)
    int_vect_ey_entry.insert(0, vector_ey_str)
    int_vect_ey_entry.grid(row=3, column=0)
    # define a button that will plot these
    # for a 2 form only
    if type_form == 2:
        int_vector_load_btn = tk.Button(int_vector_window, text='PLOT', padx=20, pady=10, command=Int_deriv_22_form)
        int_vector_load_btn.grid(row=4, column=0, pady=10)
    # for 2 form with 1 form included:
    elif type_form == 1:
        int_vector_load_btn = tk.Button(int_vector_window, text='PLOT', padx=20, pady=10, command=Int_deriv_21_form)
        int_vector_load_btn.grid(row=4, column=0, pady=10)


# define a function that will respond to the made choice reg. int deriv.
def int_deriv_choice(var):
    if var == 0:
        # only use 2 form, therefore call previous functions for this
        Int_deriv_2_form_response(2)
    elif var == 1:
        # call the function, but make it deal with 2 form together with the
        # 1 form
        Int_deriv_2_form_response(1)
    # close the previous window
    int_option_window.destroy()


# define response to int deriv button
# shows window to select if should just use 2 form
# or a comnination of 2 form and 1 form
def Int_deriv_response():
    global int_option_window
    # show new window with label and two options.
    int_option_window = tk.Toplevel()
    int_option_window.title('input a vector for the interior derivative')
    int_option_window.geometry('425x200')
    # define Label
    tk.Label(int_option_window, text='Perform w.r.t 2 form only or combine given 2 form and 1 form ?').grid(row=0, column=0, columnspan=2)
    # define response buttons to the stated question
    form_2_only_btn = tk.Button(int_option_window, text='2 form', padx=30, pady=30, command=lambda: int_deriv_choice(0))
    from_2_and_1_btn = tk.Button(int_option_window, text='Both', padx=30, pady=30, command=lambda: int_deriv_choice(1))
    form_2_only_btn.grid(row=1, column=0, pady=20)
    from_2_and_1_btn.grid(row=1, column=1, pady=20)


# perform ext deriv on the result of int_deriv and plots it as stacks
def Ext_deriv_response():
    global form_2, form_2_str, form_2_sgn
    # celar current axis
    ax.clear()
    # get the inpus from fields of x and u components
    x_comp_str = str(simplify(form_1_x_entry.get()))
    y_comp_str = str(simplify(form_1_y_entry.get()))
    # from found u and v in the interior derivative, set up sympy components
    sympy_expr_x = parse_expr(x_comp_str, evaluate=False)
    sympy_expr_y = parse_expr(y_comp_str, evaluate=False)
    # combine the 2 into a list:
    expressions = np.array([sympy_expr_x, sympy_expr_y])
    # set up an array of coordinates that need to be used (in standard order)
    coords = ['x', 'y']
    # from these, use the find_2_form function to get the 2 form
    form_2 = find_2_form(expressions, coords, m)[0]
    # get the signs of the 2 form
    form_2_sgn = np.sign(form_2)
    # get the string of this new 2 form to use it in int deriv
    # also put it into the entry
    form_2_str = str(simplify(str(unformat(result[0][0]))))
    # unformat it to display in the entry box, this way it does not
    # format twice if int deriv runs again
    form_2_str = unformat(form_2_str)
    form_2_entry.delete(0, 'end')
    form_2_entry.insert(0, form_2_str)
    # plot now depends on the chosen option - stacks or blocks
    if stack_block_int == 1:
        # first need to split up result into equal parts, as before
        # AGAIN - NOT UNIQUE - therefore show the message here
        # split this equation into two ARBITRARY parts
        # supposing that the 2 form came from those via the exterior derivative
        # for the example, split it into two EQUAL components
        eq_1 = str(simplify('(' + form_2_str + ')' + '*(' + x_split_str + ')'))
        eq_2 = str(simplify('(' + form_2_str + ')' + '*(' + y_split_str + ')'))
        # turn these equatiosn to components to use in plotting the 2 form from stacks
        u, v = eq_to_comps(eq_1, eq_2, xg, yg)
        # clear the current plot
        ax.clear()
        # use plotting stacks to display these
        form_2_components_plot(xg, yg, u, zero_field, form_2_sgn, s_max, L, fract, colour_str)
        form_2_components_plot(xg, yg, zero_field, v, form_2_sgn, s_max, L, fract, colour_str)
        # display the new plot
        canvas.draw()
        # define as string message to show in the message box
        uniqueness_message = '''
        THIS RESULT IN NOT UNIQUE!!! \n
        Cancellations which are not accounted for
        may occur when the two derivatives (dx^dy and dy^dx) are
        merged
        To avoid this issue, the plot is found by splitting
        the resulting 2 form into components \n
        HOWEVER \n
        Depeding on the way that the 2 form
        is split into dx^dy and dy^dx, the result will be
        different. \n
        All of them are quivalent in giving the same 2 form
        in the end, but local details are not consistant. \n
        This result is one possible representation, out of all which
        belong to a rotational group. The superposition of all of of its
        elements gives the total representation - blocks, not local.
        For a consistant expression, need to use the blocks method,
        lossing local detail
        '''
        
        # display it
        tk.messagebox.showinfo('Uniqueness of this result', uniqueness_message)
    elif stack_block_int == 0:
        # plot the new form using the previously define funtion
        plot_form(form_2, fract_s)
        canvas.draw()
    # display a background green on the 2 form entry to show that
    # this entry is being displayed now.
    form_2_entry.configure(bg='#C0F6BB')
    # undo it for 1 forms
    form_1_x_entry.configure(bg='#FFFFFF')
    form_1_y_entry.configure(bg='#FFFFFF')


# define a function that will wedge two 1 forms and plot them
def wedge_product():
    global to_wedge_x_1_str, to_wedge_y_1_str, to_wedge_x_2_str, to_wedge_y_2_str, form_2_str, form_2_eq, form_2, form_2_sgn
    # first, get all entries out, save as string for these to display when
    # window is opened again
    to_wedge_x_1_str = str(simplify(str(to_wedge_x_1_entry.get())))
    to_wedge_y_1_str = str(simplify(str(to_wedge_y_1_entry.get())))
    to_wedge_x_2_str = str(simplify(str(to_wedge_x_2_entry.get())))
    to_wedge_y_2_str = str(simplify(str(to_wedge_y_2_entry.get())))
    u_1, v_1 = eq_to_comps(to_wedge_x_1_str, to_wedge_y_1_str, xg, yg)
    u_2, v_2 = eq_to_comps(to_wedge_x_2_str, to_wedge_y_2_str, xg, yg)
    # clear the axis:
    ax.clear()
    # first, find the result of the 2 form
    # this if, in terms of the above commented fields:
    # 2 form = f*m - g*h
    # get it mathematically, as a string
    form_2_str = str(simplify( '(' + to_wedge_x_1_str + ')*(' +  to_wedge_y_2_str + ')' + ' - (' + to_wedge_y_1_str + ')*(' +  to_wedge_x_2_str + ')' ))
    # put it into the entry box for 2 forms
    form_2_entry.delete(0, 'end')
    form_2_entry.insert(0, form_2_str)
    
    # complete the process depending on the selected radiobutton for stacks
    # or blocks
    if stack_block_int == 1:
        # plot these as stacks, with no arrowheads, on top of one another.
        
        # Given w_1 = fdx+gdy  and w_2=hdx+mdy. The graphical representation must be:
        # fdx /\ mdy "and" -hdx/\gdy {equivalend to [(u_1,0)+(0,v2)] and [(-u_2,0)+(0,v1)]},
        # which is executed when the if condition below is satisfited. Basically two rectangular
        # stacks with the scaling factor are accounted to via giving similar colors to the vertical
        # stacks (red) and green to the horizantal ones. After the first rectagular stacks are added
        # the second group will either sit on top of the first (in which case scaling contibution is zero)
        # or sit in some gaps and hence increasing the denisty as result of its scaling function.
        # If any of the coefficients (f,g,h and m)  is zero, the stacking reduces to one "function*dx/\dy", these are executed in the elif options.
        # One situation to be added when (u_1*v_2-u2*v_1).all() = 0, the scaling function here is zero and hence no 2-form should be produced/ or produced in faded color.
        # !!! ISSUES:
        # 1- The max number of stacks possible will hinder good visualization when
        # the stacks are dense. It's a general issue, but more clear here than other cases due to the nature of 2-forms.
        # 2- Would be nice to uniformly distribute the stacks in each direction after finishing the double stacking.
        if to_wedge_x_1_str != '0' and to_wedge_y_1_str != '0' and to_wedge_x_2_str != '0' and to_wedge_y_2_str != '0':
            stack_plot(xg, yg, ax, u_1, 0, s_max, L, fract, arrowheads=False, colour='green')
            stack_plot(xg, yg, ax, 0, v_2, s_max, L, fract, arrowheads=False, colour='green')
            stack_plot(xg, yg, ax, 0, v_1, s_max, L, fract, arrowheads=False, colour='red')
            stack_plot(xg, yg, ax, -u_2, 0, s_max, L, fract, arrowheads=False, colour='red')
        elif to_wedge_x_1_str != '0' and to_wedge_y_2_str != '0':
            stack_plot(xg, yg, ax, u_1, 0, s_max, L, fract, arrowheads=False, colour='red')
            stack_plot(xg, yg, ax, 0, v_2, s_max, L, fract, arrowheads=False, colour='green')
        elif to_wedge_y_1_str != '0' and to_wedge_x_2_str != '0':
            stack_plot(xg, yg, ax, 0, v_1, s_max, L, fract, arrowheads=False, colour='green')
            stack_plot(xg, yg, ax, -u_2, 0, s_max, L, fract, arrowheads=False, colour='red')
    elif stack_block_int == 0:
        # to do it with blocks, first, get the numercial 2 from
        # from previously found string
        # format it to be python understood
        form_2_eq = format_eq(form_2_str)
        # check against constant and zero 2 forms being supplied
        form_2_eq = form_2_constant_correction(form_2_eq)
        # get the numerical evaluation of it
        form_2 = eval(form_2_eq)
        # get the signs of thsi new 2 form
        form_2_sgn = np.sign(form_2)
        # plot the new form using the previously define funtion
        plot_form(form_2, fract_s)
    # first check against cancellation effects
    if str(simplify( '(' + to_wedge_x_1_str + ')*(' +  to_wedge_y_2_str + ')' + ' - (' + to_wedge_y_1_str + ')*(' +  to_wedge_x_2_str + ')' )) == '0':
        # delete whatever plot was created
        ax.clear()
        # update axis limits, becuase now they were cleared
        ax.set_aspect('equal')
        ax_L = L + L/delta_factor
        ax.set_xlim(-ax_L, ax_L)
        ax.set_ylim(-ax_L, ax_L)
        # display a window info that the wedge product is zero
        tk.messagebox.showinfo('Nothing to show', "The wedge product is zero")
    # put these onto the canvas
    canvas.draw()
    # close the extra window
    wedge_2_window.destroy()
    # neither the entered single 1 form or the single 2 form are being plotted
    # but the resulting 2 form is input into the entry box
    # therefore make 2 form green.
    form_1_x_entry.configure(bg='#FFFFFF')
    form_1_y_entry.configure(bg='#FFFFFF')
    form_2_entry.configure(bg='#C0F6BB')


# define a reponse function, opens new window where two 1 forms to be wedged can be entered
def wedge_2_response():
    global wedge_2_window, to_wedge_x_1_entry, to_wedge_y_1_entry, to_wedge_x_2_entry, to_wedge_y_2_entry
    # open a titled new window
    wedge_2_window = tk.Toplevel()
    wedge_2_window.title('input two 1 forms to wedge')
    # define all entry boxes
    tk.Label(wedge_2_window, text='first 1 form x component :').grid(row=0, column=0)
    to_wedge_x_1_entry = tk.Entry(wedge_2_window, width=30, borderwidth=1)
    to_wedge_x_1_entry.insert(0, to_wedge_x_1_str)
    to_wedge_x_1_entry.grid(row=1, column=0)
    tk.Label(wedge_2_window, text='first 1 from y component:').grid(row=2, column=0)
    to_wedge_y_1_entry = tk.Entry(wedge_2_window, width=30, borderwidth=1)
    to_wedge_y_1_entry.insert(0, to_wedge_y_1_str)
    to_wedge_y_1_entry.grid(row=3, column=0)
    tk.Label(wedge_2_window, text='second 1 form x component :').grid(row=4, column=0)
    to_wedge_x_2_entry = tk.Entry(wedge_2_window, width=30, borderwidth=1)
    to_wedge_x_2_entry.insert(0, to_wedge_x_2_str)
    to_wedge_x_2_entry.grid(row=5, column=0)
    tk.Label(wedge_2_window, text='second 1 form y component :').grid(row=6, column=0)
    to_wedge_y_2_entry = tk.Entry(wedge_2_window, width=30, borderwidth=1)
    to_wedge_y_2_entry.insert(0, to_wedge_y_2_str)
    to_wedge_y_2_entry.grid(row=7, column=0)
    # define a button that will plot these
    plot_wedge_btn = tk.Button(wedge_2_window, text='PLOT', padx=20, pady=10, command=wedge_product)
    plot_wedge_btn.grid(row=8, column=0, pady=10)


# define a fucntion that will respond to finding the Hodge dual of the given froms
# 1 form or 2 form depending on chosen parameter  - to be implemented later
def Hodge_response():
    tk.messagebox.showwarning('Hodge error', 'This has not yet been implemented, await further updates')


# define a function to respond to Radiobuttons for stacks vs blocks
def plot_type_response(type_stack_block):
    global stack_block_int
    stack_block_int = int(type_stack_block)


# response to saving customisations
def customisation_save():
    global fract_s, fract, pt_den, x, y, xg, yg, field_unit, zero_field, vector_ex, vector_ey, zero_form, x_comp, y_comp
    # get the inputs
    fract_s = float(fract_s_entry.get())
    fract = float(fract_entry.get())
    pt_den = int(pt_den_entry.get())
    # for these to correctly apply
    # aslo need to change what depends on them
    # therefore, here: grids:
    # define x and y values
    x = np.linspace(-L, L, pt_den)
    y = np.linspace(-L, L, pt_den)
    # create a grid on x-y plane
    xg, yg = np.meshgrid(x, y)
    # the field unit needs to update:
    field_unit = np.ones(np.shape(xg))
    # so does zero_field
    zero_field = np.zeros(np.shape(xg))
    # next ones may, or may not be present:
    try:
        # so do vector_ex and vector_ey
        vector_ex = eval(vector_ex_eq)
        vector_ey = eval(vector_ey_eq)
        # and the zero_form
        x_comp = eval(x_comp_eq)
        y_comp = eval(y_comp_eq)
        zero_form = x_comp*vector_ex + y_comp*vector_ey
    except NameError:
        pass
    # destroy the window
    customise_window.destroy()

# response button to customise features
def customise_calc():
    global customise_window, fract_s_entry, fract_entry, pt_den_entry
    # show a new window
    customise_window = tk.Toplevel()
    customise_window.title('custmoise stack size, block size and number of points')
    # define entry boxes
    # block size as fraction of graph
    tk.Label(customise_window, text='block size as fraction of graph size').grid(row=0, column=0)
    fract_s_entry = tk.Entry(customise_window, width=20, borderwidth=1)
    fract_s_entry.insert(0, fract_s)
    fract_s_entry.grid(row=1, column=0)
    # stack size as fraction of graph
    tk.Label(customise_window, text='stack size as fraction of graph size').grid(row=2, column=0)
    fract_entry = tk.Entry(customise_window, width=20, borderwidth=1)
    fract_entry.insert(0, fract)
    fract_entry.grid(row=3, column=0)
    # number of points along each axis
    tk.Label(customise_window, text='number of points along axis').grid(row=4, column=0)
    pt_den_entry = tk.Entry(customise_window, width=20, borderwidth=1)
    pt_den_entry.insert(0, pt_den)
    pt_den_entry.grid(row=5, column=0)
    # get a submittion button for these
    Custom_submit_btn = tk.Button(customise_window, text='SAVE', padx=20, pady=10, command=customisation_save)
    Custom_submit_btn.grid(row=6, column=0, pady=10)


# define a function to save splittings
def splitting_choice_save():
    global x_split_str, y_split_str
    # get the component split factors from inputs
    x_split_str = str(simplify(x_split_entry.get()))
    y_split_str = str(simplify(y_split_entry.get()))
    # close the extra window
    splitting_window.destroy()


# define a function to let user input splitting options
def splitting_choice():
    global splitting_window, x_split_entry, y_split_entry
    # start a new window 
    splitting_window = tk.Toplevel()
    splitting_window.title('Choose splitting in 2 form stacks')
    splitting_window.geometry('400x300')
    # info, that factors must add up to 1 to return same 2 form
    text_warning_split = '''Both MUST add to 1 to retain 2 form \n
    or must combine with 2 form to give the same expression \n
    otherwise, the plotted 2 form will be different!
                        '''
    tk.Label(splitting_window, text=text_warning_split).grid(row=0, column=0, pady=10)
    # define entry box for splitting factor that goes to x:
    tk.Label(splitting_window, text='factor going toward x').grid(row=2, column=0)
    x_split_entry = tk.Entry(splitting_window, width=40, borderwidth=1)
    x_split_entry.insert(0, x_split_str)
    x_split_entry.grid(row=3, column=0)
    # define entry box for splitting factor that goes to y:
    tk.Label(splitting_window, text='factor going toward y').grid(row=4, column=0)
    y_split_entry = tk.Entry(splitting_window, width=40, borderwidth=1)
    y_split_entry.insert(0, y_split_str)
    y_split_entry.grid(row=5, column=0)
    # define a button to SAVE these
    Splitting_btn = tk.Button(splitting_window, text='SAVE', padx=20, pady=10, command=splitting_choice_save)
    Splitting_btn.grid(row=6, column=0, pady=10)
    

'''

Define GUI interactions

'''

# define a window to supply the 2 form
tk.Label(bot_frame, text='2 form on R2').grid(row=0, column=1)
form_2_entry = tk.Entry(bot_frame, width=20, borderwidth=2)
form_2_entry.grid(row=0, column=0)
form_2_entry.insert(0, form_2_str)
# begin  displaying it on green colour to show that this is ebing displayed to
# beign with
form_2_entry.configure(bg='#C0F6BB')

# define entries for a 1 form
tk.Label(bot_frame, text='x component of 1 form').grid(row=1, column=1)
form_1_x_entry = tk.Entry(bot_frame, width=20, borderwidth=2)
form_1_x_entry.grid(row=1, column=0)
form_1_x_entry.insert(0, string_x)

tk.Label(bot_frame, text='y component of 1 form').grid(row=2, column=1)
form_1_y_entry = tk.Entry(bot_frame, width=20, borderwidth=2)
form_1_y_entry.grid(row=2, column=0)
form_1_y_entry.insert(0, string_y)

# extra label for 0 form.
# for now not an entry box because it doesn't get used anywhere
# I make it red to remember to come back to it later and use it.
Label_zero_form = tk.Label(bot_frame, text='', fg='red')
Label_zero_form.grid(row=3, column=0)

# define a button to submit the supplied 2 form and plot it as blocks
form_2_btn = tk.Button(small_frame, text='2 form plot', padx=3, pady=5, command=form_2_response)
form_2_btn.grid(row=0, column=0)

# define radio buttons to chose stacks or blocks
stack_block = tk.IntVar()  # Tkinter variable for Radiobuttons
stack_block.set(0)
stack_block_int = stack_block.get()
blocks_rb = tk.Radiobutton(right_frame, text='blocks', variable=stack_block, value=0, command=lambda: plot_type_response(stack_block.get()))
blocks_rb.grid(row=0, column=1)
stacks_rb = tk.Radiobutton(right_frame, text='stacks', variable=stack_block, value=1, command=lambda: plot_type_response(stack_block.get()))
stacks_rb.grid(row=0, column=2)


# define a button that will let the user chose the splitting option
# for 2 forms plotted as stacks.
splitting_opts_btn = tk.Button(right_frame, text='choose splitting', command=splitting_choice)
splitting_opts_btn.grid(row=0, column=3, padx=10)

# define a button that will just plot the 1 form
# this will not be needed when its it merged with main GUI
# as there will already be a plot button there
form_1_stacks_btn = tk.Button(small_frame, text='1 form plot', padx=3, pady=5, command=form_1_stacks_response)
form_1_stacks_btn.grid(row=2, column=0)

# add a button to plot the interior derivative as superposing stack fields
INT_btn = tk.Button(right_frame, text='Int Deriv', padx=63, pady=10, command=Int_deriv_response)
INT_btn.grid(row=1, column=0)

# define a button to plot the exterior derivative from given u and v
# Note, it will get the 2 form first, then return back down
# to a one form to avoid concellations
# therefore, it will also just be one possible representation
# not the only possible one
EXT_int_btn = tk.Button(right_frame, text='Ext Deriv', padx=62, pady=10, command=Ext_deriv_response)
EXT_int_btn.grid(row=2, column=0)

# define a wedge product button that will let the user input TWO 1 forms
# in a new window to be wedged to gice a 2 form
wedge_btn = tk.Button(right_frame, text='wedge two 1 forms', padx=27, pady=10, command=wedge_2_response)
wedge_btn.grid(row=3, column=0)

# define ab utton that will Find the Hodge dual
Hodge_btn = tk.Button(right_frame, text='Hodge', padx=67, pady=10, command=Hodge_response)
Hodge_btn.grid(row=4, column=0)



# SOME customisations
# extra one sthat are not yet in main GUI (plotting_stacks_GUI)
# therefore will need to be added when the 2 are merged.

# Include ocustomisation buttonsplay options to change
# will show a window and di
# stack size and block size, stack size is already avaliable in main GUI
# but ont blocvk size
customise_calc_btn = tk.Button(small_frame, text='customise', command=customise_calc)
customise_calc_btn.grid(row=3, column=1)

'''

NOTE: I am not putting in customisations here as this will later (hopefully)
Join the main GUI which already has all that in it.

'''


'''

To integrate the components after splitting and get the original 1 forms
 need to inegrate each, then, need this:

# format it to be of correct shape for the other needed operations
form_2_eq = format_eq(form_2_str)
# check against constant and zero 2 forms being supplied
form_2_eq = form_2_constant_correction(form_2_eq)
# find its numerical evaulation
form_2 = eval(form_2_eq)
# get the signs of the 2 form
form_2_sgn = np.sign(form_2)
# if these came from exterior derivative. One was differentiated
# w.r.t x and one with y, undo this by integrating such.
# This is not unique as any split may be possible to original 2 form
# and cancellations may have occured after differentating,
# when merging dx^dy and dy^dx components.
# for this need to define symbols for sympy to use
x_symbol = Symbol('x')
y_symbol = Symbol('y')
comp_x = str(simplify( '-(' + str(integrate(eq_1, y_symbol)) + ')'))
comp_y = str(simplify(str(integrate(eq_2, x_symbol))))

'''

# return time to run
stop = timeit.default_timer()
print('Time: ', stop - start)

tk.mainloop()
