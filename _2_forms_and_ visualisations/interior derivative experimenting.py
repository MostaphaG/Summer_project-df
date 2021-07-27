# Attempt 1 - interior derivative calculation and plot on R2

# import modules
import timeit
import numpy as np
import matplotlib.pyplot as plt
from sympy import diff
from sympy.parsing.sympy_parser import parse_expr
from matplotlib.backend_bases import key_press_handler
from matplotlib.lines import Line2D
import tkinter as tk
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib import patches as patch


# %%

# start the timer
start = timeit.default_timer()

'''

Define initially needed functions

'''


# function make variables python understood
def format_eq(string):
    # replace all the x and y with xg and yg:
    string = string.replace('x', 'xg')
    string = string.replace('y', 'yg')
    string = string.replace('z', 'zg')
    string = string.replace('R', 'rg')  # otherwise: sqrt won't work, becaue of the r in it &arctan won't work because of tan in it and the r in it.
    string = string.replace('theta', 'thetag')
    # where there are special functions, replace them with library directions
    string = string.replace('pi', 'np.pi')
    string = string.replace('sqrt', 'np.sqrt')
    string = string.replace('sin', 'np.sin')
    string = string.replace('cos', 'np.cos')
    string = string.replace('tan', 'np.tan')
    string = string.replace('arcta', 'np.arctan2')
    string = string.replace('^', '**')
    string = string.replace('ln', 'np.log')
    string = string.replace('e^', 'np.exp')
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
        u = float(equation_x)*np.ones(np.shape(xg))
    if equation_y.find('x') & equation_y.find('y') == -1:
        v = float(equation_y)*np.ones(np.shape(yg))
    # return these
    return u, v


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
def form_2_components_plot(xg, yg, u, v, s_max, L, pt_den, fract, colour_str):
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
                 
                    ax.add_line(Line2D((Ax1,Bx1),(Ay1,By1), linewidth=0.7, color=colour_str[color_index]))
                    ax.add_line(Line2D((Ax2,Bx2),(Ay2,By2), linewidth=0.7, color=colour_str[color_index]))
                    
                    # change the parameter to loop over all changes in displacement for current magnitude
                    s += 1


# define a function that will complete all stack plotting:
# needed after using interior derivative on the 2 form
# as that returns a 1 form - vector field.
def stack_plot(xg, yg, u, v, s_max, L, pt_den, fract, arrows='True', orientation='mid', scale=1, w_head=1/8, h_head=1/4):
    # get the lengths of x and y from their grids
    x_len = len(xg[:, 0])
    y_len = len(yg[0, :])
    
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
                
    plt.close()


# define a function to plot the simplified 2 forms, with coloured squares
def plot_form(form_2):
    # celar the currently present plot
    ax.clear()
    # redefine the axis limits
    ax.set_xlim(-ax_L, ax_L)
    ax.set_ylim(-ax_L, ax_L)
    # find the maximum of the 2-form, over the grid
    max_val2 = np.max(form_2)
    # get an array of relative magnitudes:
    Mag = form_2/max_val2
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


# define a function that will find the 2 form from given expressions
# in a given number of dimensions and in terms of given coordinate symbols
# here it is only needed on R2 in terms of 'x' and 'y'.
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

Start the Tkinter window with basic layout

'''

# start a tkinter window
root = tk.Tk()

# set its title
root.title('form derivatives')

# set window size from screen resolution
width = root.winfo_screenwidth()
height = root.winfo_screenheight()
root.geometry(str(width) + 'x' + str(height))

# define a frame to put the plot into:
plot_frame = tk.LabelFrame(root, text='plot Frame', padx=5, pady=5)
plot_frame.grid(row=0, column=0)

# define a frame where all options will lie
right_frame = tk.LabelFrame(root, text='Options Frame', padx=20, pady=100)
right_frame.grid(row=0, column=1)

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


def on_key_press(event):
    print("you pressed {}".format(event.key))
    key_press_handler(event, canvas, toolbar)


# define a window to supply the 2 form
tk.Label(root, text='2 form on R2').grid(row=1, column=0)
form_2_entry = tk.Entry(root, width=20, borderwidth=2)
form_2_entry.grid(row=2, column=0)
form_2_entry.insert(0, 'sin(x*y)')

'''

Define needed parameters and expressions

'''


# define scale of the graph
L = 5
pt_den = 21   # number of points on each axis

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

# define the 2 form components to each elementary 2 form on R^m
string_1 = 'sin(x*y)'  # dx^dy component

# define the initial 2 form, in terms of a string, equation and numerically
form_2_str = string_1
form_2_eq = format_eq(form_2_str)
form_2 = eval(form_2_eq)

# get the signs of the 2 form
form_2_sgn = np.sign(form_2)

# put the initial plot onto the canvas
plot_form(form_2)
canvas.draw()

# #############################################################################
# define some extra parameters for stack plots
# #############################################################################

# define colours to use for 2 form plots with components
# first string defines colour for positive (ccw), second for negative (cw)
# last one is an in case, for when the magnitude is exactly zero.
colour_str = ['red', 'blue', 'grey']

# define fract of graph to set as size of stack
fract = 0.05

# define max number of sheets to start with
s_max = 5

# set up a zero vector filed to plot x and y components as 2 separate fields:
zero_field = np.zeros(np.shape(xg))

'''

define needed functions for responses

'''


# define a function that will respond to submitting the 2 form:
def Submit_form():
    global form_2_str, form_2_eq, form_2, form_2_sgn
    form_2_str = str(form_2_entry.get())
    # format it to be python understood
    form_2_eq = format_eq(form_2_str)
    # get the numerical evaluation of it
    form_2 = eval(form_2_eq)
    # get the signs of thsi new 2 form
    form_2_sgn = np.sign(form_2)
    # plot the new form using the previously define funtion
    plot_form(form_2)
    canvas.draw()


def Int_deriv_R2():
    global u, v, u_str, v_str
    # clear the already present axis
    ax.clear()
    # using interior product, get the u and v (dx and dy) components
    # of the resulting 1 form
    u = -form_2
    v = form_2
    # to be usable in ext_deriv, define strings of these variables
    u_str = '-' + form_2_str
    v_str = form_2_str
    # use the stacks plotter to present this
    stack_plot(xg, yg, u, v, s_max, L, pt_den, fract, False)
    canvas.draw()


def Ext_deriv_R2():
    global form_2_sgn
    # celar current axis
    ax.clear()
    # from found u and v in the interior derivative, set up sympy components
    sympy_expr_x = parse_expr(u_str, evaluate=False)
    sympy_expr_y = parse_expr(v_str, evaluate=False)
    # combine the 2 into a list:
    expressions = np.array([sympy_expr_x, sympy_expr_y])
    # set up an array of coordinates that need to be used (in standard order)
    coords = ['x', 'y']
    # from these, use the find_2_form function to get the 2 form
    form_2 = find_2_form(expressions, coords, m)
    # get the signs of the 2 form
    form_2_sgn = np.sign(form_2[0])
    # evaluate the u and v given previously, with formating for them to be
    # python understood
    u, v = eq_to_comps(str(expressions[0]), str(expressions[1]), xg, yg)
    # plot the fields with desired parameters as specidfied above
    # arrow params not needed as arrows arent plotted
    form_2_components_plot(xg, yg, u, zero_field, s_max, L, pt_den, fract, colour_str)
    form_2_components_plot(xg, yg, zero_field, v, s_max, L, pt_den, fract, colour_str)
    # these fields are the original components
    canvas.draw()


'''

More Tkinter options, ones that need responses to functions

'''


# define a button to submit the supplied 2 form
SUBMIT_btn = tk.Button(root, text='SUBMIT', padx=60, pady=30, command=Submit_form)
SUBMIT_btn.grid(row=3, column=0)

# add a button to plot the interior derivative as superposing stack fields
INT_btn = tk.Button(right_frame, text='Int Deriv', padx=10, pady=10, command=Int_deriv_R2)
INT_btn.grid(row=0, column=0)

# define a button to plot the exterior derivative from given u and v
# in this case, given by first evaluating the interior derivative
# of the supplied 2 form
# it will do nothing when pressed first in this code
# NOTE - not the case in the general GUI code
EXT_btn = tk.Button(right_frame, text='Ext Deriv', padx=10, pady=10, command=Ext_deriv_R2)
EXT_btn.grid(row=1, column=0)


# return time to run
stop = timeit.default_timer()
print('Time to run: ', stop - start)

# keep the program running until otherwise specified by user
tk.mainloop()
