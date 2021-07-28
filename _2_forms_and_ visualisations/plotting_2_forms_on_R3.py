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

# %%

# start the timer
start = timeit.default_timer()

# start a tkinter window
root = tk.Tk()

# set its title
root.title('2 forms on R3')

# set window size
root.geometry('1000x800')

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
    return string


# define a function that takes input string that is python understood and turn into vector components:
def eq_to_comps(equation_x, equation_y, equation_z, xg, yg, zg):
    # use this fucntion to replace given string to python understood equations:
    equation_x = format_eq(string_x)
    equation_y = format_eq(string_y)
    equation_z = format_eq(string_z)
    # use these to define the field:
    # also: checking if equation equals zero, to then replace it with an array and not just 0:
    F_x = eval(equation_x)
    F_y = eval(equation_y)
    F_z = eval(equation_z)
    if equation_x.find('x') & equation_x.find('y') & equation_x.find('z') == -1:
        F_x = float(equation_x)*np.ones(np.shape(xg))
    if equation_y.find('x') & equation_y.find('y') & equation_y.find('z') == -1:
        F_y = float(equation_y)*np.ones(np.shape(yg))
    if equation_z.find('x') & equation_z.find('y') & equation_z.find('z') == -1:
        F_z = float(equation_z)*np.ones(np.shape(zg))
    # return these
    return F_x, F_y, F_z


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
            if (temp == ' + 0') or (temp == ' - (0)') or (temp == '0*x'):
                pass
            else:
                result[pair, 0] += temp
            if (temp1 == ' + 0') or (temp1 == ' - (0)') or (temp1 == '0*x'):
                pass
            else:
                result[pair, 0] += temp1
            # update the result row counter
            pair += 1
    # format string in each result row
    # making sure to format it correctly even if it contains constants or '0'
    # this is done if result is to be directly used later for any reason
    # instead of form_2 numerically.
    for d in range(pair):
        if result[d, 0].find('x') & result[d, 0].find('y') & result[d, 0].find('z') == -1:
            if result[d, 0] == '':
                result[d, 0] = '0*x'
                result[d, 0] = format_eq(result[d, 0])
            else:
                result[d, 0] = '(' + result[d, 0] + ')*x'
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
        u = u[:, h_index, :]
        v = v[:, h_index, :]
    elif axis_view == 'x':
        grid_x = grid_x[:, h_index, :]
        grid_y = grid_y[:, h_index, :]
        u = u[h_index, :, :]
        v = v[h_index, :, :]
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

'''

Set up all needed parameters, plots etc

'''

# define parameters of the axis
L = 5
pt_den = 17  # number of points on each axis

# define x, y and z values
x = np.linspace(-L, L, pt_den)
y = np.linspace(-L, L, pt_den)
z = np.linspace(-L, L, pt_den)

# create the grids
xg, yg, zg = np.meshgrid(x, y, z)

# define the wanted 1 form on R3 in terms of each component:
string_x = 'x*sin(y)'  # x component
string_y = 'y*cos(x)'  # y component
string_z = '0'  # z component

# to start with, set as viewing aling z axis onto x-y plane
axis_view = 'z'

# set the initial index of height of the viewed plane along viewing axis
h_index = 11

# set the dimensionality
m = 3

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

# because on R3, the x, y and z components are made up of
# dy^dz, dx^dz and dx^dy
# need to superpose stacks from the x y and z fields (depending on viewing)
# the filed will point in the extra dimension, or by polar sense
# will be defined in plane based on orientation of stacks. Indicated by colour
F_x, F_y, F_z = eq_to_comps(str(expressions[0]), str(expressions[1]), str(expressions[2]), xg, yg, zg)

# Note that that order on RHS is not standard because result gives
# elemental 2 forms in order: dx^dy, dx^dz and dy^dz.
# so the z component is first, then y and then x.

# set up a zero vector filed to plot x and y components as 2 separate fields:
zero_field = np.zeros(np.shape(xg))

# set up the delta_factor of additional axis space L/delta_factor gives extra space on axis
delta_factor = 10

# fraction of sheet length to graph length
fract = 0.05

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

# plot the starting field with desired parameters as specidfied above
# arrow params not needed as arrows arent plotted
form_2_components_plot_3(xg, yg, h_index, axis_view, F_x, zero_field, s_max, L, pt_den, fract, colour_str)
form_2_components_plot_3(xg, yg, h_index, axis_view, zero_field, F_y, s_max, L, pt_den, fract, colour_str)

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


# define a function to update z Index and redraw the plot based on the slider
def slide(var):
    global h_index
    # extract current value from slider
    h_index = slider_z.get()
    # remove the currently displayed plot
    ax.clear()
    # replot the graph with that new h_index:
    if axis_view == 'z':
        form_2_components_plot_3(xg, yg, h_index, axis_view, F_x, zero_field, s_max, L, pt_den, fract, colour_str)
        form_2_components_plot_3(xg, yg, h_index, axis_view, zero_field, F_y, s_max, L, pt_den, fract, colour_str)
        ax.set_xlabel('$x$')
        ax.set_ylabel('$y$')
    elif axis_view == 'y':
        form_2_components_plot_3(xg, zg, h_index, axis_view, F_x, zero_field, s_max, L, pt_den, fract, colour_str)
        form_2_components_plot_3(xg, zg, h_index, axis_view, zero_field, F_z, s_max, L, pt_den, fract, colour_str)
        ax.set_xlabel('$x$')
        ax.set_ylabel('$z$')
    elif axis_view == 'x':
        form_2_components_plot_3(yg, zg, h_index, axis_view, F_y, zero_field, s_max, L, pt_den, fract, colour_str)
        form_2_components_plot_3(yg, zg, h_index, axis_view, zero_field, F_z, s_max, L, pt_den, fract, colour_str)
        ax.set_xlabel('$y$')
        ax.set_ylabel('$z$')
    # draw that onto the screen
    canvas.draw()


# Label a slider
tk.Label(root, text='index of plane along the viewing axis').grid(row=1, column=0)
# define a slider to update h_index
slider_z = tk.Scale(root, from_ = 0, to=len(z)-1, orient=tk.HORIZONTAL)

# bind the button to an event of releasing the mouse
slider_z.bind("<ButtonRelease-1>", slide)
# updating in real time is too slow, but a button is clumsy
# now the slider will update to a value that the mouse was released at

# set the initial value and put the slider on the screen
slider_z.set(z[0])
slider_z.grid(row=2, column=0)


# define a function that will repond to changing axis view with radiobuttons
def view_response(view_var):
    # get and make global the axis view variable
    global axis_view
    axis_view = view_tk.get()
    # clear the current plot
    ax.clear()
    # draw the new plots, depending on the chosen viewing direction
    if axis_view == 'z':
        form_2_components_plot_3(xg, yg, h_index, axis_view, F_x, zero_field, s_max, L, pt_den, fract, colour_str)
        form_2_components_plot_3(xg, yg, h_index, axis_view, zero_field, F_y, s_max, L, pt_den, fract, colour_str)
        ax.set_xlabel('$x$')
        ax.set_ylabel('$y$')
    elif axis_view == 'y':
        form_2_components_plot_3(xg, zg, h_index, axis_view, F_x, zero_field, s_max, L, pt_den, fract, colour_str)
        form_2_components_plot_3(xg, zg, h_index, axis_view, zero_field, F_z, s_max, L, pt_den, fract, colour_str)
        ax.set_xlabel('$x$')
        ax.set_ylabel('$z$')
    elif axis_view == 'x':
        form_2_components_plot_3(yg, zg, h_index, axis_view, F_y, zero_field, s_max, L, pt_den, fract, colour_str)
        form_2_components_plot_3(yg, zg, h_index, axis_view, zero_field, F_z, s_max, L, pt_den, fract, colour_str)
        ax.set_xlabel('$y$')
        ax.set_ylabel('$z$')
    # draw that onto the screen
    canvas.draw()


# define rediobuttons to chose from which axis the user is looking:
view_tk = tk.StringVar()
view_tk.set('z')
view_z_btn = tk.Radiobutton(root, text='z', variable=view_tk, value='z', command=lambda: view_response(view_tk.get())).grid(row=1, column=1)
view_z_btn = tk.Radiobutton(root, text='y', variable=view_tk, value='y', command=lambda: view_response(view_tk.get())).grid(row=2, column=1)
view_z_btn = tk.Radiobutton(root, text='x', variable=view_tk, value='x', command=lambda: view_response(view_tk.get())).grid(row=3, column=1)


# return time to run
stop = timeit.default_timer()
print('Time to run: ', stop - start)

# keep the program running until otherwise specified by user
tk.mainloop()
