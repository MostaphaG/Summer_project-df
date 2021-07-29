
#%% Import Modules
import timeit
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
import matplotlib as mpl
from sympy import diff
from sympy.parsing.sympy_parser import parse_expr

# %% VFA GUI

# start the timer
start = timeit.default_timer()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define needed functions for the initial plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
def stack_plot(xg, yg, axis, u, v, s_max, L, pt_den, fract, arrows=False, orientation='mid', scale=1, w_head=1/8, h_head=1/4, axis_check=0):
    # get the lengths of x and y from their grids
    
    x_len = len(xg[:, 0])
    y_len = len(yg[0, :])
    
    # Scaling of axes and setting equal proportions circles look like circles
    # axis.set_aspect('equal')
    # ax_L = L + L/delta_factor
    # axis.set_xlim(-ax_L, ax_L)
    # axis.set_ylim(-ax_L, ax_L)
    
    # Account for change to grid centre for divergence plot
    if axis_check == 1:
        if click_opt_int > 2:       
            axis.set_xlim(-L-L/5, L+L/5)
            axis.set_ylim(-L-L/5, L+L/5)
        else:
            axis.set_xlim(-L+x_m-L/5, L+x_m+L/5)
            axis.set_ylim(-L+y_m-L/5, L+y_m+L/5)
    else:
        # redefine axis limits here, as: before plotting a new figure, axis
        # are cleared, it will rescale itself back to (0, 0), (1, 1) otherwise
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
            
            if axis_check == 1 and click_opt_int > 1 and i == i_m and j == j_m:
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
                    axis.add_line(Line2D((Ax1, Bx1), (Ay1, By1), linewidth=1, color='green'))
                    axis.add_line(Line2D((Ax2, Bx2), (Ay2, By2), linewidth=1, color='green'))
                    
                    # update parameter to reapet and draw all needed arrows
                    s += 1
            # deal with the odd number of stacks:
            elif parity(n) is False:
                # Add the centre line for odd numbers of stacks
                axis.add_line(Line2D((A_x[i, j], B_x[i, j]), (A_y[i, j], B_y[i, j]), linewidth=1, color='green'))
                
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
                    axis.add_line(Line2D((Ax1,Bx1),(Ay1,By1), linewidth=1, color='green'))
                    axis.add_line(Line2D((Ax2,Bx2),(Ay2,By2), linewidth=1, color='green'))
                    
                    # change the parameter to loop over all changes in displacement for current magnitude
                    s += 1
            # plot lines of arrowheads from central sheet for n = 1 or on top sheet for n>1 
            if n > 1:   # for all lines ubt the single sheet one
                axis.add_line(Line2D((p_sh1x[i, j],p_sh3x[i, j]),(p_sh1y[i, j],p_sh3y[i, j]), linewidth=1, color='green'))
                axis.add_line(Line2D((p_sh2x[i, j],p_sh3x[i, j]),((p_sh2y[i, j],p_sh3y[i, j])), linewidth=1, color='green'))
            # then define it for the stacks with only 1 sheet:
            else:
                axis.add_line(Line2D((P_sh1x[i, j], P_sh3x[i, j]), (P_sh1y[i, j], P_sh3y[i, j]), linewidth=1, color='green'))
                axis.add_line(Line2D((P_sh2x[i, j], P_sh3x[i, j]), ((P_sh2y[i, j], P_sh3y[i, j])), linewidth=1, color='green'))
                
    plt.close()


# define a function to replace input string to be 'python understood'
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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set up basic layout of the window
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# define an object tracker - for GUI set up
root = tk.Tk()

# set its title
root.title('Vector field analyser - differential forms')

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
right_frame = tk.LabelFrame(root, text='Options Frame', padx=5, pady=5)
right_frame.grid(row=1, column=1)

# bot frame:
bot_frame = tk.LabelFrame(root, text='Field Input Frame', padx=100, pady=5)
bot_frame.grid(row=2, column=0)

# plot frame:
plot_frame = tk.LabelFrame(root, text='Vector Field Frame', padx=5, pady=5)
plot_frame.grid(row=1, column=0)

# plot characteristics frame and plot button
small_frame = tk.LabelFrame(root, text='Plot Customisation Frame', padx=35, pady=5)
small_frame.grid(row=2, column=1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set up initial parameters and plot the initial graph, put it in plot frame
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# define scale of the graph
L = 5
pt_den = 10   # number of points on each axis

# define x and y values
x = np.linspace(-L, L, pt_den)
y = np.linspace(-L, L, pt_den)

# create a grid on x-y plane
xg, yg = np.meshgrid(x, y)

# define an example vector field
u = yg*np.sin(xg)  # x component
v = -xg*np.cos(yg)  # y component
# for no dependance in any initial component, use : np.zeros(np.shape(xg)) or yg

'''
SET UP A LIST OF DEFAULT VECTOR FIELDS TO DISPLAY IN DROPDOWN MENU
'''
# list of names of fields to display
field_name_list = ['Default: y*sin(x)i - x*cos(y)j',
              'Simple pendulum: yi  - sin(x)j',
              'Harmonic oscillator: yi -xj',
              'Linear field example 1: (14*x - 4*y)i + (-1*x + 4*y)j',
              'Linear field example 2: xi',
              'Constant field: 6i + 3j',
              'Falling cat field (Planar 3 link robot)',
              'Gravitational/Electric Point Charge: -x/(x**2+y**2)i + -y/(x**2+y**2)j',
              'Magnetic Field of Current Carrying Wire: -y/(x**2+y**2)i + x/(x**2+y**2)j'
              ]


# list of x components, in order of field_name_list
field_x_list = ['y*sin(x)',
                'y',
                'y',
                '14*x - 4*y',
                'x',
                '6',
                '(3*cos(y) + 4)/(15 + 6*cos(x) + 6*cos(y))',
                '-x/(x**2+y**2)*step(4*(x**2+y**2))',
                '-y/(x**2+y**2)*step(4*(x**2+y**2))'
                ]


# list of y components, in order of field_name_list
field_y_list = ['- x*cos(y)',
                '-sin(x)',
                '-x',
                '(-1*x + 4*y)',
                '0',
                '3',
                '-(3*cos(x) + 4)/(15 + 6*cos(x) + 6*cos(y))',
                '-y/(x**2+y**2)*step(4*(x**2+y**2))',
                'x/(x**2+y**2)*step(4*(x**2+y**2))'
                ]


# set up quiver factors
arrows = False  # set up if arrows should be plotted on stacks or not.
orientation = 'mid'  # how the arrow rotates about its assigned grid point - options: tail, mid and tip as string
scale = 10  # the scale reduction factor, if None (as None-type, not str), automatically computed by average, if 1 = mag

# set up the delta_factor of additional axis space L/delta_factor gives extra space on axis
delta_factor = 10

# fraction of sheet length to graph length
fract = 0.05

# define the maximum number of stack to plot, dep. on magnitude (initialy)
s_max = 2

# set screen dpi
my_dpi = 100

# define denominator of fractional height and width of arrowhead based on stack size
w_head = 1/8
h_head = 1/4

# create a figure, use dpi to fit it more precisely to size of the frame
fig = plt.figure(figsize=(855/my_dpi, 573/my_dpi), dpi=my_dpi)

# set up axis
main_axis = fig.gca()

# set up visuals - axis labels
main_axis.set_xlabel('$x$')
main_axis.set_ylabel('$y$')

# Scaling of axes and setting equal proportions circles look like circles
main_axis.set_aspect('equal')
ax_L = L + L/delta_factor
main_axis.set_xlim(-ax_L, ax_L)
main_axis.set_ylim(-ax_L, ax_L)


'''define the initial polar plot also. Despite it not being plotted to start
with needed for when Polar plot option is used'''

r_max = 5
r_den = 10
theta_den = 20

# define the axis
r = np.linspace(0.2, r_max, r_den)
theta = np.linspace(0, 360, theta_den) * np.pi/180

# define a polar grid
rg, thetag = np.meshgrid(r, theta)

# define an initial linear scaling of the poalr field to be able to scale in
# poalr gird without converting to cartesian grid first
a_polar = 1

# define an initial polar field (same as initial cartesian field.)
F_r_str_initial = 'r*sin(theta)*sin(r*cos(theta))'
F_theta_str_initial = '-r*cos(theta)*cos(r*sin(theta))'

# from the unformated string called 'initial', set up a python understood string
F_r_str = format_eq(F_r_str_initial)
F_theta_str = format_eq(F_theta_str_initial)

# using these, evaluate the initial radial and angular components of the field.
F_r = eval(F_r_str)
F_theta = eval(F_theta_str)

# convert the field back to cartesian
u_p = F_r*np.cos(thetag) - F_theta*np.sin(thetag)  # x component
v_p = F_r*np.sin(thetag) + F_theta*np.cos(thetag)  # y component

''' end of polar setting up'''

# plot the cartessian field with desired parameters as specidfied above
plottedfield = stack_plot(xg, yg, main_axis, u, v, s_max, L, pt_den, fract, arrows, orientation, scale, w_head, h_head, 0)

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
    # respond with only toolbar actions when only tools are to be used
    if click_opt_int == 0:
        print("you pressed {}".format(event.key))
        key_press_handler(event, canvas, toolbar)
    # when the derivative option is selected, cerry out the derivative when clicked
    else:
        global x_pix, y_pix, x_m, y_m
        # get cartesian Coordinates of click
        ix_plot, iy_plot = event.xdata, event.ydata
        # from those, get pixel Coordinates of click
        x_pix, y_pix = event.x , event.y
        x_m = float(ix_plot)
        y_m = float(iy_plot)
        mpl.rcParams['toolbar'] = 'None'  # this does not do what it should yet
        deriv_calc(x_m, y_m)
        
# connect figure event to a function that responds to clicks, defined above
fig.canvas.mpl_connect("button_press_event", on_key_press)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define other needed functions, for input reponses
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# define a function that will change r and theta components to equivalents in x and y
# this will be needed to display the filed in cartesian after submitting a polar field
# so that the PLOT button does not break.
def p_cart(string):
    string= string.replace('R', 'r')
    string = string.replace('r', 'sqrt(x**2 + y**2)')
    string = string.replace('theta', 'arcta(y, x)')
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


# define a function that will respond to radio buttons behind choosing vector types:
def vect_type_response(tensor):
    # clear the plot that is already there:
    main_axis.clear()
    # use the tensor to determine what to plot:
    # 0 is just stacks, 1 is for only arrows and 2 is for both
    if tensor == 0:
        arrows = False
        stack_plot(xg, yg, main_axis, u, v, s_max, L, pt_den, fract, arrows, orientation, scale, w_head, h_head, 0)
        canvas.draw()
    elif tensor == 1:
        main_axis.quiver(xg, yg, u, v, pivot=orientation, scale=scale, scale_units='xy')
        # repeat the displaying of the figure so that it updates in GUI
        canvas.draw()
    elif tensor == 2:
        arrows = True
        stack_plot(xg, yg, main_axis, u, v, s_max, L, pt_den, fract, arrows, orientation, scale, w_head, h_head, 0)
        # repeat the displaying of the figure so that it updates in GUI
        canvas.draw()


# define the PLOT button response function
def PLOT_response():
    # first, take from entry boxes, wanted parameters and make them global:
    # these must be changed globally for other functions to work with the new field.
    global L, pt_den, s_max, x, y, xg, yg, u, v, tensor, main_axis, string_x, string_y
    # clear the current axis
    main_axis.clear()
    # take the new axis parameters and field definitions out of the boxes
    L = float(L_entry.get())
    pt_den = int(pt_den_entry.get())
    s_max = int(s_max_entry.get())
    string_x = str(x_comp_entry.get())
    string_y = str(y_comp_entry.get())
    # from L redefine the axis
    ax_L = L + L/delta_factor
    main_axis.set_xlim(-ax_L, ax_L)
    main_axis.set_ylim(-ax_L, ax_L)
    # from pt_den and L, change the axis coordinates and the grid:
    x = np.linspace(-L, L, pt_den)
    y = np.linspace(-L, L, pt_den)
    xg, yg = np.meshgrid(x, y)
    # take all these values, and the input from field component bnoxes to set up the field:
    u, v = eq_to_comps(string_x, string_y, xg, yg)
    # plot the new field
    stack_plot(xg, yg, main_axis, u, v, s_max, L, pt_den, fract, arrows, orientation, scale, w_head, h_head, 0)
    # put it onto the screen
    canvas.draw()
    # change the radio button ticks back to stack only
    tensor.set(0)


# define a function to respond to submitting arrohead changes in the new window
def custom_submission():
    # first, take from entry boxes, wanted parameters and make them global:
    global w_head, h_head, fract, scale
    w_head = float(w_entry.get())
    h_head = float(h_entry.get())
    fract = float(fract_entry.get())
    scale =  float(arr_scale_entry.get())
    # depending on the Radio buttons, replot the graph and put it onto the GUI
    vect_type_response(tensor.get())
    # then close the window
    arrowH_opt_window.destroy()


# define a reponse function to open a new window when arrowh_btn is pressed:
def custom_btn_reponse():
    global w_entry, h_entry, fract_entry, arr_scale_entry, arrowH_opt_window
    # open a titled new window
    arrowH_opt_window = tk.Toplevel()
    arrowH_opt_window.title('optimisation settings')
    # define and label and first entry, for width
    tk.Label(arrowH_opt_window, text='arrowhead base width as sheet width fraction:').grid(row=0, column=0)
    w_entry = tk.Entry(arrowH_opt_window, width=30, borderwidth=1)
    w_entry.insert(0, w_head)
    w_entry.grid(row=1, column=0)
    # define and label second entry, for height
    tk.Label(arrowH_opt_window, text='arrowhead perp. height as sheet length fraction:').grid(row=2, column=0)
    h_entry = tk.Entry(arrowH_opt_window, width=30, borderwidth=1)
    h_entry.insert(0, h_head)
    h_entry.grid(row=3, column=0)
    # define an entry for fract update, to change the size of each stack as a frac of graph size L
    tk.Label(arrowH_opt_window, text='fraction of graph to be set as the stack size:').grid(row=4, column=0)
    fract_entry = tk.Entry(arrowH_opt_window, width=30, borderwidth=1)
    fract_entry.insert(0, fract)
    fract_entry.grid(row=5, column=0)
    # define an entry for fract update, to change the size of each stack as a frac of graph size L
    tk.Label(arrowH_opt_window, text='arrow size linear scaling:').grid(row=6, column=0)
    arr_scale_entry = tk.Entry(arrowH_opt_window, width=30, borderwidth=1)
    arr_scale_entry.insert(0, scale)
    arr_scale_entry.grid(row=7, column=0)
    # define a button to submit those changes:
    submit_arr_btn = tk.Button(arrowH_opt_window, text='SUBMIT ALL', padx=20, pady=10, command=custom_submission)
    submit_arr_btn.grid(row=8, column=0, pady=10)


def polar_submit(tensorp):
    # take the input values into new variables
    global F_r, F_theta, r_max, r_den, theta_den, r, theta, rg, thetag, u_p, v_p, a_polar, L, F_r_str_initial, F_theta_str_initial
    F_r = Fr_entry.get()
    F_theta = Ftheta_entry.get()
    r_max = float(r_max_entry.get())
    r_den = int(r_den_entry.get())
    theta_den = int(theta_den_entry.get())
    L = float(L_entry1.get())
    a_polar = float(a_entry.get())
    # rescale the axis:
    ax_L = L + L/delta_factor
    main_axis.set_xlim(-ax_L, ax_L)
    main_axis.set_ylim(-ax_L, ax_L)
    # based on these, change the axis coordinates
    r = np.linspace(0.2, r_max, r_den)
    theta = np.linspace(0, 360, theta_den) * np.pi/180
    # and mesh the new ones
    rg, thetag = np.meshgrid(r, theta)
    # update intial strings to display
    F_r_str_initial = F_r
    F_theta_str_initial = F_theta
    # format radial and angular components to be python understood
    F_r_str = format_eq(F_r)
    F_theta_str = format_eq(F_theta)
    # evalueate them, bearing in mind the process to follow when user inputs 0
    if F_r_str == '0':
        F_r = np.zeros(np.shape(rg))
        F_theta = eval(F_theta_str)
    elif F_theta_str == '0':
        F_r = eval(F_r_str)
        F_theta = np.zeros(np.shape(thetag))
    elif F_r_str and F_theta_str == '0':
        F_r = np.zeros(np.shape(rg))
        F_theta = np.zeros(np.shape(thetag))
    else:
        F_r = eval(F_r_str)
        F_theta = eval(F_theta_str)
    # redefine the field with these (in cartesian)
    u_p = F_r*np.cos(thetag) - F_theta*np.sin(thetag)  # x component
    v_p = F_r*np.sin(thetag) + F_theta*np.cos(thetag)  # y component
    # scale the field with given a:
    u_p *= a_polar
    v_p *= a_polar
    # convert the cooridnates rg and thetag back to cartesian, needed for the plotting functions
    xg = rg*np.cos(thetag)
    yg = rg*np.sin(thetag)
    # clear the plot that is already there:
    main_axis.clear()
    # use the selected tensor to determine what to plot:
    # 0 is just stacks, 1 is for only arrows and 2 is for both
    if tensorp == 0:
        arrows = False  # set correct variable as asked by user
        stack_plot(xg, yg, main_axis, u_p, v_p, s_max, L, pt_den, fract, arrows, orientation, scale, w_head, h_head, 0)  # plot
        # display the figure so that it updates in GUI
        canvas.draw()
    elif tensorp == 1:
        main_axis.quiver(xg, yg, u_p, v_p, pivot=orientation, scale=scale, scale_units='xy')
        canvas.draw()
    elif tensorp == 2:
        arrows = True
        stack_plot(xg, yg, main_axis, u_p, v_p, s_max, L, pt_den, fract, arrows, orientation, scale, w_head, h_head, 0)
        # repeat the displaying of the figure so that it updates in GUI
        canvas.draw()
    # display in x and y components, in a somewhat messy form, straight from poalr to cart
    x_comp_entry.delete(0, 'end')
    y_comp_entry.delete(0, 'end')
    x_comp_entry.insert(0, p_cart(F_r_str_initial) + p_cart(' * cos(theta)') + ' - ' + p_cart(F_theta_str_initial) + p_cart(' * sin(theta)'))
    y_comp_entry.insert(0, p_cart(F_r_str_initial) + p_cart(' * sin(theta)') + ' + ' + p_cart(F_theta_str_initial) + p_cart(' * cos(theta)'))
    ''' NOTE - the cartesian field will not be quite the same because range issues in arctan '''
    # then close the window
    polar_fld_window.destroy()


# define a function that responds to the polar button, to allow user to input
# details about the polar field they fish to plot
def Polar_btn_response():
    # need these global to pass them onto the function that responds to submission
    global polar_fld_window, Fr_entry, Ftheta_entry, r_max_entry, r_den_entry, theta_den_entry, L_entry1, a_entry
    global p_arr_btn, p_both_btn, p_stack_btn
    # open a window with input fields to enter polar components
    polar_fld_window = tk.Toplevel()
    polar_fld_window.title('polar field input')
    # define and label and first entry, for radial
    tk.Label(polar_fld_window, text='radial component in terms of \'R\':').grid(row=0, column=0)
    Fr_entry = tk.Entry(polar_fld_window, width=30, borderwidth=1)
    Fr_entry.insert(0, str(F_r_str_initial))
    Fr_entry.grid(row=1, column=0)
    # define and label second entry, for height
    tk.Label(polar_fld_window, text='angular component in terms of \'theta\' :').grid(row=2, column=0)
    Ftheta_entry = tk.Entry(polar_fld_window, width=30, borderwidth=1)
    Ftheta_entry.insert(0, str(F_theta_str_initial))
    Ftheta_entry.grid(row=3, column=0)
    # define an entry for radius limit r_max
    tk.Label(polar_fld_window, text='maximum radius:').grid(row=4, column=0)
    r_max_entry = tk.Entry(polar_fld_window, width=30, borderwidth=1)
    r_max_entry.insert(0, r_max)
    r_max_entry.grid(row=5, column=0)
    # define an entry for number of points along r
    tk.Label(polar_fld_window, text='number of points along r:').grid(row=6, column=0)
    r_den_entry = tk.Entry(polar_fld_window, width=30, borderwidth=1)
    r_den_entry.insert(0, r_den)
    r_den_entry.grid(row=7, column=0)
    # define an entry for number of points along theta
    tk.Label(polar_fld_window, text='number of points along theta:').grid(row=8, column=0)
    theta_den_entry = tk.Entry(polar_fld_window, width=30, borderwidth=1)
    theta_den_entry.insert(0, theta_den)
    theta_den_entry.grid(row=9, column=0)
    # define linear axis size
    tk.Label(polar_fld_window, text='graph size').grid(row=10, column=0)
    L_entry1 = tk.Entry(polar_fld_window, width=30, borderwidth=1)
    L_entry1.grid(row=11, column=0, padx = 2)
    L_entry1.insert(0, L)
    # because the arrow scaling only works in 'optimisation' button response
    # for which, the PLOT must be used first, chaning the plot to a
    # cartesian grid, need to also define a linear scaling in this window !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    tk.Label(polar_fld_window, text='linear scaling (used in plotting arrows)').grid(row=12, column=0)
    a_entry = tk.Entry(polar_fld_window, width=30, borderwidth=1)
    a_entry.grid(row=13, column=0, padx = 2)
    a_entry.insert(0, a_polar)
    # define a number that will tarck which vector field is wanted
    tensorp = tk.IntVar()
    # define buttons to chose between different fields
    # when one is selected, go to the submission function.
    # treat chosing the field as submitting all given variables.
    p_arr_btn = tk.Radiobutton(polar_fld_window, text='arrow', variable=tensorp, value=1, command=lambda: polar_submit(tensorp.get())).grid(row=7, column=1)
    p_both_btn = tk.Radiobutton(polar_fld_window, text='both', variable=tensorp, value=2, command=lambda: polar_submit(tensorp.get())).grid(row=7, column=2)
    p_stack_btn = tk.Radiobutton(polar_fld_window, text='stack', variable=tensorp, value=0, command=lambda: polar_submit(tensorp.get())).grid(row=7, column=3)
    # make sure that no Radiobutton is selected to begin with
    tensorp.set(None)  # as this one is showed as selected when opening the window


# define a function that will respons to field selection in the drop down menu
def field_selection_response(event):
    # clear the x and y component boxes
    x_comp_entry.delete(0, 'end')
    y_comp_entry.delete(0, 'end')
    # get the index at which this entry is
    selected_index = field_name_list.index(str(field_select_drop.get()))
    # using that index, get the x and y components from their lists
    # and insert these into x and y comp. entry boxes
    x_comp_selected = field_x_list[selected_index]
    y_comp_selected = field_y_list[selected_index]
    x_comp_entry.insert(0, x_comp_selected)
    y_comp_entry.insert(0, y_comp_selected)
    # now call the plot function to finalise all these onto the plot
    PLOT_response()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define wanted standard buttons
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# define the PLOT button
PLOT_btn = tk.Button(small_frame, text='PLOT', padx=60, pady=30, command=PLOT_response)
PLOT_btn.grid(row=0, column=0, columnspan=2, rowspan=2)

# define a button in small frame that will open new window to adjust arrowheads
custom_btn = tk.Button(small_frame, text='customise', padx=1, pady=1, command=custom_btn_reponse)
custom_btn.grid(row=0, column=3)

# define a button to open a window to input a polar field
polar_btn = tk.Button(bot_frame, text='Polar', padx=50, pady=20, command=Polar_btn_response)
polar_btn.grid(row=0, column=0)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define wanted entry boxes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
''' define entry boxes to update: L, pt_den, s_max and a  and a PLOT button '''
# define entry boxes for each (in order): L, pt_den, s_max and a ; and info txt
# Also input into them the initial values
tk.Label(small_frame, text='Size').grid(row=2, column=0)
L_entry = tk.Entry(small_frame, width=11, borderwidth=1)
L_entry.grid(row=3, column=0, padx = 2)
L_entry.insert(0, L)

tk.Label(small_frame, text='grid').grid(row=2, column=1)
pt_den_entry = tk.Entry(small_frame, width=11, borderwidth=1)
pt_den_entry.grid(row=3, column=1, padx = 2)
pt_den_entry.insert(0, pt_den)

tk.Label(small_frame, text='max sheets').grid(row=2, column=2)
s_max_entry = tk.Entry(small_frame, width=11, borderwidth=1)
s_max_entry.grid(row=3, column=2, padx = 2)
s_max_entry.insert(0, s_max)

# define entry boxes for the field equations in x and y
tk.Label(bot_frame, text='x component').grid(row=1, column=0)
x_comp_entry = tk.Entry(bot_frame, width=20, borderwidth=2)
x_comp_entry.grid(row=2, column=0)
x_comp_entry.insert(0, 'y*sin(x)')

tk.Label(bot_frame, text='y component').grid(row=1, column=1)
y_comp_entry = tk.Entry(bot_frame, width=20, borderwidth=2)
y_comp_entry.grid(row=2, column=1)
y_comp_entry.insert(0, '-x*cos(y)')

# define strings from initial components 
# these are needed by the derivative function, therefore for the derivative
# to work on the initial field, need to initially define them
string_x = str(x_comp_entry.get())
string_y = str(y_comp_entry.get())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define wanted Radio buttons
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

'''define buttons to choose quiver, quiver and stack or just stack plots'''
# define a number that will tarck which vector field is wanted
tensor = tk.IntVar()

# define each button and put them on the screen, in the right_frame
arrow_btn = tk.Radiobutton(right_frame, text='arrow', variable=tensor, value=1, command=lambda: vect_type_response(tensor.get())).grid(row=10, column=0)
arrow_stack_btn = tk.Radiobutton(right_frame, text='both', variable=tensor, value=2, command=lambda: vect_type_response(tensor.get())).grid(row=10, column=1)
stack_btn = tk.Radiobutton(right_frame, text='stack', variable=tensor, value=0, command=lambda: vect_type_response(tensor.get())).grid(row=10, column=2)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Define wanted dropdown menus
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

field_select = tk.StringVar()
field_select.set(field_name_list[0])

field_select_drop_label = tk.Label(bot_frame, text='Select Pre-Defined Field:')
field_select_drop_label.grid(row=3, column=0)

field_select_drop = ttk.Combobox(bot_frame, value = field_name_list, width=40)
field_select_drop.current(0)
field_select_drop.grid(row=4, column=0)
field_select_drop.bind("<<ComboboxSelected>>", field_selection_response)

# Testing matplotlib objects in functions

# def fun_test(axis):
#     ax_test = axis.inset_axes([0.1,0.1,0.2,0.2])
#     fig.canvas.draw()
#     deriv_inset_ax.clear()
#     deriv_inset_ax.remove()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Derivative Plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# define initial pixes coordinates for the function to use
# when mouse is clicked, these are changed to the pressed position
x_m = float(0)
y_m = float(0)


# define a function that will calculate the local, geometrical derivative
def deriv_calc(x_m, y_m):
    global i_m, j_m  # deriv_inset_ax
    
    # Range and point density of the derivative plot
    d_range = 0.33*L/(zoom_slider.get())
    d_length = d_length_select.get()
    dpd = dpd_select.get()
    d_scale = scale*(zoom_slider.get())
    
    # Select index for the middle of the new derivative axis
    i_m = int(round((dpd/2)-0.1))
    j_m = int(round((dpd/2)-0.1))
    
    # Divergence Plot
    if click_opt_int > 2:
        dx = np.linspace(-d_range, d_range, dpd)
        dy = np.linspace(-d_range, d_range, dpd)
        J = jacobian(2, string_x, string_y)
        # Evaluate the Jacobian elements at (x_m,y_m) click location 
        du_dx = eval(format_eq_div(format_eq(J[0, 0])))
        du_dy = eval(format_eq_div(format_eq(J[0, 1])))
        dv_dx = eval(format_eq_div(format_eq(J[1, 0])))
        dv_dy = eval(format_eq_div(format_eq(J[1, 1])))
        dxg, dyg = np.meshgrid(dx, dy)
        
        u_div = (du_dx + dv_dy)*dxg
        v_div = (du_dx + dv_dy)*dyg
        
        u_curl = (du_dy - dv_dx)*dyg
        v_curl = -(du_dy - dv_dx)*dxg
        
    # Zoom/Derivative Plot
    else:
        # define new axis in the derivative plot
        dx = np.linspace(-d_range+x_m, d_range+x_m, dpd)
        dy = np.linspace(-d_range+y_m, d_range+y_m, dpd)
        dxg, dyg = np.meshgrid(dx, dy)
        # define the vector field in these new axis
        u1, v1 = eq_to_comps(string_x, string_y, dxg, dyg)
        # Calculate derivative field components
        # This is done geometrically by subracting thevector at the centre of the
        # grid, from vectors at other grid positions in the derivative axis.
        du1 = u1 - u1[i_m, j_m]
        dv1 = v1 - v1[i_m, j_m]
  
    # Create axes at clicked position from supplied position and given axis sizes
    deriv_inset_ax = main_axis.inset_axes([(x_pix-178)/500 - (d_length/2), (y_pix-59)/500 - (d_length/2), d_length, d_length])
    
    # Check radiobutton selection
    if click_opt_int == 1:
        u_s = u1
        v_s = v1
        scale_s = d_scale
        
    elif click_opt_int == 2:
        u_s = du1
        v_s = dv1
        scale_s = scale
        
    elif click_opt_int == 3:
        u_s = u_div
        v_s = v_div
        scale_s = scale
        
    elif click_opt_int == 4:
        u_s = u_curl
        v_s = v_curl
        scale_s = scale
        
    # Stack or arrows plot
    if tensor.get() == 0:        
        # Stack
        arrows = False
        #stack_plot_deriv(dxg, dyg, u_s, v_s, 5, d_range, dpd, 0.1, arrows, orientation, scale_s)
        stack_plot(dxg, dyg, deriv_inset_ax, u_s, v_s, 5, d_range, dpd, 0.1, arrows, orientation, scale_s, w_head, h_head, 1)
    elif tensor.get() == 1:
        # Arrows        
        deriv_inset_ax.quiver(dxg, dyg, u_s, v_s, pivot='mid', scale=scale_s, scale_units='xy')
    elif tensor.get() == 2:
        # Arrows + Stack
        arrows = True
        #stack_plot_deriv(dxg, dyg, u_s, v_s, 5, d_range, dpd, 0.1, arrows, orientation, scale_s)
        stack_plot(dxg, dyg, deriv_inset_ax, u_s, v_s, 5, d_range, dpd, 0.1, arrows, orientation, scale_s, w_head, h_head, 1)

    # Don't display the x and y axis values
    deriv_inset_ax.set_xticks([])
    deriv_inset_ax.set_yticks([])
    
    # Redraw the figure canvas, showing the inset axis
    fig.canvas.draw()
    deriv_inset_ax.clear()
    deriv_inset_ax.remove()


# =============================================================================
# Define new function for plotting stacks in the derivative plot
# =============================================================================


# def stack_plot_deriv(xg, yg, u, v, s_max, L, pt_den, fract, arrows, orientation, scale, w_head=1/8, h_head=1/4):
#     # get the lengths of x and y from their grids
#     x_len = len(xg[:, 0])
#     y_len = len(yg[0, :])
#     # set the visuals for the derivative axis
#     deriv_inset_ax.set_aspect('equal')
    
#     # Account for change to grid centre for divergence plot
#     if click_opt_int > 2:       
#         deriv_inset_ax.set_xlim(-L-L/5, L+L/5)
#         deriv_inset_ax.set_ylim(-L-L/5, L+L/5)
#     else:
#         deriv_inset_ax.set_xlim(-L+x_m-L/5, L+x_m+L/5)
#         deriv_inset_ax.set_ylim(-L+y_m-L/5, L+y_m+L/5)
    
#     # AS BEFORE:
#     R_int = np.zeros(shape=((x_len), (y_len)))
    
#     if arrows is True:
#         deriv_inset_ax.quiver(xg, yg, u, v, pivot=orientation, scale=scale, scale_units='xy')
#     else:
#         pass

#     mag = np.sqrt(u**2 + v**2)
#     theta = np.arctan2(v, u)

#     sheet_L = L*fract
#     s_L = fract*L
#     max_size = np.max(mag)
#     R = mag/max_size

#     I_sin = np.sin(theta)
#     I_cos = np.cos(theta)
    
#     A_x = xg + (sheet_L/2)*I_sin
#     A_y = yg - (sheet_L/2)*I_cos
#     B_x = xg - (sheet_L/2)*I_sin
#     B_y = yg + (sheet_L/2)*I_cos
    
#     p_sh1x = xg + (s_L/2)*I_cos + (sheet_L*w_head)*I_sin
#     p_sh1y = yg + (s_L/2)*I_sin - (sheet_L*w_head)*I_cos
#     p_sh2x = xg + (s_L/2)*I_cos - (sheet_L*w_head)*I_sin
#     p_sh2y = yg + (s_L/2)*I_sin + (sheet_L*w_head)*I_cos
#     p_sh3x = xg + (s_L*0.5 + s_L*h_head)*I_cos
#     p_sh3y = yg + (s_L*0.5 + s_L*h_head)*I_sin
    
#     P_sh1x = xg + (sheet_L*w_head)*I_sin
#     P_sh1y = yg - (sheet_L*w_head)*I_cos
#     P_sh2x = xg - (sheet_L*w_head)*I_sin
#     P_sh2y = yg + (sheet_L*w_head)*I_cos
#     P_sh3x = xg + (s_L*h_head)*I_cos
#     P_sh3y = yg + (s_L*h_head)*I_sin
    
#     for i in range(x_len):
#         for j in range(y_len):
#             for t in range(1, s_max+1):
#                 if (t-1)/s_max <= R[i, j] <= t/s_max:
#                     R_int[i, j] = t
            
#             n = R_int[i, j]
            
#             # Prevent stack plotting in centre point of the derivative and div plot
#             if click_opt_int > 1 and i == i_m and j == j_m:
#                 continue
            
#             if parity(n) is True:
#                 s = 0
                
#                 while s <= 0.5*(n-2):
#                     Ax1 = A_x[i, j] + G(s, n, 0)*s_L*I_cos[i, j]
#                     Ay1 = A_y[i, j] + G(s, n, 0)*s_L*I_sin[i, j]
#                     Bx1 = B_x[i, j] + G(s, n, 0)*s_L*I_cos[i, j]
#                     By1 = B_y[i, j] + G(s, n, 0)*s_L*I_sin[i, j]
#                     Ax2 = A_x[i, j] - G(s, n, 0)*s_L*I_cos[i, j]
#                     Ay2 = A_y[i, j] - G(s, n, 0)*s_L*I_sin[i, j]
#                     Bx2 = B_x[i, j] - G(s, n, 0)*s_L*I_cos[i, j]
#                     By2 = B_y[i, j] - G(s, n, 0)*s_L*I_sin[i, j]
#                     deriv_inset_ax.add_line(Line2D((Ax1, Bx1), (Ay1, By1), linewidth=1, color='green'))
#                     deriv_inset_ax.add_line(Line2D((Ax2, Bx2), (Ay2, By2), linewidth=1, color='green'))
#                     s += 1
                    
#             elif parity(n) is False:
#                 deriv_inset_ax.add_line(Line2D((A_x[i, j], B_x[i, j]), (A_y[i, j], B_y[i, j]), linewidth=1, color='green'))
#                 s = 1 
                
#                 while s <= 0.5*(n-1):
#                     Ax1 = A_x[i, j] + G(s, n, 1)*s_L*I_cos[i, j]
#                     Ay1 = A_y[i, j] + G(s, n, 1)*s_L*I_sin[i, j]
#                     Bx1 = B_x[i, j] + G(s, n, 1)*s_L*I_cos[i, j]
#                     By1 = B_y[i, j] + G(s, n, 1)*s_L*I_sin[i, j]
#                     Ax2 = A_x[i, j] - G(s, n, 1)*s_L*I_cos[i, j]
#                     Ay2 = A_y[i, j] - G(s, n, 1)*s_L*I_sin[i, j]
#                     Bx2 = B_x[i, j] - G(s, n, 1)*s_L*I_cos[i, j]
#                     By2 = B_y[i, j] - G(s, n, 1)*s_L*I_sin[i, j]
#                     deriv_inset_ax.add_line(Line2D((Ax1, Bx1), (Ay1, By1), linewidth=1, color='green'))
#                     deriv_inset_ax.add_line(Line2D((Ax2, Bx2), (Ay2, By2), linewidth=1, color='green'))
#                     s += 1
                    
#             if n > 1:
#                 deriv_inset_ax.add_line(Line2D((p_sh1x[i, j], p_sh3x[i, j]), (p_sh1y[i, j], p_sh3y[i, j]), linewidth=1, color='green'))
#                 deriv_inset_ax.add_line(Line2D((p_sh2x[i, j], p_sh3x[i, j]), ((p_sh2y[i, j], p_sh3y[i, j])), linewidth=1, color='green'))
#             else:
#                 deriv_inset_ax.add_line(Line2D((P_sh1x[i, j], P_sh3x[i, j]), (P_sh1y[i, j], P_sh3y[i, j]), linewidth=1, color='green'))
#                 deriv_inset_ax.add_line(Line2D((P_sh2x[i, j], P_sh3x[i, j]), ((P_sh2y[i, j], P_sh3y[i, j])), linewidth=1, color='green'))


# set up the initial variable (code starts in option to use matplotlib tools)
click_opt_int = 0


# define a function that will update the variable that defines click action
def click_option_handler(click_option):
    global click_opt_int
    click_opt_int = click_option
    #fun_test(ax)
    # and for the the initial plot:
    if click_opt_int == 0:
        fig.canvas.draw()


# =============================================================================
# Calculate the Jacobian matrix of the defined vector field
# =============================================================================


def jacobian(m, u_str, v_str):
    # take the input strings and turn them into sympy expressions to be able to
    # use sympy's partial differentiation
    
    u_str = u_str.replace('^','**')
    v_str = v_str.replace('^','**')
    u_str = u_str.replace('step', '+0*')
    v_str = v_str.replace('step', '+0*')

    
    sympy_expr_x = parse_expr(u_str, evaluate=False)
    sympy_expr_y = parse_expr(v_str, evaluate=False)
    
    # define a sympy expression for string 0
    # sympy_expr_zero = parse_expr('0*x', evaluate=False)
    
    # combine the 2 intoa list:
    expressions = np.array([sympy_expr_x, sympy_expr_y])
    
    # set up an array to store derrivatives.
    J = np.empty((m, m), dtype='object')
    
    # set up an array of coordinates that need to be used (in standard order)
    coords = ['x', 'y']
    
    # set up an array to store the results
    # result = np.empty((int((m-1)*m/2), 1), dtype='object')
    # for i in range(int((m-1)*m/2)):
    #     result[i] = str(result[i])
    
    # loop over differentiating each, when differentiating w.r.t its coord, set to 0
    for coord_index in range(len(coords)):
        # loop over differentiating each component:
        for comp_index in range(len(expressions)):
            J[comp_index, coord_index] = str(diff(expressions[comp_index], coords[coord_index]))
            
    return J

# =============================================================================
# Additional formatting function used in divergence plots
# =============================================================================


def format_eq_div(string):
    string = string.replace('xg', 'x_m')
    string = string.replace('yg', 'y_m')
    return string


# =============================================================================
# Radiobutton to select what happens when clicking the plot
# =============================================================================

click_option = tk.IntVar()
click_option.set(0)

click_option_Tools_btn = tk.Radiobutton(right_frame, text='Tools', variable=click_option, value=0, command=lambda: click_option_handler(click_option.get()))
click_option_Zoom_btn = tk.Radiobutton(right_frame, text='Zoom', variable=click_option, value=1, command=lambda: click_option_handler(click_option.get()))
click_option_Deriv_btn = tk.Radiobutton(right_frame, text='Deriv.', variable=click_option, value=2, command=lambda: click_option_handler(click_option.get()))
click_option_Div_btn = tk.Radiobutton(right_frame, text='Div.', variable=click_option, value=3, command=lambda: click_option_handler(click_option.get()))
click_option_Curl_btn = tk.Radiobutton(right_frame, text='Curl', variable=click_option, value=4, command=lambda: click_option_handler(click_option.get()))

click_option_Tools_btn.grid(row=0, column=0)
click_option_Zoom_btn.grid(row=0, column=1)
click_option_Deriv_btn.grid(row=0, column=2)
click_option_Div_btn.grid(row=1, column=0)
click_option_Curl_btn.grid(row=1, column=1)

# =============================================================================
# Zooming window zoom slider
# =============================================================================

tk.Label(right_frame, text='Zoom').grid(row=2, column=0)
zoom_slider = tk.Scale(right_frame, from_=1, to=50, orient=tk.HORIZONTAL)
zoom_slider.grid(row=2, column=1)

# =============================================================================
# Drop down to select the derivative plot point density (dpd)
# =============================================================================

dpd_select = tk.IntVar()
dpd_select.set(5)
dpd_list = [5, 7, 9]

tk.Label(right_frame, text='Select Inset Plot Point Density:').grid(row=3, column=0)
dpd_drop = tk.OptionMenu(right_frame, dpd_select, *dpd_list)
dpd_drop.grid(row=3, column=1)

# =============================================================================
# Drop down to select inset axis size (d_length)
# =============================================================================

d_length_select = tk.DoubleVar()
d_length_select.set(0.3)
d_length_list = [0.2, 0.25, 0.3, 0.35, 0.4]

tk.Label(right_frame, text='Select Inset Plot Size (units?):').grid(row=4, column=0)
d_length_drop = tk.OptionMenu(right_frame, d_length_select, *d_length_list)
d_length_drop.grid(row=4, column=1)

# =============================================================================
# Step function - for singularities
# =============================================================================

# takes a matrix a and sorts through elements setting to zero if condition is met and one otherwise  
# used to remove singularity in Grav&Mag predefined field
def step(a):
    rows = len(a[:,0])
    columns = len(a[0,:])
    for i in range(rows):  
        for j in range(columns):
            if -1 < a[i,j] < 1:
                a[i,j] = 0
            else:
                a[i,j] = 1
    return a

# return time to run
stop = timeit.default_timer()
print('Time: ', stop - start)

# keep the program running until otherwise specified by user
tk.mainloop()
