
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
import os
from PIL import Image, ImageTk
from math import isnan
from matplotlib import patches as patch

# input many numpy functions to deal with user input
from numpy import sin, cos, tan, sqrt, log, arctan, arcsin, arccos, tanh
from numpy import sinh, cosh, arcsinh, arccosh, arctanh, exp, pi

# %% VFA GUI

# start the timer
start = timeit.default_timer()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define needed functions for the initial plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# define a function that will search for singularities and mark them on a very
# fine grid, but looking point by point without saving
def singularity_fine(u_str, v_str, dist_points, N=500):
    global points
    # set dot size as fraction of stack
    dot_size = 10
    # set grey square size as fract of point separation
    square_size = 10
    # set max accepted value before inf classed
    max_inf = 100
    # set up an array of NxN points
    points = np.zeros(shape=(2, N+1))
    # get interval_step size
    interval_size = (2*L)/N
    for i in range(N+1):
        points[0, i] = -L + i*interval_size
        points[1, i] = -L + i*interval_size
    # loop over an array of N by N to check for singularities
    for i in range(N+1):
        for j in range(N+1):
            # find the value of magnitudes at this point
            value_u = eval(format_eq(u_str, 0, 1, i, j))
            value_v = eval(format_eq(v_str, 0, 1, i, j))
            if isnan(value_u) is True or isnan(value_v) is True:
                #colour this region as a shaded square
                rect = patch.Rectangle((points[0, i] - dist_points/square_size, points[1, j]  - dist_points/square_size), 2*dist_points/square_size, 2*dist_points/square_size, color='#B5B5B5')
                main_axis.add_patch(rect)
            if abs(value_u) == np.inf  or abs(value_u) > max_inf or abs(value_v) == np.inf  or abs(value_v) > max_inf:
                # colour this point as a red dot
                circ = patch.Circle((points[0, i], points[1, j]), L*fract/dot_size, color='red')
                main_axis.add_patch(circ)


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
def stack_plot(xg, yg, axis, F_x, F_y, s_max, L, pt_den, fract, arrows=False, stacks=True, orientation='mid', scale=1, w_head=1/8, h_head=1/4, axis_check=0, arrowheads=True, colour='green'):
    global s_L, mag, bool_array, test_mag
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
    
    # find the distance between neightbouring points on the grid
    dist_points = xg[0, 1] - xg[0, 0]
    
    # define an empty array of magnitudes, to then fill with integer rel. mags
    R_int = np.zeros(shape=((x_len), (y_len)))

    # #########################################################################
    # get variables needed for the initial, simplified stack plot
    # #########################################################################
    
    # find the arrow length corresponding to each point and store in mag array
    mag = np.sqrt(F_x**2 + F_y**2)
    
    test_mag = mag * 1
    
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
    
    # Define scaling factor
    ScaleFactor = max_size/(0.9*(2*L/pt_den))

    # find the relative magnitude of vectors to maximum, as an array
    R = mag/max_size    
    
    if ascale.get() == 0:
        ScaleFactor = scale
    elif ascale.get() == 1:
        ScaleFactor = max_size/(0.9*(2*L/pt_den))
    
    # for arrows to work, with nan and infs
    # make a local variable of F_x and F_y
    # so that thye don't alter globally
    F_x_local = F_x * 1
    F_y_local = F_y * 1
    
    # plot the quiver plot on grid points if chosen in original function
    if arrows is True:
        # prevent any magnitudes from being inf or nan
        # only here, need to do it to u and v not just mag
        for i in range(x_len):
            for j in range(y_len):
                if isnan(F_x_local[i,j]) == True or isnan(F_y_local[i,j]) == True or abs(F_x_local[i, j]) == np.inf or abs(F_y_local[i, j]) == np.inf or abs(F_y_local[i, j]) > 1e15 or abs(F_x_local[i, j]) > 1e15:
                    F_x_local[i,j] = F_y_local[i,j] = 0
        axis.quiver(xg, yg, F_x_local, F_y_local, pivot=orientation, scale=ScaleFactor, scale_units='xy') 
    else:
        pass
    
    if stacks is True:
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


# define a function to replace input string to be 'python understood'
def format_eq(string, LI=0, singular_ty=0, i=0, j=0):
    # replace all the x and y with xg and yg:
    if LI == 1 :
        string = string.replace('x', 'intervals[0,:]')
        string = string.replace('y', 'intervals[1,:]')
    elif singular_ty == 1:
        string = string.replace('x', 'points[' + str(0) + ', ' + str(i) + ']')
        string = string.replace('y', 'points[' + str(1) + ', ' + str(j) + ']')
    else:
        string = string.replace('x', 'xg')
        string = string.replace('y', 'yg')
        string = string.replace('z', 'zg')
    # where there are special functions, replace them with library directions
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
    string = string.replace('**', '^')
    string = string.replace('np.log', 'ln')
    string = string.replace('np.exp', 'e**')
    string = string.replace('field_unit', '1')
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
    # return these
    return u, v


# define a function to take care of tab changes
def tab_selection(event):
    global click_opt_int, toolbar, LI_coord, LI_total_label, LI_total, LI_instruction_label, LI_restart_btn, LI_use_var, LI_shape_instruction, LI_shape_drop, LI_shape_select
    global x_m, y_m
    selected_tab = event.widget.select()
    tab_text = event.widget.tab(selected_tab, "text")
    print('you changed tab to ', tab_text)
    if tab_text == 'Line Integrals':
        x_m = None
        y_m = None
        fig.canvas.draw()
        # Initialise a global variable for storing click coordinates
        # and total line integral
        # home the main screen
        toolbar.home()
        # set the variable to show that LI was used
        LI_use_var = 1
        # unclick any chosen options from toolbar
        state = fig.canvas.toolbar.mode
        if state == 'zoom rect':
            toolbar.zoom()
        if state == 'pan/zoom':
            toolbar.pan()
        # get rid of the 2 buttons we don't want
        toolbar.children['!button4'].pack_forget()
        toolbar.children['!button5'].pack_forget()
        # set up an empty list to later store click coords
        LI_coord = []
        # initialise a variable that will keep track of total LI
        LI_total = 0
        # set variable for mouse to respond to integrals in new tab
        click_opt_int = 5
        # restart the plot, in case calculus code was used last
        PLOT_response()
        # set up grid initally as polygon is selected
        main_axis.grid(True)
        # draw it on
        canvas.draw()
    elif tab_text == 'Main':
        LI_restart()
        # by default return to initial 'tools'
        click_opt_int = 0
        click_option.set(0)
        click_option_handler(click_option.get())
        LI_shape_select.set(LI_shape_list[0])
        # restart the plot, in case calculus code was used last
        PLOT_response()
        # get rid of the grid
        main_axis.grid(False)
        # draw it on
        canvas.draw()
    elif tab_text == 'Calculus':
        # show an empty plot as nothing has been implemented yet
        # and no wedge has been done
        # TEMPORARY
        main_axis.clear()
        canvas.draw()
        tk.messagebox.showinfo('EMPTY PLOT INFO', 'Calculus code has not been merged here yet and no wedge was selected yet')

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

# get a toggle on and off switch images and scale the images to wanted size
toggle_image_on = Image.open('toggle_on_image.png')
toggle_image_off = Image.open('toggle_off_image.png')

toggle_image_on = toggle_image_on.resize((65, 25))
toggle_image_off = toggle_image_off.resize((65, 25))

toggle_image_on = ImageTk.PhotoImage(toggle_image_on)
toggle_image_off = ImageTk.PhotoImage(toggle_image_off)

# set up frames for each of:
# bottom side (field, scalings etc) and the right side (with detailed options)
# and top left for plot

# right frame:
right_frame_frame = tk.LabelFrame(root, text='Options Frame', padx=5, pady=5)
right_frame_frame.grid(row=1, column=1)

# bot frame:
bot_frame = tk.LabelFrame(root, text='Field Input Frame', padx=100, pady=5)
bot_frame.grid(row=2, column=0)

# plot frame:
plot_frame = tk.LabelFrame(root, text='Vector Field Frame', padx=5, pady=5)
plot_frame.grid(row=1, column=0)

# plot characteristics frame and plot button
small_frame = tk.LabelFrame(root, text='Plot Customisation Frame', padx=35, pady=5)
small_frame.grid(row=2, column=1)

# define notebook for tabs
notebook = ttk.Notebook(right_frame_frame)
notebook.grid(row=0, column=0)

# main options:
right_frame = tk.LabelFrame(notebook)
right_frame.grid(row=1, column=1)
# Line integrals
LI_frame = tk.Frame(notebook)
LI_frame.grid(row=0, column=2)
# calculus
calculus_frame = tk.Frame(notebook)
calculus_frame.grid(row=0, column=2)

# finsalise them
notebook.add(right_frame, text='Main')
notebook.add(LI_frame, text='Line Integrals')
notebook.add(calculus_frame, text='Calculus')

# bind the clicks on tabs to a function
notebook.bind_all('<<NotebookTabChanged>>', tab_selection)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set up initial parameters and plot the initial graph, put it in plot frame
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# define scale of the graph
L = 5
pt_den = 11   # number of points on each axis

# Initialise auto scaling variable
ascale = tk.IntVar()
ascale.set(0)

# define x and y values
x = np.linspace(-L, L, pt_den)
y = np.linspace(-L, L, pt_den)

# create a grid on x-y plane
xg, yg = np.meshgrid(x, y)

# define an example vector field
u = yg*np.sin(xg)  # x component
v = -xg*np.cos(yg)  # y component
# for no dependance in any initial component, use : np.zeros(np.shape(xg)) or yg

'''SET UP A LIST OF DEFAULT VECTOR FIELDS TO DISPLAY IN DROPDOWN MENU'''
# list of names of fields to display
field_name_list = ['Default: y*sin(x)dx - x*cos(y)dy',
                   'Simple pendulum: ydx  - sin(x)dy',
                   'Harmonic oscillator: ydx -xdy',
                   'Linear field example 1: (14*x - 4*y)dx + (-1*x + 4*y)dy',
                   'Linear field example 2: xdx',
                   'Constant field: 6dx + 3dy',
                   'Falling cat field (Planar 3 link robot)',
                   'Gravitational/Electric Point Charge: -x/(x**2+y**2)dx + -y/(x**2+y**2)dy',
                   'Magnetic Field of Current Carrying Wire: -y/(x**2+y**2)dx + x/(x**2+y**2)dy',
                   'Flamms paraboloid',
                   'BLACK HOLE!'
                   ]
# NOTE:
# Flamm's paraboloid ( https://rreusser.github.io/flamms-paraboloid/, https://en.wikipedia.org/wiki/Schwarzschild_metric")
# Black hole field analogue "taking the Schwarzchild contraction factor, at \theta = pi/2, g = (1-(r_s/r))^(-1) and defining the one form w = \del(g)/\del(x) dx + \del(g)/\del(y) dy    

# list of x components, in order of field_name_list
field_x_list = ['y*sin(x)',
                'y',
                'y',
                '14*x - 4*y',
                'x',
                '6',
                '(3*cos(y) + 4)/(15 + 6*cos(x) + 6*cos(y))',
                '-x/(x**2+y**2)',
                '-y/(x**2+y**2)',
                'x/(sqrt(x**2 + y**2)*(1-2/(sqrt(x**2 + y**2)))) - y',
                '-2*x*((x^2+y^2)^(-1.5))*(1-(2/sqrt(x^2+y^2)))^(-2)'
                ]

# list of y components, in order of field_name_list
field_y_list = ['- x*cos(y)',
                '-sin(x)',
                '-x',
                '(-1*x + 4*y)',
                '0',
                '3',
                '-(3*cos(x) + 4)/(15 + 6*cos(x) + 6*cos(y))',
                '-y/(x**2+y**2)',
                'x/(x**2+y**2)',
                'y/(sqrt(x**2 + y**2)*(1-2/(sqrt(x**2 + y**2)))) + x',
                '-2*y*((x^2+y^2)^(-1.5))*(1-(2/sqrt(x^2+y^2)))^(-2)'
                ]

# set up quiver factors
arrows = False  # set up if arrows should be plotted on stacks or not.
orientation = 'mid'  # how the arrow rotates about its assigned grid point - options: tail, mid and tip as string
scale = 5  # the scale reduction factor, if None (as None-type, not str), automatically computed by average, if 1 = mag

# set up the delta_factor of additional axis space L/delta_factor gives extra space on axis
delta_factor = 10

# fraction of sheet length to graph length
fract = 0.05

# define the maximum number of stack to plot, dep. on magnitude (initialy)
s_max = 4

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


''' define the initial polar plot also '''
# Despite it not being plotted to start with needed for when
# Polar plot option is used

r_min = 0
r_den = 11
theta_den = 25

# define the axis
r = np.linspace(r_min, L, r_den)
theta = np.linspace(360/(theta_den-1), 360, theta_den) * np.pi/180

# define a polar grid
rg, thetag = np.meshgrid(r, theta)

# define an initial linear scaling of the poalr field to be able to scale in
# poalr gird without converting to cartesian grid first
a_polar = 1

# set up a variable to keep track if LI was used recently
LI_use_var = 0

# set up line intergal enpty variables
LI_total = 0
LI_coord = []
shape_area = 0
ratio1 = 0

''' end of polar setting up'''

# set up initial strings for 2 forms window to display, for it to save properly after
to_wedge_x_1_str = ''
to_wedge_y_1_str = ''
to_wedge_x_2_str = ''
to_wedge_y_2_str = ''

# Initialise the click button selection
click_opt_int = 0

# define initial pixes coordinates for the function to use
# when mouse is clicked, these are changed to the pressed position
x_m = float(0)
y_m = float(0)

# initial orientation for circle integration:
orient_int = 'cw'

# define initial stack bool
stacks = True

# =============================================================================
# set up initial plot and canvas
# =============================================================================

# plot the cartessian field with desired parameters as specidfied above
plottedfield = stack_plot(xg, yg, main_axis, u, v, s_max, L, pt_den, fract, arrows, stacks, orientation, scale, w_head, h_head, 0)

# reduce white space from the figure in the plot frame
fig.tight_layout()

# set up the space for the plot to be put into AKA the plot frame
canvas = FigureCanvasTkAgg(fig, master=plot_frame)
canvas.draw()
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

# put the default matplotlib toolbar, below the figure
# NavigationToolbar2Tk.toolitems = [t for t in NavigationToolbar2Tk.toolitems if t[0] in ('Home','Back','Forward','Pan', 'Zoom','Save',)]

toolbar = NavigationToolbar2Tk(canvas, plot_frame)
toolbar.update()  # allow the plot to update based on the toolbar
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)


# track the mouse presses for the toolbar to respond to
def on_key_press(event):
    global x_pix, y_pix, x_m, y_m
    # get cartesian Coordinates of click
    # from those, get pixel Coordinates of click
    x_pix, y_pix = event.x , event.y
    x_m = float(event.xdata)
    y_m = float(event.ydata)
    # respond with only toolbar actions when only tools are to be used
    if click_opt_int == 0:
        key_press_handler(event, canvas, toolbar)
    # when the derivative option is selected, cerry out the derivative when clicked
    elif 0 < click_opt_int < 5:
        deriv_calc(x_m,y_m)
    elif click_opt_int == 5:
        if LI_shape_select.get() == 'Polygon':
            # Store the coordinates of the click in list
            LI_coord.append([x_m, y_m])
            line_int_poly(100000, string_x, string_y)
        # elif LI_shape_select.get() == 'square':
        #     print('not implemented')   # insert function later
        elif LI_shape_select.get() == 'Circle':
            # get the radius and call approperiate function
            Radius_LI_circ = eval(Radius_LI_circ_entry.get())
            line_int_circ([x_m,y_m], Radius_LI_circ, 100000, string_x, string_y, orient_int)


# connect figure event to a function that responds to clicks, defined above
cid = fig.canvas.mpl_connect("button_press_event", on_key_press)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define other needed functions, for input reponses
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


''' LINE INTEGRALS'''


# define a funciton to restart line integral calculation and lines
def LI_restart():
    global LI_total, LI_coord, shape_area
    # first, initialise variables again
    LI_total = 0
    LI_coord = []
    shape_area = 0
    ratio1 = 0
    # update the label
    LI_total_label.configure(text=LI_total)
    shape_area_label.configure(text=shape_area)
    ratio1_label.configure(text=ratio1)
    # NOT IDEAL BUT HOPEFULLY TEMPORARY
    # call plot-respose to redraw (so lines are deleted)
    # ideally would want a way of deleting lines without restarting the plot
    PLOT_response()
    # get the grid drawn on as needed
    if LI_shape_select.get() == 'Polygon':
        main_axis.grid(True)
        canvas.draw()
    else:
        pass


# define a button to change orientation of drawn shape
def orient_int_response():
    global orient_int
    if orient_int == 'cw':
        # change it to the other way around
        orient_int = 'ccw'
        # change the button display
        orient_int_btn.configure(text=orient_int)
    elif orient_int == 'ccw':
        # change it to the other way around
        orient_int = 'cw'
        # change the button display
        orient_int_btn.configure(text=orient_int)


# define LI shape selection response to dropdown
def LI_shape_select_response(selected_shape):
    # deal with lines
    global LI_shape_selected, Radius_LI_circ_entry, Radius_LI_label, orient_int_btn
    
    LI_shape_select.set(selected_shape)

    if selected_shape == 'Circle':
        # if circle is selected, display an entry box for the radius of it
        Radius_LI_label = tk.Label(LI_frame, text='Circle Radius:')
        Radius_LI_label.grid(row=3, column=0)
        Radius_LI_circ_entry = tk.Entry(LI_frame, width=10)
        Radius_LI_circ_entry.grid(row=3, column=1)
        Radius_LI_circ_entry.insert(0, '3')
        # also define a button to change orientation
        # define a button that will change circle orientation
        orient_int_btn = tk.Button(LI_frame, text=orient_int, command=orient_int_response)
        orient_int_btn.grid(row=4, column=1)
        # get rid of grid on plot
        main_axis.grid(False)
        # update canvas
        canvas.draw()
    else:
        # restart the plot and the integral values
        LI_restart()
        # set up a grid
        main_axis.grid(True)
        # update canvas
        canvas.draw()
        # get rid of circle options
        try:
            Radius_LI_label.destroy()
            Radius_LI_circ_entry.destroy()
            orient_int_btn.destroy()
        except UnboundLocalError:  
            pass

# Compute line integral for circles
def line_int_circ(cent, R, N, u_str, v_str, orient_int):
    global dt, LI_total_label, LI_total
    # Parametric increment (theta)
    dt = np.linspace(0, 2*np.pi, N)
    
    # Centre point
    xc = cent[0]
    yc = cent[1]
    
    # Create array to store interval point coordinates
    intervals = np.zeros(shape=(2, N))
    uv_store = np.zeros(shape=(2, N))
    dx = np.zeros(shape=(1, N))
    dy = np.zeros(shape=(1, N))
    
    # Magnitude of increment equals that of the circumference subtended by dt for large N
    A = 2*np.pi*R/(N)
    
    # Loop through to assign coordinates to interval points and plot them
    # depending on chosen orientation
    if orient_int == 'ccw':
        intervals[0, :] = xc + R*np.cos(dt)
        intervals[1, :] = yc + R*np.sin(dt)
    else:
        intervals[0, :] = xc - R*np.cos(dt)
        intervals[1, :] = yc - R*np.sin(dt)
    # get the points along the circle and save
    uv_store[0, :] = eval(format_eq(u_str, LI=1))
    uv_store[1, :] = eval(format_eq(v_str, LI=1))
    
    # Increment vector components
    dx = -A*np.sin(dt)
    dy = A*np.cos(dt)
    
    # res = np.sum(dx[:-1]*uv_store[0, :-1] + dy[:-1]*uv_store[1, :-1])
    res = np.sum(dx[:-1]*uv_store[0, :-1] + dy[:-1]*uv_store[1, :-1])
    
    # Plot the circle
    circle1 = mpl.patches.Circle(cent, R, fill=False, color='red')
    main_axis.add_artist(circle1)
    fig.canvas.draw()
    circle1.remove()
    
    # update the total
    LI_total = res
    # update its label
    LI_total_label.configure(text=str(round(LI_total, 6)))
    
    # Shape area label
    shape_area = np.pi*R**2
    shape_area_label.configure(text=str(round(shape_area, 3)))
    
    ratio1 = LI_total/shape_area
    ratio1_label.configure(text=str(round(ratio1, 3)))
    
    return res


# define a function that will complete the line integral
def line_int_poly(N, u_str, v_str):
    global LI_total, coord_array, coord_diff, LI_verts
    # set up a conuter to know how many time to complete the sum process
    c_count = len(LI_coord)
    
    # Tolerance for auto-joining lines together i.e. distance below which lines will join
    ctol = 0.1
    
    # array of coordinates
    coord_array = np.array(LI_coord)
    
    if c_count > 1:
        # Auto-join lines (check previous coords for if any are very close to eachother)
        coord_diff = coord_array - coord_array[c_count-1,:]
        for i in range(c_count):
            if sqrt((coord_diff[i,0])**2 + (coord_diff[i,1])**2) < ctol:
                LI_coord[c_count-1] = LI_coord[i]
                LI_verts = LI_coord[i:]
                break
        # get coordinates from mouse clicks
        a = LI_coord[c_count - 2]
        b = LI_coord[c_count - 1]
        # Plot line between points a and b
        main_axis.add_line(Line2D((a[0], b[0]), (a[1], b[1]), linewidth=2, color='red'))
        
        # linegrad = (b[1]-a[1])/(b[0]-a[0])  # keep as comment for now
        # find line length
        # linelength = np.sqrt((b[0]-a[0])**2 + (b[1]-a[1])**2)
        
        # find the length of interval
        # interval_len = linelength/N
        
        # find steps along the line, set by accuracy N
        dx = (b[0]-a[0])/N
        dy = (b[1]-a[1])/N
        
        # Create array to store interval point coordinates
        intervals = np.zeros(shape=(2, N))
        uv_store = np.zeros(shape=(2, N))
        
        # Loop through to assign coordinates to interval points and plot them
        for i in range(N):
            intervals[0, i] = a[0] + i*dx
            intervals[1, i] = a[1] + i*dy
            # ax.plot(intervals[0,i], intervals[1,i], 'bo', markersize=100/N)
        
        # Evaluate the vector components at the interval points
        uv_store[0, :] = eval(format_eq(u_str, 1))
        uv_store[1, :] = eval(format_eq(v_str, 1))
        
        # Evaluate line integral as sum of vector components multiplied by the small x and y displacements 
        res = dx*np.sum(uv_store[0, :]) + dy*np.sum(uv_store[1, :])
        
        # display the drawn on lines
        canvas.draw()
        
        # update the total
        LI_total += res
        # update its label
        LI_total_label.configure(text=str(round(LI_total, 6)))
        
        if len(LI_verts) > 3:
            shape_area = calc_area(LI_verts)
            shape_area_label.configure(text=str(round(shape_area, 3)))
            ratio1 = LI_total/shape_area
            ratio1_label.configure(text=str(round(ratio1, 3)))
            
        return res
    
    else:
        # If only clicked once, plot small red circle
        circle = patch.Circle(LI_coord[0], L*fract/6, color='red')
        main_axis.add_patch(circle)
        canvas.draw()


# Calc polygon area when user completes a shape
def calc_area(vert_list):
    # gte number of verticies
    n = len(vert_list)
    # set up array of vertiies from list
    M = np.array(vert_list)
    # initialise variables
    S1 = 0
    S2 = 0
    # loop over the verticies
    for i in range(n-1):
        # find total side lengths
        S1 = S1 + (M[i, 0] * M[i+1, 1])
        S2 = S2 + (M[i, 1] * M[i+1, 0])
    # get area from these
    A = 0.5*abs(S1-S2)
    return A


''' RESPONSE FUNCTIONS TO PLOT '''


# define a function that will respond to radio buttons behind choosing vector types:
# Note, tensor is a local variable. So to get value in other functions, use tensor.get() to extract current radiobutton value
def vect_type_response(tensor):
    # clear the plot that is already there:
    main_axis.clear()
    # use the tensor to determine what to plot:
    # 0 is just stacks, 1 is for only arrows and 2 is for both
    if tensor == 0:
        arrows = False
        stacks = True
    elif tensor == 1:
        arrows = True
        stacks = False
    elif tensor == 2:
        arrows = True
        stacks = True 
    stack_plot(xg, yg, main_axis, u, v, s_max, L, pt_den, fract, arrows, stacks, orientation, scale, w_head, h_head, 0)
    canvas.draw()

# define the PLOT button response function
def PLOT_response():
    # first, take from entry boxes, wanted parameters and make them global:
    # these must be changed globally for other functions to work with the new field.
    global L, pt_den, s_max, x, y, xg, yg, u, v, tensor, main_axis, string_x, string_y, arrows, stacks
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
    # plot depending on chosen type of vector
    if tensor.get() == 0:
        arrows = False  
        stacks = True
    elif tensor.get() == 1:
        arrows = True  
        stacks = False
    elif tensor.get() == 2:
        arrows = True
        stacks = True
    # create a figure and display it
    stack_plot(xg, yg, main_axis, u, v, s_max, L, pt_den, fract, arrows, stacks, orientation, scale, w_head, h_head, 0)
    canvas.draw()
    # recolour pt_den to white, if it was red from polar plots
    pt_den_entry.configure(bg='white')


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


''' CUSTOMISATIONS '''


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
    # recolour pt_den to white, if it was red from polar plots
    pt_den_entry.configure(bg='white')


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


# define a response funcction to autoscale toggle button
def scale_toggle_response():
    global ascale
    if ascale.get() == 0:
        # the burron is off, and has been clicked therefore change the
        # variable to an and the image to on
        ascale.set(1)
        ascale_toggle.configure(image=toggle_image_on)
        # for it to update, reclick whatever radiobutton is selected
        # or, if stacks only is chosen, change it to both, to show some change
        vect_type_response(tensor.get())
    else:
        # the button is on and has been clicked
        # set it to off and change image
        ascale.set(0)
        ascale_toggle.configure(image=toggle_image_off)
        # for it to update, reclick whatever radiobutton is selected
        # or, if stacks only is chosen, change it to both, to show some change
        vect_type_response(tensor.get())


# deifne a response to the SAVE button in the polar grid customisation window
def save_polar_grid():
    global r_min, r_den, theta_den, r, theta, rg, thetag
    r_min = float(r_min_entry.get())
    r_den = int(r_den_entry.get())
    theta_den = int(theta_den_entry.get())
    # using these redefine the new polar grids
    r = np.linspace(r_min, L, r_den)
    theta = np.linspace(360/(theta_den-1), 360, theta_den) * np.pi/180
    rg, thetag = np.meshgrid(r, theta)
    Polar_grid_plot_response(tensor.get())
    # once these are saved (made global), destroy the new window
    polar_grid_window.destroy()


# define a button that will open a new window where the user can
# customise the polar grid parameters
def polar_grid_custom_reponse():
    global r_min_entry, r_den_entry, theta_den_entry, polar_grid_window
    # open a titled new window
    polar_grid_window = tk.Toplevel()
    polar_grid_window.title('optimisation settings for the polar grid')
    # define an entry for minumum radius
    tk.Label(polar_grid_window, text='minimum radius:').grid(row=0, column=0)
    r_min_entry = tk.Entry(polar_grid_window, width=30, borderwidth=1)
    r_min_entry.insert(0, r_min)
    r_min_entry.grid(row=1, column=0)
    # define an entry for number of points along r
    tk.Label(polar_grid_window, text='number of points along r:').grid(row=2, column=0)
    r_den_entry = tk.Entry(polar_grid_window, width=30, borderwidth=1)
    r_den_entry.insert(0, r_den)
    r_den_entry.grid(row=3, column=0)
    # define an entry for number of points along theta
    tk.Label(polar_grid_window, text='number of points along theta:').grid(row=4, column=0)
    theta_den_entry = tk.Entry(polar_grid_window, width=30, borderwidth=1)
    theta_den_entry.insert(0, theta_den)
    theta_den_entry.grid(row=5, column=0)
    # define a button that will allow the user to save these inputs
    save_polar_grid_btn = tk.Button(polar_grid_window, text='SAVE', padx=20, pady=10, command=save_polar_grid)
    save_polar_grid_btn.grid(row=6, column=0, pady=10)


''' POLAR PLOTS '''

# define a function to repond to plotting apolar grid
# takes the same field, but plots it on a polar grid
def Polar_grid_plot_response(tensor):
    global xg, yg, u, v, s_max, pt_den_entry
    # set the number of sheets to use from input box
    s_max = int(s_max_entry.get())
    # the polar grid comes from global already defined
    # to change it, change it in the poalr field window
    # apart from size, this should be based on L
    # therefore use it to redefine it with that.
    L = float(L_entry.get())
    # using these redefine the new polar grids
    r = np.linspace(r_min, L, r_den)
    theta = np.linspace(360/(theta_den-1), 360, theta_den) * np.pi/180
    # mesh into a grid
    rg, thetag = np.meshgrid(r, theta)
    # convert grid to cartesian
    xg = rg*np.cos(thetag)
    yg = rg*np.sin(thetag)
    # reevaluate the given fields with these new grids:
    string_x = str(x_comp_entry.get())
    string_y = str(y_comp_entry.get())
    u, v = eq_to_comps(string_x, string_y, xg, yg)
    # clear the plot that is already there:
    main_axis.clear()
    # use the selected tensor to determine what to plot:
    # 0 is just stacks, 1 is for only arrows and 2 is for both
    if tensor == 0:
        arrows = False  
        stacks = True
    elif tensor == 1:
        arrows = True  
        stacks = False
    elif tensor == 2:
        arrows = True
        stacks = True
    # using those, create the plot and display it
    stack_plot(xg, yg, main_axis, u, v, s_max, L, pt_den, fract, arrows, stacks, orientation, scale, w_head, h_head, 0)
    canvas.draw()
    # colour pt_den red to show that it is not approperiate to use it now
    # need to def # of points along r and theta, in the additional window
    pt_den_entry.configure(bg='red')


'''

DIFFERENTIAL CALCULUS FUNCTIONS

'''

# define a function that will wedge two 1 forms and plot them
def wedge_product():
    global to_wedge_x_1_str, to_wedge_y_1_str, to_wedge_x_2_str, to_wedge_y_2_str
    # first, get all entries out, save as string for these to display when
    # window is opened again
    to_wedge_x_1_str = str(to_wedge_x_1_entry.get())
    to_wedge_y_1_str = str(to_wedge_y_1_entry.get())
    to_wedge_x_2_str = str(to_wedge_x_2_entry.get())
    to_wedge_y_2_str = str(to_wedge_y_2_entry.get())
    u_1, v_1 = eq_to_comps(to_wedge_x_1_str, to_wedge_y_1_str, xg, yg)
    u_2, v_2 = eq_to_comps(to_wedge_x_2_str, to_wedge_y_2_str, xg, yg)
    # clear the axis:
    main_axis.clear()
    # plot these as stacks, with no arrowheads, on top of one another.
    arrows = False
    stacks = True
    '''
    COMMENT
    # Given w_1 = fdx+gdy  and w_2=hdx+mdy. The graphical representation must be:
    # fdx /\ mdy "and" -hdx/\gdy {equivalend to [(u_1,0)+(0,v2)] and [(-u_2,0)+(0,v1)]},
    # which is executed when the if condition below is satisfited. Basically two rectangular
    # stacks with the scaling factor are accounted to via giving similar colors to the vertical
    # stacks (red) and green to the horizantal ones. After the first rectagular stacks are added
    # the second group will either sit on top of the first (in which case scaling contibution is zero)
    # or sit in some gaps and hence increasing the denisty as result of its scaling function.
    # If any of the coefficients (f,g,h and m)  is zero, the stacking reduces to one "function*dx/\dy", these are executed in the elif options.
    # One situation to be added when (u_1*v_2-u2*v_1).all() = 0, the scaling function here is zero and hence no 2-form should be produced/ or produced in faded color.
     
    # Issues:
    # 2- Would be nice to uniformly distribute the stacks in each direction after finishing the double stacking.
    '''
    if u_1.any() != 0 and v_1.any() != 0 and u_2.any() != 0 and v_2.any() != 0:
        stack_plot(xg, yg, main_axis, u_1, 0, s_max, L, pt_den, fract, arrows, stacks, orientation, scale, w_head, h_head, 0, arrowheads=False, colour='red')
        stack_plot(xg, yg, main_axis, 0, v_2, s_max, L, pt_den, fract, arrows, stacks, orientation, scale, w_head, h_head, 0, arrowheads=False, colour='green')
        stack_plot(xg, yg, main_axis, 0, v_1, s_max, L, pt_den, fract, arrows, stacks, orientation, scale, w_head, h_head, 0, arrowheads=False, colour='green')
        stack_plot(xg, yg, main_axis, -u_2, 0, s_max, L, pt_den, fract, arrows, stacks, orientation, scale, w_head, h_head, 0, arrowheads=False, colour='red')
    elif u_1.any() != 0 and v_2.any() != 0:
        stack_plot(xg, yg, main_axis, u_1, 0, s_max, L, pt_den, fract, arrows, stacks, orientation, scale, w_head, h_head, 0, arrowheads=False, colour='red')
        stack_plot(xg, yg, main_axis, 0, v_2, s_max, L, pt_den, fract, arrows, stacks, orientation, scale, w_head, h_head, 0, arrowheads=False, colour='green')
    elif v_1.any() != 0 and u_2.any() != 0:
        stack_plot(xg, yg, main_axis, 0, v_1, s_max, L, pt_den, fract, arrows, stacks, orientation, scale, w_head, h_head, 0, arrowheads=False, colour='green')
        stack_plot(xg, yg, main_axis, -u_2, 0, s_max, L, pt_den, fract, arrows, stacks, orientation, scale, w_head, h_head, 0, arrowheads=False, colour='red')
    else:
        print("The wedge product is zero")
        
    # put these onto the canvas
    canvas.draw()
    # close the extra window
    wedge_2_window.destroy()


# define a function that will repond to wedge_2_btn
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


# define a function that will calculate the local, geometrical derivative
def deriv_calc(x_m, y_m):
    global i_m, j_m  # deriv_inset_ax
    
    # Try and except to account for the error caused when zoom window selected without first clicking
    # (i.e. undefined x_pix and y_pix)
    
    try:
    
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
            du_dx = round(eval(format_eq_div(format_eq(J[0, 0]))),8)
            du_dy = round(eval(format_eq_div(format_eq(J[0, 1]))),8)
            dv_dx = round(eval(format_eq_div(format_eq(J[1, 0]))),8)
            dv_dy = round(eval(format_eq_div(format_eq(J[1, 1]))),8)
            dxg, dyg = np.meshgrid(dx, dy)
            
            # Div --> Trace of the Jacobian Matrix
            u_div = (du_dx + dv_dy)*dxg
            v_div = (du_dx + dv_dy)*dyg
            
            #Curl --> Skew Symmetric Part of Jacobian Matrix 0.5*(A-A^T)
            u_curl = -0.5*(du_dy - dv_dx)*dyg
            v_curl = 0.5*(du_dy - dv_dx)*dxg
            
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
        deriv_inset_ax = main_axis.inset_axes([(x_pix-178)/500 - (0.931*d_length/(2*L)), (y_pix-59)/500 - (0.931*d_length/(2*L)), 0.931*d_length/L, 0.931*d_length/L])
        
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
    
        if tensor.get() == 0:
            arrows = False
            stacks = True
        if tensor.get() == 1:
            arrows = True
            stacks = False
        if tensor.get() == 2:
            arrows = True
            stacks = True            
        
        stack_plot(dxg, dyg, deriv_inset_ax, u_s, v_s, 5, d_range, dpd, 0.1, arrows, stacks, orientation, scale_s, w_head, h_head, 1) 
    
        # Don't display the x and y axis values
        if click_opt_int > 2:  
            deriv_inset_ax.set_xticks([])
            deriv_inset_ax.set_yticks([])
        
        # Redraw the figure canvas, showing the inset axis
        fig.canvas.draw()
        deriv_inset_ax.clear()
        deriv_inset_ax.remove()
        
    # i.e. if click coordinates are undefined, do nothing
    except (NameError, UnboundLocalError):
        pass


# Calculate the Jacobian matrix of the defined vector field
def jacobian(m, u_str, v_str):
    # take the input strings and turn them into sympy expressions to be able to
    # use sympy's partial differentiation
    u_str = u_str.replace('^','**')
    v_str = v_str.replace('^','**')
    sympy_expr_x = parse_expr(u_str, evaluate=False)
    sympy_expr_y = parse_expr(v_str, evaluate=False)
    # define a sympy expression for string 0
    # sympy_expr_zero = parse_expr('0*x', evaluate=False)
    
    # combine the 2 into a list:
    expressions = np.array([sympy_expr_x, sympy_expr_y])
    # set up an array to store derrivatives.
    J = np.empty((m, m), dtype='object')
    # set up an array of coordinates that need to be used (in standard order)
    coords = ['x', 'y']
    # loop over differentiating each, when differentiating w.r.t its coord, set to 0
    for coord_index in range(len(coords)):
        # loop over differentiating each component:
        for comp_index in range(len(expressions)):
            J[comp_index, coord_index] = str(diff(expressions[comp_index], coords[coord_index]))
    return J


# define a function that will update the variable that defines click action
# and deals with setting up responses to click option changes
def click_option_handler(click_option):
    global click_opt_int, toolbar, LI_coord, LI_total_label, LI_total, LI_instruction_label, LI_restart_btn, LI_use_var, LI_shape_instruction, LI_shape_drop, LI_shape_select
    click_opt_int = click_option
    # tools being selected
    if click_opt_int == 0:
        x_m = None
        y_m = None
        fig.canvas.draw()
        # if the tools is selected again, add the zoom and pan buttons
        # get rid of the modified toolbar:
        toolbar.destroy()
        # put the default matplotlib toolbar, back on:
        toolbar = NavigationToolbar2Tk(canvas, plot_frame)
        toolbar.update()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    # when 'tool' is not selected, disable the pan and zoom:
    elif 0 < click_opt_int < 5:
        fig.canvas.draw()
        toolbar.home()
        # close the selected mouse options
        state = fig.canvas.toolbar.mode
        if state == 'zoom rect':
            toolbar.zoom()
        if state == 'pan/zoom':
            toolbar.pan()
        # get rid of the 2 buttons we don't want
        toolbar.children['!button4'].pack_forget()
        toolbar.children['!button5'].pack_forget()
        # run deriv calc as per its calling
        deriv_calc(x_m, y_m)


# Additional formatting function used in divergence plots
def format_eq_div(string):
    string = string.replace('xg', 'x_m')
    string = string.replace('yg', 'y_m')
    return string


def update_deriv(self):
    deriv_calc(x_m,y_m)


''' SINGULARITIES '''


# define response to button to show singularities
def show_singularities():
    # get the number of points to use:
    fine_grid_N = int(fine_grid_N_entry.get())
    # warn user
    if fine_grid_N > 200:
        tk.messagebox.showwarning('WARNING', 'This will run for a long time if there exist many singulairites and if the fine grid has too many points')
    # ask for how many points to show
    dist_points = xg[0, 1] - xg[0, 0]
    singularity_fine(string_x, string_y, dist_points, N=fine_grid_N)
    canvas.draw()


# define a function that will respond to plotting a known singularity
def known_singularity_response():
    global magnitude_check
    # get the variable
    singularity_eq_type = singular_var.get()
    # set up a finer grid than normal to show clearly
    x_sing, y_sing = np.linspace(-L, L, pt_den*5), np.linspace(-L, L, pt_den*5)
    known_singularity = known_singularity_entry.get()
    known_singularity = format_eq(known_singularity)
    if singularity_eq_type == singular_list[0]:
        # get each input separately
        inputs = known_singularity.split('; ')
        # find out how many of them have been put in:
        number = len(inputs)
        # loop over plotting all
        for i in range(number):
            known_singularity = inputs[i]
            # from the input, define y values
            if known_singularity.find('xg') == -1 & known_singularity.find('yg') == -1:
                y_vals_singular = eval('(' + known_singularity + ') + 0*y_sing')
            else:
                # replace variables
                known_singularity = known_singularity.replace('xg', 'x_sing')
                known_singularity = known_singularity.replace('yg', 'y_sing')
                y_vals_singular = eval(known_singularity)
            # plot
            main_axis.plot(x_sing, y_vals_singular, 'r-.')
            # check if these actually are singularities
            checker_singularities = 0
            string_check_x = string_x.replace('x', 'x_sing')
            string_check_x = string_check_x.replace('y', 'y_vals_singular')
            string_check_y = string_y.replace('x', 'x_sing')
            string_check_y = string_check_y.replace('y', 'y_vals_singular')
            # get the magnitude
            magnitude_check = np.sqrt(eval(string_check_x)**2 + eval(string_check_y)**2)
            # check if this is a singularity
            for i in range(len(magnitude_check)):
                if isnan(magnitude_check[i]) is True or abs(magnitude_check[i]) == np.inf or magnitude_check[i] > 1e3:
                    checker_singularities = 1
                else:
                    pass
            if checker_singularities == 0:
                tk.messagebox.showwarning('WARNING', 'The point you have input does not register as a singularity')
    elif singularity_eq_type == singular_list[1]:
        # get each input separately
        inputs = known_singularity.split('; ')
        # find out how many of them have been put in:
        number = len(inputs)
        # loop over plotting all
        for i in range(number):
            known_singularity = inputs[i]
            # as above but the other way around
            if known_singularity.find('xg') == -1 & known_singularity.find('yg') == -1:
                x_vals_singular = eval('(' + known_singularity + ') + 0*x_sing')
            else:
                # replace variables
                known_singularity = known_singularity.replace('xg', 'x_sing')
                known_singularity = known_singularity.replace('yg', 'y_sing')
                x_vals_singular = eval(known_singularity)
            # plot
            main_axis.plot(x_vals_singular, y_sing, 'r-.')
            # check if these actually are singularities
            checker_singularities = 0
            string_check_x = string_x.replace('x', 'x_vals_singular')
            string_check_x = string_check_x.replace('y', 'y_sing')
            string_check_y = string_y.replace('x', 'x_vals_singular')
            string_check_y = string_check_y.replace('y', 'y_sing')
            # get the magnitude
            magnitude_check = np.sqrt(eval(string_check_x)**2 + eval(string_check_y)**2)
            # check if this is a singularity
            for i in range(len(magnitude_check)):
                if isnan(magnitude_check[i]) is True or abs(magnitude_check[i]) == np.inf or magnitude_check[i] > 1e3:
                    checker_singularities = 1
                else:
                    pass
            if checker_singularities == 0:
                tk.messagebox.showwarning('WARNING', 'The point you have input does not register as a singularity')
    elif singularity_eq_type == singular_list[2]:
        # get each input separately
        inputs = known_singularity.split('; ')
        # find out how many of them have been put in:
        number = len(inputs)
        # loop over plotting all
        for i in range(number):
            known_singularity = inputs[i]
            # split the string into the 2 coordinates
            known_singularity = known_singularity.split(',')
            point_x = eval(known_singularity[0])
            point_y = eval(known_singularity[1])
            circ = patch.Circle((point_x, point_y), L*fract/3, color='red')
            main_axis.add_patch(circ)
            # check if this is actually a singularity and
            # warn user if it is not.
            string_check_x = string_x.replace('x', 'point_x')
            string_check_x = string_check_x.replace('y', 'point_y')
            string_check_y = string_y.replace('x', 'point_x')
            string_check_y = string_check_y.replace('y', 'point_y')
            # get the magnitude
            magnitude_check = np.sqrt(eval(string_check_x)**2 + eval(string_check_y)**2)
            # check if this is a singularity
            if isnan(magnitude_check) is True or abs(magnitude_check) == np.inf or magnitude_check > 1e15:
                pass
            else:
                tk.messagebox.showwarning('WARNING', 'The point you have input does not register as a singularity')
    canvas.draw()


# define a function that will respond to dropdown for singularities
# doesn't need to do anything.
def singular_drop_response(var):
    return



# =============================================================================
# DEFINE ALL NEEDED WIDGETS
# =============================================================================

'''

DEFINE ALL WIDGETS IN MAIN TAB

'''

# define a number that will tarck which vector field is wanted
tensor = tk.IntVar()
tensor.set(0)

tensor_label = tk.Label(right_frame, text='Arrows/Stacks:')
tensor_label.grid(row=7, column=0)

# define each button and put them on the screen, in the right_frame
arrow_btn = tk.Radiobutton(right_frame, text='arrow', variable=tensor, value=1, command=lambda: vect_type_response(tensor.get())).grid(row=7, column=1)
arrow_stack_btn = tk.Radiobutton(right_frame, text='both', variable=tensor, value=2, command=lambda: vect_type_response(tensor.get())).grid(row=7, column=3)
stack_btn = tk.Radiobutton(right_frame, text='stack', variable=tensor, value=0, command=lambda: vect_type_response(tensor.get())).grid(row=7, column=2)

# define a button for 2 1 forms to wedge
# this will open a new window where the uder can input 2 2 forms that will be wedged
# and the result will be plotted
wedge_2_btn = tk.Button(calculus_frame, text='wedge', command= wedge_2_response)
wedge_2_btn.grid(row=0, column=0)

# get a button to draw on singularities
singularity_button = tk.Button(right_frame, text='search singularities', command=show_singularities)
singularity_button.grid(row=8, column=0)
# entry for N
tk.Label(right_frame, text='<-- sampling points').grid(row=8, column=2, columnspan=2)
fine_grid_N_entry = tk.Entry(right_frame, width=5)
fine_grid_N_entry.grid(row=8, column=1)
fine_grid_N_entry.insert(0, 10)

# define an entry where the user can inpu known singularity equation
# this will be taken and plotted as a red, dotted line
tk.Label(right_frame, text='equation of known singularity :').grid(row=9, column=0, columnspan=2)

# define a dropdown to select y= or x=
singular_var = tk.StringVar()
singular_list = ['y=', 'x=', 'point']
singular_var.set(singular_list[0])
dpd_drop = tk.OptionMenu(right_frame, singular_var, *singular_list, command=singular_drop_response)
dpd_drop.grid(row=10, column=0)
# equation entry box
known_singularity_entry = tk.Entry(right_frame, width=20)
known_singularity_entry.grid(row=10, column=1)
known_singularity_entry.insert(0, '')

# define asubmit button to that entry
submit_known_singularity_btn = tk.Button(right_frame, text='show expression', command=known_singularity_response)
submit_known_singularity_btn.grid(row=11, column=0)

# DERIVATIVE FUNCTIONS

# Radiobuttons to select what happens when clicking the plot
click_option = tk.IntVar()
click_option.set(0)
click_option_Tools_btn = tk.Radiobutton(right_frame, text='Tools', variable=click_option, value=0, command=lambda: click_option_handler(click_option.get()))
click_option_Zoom_btn = tk.Radiobutton(right_frame, text='Zoom', variable=click_option, value=1, command=lambda: click_option_handler(click_option.get()))
click_option_Deriv_btn = tk.Radiobutton(right_frame, text='Deriv.', variable=click_option, value=2, command=lambda: click_option_handler(click_option.get()))
click_option_Div_btn = tk.Radiobutton(right_frame, text='Div.', variable=click_option, value=3, command=lambda: click_option_handler(click_option.get()))
click_option_Curl_btn = tk.Radiobutton(right_frame, text='Curl', variable=click_option, value=4, command=lambda: click_option_handler(click_option.get()))
click_option_Tools_btn.grid(row=1, column=0)
click_option_Zoom_btn.grid(row=1, column=1)
click_option_Deriv_btn.grid(row=1, column=2)
click_option_Div_btn.grid(row=2, column=0)
click_option_Curl_btn.grid(row=2, column=1)

# Zooming window zoom slider
tk.Label(right_frame, text='Zoom').grid(row=3, column=0)
zoom_slider = tk.Scale(right_frame, from_=1, to=50, orient=tk.HORIZONTAL)
zoom_slider.bind("<ButtonRelease-1>", update_deriv)
zoom_slider.grid(row=3, column=1)

# Drop down to select the derivative plot point density (dpd)
dpd_select = tk.IntVar()
dpd_select.set(5)
dpd_list = [5, 7, 9]

tk.Label(right_frame, text='Select Inset Plot Point Density:').grid(row=4, column=0)
dpd_drop = tk.OptionMenu(right_frame, dpd_select, *dpd_list, command = update_deriv)
dpd_drop.grid(row=4, column=1)

# Drop down to select inset axis size (d_length)
d_length_select = tk.DoubleVar()
d_length_list = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]
d_length_select.set(d_length_list[2])
tk.Label(right_frame, text='Select Inset Plot Size :').grid(row=5, column=0)
d_length_drop = tk.OptionMenu(right_frame, d_length_select, *d_length_list, command = update_deriv)
d_length_drop.grid(row=5, column=1)

# Autoscale Toggle
ascale_label = tk.Label(right_frame, text='Toggle Autoscaling:')
ascale_label.grid(row=6, column=0)
ascale_toggle = tk.Button(right_frame, image=toggle_image_off, bd=0, command=scale_toggle_response)
ascale_toggle.grid(row=6, column=1, pady=5)

# Step function - for singularities
#def step(a):
#    rows = len(a[:,0])
#    columns = len(a[0,:])
#    for i in range(rows):  
#        for j in range(columns):
#            if -1 < a[i,j] < 1:
#                a[i,j] = 0
#            else:
#                a[i,j] = 1
#    return a


'''

set up all in SMALL FRAME

'''


# define the PLOT button
PLOT_btn = tk.Button(small_frame, text='PLOT', padx=60, pady=30, command=PLOT_response)
PLOT_btn.grid(row=0, column=0, columnspan=2, rowspan=1)

# define a button in small frame that will open new window to adjust arrowheads
custom_btn = tk.Button(small_frame, text='customise visuals', padx=1, pady=1, command=custom_btn_reponse)
custom_btn.grid(row=0, column=3)

# define a button to customise the polar grids
polar_grid_custom_btn = tk.Button(small_frame, text='customise polar grid', padx=1, pady=1, command=polar_grid_custom_reponse)
polar_grid_custom_btn.grid(row=1, column=3)

# define a button that will just plot the given cartesian field
# on a polar grid
polar_grid_plot_btn = tk.Button(small_frame, text='polar grid plot', command= lambda: Polar_grid_plot_response(tensor.get()))
polar_grid_plot_btn.grid(row=1, column=0, columnspan=2)

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


'''

set up all in BOTTOM FRAME

'''


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

field_select = tk.StringVar()
field_select.set(field_name_list[0])

field_select_drop_label = tk.Label(bot_frame, text='Select Pre-Defined Field:')
field_select_drop_label.grid(row=3, column=0)

field_select_drop = ttk.Combobox(bot_frame, value = field_name_list, width=40)
field_select_drop.current(0)
field_select_drop.grid(row=4, column=0)
field_select_drop.bind("<<ComboboxSelected>>", field_selection_response)


'''

set up all in LI tab

'''


# define a label that will display it
LI_instruction_label = tk.Label(LI_frame, text='LI Total:')
LI_instruction_label.grid(row=0, column=0, padx=10)
LI_total_label = tk.Label(LI_frame, text=LI_total)
LI_total_label.grid(row=0, column=1)

tk.Label(LI_frame, text='Shape Area:').grid(row=0, column=2)
shape_area_label = tk.Label(LI_frame, text=shape_area)
shape_area_label.grid(row=0, column=3)

tk.Label(LI_frame, text='Ratio:').grid(row=0, column=4)
ratio1_label = tk.Label(LI_frame, text=ratio1)
ratio1_label.grid(row=0, column=5)

# display a restart button that will clear the lines
# and restart the variables.
LI_restart_btn = tk.Button(LI_frame, text='LI Restart', padx=20, command=LI_restart)
LI_restart_btn.grid(row=1, column=0, columnspan=2)

# define a drop down to draw: connected lines, square or circle.
LI_shape_select = tk.StringVar()
LI_shape_list = ['Polygon', 'Circle']
LI_shape_select.set(LI_shape_list[0])
LI_shape_instruction = tk.Label(LI_frame, text='Select what to draw:')
LI_shape_instruction.grid(row=2, column=0)
LI_shape_drop = tk.OptionMenu(LI_frame, LI_shape_select, *LI_shape_list, command=LI_shape_select_response)
LI_shape_drop.grid(row=2, column=1)


# return time to run
stop = timeit.default_timer()
print('Time: ', stop - start)

# keep the program running until otherwise specified by user
tk.mainloop()
