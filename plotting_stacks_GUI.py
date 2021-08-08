
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

# %% VFA GUI

# start the timer
start = timeit.default_timer()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define needed functions for the initial plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# define a function that will find region ON GRID
# that are not defined as regions, as opposed to a lines or points
# do so by finding  values of neightbouring values
# if there is 3 or more that are also not defined, its a region
def undef_region(mag):
    # set up a same size array to store classifications
    # 1 will mean a region, zero will mean a point or a line (of inf or nan)
    bool_array = np.zeros(np.shape(u))
    # get the lengths of supplied fields
    x_len = len(mag[:, 0])
    y_len = len(mag[0, :])
    # loop over all indexes
    for i in range(x_len):
        for j in range(y_len):
            # set up a conuting variable
            counter = 0
            # for the current index, 'look' up, down, left and right
            # check how many of these are inf or nan.
            # making sure, ends of grids are taken care of
            # check up
            try:
                if abs(mag[i - 1, j]) == np.inf or isnan(mag[i - 1, j]) is True or abs(mag[i - 1, j]) > 1e15:
                    counter += 1
            except IndexError:
                pass
            # check down
            try:
                if abs(mag[i + 1, j]) == np.inf or isnan(mag[i + 1, j]) is True or abs(mag[i + 1, j]) > 1e15:
                    counter += 1
            except IndexError:
                pass
            # check left
            try:
                if abs(mag[i, j - 1]) == np.inf or isnan(mag[i, j - 1]) is True or abs(mag[i, j - 1]) > 1e15:
                    counter += 1
            except IndexError:
                pass
            # check right
            try:
                if abs(mag[i, j + 1]) == np.inf or isnan(mag[i, j + 1]) is True or abs(mag[i, j + 1]) > 1e15:
                    counter += 1
            except IndexError:
                pass
            # now check corners:
            try:
                if abs(mag[i - 1, j - 1]) == np.inf or isnan(mag[i - 1, j - 1]) is True or abs(mag[i -1, j - 1]) > 1e15:
                    counter += 1
            except IndexError:
                pass
            try:
                if abs(mag[i + 1, j + 1]) == np.inf or isnan(mag[i + 1, j + 1]) is True or abs(mag[i + 1, j + 1]) > 1e15:
                    counter += 1
            except IndexError:
                pass
            try:
                if abs(mag[i + 1, j - 1]) == np.inf or isnan(mag[i + 1, j - 1]) is True or abs(mag[i + 1, j - 1]) > 1e15:
                    counter += 1
            except IndexError:
                pass
            try:
                if abs(mag[i - 1, j + 1]) == np.inf or isnan(mag[i - 1, j + 1]) is True or abs(mag[i - 1, j + 1]) > 1e15:
                    counter += 1
            except IndexError:
                pass
            # CRUCIAL
            # check its own point too!!
            try:
                if abs(mag[i, j]) == np.inf or isnan(mag[i, j]) is True or abs(mag[i, j]) > 1e15:
                    counter += 1
            except IndexError:
                pass
            # now, depending on couner, define truth value in bool array
            # True if counter >= 3, then its a region
            if counter > 3:
                bool_array[i, j] = 1
            elif counter <= 3:
                bool_array[i, j] = 0
    return bool_array
            


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
def stack_plot(xg, yg, axis, u, v, s_max, L, pt_den, fract, arrows=False, stacks=True, orientation='mid', scale=1, w_head=1/8, h_head=1/4, axis_check=0, arrowheads=True, colour='green'):
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
    mag = np.sqrt(u**2 + v**2)
    
    test_mag = mag
    
    # find direction of each arrow
    angles = np.arctan2(v, u)   # theta defined from positive x axis ccw
    
    # find regions ON GRID that are nan or inf as a bool array
    bool_array = undef_region(mag)
    
    # deal with infs and nans in mag
    for i in range(x_len):
        for j in range(y_len):
            # set to zero points that are not defined or inf
            if isnan(mag[i, j]) is True or abs(mag[i, j]) == np.inf  or abs(mag[i, j]) > 1e15:
                # depending on bool_array, shade points on grid that are in undefined
                # region
                if bool_array[i, j] == 1:
                    # colour this region as a shaded square
                    rect = patch.Rectangle((xg[i, j] - dist_points/2, yg[i, j]  - dist_points/2), dist_points, dist_points, color='#B5B5B5')
                    axis.add_patch(rect)
                if bool_array[i, j] == 0:
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
    
    # plot the quiver plot on grid points if chosen in original function
    if arrows is True:
        # prevent any magnitudes from being inf or nan
        # only here, need to do it to u and v not just mag
        for i in range(x_len):
            for j in range(y_len):
                if isnan(u[i,j]) == True or isnan(v[i,j]) == True or abs(u[i, j]) == np.inf or abs(v[i, j]) == np.inf or abs(v[i, j]) > 1e10 or abs(u[i, j]) > 1e10:
                    u[i,j] = v[i,j] = 0
        axis.quiver(xg, yg, u, v, pivot=orientation, scale=ScaleFactor, scale_units='xy') 
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
                '-x/(x**2+y**2)',
                '-y/(x**2+y**2)'
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
                'x/(x**2+y**2)'
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


'''define the initial polar plot also. Despite it not being plotted to start
with needed for when Polar plot option is used'''

r_min = 0.2
r_max = 5
r_den = 10
theta_den = 20

# define the axis
r = np.linspace(r_min, r_max, r_den)
theta = np.linspace(0, 360, theta_den) * np.pi/180

# define a polar grid
rg, thetag = np.meshgrid(r, theta)

# define an initial linear scaling of the poalr field to be able to scale in
# poalr gird without converting to cartesian grid first
a_polar = 1

''' end of polar setting up'''

# set up initial strings for 2 forms window to display, for it to save properly after
to_wedge_x_1_str = ''
to_wedge_y_1_str = ''
to_wedge_x_2_str = ''
to_wedge_y_2_str = ''

stacks = True

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
    # respond with only toolbar actions when only tools are to be used
    if click_opt_int == 0:
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
        deriv_calc(x_m,y_m)

# connect figure event to a function that responds to clicks, defined above
cid = fig.canvas.mpl_connect("button_press_event", on_key_press)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define other needed functions, for input reponses
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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

    
    if tensor.get() == 0:
        arrows = False  
        stacks = True
    elif tensor.get() == 1:
        arrows = True  
        stacks = False
    elif tensor.get() == 2:
        arrows = True
        stacks = True

    stack_plot(xg, yg, main_axis, u, v, s_max, L, pt_den, fract, arrows, stacks, orientation, scale, w_head, h_head, 0)

    canvas.draw()


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


# define a function to repond to plotting apolar grid
# takes the same field, but plots it on a polar grid
def Polar_grid_plot_response(tensor):
    global xg, yg, u, v, s_max
    # set the number of sheets to use from input box
    s_max = int(s_max_entry.get())
    # the polar grid comes from global already defined
    # to change it, change it in the poalr field window
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
    
    stack_plot(xg, yg, main_axis, u, v, s_max, L, pt_den, fract, arrows, stacks, orientation, scale, w_head, h_head, 0)
    
    canvas.draw()


# deifne a response to the SAVE button in the polar grid customisation window
def save_polar_grid():
    global r_min, r_max, r_den, theta_den, r, theta, rg, thetag
    r_min = float(r_min_entry.get())
    r_max = float(r_max_entry.get())
    r_den = int(r_den_entry.get())
    theta_den = int(theta_den_entry.get())
    # using these redefine the new polar grids
    r = np.linspace(r_min, r_max, r_den)
    theta = np.linspace(0, 360, theta_den) * np.pi/180
    rg, thetag = np.meshgrid(r, theta)
    Polar_grid_plot_response(tensor.get())
    # once these are saved (made global), destroy the new window
    polar_grid_window.destroy()


# define a button that will open a new window where the user can
# customise the polar grid parameters
def polar_grid_custom_reponse():
    global r_min_entry, r_max_entry, r_den_entry, theta_den_entry, polar_grid_window
    # open a titled new window
    polar_grid_window = tk.Toplevel()
    polar_grid_window.title('optimisation settings for the polar grid')
    # define an entry for minumum radius
    tk.Label(polar_grid_window, text='minimum radius:').grid(row=0, column=0)
    r_min_entry = tk.Entry(polar_grid_window, width=30, borderwidth=1)
    r_min_entry.insert(0, r_min)
    r_min_entry.grid(row=1, column=0)
    # define an entry for radius limit r_max
    tk.Label(polar_grid_window, text='maximum radius:').grid(row=2, column=0)
    r_max_entry = tk.Entry(polar_grid_window, width=30, borderwidth=1)
    r_max_entry.insert(0, r_max)
    r_max_entry.grid(row=3, column=0)
    # define an entry for number of points along r
    tk.Label(polar_grid_window, text='number of points along r:').grid(row=4, column=0)
    r_den_entry = tk.Entry(polar_grid_window, width=30, borderwidth=1)
    r_den_entry.insert(0, r_den)
    r_den_entry.grid(row=5, column=0)
    # define an entry for number of points along theta
    tk.Label(polar_grid_window, text='number of points along theta:').grid(row=6, column=0)
    theta_den_entry = tk.Entry(polar_grid_window, width=30, borderwidth=1)
    theta_den_entry.insert(0, theta_den)
    theta_den_entry.grid(row=7, column=0)
    # define a button that will allow the user to save these inputs
    save_polar_grid_btn = tk.Button(polar_grid_window, text='SAVE', padx=20, pady=10, command=save_polar_grid)
    save_polar_grid_btn.grid(row=8, column=0, pady=10)


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
    # 1- The max number of stacks possible will hinder good visualization when
    # the stacks are dense. It's a general issue, but more clear here than other cases due to the nature of 2-forms.
    # 2- Would be nice to uniformly distribute the stacks in each direction after finishing the double stacking.
    
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


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define wanted Radio buttons
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

'''define buttons to choose quiver, quiver and stack or just stack plots'''
# define a number that will tarck which vector field is wanted
tensor = tk.IntVar()
tensor.set(0)

tensor_label = tk.Label(right_frame, text='Toggle Arrows/Stacks:')
tensor_label.grid(row = 10, column = 0)

# define each button and put them on the screen, in the right_frame
arrow_btn = tk.Radiobutton(right_frame, text='arrow', variable=tensor, value=1, command=lambda: vect_type_response(tensor.get())).grid(row=10, column=1)
arrow_stack_btn = tk.Radiobutton(right_frame, text='both', variable=tensor, value=2, command=lambda: vect_type_response(tensor.get())).grid(row=10, column=2)
stack_btn = tk.Radiobutton(right_frame, text='stack', variable=tensor, value=0, command=lambda: vect_type_response(tensor.get())).grid(row=10, column=3)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define wanted standard buttons
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

# define a button for 2 1 forms to wedge
# this will open a new window where the uder can input 2 2 forms that will be wedged
# and the result will be plotted
wedge_2_btn = tk.Button(right_frame, text='wedge', command= wedge_2_response)
wedge_2_btn.grid()

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
            du_dx = eval(format_eq_div(format_eq(J[0, 0])))
            du_dy = eval(format_eq_div(format_eq(J[0, 1])))
            dv_dx = eval(format_eq_div(format_eq(J[1, 0])))
            dv_dy = eval(format_eq_div(format_eq(J[1, 1])))
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
        deriv_inset_ax.set_xticks([])
        deriv_inset_ax.set_yticks([])
        
        # Redraw the figure canvas, showing the inset axis
        fig.canvas.draw()
        deriv_inset_ax.clear()
        deriv_inset_ax.remove()
        
    # i.e. if click coordinates are undefined, do nothing
    except: NameError

# Initialise the click button selection
click_opt_int = 0

# define a function that will update the variable that defines click action
def click_option_handler(click_option):
    global click_opt_int, toolbar
    click_opt_int = click_option
    
    if click_opt_int == 0:
        fig.canvas.draw()
        # if the tools is selected again, add the zoom and pan buttons
        # get rid of the modified toolbar:
        toolbar.destroy()
        # put the default matplotlib toolbar, back on:
        toolbar = NavigationToolbar2Tk(canvas, plot_frame)
        toolbar.update()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        
    # when 'tool' is not selected, disable the pan and zoom:
    elif click_opt_int > 0:
        fig.canvas.draw()
        
        toolbar.home()
        #toolbar.destroy()
        
        # close the selected mouse options
        state = fig.canvas.toolbar.mode
        if state == 'zoom rect':
            toolbar.zoom()
        if state == 'pan/zoom':
            toolbar.pan()
        
        # get rid of the 2 buttons we don't want
        
        # toolbar.destroy()
        # NavigationToolbar2Tk.toolitems = [t for t in NavigationToolbar2Tk.toolitems if t[0] not in ('Pan', 'Zoom',)]
        # toolbar = NavigationToolbar2Tk(canvas, plot_frame)
        # toolbar.update()
        
        toolbar.children['!button4'].pack_forget()
        toolbar.children['!button5'].pack_forget()
        
        deriv_calc(x_m, y_m)

# =============================================================================
# Calculate the Jacobian matrix of the defined vector field
# =============================================================================

def jacobian(m, u_str, v_str):
    # take the input strings and turn them into sympy expressions to be able to
    # use sympy's partial differentiation
    
    u_str = u_str.replace('^','**')
    v_str = v_str.replace('^','**')
    
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

def update_deriv(self):
    deriv_calc(x_m,y_m)
    
# =============================================================================
# Radiobutton to select what happens when clicking the plot
# =============================================================================

# click_option_label = tk.Label(right_frame, text='Select Toolbar and Zooming Windows:')
# click_option_label.grid(row = -1, column = 0)

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
zoom_slider.bind("<ButtonRelease-1>", update_deriv)
zoom_slider.grid(row=2, column=1)

# =============================================================================
# Drop down to select the derivative plot point density (dpd)
# =============================================================================

dpd_select = tk.IntVar()
dpd_select.set(5)
dpd_list = [5, 7, 9]

tk.Label(right_frame, text='Select Inset Plot Point Density:').grid(row=3, column=0)
dpd_drop = tk.OptionMenu(right_frame, dpd_select, *dpd_list, command = update_deriv)
dpd_drop.grid(row=3, column=1)

# =============================================================================
# Drop down to select inset axis size (d_length)
# =============================================================================

d_length_select = tk.DoubleVar()
d_length_list = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0]
d_length_select.set(d_length_list[2])

tk.Label(right_frame, text='Select Inset Plot Size :').grid(row=4, column=0)
d_length_drop = tk.OptionMenu(right_frame, d_length_select, *d_length_list, command = update_deriv)
d_length_drop.grid(row=4, column=1)

# =============================================================================
# Autoscale Toggle
# =============================================================================

ascale_label = tk.Label(right_frame, text='Toggle Autoscaling:')
ascale_label.grid(row=5, column=0)

ascale_toggle = tk.Button(right_frame, image=toggle_image_off, bd=0, command=scale_toggle_response)
ascale_toggle.grid(row=5, column=1, pady=5)

# =============================================================================
# Step function - for singularities
# =============================================================================

# takes a matrix a and sorts through elements setting to zero if condition is met and one otherwise  
# used to remove singularity in Grav&Mag predefined field

# Problem with this strategy is that 0*nan = nan so this was not eliminating the singularity
# and causing error when automatically scaling. Still could be useful if used to define the zero region
# in a plot with a central singularity. 

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

# return time to run
stop = timeit.default_timer()
print('Time: ', stop - start)

# keep the program running until otherwise specified by user
tk.mainloop()
