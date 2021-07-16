# Attempt 1 - zooming tool on the stack plot

# import modules
import timeit
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import tkinter as tk
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler

# %%

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
def stack_plot(xg, yg, ax, u, v, s_max, L, pt_den, fract, arrows='True', orientation='mid', scale=1, w_head=1/8, h_head=1/4):
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

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set up basic layout of the window
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# define an object tracker - for GUI set up
root = tk.Tk()

# set its title
root.title('Vector field analyser - differential forms')

# set a window size for it all to initially appear in
root.geometry("1400x920")

# set up frames for each of:
# bottom side (field, scalings etc) and the right side (with detailed options)
# and top left for plot

# right frame:
right_frame = tk.LabelFrame(root, text='Options Frame', padx=120, pady=228)
right_frame.grid(row=0, column=1)

# bot frame:
bot_frame = tk.LabelFrame(root, text='Field Input Frame', padx=192, pady=87)
bot_frame.grid(row=1, column=0)

# plot frame:
plot_frame = tk.LabelFrame(root, text='Vector Field Frame', padx=10, pady=10)
plot_frame.grid(row=0, column=0)

# plot characteristics frame and plot button
small_frame = tk.LabelFrame(root, text='Plot Customisation Frame', padx=29, pady=41)
small_frame.grid(row=1, column=1)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set up initial parameters and plot the initial graph, put it in plot frame
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# define scale of the graph
L = 5
pt_den = 20   # number of points on each axis

# define x and y values
x = np.linspace(-L, L, pt_den)
y = np.linspace(-L, L, pt_den)

# create a grid on x-y plane
xg, yg = np.meshgrid(x, y)

# define an example vector field
a = 0.2
u = a*yg  # x component
v = -a*xg  # y component
# for no dependance, use : np.zeros(np.shape(xg))  ----- or yg

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

# set up quiver factors
arrows = False  # set up if arrows should be plotted on stacks or not.
orientation = 'mid'  # how the arrow rotates about its assigned grid point - options: tail, mid and tip as string
scale = 1  # the scale reduction factor, if None (as None-type, not str), automatically computed by average, if 1 = mag

# set up the delta_factor of additional axis space L/delta_factor gives extra space on axis
delta_factor = 10

# fraction of sheet length to graph length
fract = 0.05

# define the maximum number of stack to plot, dep. on magnitude
s_max = 6

# set screen dpi
my_dpi = 100

# define denominator of fractional height and width of arrowhead based on stack size
w_head = 1/8
h_head = 1/4
    
# create a figure
fig = plt.figure(figsize=(855/my_dpi, 573/my_dpi), dpi=my_dpi)

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

# plot the field with desired parameters as specidfied above
plottedfield = stack_plot(xg, yg, ax, u, v, s_max, L, pt_den, fract, arrows, orientation, scale, w_head, h_head)

# reduce white space on the plot figure
fig.tight_layout()

# set up the space foe the plot to be put into AKA the plot frame
canvas = FigureCanvasTkAgg(fig, master=plot_frame)
canvas.draw()
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

# put the default matplotlib toolbar, below the figure
toolbar = NavigationToolbar2Tk(canvas, plot_frame)
toolbar.update()
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

# tack the mouse presses for the toolbar to respond to
def on_key_press(event):
    print("you pressed {}".format(event.key))
    key_press_handler(event, canvas, toolbar)

# connect the space to function that records clicks
canvas.mpl_connect("key_press_event", on_key_press)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define other needed functions, for input reponses
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# define a function to replace input string to be 'python understood'
def format_eq(string):
    # replace all the x and y with xg and yg:
    string = string.replace('x', 'xg')
    string = string.replace('y', 'yg')
    # where there are special functions, replace them with library directions
    string = string.replace('pi', 'np.pi')
    string = string.replace('sqrt', 'np.sqrt')
    string = string.replace('sin', 'np.sin')
    string = string.replace('cos', 'np.cos')
    string = string.replace('tan', 'np.tan')
    string = string.replace('^', '**')
    string = string.replace('ln', 'np.log') 
    string = string.replace('exp', 'np.exp')       
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
    # scale with given a:
    u *= a
    v *= a
    # return these
    return u, v


# defina a function that will respond to radio buttons behind choosing vector types:
def vect_type_response(tensor):
    # clear the plot that is already there:
    ax.clear()
    # use the tensor to determine what to plot:
    # 0 is just stacks, 1 is for only arrows and 2 is for both
    if tensor == 0:
        arrows = False
        stack_plot(xg, yg, ax, u, v, s_max, L, pt_den, fract, arrows, orientation, scale, w_head, h_head)
        canvas.draw()
    elif tensor == 1:
        ax.quiver(xg, yg, u, v, pivot=orientation, scale=scale, scale_units='xy')
        # repeat the displaying of the figure so that it updates in GUI
        canvas.draw()
    elif tensor == 2:
        arrows = True
        stack_plot(xg, yg, ax, u, v, s_max, L, pt_den, fract, arrows, orientation, scale, w_head, h_head)
        # repeat the displaying of the figure so that it updates in GUI
        canvas.draw()


# define the PLOT button response function
def PLOT_response():
    global L, pt_den, s_max, a, x, y, xg, yg, u, v, tensor, ax
    # clear the current axis
    ax.clear()
    # take the new axis parameters and field definitions out of the boxes
    L = float(L_entry.get())
    pt_den = int(pt_den_entry.get())
    s_max = int(s_max_entry.get())
    a = float(a_entry.get())
    string_x = str(x_comp_entry.get())
    string_y = str(y_comp_entry.get())
    # from L redefine the axis
    ax_L = L + L/delta_factor
    ax.set_xlim(-ax_L, ax_L)
    ax.set_ylim(-ax_L, ax_L)
    # from pt_den and L, change the axis coordinates and the grid:
    x = np.linspace(-L, L, pt_den)
    y = np.linspace(-L, L, pt_den)
    xg, yg = np.meshgrid(x, y)
    # take all these values, and the input from field component bnoxes to set up the field:
    u, v = eq_to_comps(string_x, string_y, xg, yg, u, v)
    # plot the new field
    stack_plot(xg, yg, ax, u, v, s_max, L, pt_den, fract, arrows, orientation, scale, w_head, h_head)
    # put it onto the screen
    canvas.draw()
    # change the radio button ticks back to stack only
    tensor.set(0)


# define a function to respond to submitting arrohead changes in the new window
def custom_submission():
    # first, take from entry boxes, wanted parameters and make them global:
    global w_head, h_head, fract
    w_head = float(w_entry.get())
    h_head = float(h_entry.get())
    fract = float(fract_entry.get())
    # clear the axis
    ax.clear()
    # replot the graph with new arrows:
    stack_plot(xg, yg, ax, u, v, s_max, L, pt_den, fract, arrows, orientation, scale, w_head, h_head)
    # put it onto the screen
    canvas.draw()
    # change the radio button ticks back to stack only
    tensor.set(0)
    # then close the window
    arrowH_opt_window.destroy()


# define a reponse function to open a new window when arrowh_btn is pressed:
def custom_btn_reponse():
    global w_entry, h_entry, fract_entry, arrowH_opt_window
    # open a titled new window
    arrowH_opt_window = tk.Toplevel()
    arrowH_opt_window.title('optimisation settings')
    # define and label a first entry, for width
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
    # define a button to submit those changes:
    submit_arr_btn = tk.Button(arrowH_opt_window, text='SUBMIT ALL', padx=20, pady=10, command=custom_submission)
    submit_arr_btn.grid(row=6, column=0, pady=10)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define wanted standard buttons
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# define the PLOT button
PLOT_btn = tk.Button(small_frame, text='PLOT', padx=60, pady=30, command=PLOT_response)
PLOT_btn.grid(row=0, column=0, columnspan=2, rowspan=2)

# define a small button in small frame that will open new window to adjust arrowheads
arrowh_btn = tk.Button(small_frame, text='customise', padx=1, pady=1, command=custom_btn_reponse)
arrowh_btn.grid(row=0, column=3)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define wanted entry boxes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
''' define entry boxes to update: L, pt_den, s_max and a  and a PLOT button '''
# define entry boxes for each (in order): L, pt_den, s_max and a ; and info txt
# Also input into them the initial values
L_label = tk.Label(small_frame, text='Size').grid(row=2, column=0)
L_entry = tk.Entry(small_frame, width=11, borderwidth=1)
L_entry.grid(row=3, column=0, padx = 2)
L_entry.insert(0, L)

pt_den_label = tk.Label(small_frame, text='grid').grid(row=2, column=1)
pt_den_entry = tk.Entry(small_frame, width=11, borderwidth=1)
pt_den_entry.grid(row=3, column=1, padx = 2)
pt_den_entry.insert(0, pt_den)

s_max_label = tk.Label(small_frame, text='max sheets').grid(row=2, column=2)
s_max_entry = tk.Entry(small_frame, width=11, borderwidth=1)
s_max_entry.grid(row=3, column=2, padx = 2)
s_max_entry.insert(0, s_max)

a_label = tk.Label(small_frame, text='field scaling  \'a\'').grid(row=2, column=3)
a_entry = tk.Entry(small_frame, width=11, borderwidth=1)
a_entry.grid(row=3, column=3, padx = 2)
a_entry.insert(0, a)

# define entry boxes for the field equations in x and y
x_comp_label = tk.Label(bot_frame, text='x component').grid(row=0, column=0)
x_comp_entry = tk.Entry(bot_frame, width=20, borderwidth=2)
x_comp_entry.grid(row=1, column=0)
x_comp_entry.insert(0, 'a*y')

y_comp_label = tk.Label(bot_frame, text='y component').grid(row=0, column=1)
y_comp_entry = tk.Entry(bot_frame, width=20, borderwidth=2)
y_comp_entry.grid(row=1, column=1)
y_comp_entry.insert(0, '-a*x')

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
# define wanted sliders
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define wanted checkboxes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define wanted checkboxes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mouse Click Plotting for the Derivative Field
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Derivative Plot - works well for linear fields!
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Comment out this section if not working on derivatives! 

# Open the window for the derivative plot
x_m = float(0)
y_m = float(0)

def deriv_calc(x_m,y_m):
    # index of point where derivative is taken
    # Set index from the x_mouse and y_mouse coords    
    global i_m, j_m
    i_m = int(round(((y_m/(2*L))*(pt_den))+(pt_den/2),0))
    j_m = int(round(((x_m/(2*L))*(pt_den))+(pt_den/2),0))
    # constant to change the size of the arrows
    b = 0.2
    # subtract the components of the point i,j from the rest of the vectors
    du = u[i_m, j_m]
    dv = v[i_m, j_m]
    DF_u = u - du
    DF_v = v - dv
    # create arrays local to the point. P is the point density in the plot
    P = 5
    K = int((P-1)/2)
    
    z_xg = xg[(i_m-K):(i_m+K)+1, (j_m-K):(j_m+K+1)]
    z_yg = yg[(i_m-K):(i_m+K)+1, (j_m-K):(j_m+K+1)]
    
    DFz_u = b*DF_u[(i_m-K):(i_m+K)+1, (j_m-K):(j_m+K+1)]
    DFz_v = b*DF_v[(i_m-K):(i_m+K)+1, (j_m-K):(j_m+K+1)]
    
    # create new window for the plot - edit size?
    deriv_window = tk.Toplevel()
    # to do - give coord rather than the index - more useful for user
    deriv_window.title('(Approximate) Derivative Field at coord: (x=' + str(x_m) + ', y=' + str(y_m) + ')')
    d_fig = plt.figure()
    d_ax = d_fig.gca()
    # local quiver plot
    d_ax.quiver(z_xg, z_yg, DFz_u, DFz_v, pivot='mid', scale_units='xy')
    # draw the plot on new window
    d_canvas = FigureCanvasTkAgg(d_fig, master=deriv_window)
    d_canvas.draw()
    d_canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    plt.close()
    
    
# Define button to open the derivative plot
deriv_button = tk.Button(right_frame, pady=10, text='Local Derivative', command=deriv_calc(x_m,y_m)).grid(row=1, column=1)

def onclick(event):
    global ix, iy, coords, x_m, y_m

    ix_plot,iy_plot = event.xdata,event.ydata
    print (ix_plot,iy_plot)
    
    x_m = float(ix_plot)
    y_m = float(iy_plot)
    
    deriv_calc(x_m,y_m)
        
cid = fig.canvas.mpl_connect('button_press_event', onclick)

# return time to run
stop = timeit.default_timer()
print('Time: ', stop - start)

# keep the program running until otherwise specified by user
tk.mainloop()
