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
# define needed functions
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
def stack_plot(xg, yg, ax, u, v, s_max, L, pt_den, fract, arrows='True', orientation='mid', scale=1, a=1):
    # get axis lengths:
    x_len = len(xg[:, 0])
    y_len = len(yg[:, 0])
    # set us parameters that need to be global:
    global ax_L
    
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
    
    # define points of stack arrowheads as arrays for all stacks
    p_sh1x = xg + (s_L/2)*np.cos(theta) + (sheet_L/8)*np.sin(theta)
    p_sh1y = yg + (s_L/2)*np.sin(theta) - (sheet_L/8)*np.cos(theta)
    p_sh2x = xg + (s_L/2)*np.cos(theta) - (sheet_L/8)*np.sin(theta)
    p_sh2y = yg + (s_L/2)*np.sin(theta) + (sheet_L/8)*np.cos(theta)
    p_sh3x = xg + (3*s_L/4)*np.cos(theta)
    p_sh3y = yg + (3*s_L/4)*np.sin(theta)
    
    # define these for when there is only 1 line in the stack plot:
    P_sh1x = xg + (sheet_L/8)*np.sin(theta)
    P_sh1y = yg - (sheet_L/8)*np.cos(theta)
    P_sh2x = xg - (sheet_L/8)*np.sin(theta)
    P_sh2y = yg + (sheet_L/8)*np.cos(theta)
    P_sh3x = xg + (s_L/4)*np.cos(theta)
    P_sh3y = yg + (s_L/4)*np.sin(theta)
    
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


# defina a function that will respond to radio buttons behind choosing vector types:
# def vect_type_response(var):
    

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set up basic layout of the window
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# define an object tracker - for GUI set up
root = tk.Tk()

# set the icon
#root.iconbitmap('C:\\Users\\macus\\Desktop\\Uni\\summer internships\\Moustafa - Differential Forms\\Greek-omega icon.ico')

# set its title
root.title('Vector field analyser - differential forms')

# set a window size for it all to initially appear in
root.geometry("1300x920")

# set up frames for each of:
# bottom side (field, scalings etc) and the right side (with detailed options)
# and top left for plot

# right frame:
right_frame = tk.LabelFrame(root, text='options frame', padx=150, pady=300)
right_frame.grid(row=0, column=1)

# put something in it for now
right_frame_dummy_label = tk.Label(right_frame, text='options frame').pack()

# bot frame:
bot_frame = tk.LabelFrame(root, text='field frame', padx=400, pady=100)
bot_frame.grid(row=1, column=0)

# put something in it for now
bot_frame_dummy_label = tk.Label(bot_frame, text='field frame').pack()

# plot frame:
plot_frame = tk.LabelFrame(root, text='plot frame', padx=10, pady=10)
plot_frame.grid(row=0, column=0)

# plot characteristics frame and plot button
small_frame = tk.LabelFrame(root, text='plot frame', padx=136, pady=100)
small_frame.grid(row=1, column=1)

# put something in it for now
small_frame_dummy_label = tk.Label(small_frame, text='axis options frame').pack()


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define wanted standard buttons
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define wanted Radio buttons
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

'''define buttons to choose quiver, quiver and stack or just stack plots'''
# define a number that will tarck which vector field is wanted
tensor = tk.IntVar()

# define each button and put them on the screen, in the right_frame
#vector_type_btn1 = tk.Radiobutton(right_frame, text='arrow', variable=tensor, value=1, command=lambda:vect_type_response(tensor.get())).grid(row=10, column=0)
#vector_type_btn2 = tk.Radiobutton(right_frame, text='both', variable=tensor, value=2, command=lambda:vect_type_response(tensor.get())).grid(row=10, column=1)
#vector_type_btn2 = tk.Radiobutton(right_frame, text='stack', variable=tensor, value=3, command=lambda:vect_type_response(tensor.get())).grid(row=10, column=2)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define wanted sliders
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define wanted checkboxes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define the input boxes for fileds in the field frame
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  plot the graph and put down the plot in its plot frame
# do this with the initialfield to be plotted
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# define scale of the graph
L = 5
pt_den = 26   # number of points on each axis

# define x and y values
x = np.linspace(-L, L, pt_den)
y = np.linspace(-L, L, pt_den)

# create a grid on x-y plane
xg, yg = np.meshgrid(x, y)

# define an example vector field
a = 0.01
u = 2*xg  # x component
v = np.zeros(np.shape(xg))  # y component
# for no dependance, use : np.zeros(np.shape(grid))

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
s_max = 10

# set screen dpi
my_dpi = 100

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
plottedfield = stack_plot(xg, yg, ax, u, v, s_max, L, pt_den, fract, arrows, orientation, scale, a)

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
# define the axis options in the axis frame
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# input to set the graph L

# input to change the s_max

# input to change pt_den

# input to change scale

# input to change fract









# return time to run
stop = timeit.default_timer()
print('Time: ', stop - start)

# keep the program running until otherwise specified by user
tk.mainloop()

