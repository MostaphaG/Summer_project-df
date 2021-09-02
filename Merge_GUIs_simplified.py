# Merging all files with simplifications, no blocks, no splitting options
# for 2-forms.

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
from sympy import diff, simplify
from sympy.parsing.sympy_parser import parse_expr
from PIL import Image, ImageTk
from math import isnan
from matplotlib import patches as patch
import matplotlib.path as mplPath
from matplotlib import animation
from scipy.integrate import odeint
import matplotlib.path as mPath

# input many numpy functions to deal with user input
from numpy import sin, cos, tan, sqrt, log, arctan, arcsin, arccos, tanh
from numpy import sinh, cosh, arcsinh, arccosh, arctanh, exp, pi, e

# %% VFA GUI

# start the timer
start = timeit.default_timer()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define needed functions for the initial plot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# define a function that will search for singularities and mark them on a very
# fine grid, but looking point by point without saving
def singularity_fine(u_str, v_str, dist_points, N=500):
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


# define a function that will take inputs from under the plot button
# and upadte them globally, as well as thier dependant variables
def update_variables(r3_check=0):
    global L, pt_den, s_max, x, y, z, xg, yg, zg, u, v, string_x, string_y, zero_field, field_unit
    global ax_L, main_axis, form_2
    # take the new axis parameters and field definitions out of the boxes
    L = eval(L_entry.get())
    pt_den = int(pt_den_entry.get())
    s_max = int(s_max_entry.get())
    string_x = str(x_comp_entry.get())
    string_y = str(y_comp_entry.get())
    # from L redefine the axis
    ax_L = L + L/delta_factor
    main_axis.set_xlim(-ax_L, ax_L)
    main_axis.set_ylim(-ax_L, ax_L)
    # from pt_den and L, change the axis coordinates and the grid:
    if r3_check == 0:
        x = np.linspace(-L, L, pt_den)
        y = np.linspace(-L, L, pt_den)
        xg, yg = np.meshgrid(x, y)
        # define zero field
        zero_field = np.zeros(np.shape(xg))
        # define the unit field
        field_unit = np.ones(np.shape(xg))
        # get 2-form to update too
        form_2 = eval(form_2_eq)
    else:
        x = np.linspace(-L, L, pt_den)
        y = np.linspace(-L, L, pt_den)
        z = np.linspace(-L, L, pt_den)
        xg, yg, zg = np.meshgrid(x, y, z)
        # define zero field
        zero_field = np.zeros(np.shape(xg))
        # define the unit field
        field_unit = np.ones(np.shape(xg))


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
def stack_plot(xg, yg, axis, F_x, F_y, s_max, L, pt_den, fract, arrows=False, stacks=True, orientation='mid', scale=1, w_head=1/8, h_head=1/4, axis_check=0, arrowheads=True, colour='green', check_2_frm=0, s_min=2):
    global s_L
    # get the lengths of x and y from their grids
    
    x_len = len(xg[:, 0])
    y_len = len(yg[0, :])
    
    # Account for change to grid centre for divergence plot
    if axis_check == 1:
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
                if check_2_frm == 0:
                    # define it for all magnitudes. Separately for odd and even corr. number of sheets:
                    # Label each element with the number of stacks required: linear scaling
                    for t in range(1, s_max+1):
                        if (t-1)/s_max <= R[i, j] <= t/s_max:
                            R_int[i, j] = t
                else:
                    for t in range(s_min, s_max+1):
                        if (t-2)/s_max <= R[i, j] <= (t-1)/s_max:
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
    string = string.replace('ln', 'log')
    return string


# function to unformat equation to display back to user in same form
def unformat(string):
    # replace all the x and y with xg and yg:
    string = string.replace('xg', 'x')
    string = string.replace('yg', 'y')
    string = string.replace('zg', 'z')
    # where there are special functions, replace them with library directions
    string = string.replace('**', '^')
    string = string.replace('log', 'ln')
    string = string.replace('exp', 'e**')
    string = string.replace('field_unit', '1')
    return string


# define a function that takes input string that is python understood and turn into vector components:
def eq_to_comps(string_x, string_y, xg, yg, check_2_frm=0):
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
    if check_2_frm == 1:
        if equation_x == '0' and equation_y != '0':
            u = np.ones(np.shape(xg))
        if equation_y == '0' and equation_x != '0':
            v = np.ones(np.shape(yg))
    else:
        pass
    # return these
    return u, v


# define a function that will take care of constant and zero 2-forms
def form_2_constant_correction(form_2_eq):
    # want to check if it contains x or y and if not, make the shape correct
    # without using eq_to_comps becuase we do not change it to any components
    # here, although the process is very similar
    if form_2_eq.find('xg') & form_2_eq.find('yg') == -1:
        form_2_eq = str(form_2_eq) + '* np.ones(np.shape(xg))'
    else:
        pass
    return form_2_eq


# set a variable thet tracks use of R3 tab
R3_use_track = 0

# dynamics tracker
dyn_use_track = 0

# define a function to take care of tab changes
def tab_selection(event):
    global x_m, y_m, R3_use_track, click_opt_int, tab_text
    global fract, form_2, notebook_bottom, m
    global expressions, coords, form_2_str, form_2_eq, toolbar
    global dyn_use_track
    # switching out of dynamics, stop the animation and reset:
    if dyn_use_track == 1:
        try:
            clear_response()
        except NameError:
            pass
    else:
        pass
    # get tab that was selected as text
    selected_tab = event.widget.select()
    tab_text = event.widget.tab(selected_tab, "text")
    # when R3 is unclicked, reinitialise 2D variables correctly
    if R3_use_track == 1:
        # first, undo it.
        R3_use_track = 0
        # then change all the variables back to what they need to be:
        update_variables()
        m = 2
        coords = ['x', 'y']
        form_2_str = 'x*y**2'  # dx^dy component, (back to R2 default)
        form_2_eq = format_eq(form_2_str)
        form_2 = eval(form_2_eq)
        # 2-form was returned to defualt on R2 so put that in the entry box
        form_2_entry.delete(0, 'end')
        form_2_entry.insert(0, form_2_str)
        # restore the fields input frame
        notebook_bottom.add(bot_frame, text='fields')
        notebook_bottom.select(0)
        notebook_singular.add(singular_frame, text='singularities')
        notebook_singular.select(0)
        # enable the plot and polar buttons again
        PLOT_btn['state'] = tk.NORMAL
        polar_grid_plot_btn['state'] = tk.NORMAL
    if tab_text == 'Line Integrals':
        global LI_coord, LI_total, tensor, click_opt_int
        global LI_selected
        x_m = None
        y_m = None
        fig.canvas.draw()
        # Initialise a global variable for storing click coordinates
        # and total line integral
        # home the main screen
        toolbar.home()
        # unclick any chosen options from toolbar
        state = fig.canvas.toolbar.mode
        if state == 'zoom rect':
            toolbar.zoom()
        if state == 'pan/zoom':
            toolbar.pan()
        # get rid of the 2 buttons we don't want
        toolbar.children['!button4'].pack_forget()
        toolbar.children['!button5'].pack_forget()
        # set variable for mouse to respond to integrals in new tab
        LI_selected = True
        # chenage fract to a better one for 1-forms
        fract = 0.05
        # deal with consistancy if both was initially selected
        if tensor.get() == 2:
            tensor.set(0)
            vect_type_response(tensor.get())
        # return the option to default
        LI_shape_select_response('Polygon')
        # draw it on
        canvas.draw()
        # enable plot buttons again
        PLOT_btn['state'] = tk.NORMAL
        polar_grid_plot_btn['state'] = tk.NORMAL
    elif tab_text == 'VF':
        LI_restart()
        # by default return to initial 'tools'
        click_opt_int = 0
        click_option.set(0)
        click_option_handler(click_option.get())
        # chenage fract to a better one for 1-forms
        fract = 0.05
        # define x, y and z values
        # restart the plot, in case calculus code was used last
        PLOT_response()
        # get rid of the grid
        main_axis.grid(False)
        # draw it on
        canvas.draw()
        # enable plot buttons again
        PLOT_btn['state'] = tk.NORMAL
        polar_grid_plot_btn['state'] = tk.NORMAL
    elif tab_text == 'Ext. Alegebra':
        global calculus_form_tracker, R2_tools_opt
        main_axis.clear()
        # get globals form entry boxes
        update_variables()
        # from these, establish the new fract, approperiate for 2-forms
        fract = 2/((pt_den-1))
        # update variable to track 1 forms and 2 forms in this tab
        calculus_form_tracker = 2
        # put the initial plot onto the canvas
        form_2_components_plot(xg, yg, form_2/2, np.pi/2, s_max, L, fract, colour_str, 2)
        form_2_components_plot(xg, yg, form_2/2, 0, s_max, L, fract, colour_str, 2)
        canvas.draw()
        # disable the PLOT button and the ploar grid button
        PLOT_btn['state'] = tk.DISABLED
        polar_grid_plot_btn['state'] = tk.DISABLED
        # clear the zero form label as plot comes back to default
        form_0_entry.configure(bg='#FFFFFF')
        # this is now never Vector fields and never arrows therefore
        # set the labels as such
        component_x_entry_label.configure(text='dx component')
        component_y_entry_label.configure(text='dy component')
        field_select_drop_label.configure(text='Select Pre-Defined 1-Form:')
        # change frame name too
        bot_frame_frame.configure(text='1-Form input frame')
        # colour the 1-form inputs white as these are not default plotin this
        # tab
        x_comp_entry.configure(bg='#FFFFFF')
        y_comp_entry.configure(bg='#FFFFFF')
        # set the tools variable as it should be
        R2_tools_opt.set(0)
        # respond to that too
        R2_tools_handler(R2_tools_opt.get())
    elif tab_text == 'R^3':
        global form_2_frame
        global F_xy_x, F_xy_y, F_xz_x, F_xz_z, F_yz_y, F_yz_z
        global max_global_dxdy, max_global_dxdz, max_global_dydz
        global h_index, hvalue_string
        # set the variable that tracks when R3 is used to 1
        R3_use_track = 1
        # clear the axis
        main_axis.clear()
        # read inpout variables and change them globally
        update_variables(1)
        # from these, establish the new fract, approperiate for 2-forms
        fract = 2/((pt_den-1))
        # set dimensionality
        m = 3
        # set strings for x and y back to default
        string_x = 'y*sin(x)'
        string_y = '- x*cos(y)'
        string_z = 'x*y*z'
        # update these in their boxes, to return to defualt each
        # time this tab opens, as the plot does that
        form_1_x_entry.delete(0, 'end')
        form_1_y_entry.delete(0, 'end')
        form_1_z_entry.delete(0, 'end')
        form_1_x_entry.insert(0, string_x)
        form_1_y_entry.insert(0, string_y)
        form_1_z_entry.insert(0, string_z)
        # put them into their correct 
        # parse them into expression
        sympy_expr_x = parse_expr(string_x, evaluate=False)
        sympy_expr_y = parse_expr(string_y, evaluate=False)
        sympy_expr_z = parse_expr(string_z, evaluate=False)
        # combine the 2 into a list:
        expressions = np.array([sympy_expr_x, sympy_expr_y, sympy_expr_z])
        # set 3D coords for sympy
        coords = ['x', 'y', 'z']
        # from these, use the find_2_form function to get the 2-form
        form_2 = find_2_form(expressions, coords, pt_den, m)
        # get it's string
        form_2_str_dxdy = str(simplify(str(unformat(result[0][0]))))
        form_2_str_dxdz = str(simplify(str(unformat(result[1][0]))))
        form_2_str_dydz = str(simplify(str(unformat(result[2][0]))))
        # get global maxima for correct tubing
        max_global_dxdy = find_global_max(form_2_str_dxdy)
        max_global_dxdz = find_global_max(form_2_str_dxdz)
        max_global_dydz = find_global_max(form_2_str_dydz)
        # define components in the x-y plane:
        eq_1_xy = str(simplify('(' + form_2_str_dxdy + ')' + '*(' + str(1/2) + ')'))
        eq_2_xy = str(simplify('(' + form_2_str_dxdy + ')' + '*(' + str(1/2) + ')'))
        F_xy_x, F_xy_y = eq_to_comps(eq_1_xy, eq_2_xy, xg, yg, check_2_frm=1)
        # define them in x-z plane
        eq_1_xz = str(simplify('(' + form_2_str_dxdz + ')' + '*(' + str(1/2) + ')'))
        eq_2_xz = str(simplify('(' + form_2_str_dxdz + ')' + '*(' + str(1/2) + ')'))
        F_xz_x, F_xz_z = eq_to_comps(eq_1_xz, eq_2_xz, xg, yg, check_2_frm=1)
        # define them in y-z plane
        eq_1_yz = str(simplify('(' + form_2_str_dydz + ')' + '*(' + str(1/2) + ')'))
        eq_2_yz = str(simplify('(' + form_2_str_dydz + ')' + '*(' + str(1/2) + ')'))
        F_yz_y, F_yz_z = eq_to_comps(eq_1_yz, eq_2_yz, xg, yg, check_2_frm=1)
        # plot the starting field with desired parameters as specidfied above
        # arrow params not needed as arrows arent plotted
        # starting with viewing axis ='z' therefore xy plane
        form_2_components_plot_3(xg, yg, h_index, axis_view, F_xy_x, zero_field, s_max, L, pt_den, fract, colour_str, 2)
        form_2_components_plot_3(xg, yg, h_index, axis_view, zero_field, F_xy_y, s_max, L, pt_den, fract, colour_str, 2)
        canvas.draw()
        # into the entry boxes for 2-forms, but the result down
        Entry_2_form_R3_dxdy.delete(0, 'end')
        Entry_2_form_R3_dxdz.delete(0, 'end')
        Entry_2_form_R3_dydz.delete(0, 'end')
        Entry_2_form_R3_dxdy.insert(0, str(form_2_str_dxdy))
        Entry_2_form_R3_dxdz.insert(0, str(form_2_str_dxdz))
        Entry_2_form_R3_dydz.insert(0, str(form_2_str_dydz))
        # colour what plots as defualt in this tab
        Entry_2_form_R3_dxdy.configure(bg='#C0F6BB')
        Entry_2_form_R3_dxdz.configure(bg='#C0F6BB')
        Entry_2_form_R3_dydz.configure(bg='#C0F6BB')
        # hide the VFA input frame for now
        notebook_bottom.hide(0)
        # hide singularity notebook
        notebook_singular.hide(0)
        # disable the PLOT button and the ploar grid button
        PLOT_btn['state'] = tk.DISABLED
        polar_grid_plot_btn['state'] = tk.DISABLED
        # get the correct h index asnd update its label
        try:
            hvalue_string = str(z[h_index])
            axis_height_txt.configure(text=str(round(eval(hvalue_string), 2)))
        except IndexError:
            # indexes run out, set the index to the first
            h_index = 0
            # redo the label update
            hvalue_string = str(z[h_index])
            axis_height_txt.configure(text=str(round(eval(hvalue_string), 2)))
            # show warning about that
            tk.messagebox.showwarning('INDEX ERROR', 'you selected a value in slider that might no longer be avalaible, it has been changed avoid errors')
        # this is now never Vector fields and never arrows therefore
        # set the labels as such
        component_x_entry_label.configure(text='dx component')
        component_y_entry_label.configure(text='dy component')
        field_select_drop_label.configure(text='Select Pre-Defined 1-Form:')
        # change frame name too
        bot_frame_frame.configure(text='1-Form input frame')
    # if the Dynamics tab is selected:
    elif tab_text == 'Dynamics':
        dyn_use_track = 1
        # No need for tools here:
        toolbar.home()
        state = fig.canvas.toolbar.mode
        if state == 'zoom rect':
            toolbar.zoom()
        if state == 'pan/zoom':
            toolbar.pan()
        toolbar.children['!button4'].pack_forget()
        toolbar.children['!button5'].pack_forget()
        # only works for 1-forms on R2:
        fract = 0.05
        PLOT_response()
        PLOT_btn['state'] = tk.NORMAL
        polar_grid_plot_btn['state'] = tk.NORMAL
    if tab_text != 'Dynamics':
        dyn_use_track = 0
    # if anything but the main window is selected, change to tools
    if tab_text != 'VF' and tab_text != 'Line Integrals' and tab_text != 'Dynamics':
        # unclick them
        click_option.set(0)
        click_option_handler(click_option.get())


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

style_notebook = ttk.Style()
style_notebook.configure('TNotebook.Tab', font=('URW Gothic L','8','bold') )

# right frame:
right_frame_frame = tk.LabelFrame(root, text='', padx=5, pady=5)
right_frame_frame.grid(row=1, column=1, rowspan=2, sticky='N')

# bot frame:
bot_frame_frame = tk.LabelFrame(root, text='', padx=5, pady=5)
bot_frame_frame.grid(row=2, column=0)

# plot frame:
plot_frame = tk.LabelFrame(root, text='', padx=5, pady=5)
plot_frame.grid(row=1, column=0)

# define notebook for tabs
notebook = ttk.Notebook(right_frame_frame)
notebook.grid(row=0, column=0)
# notebook for bottom field input, when needed to disappear.
notebook_bottom = ttk.Notebook(bot_frame_frame)
notebook_bottom.grid(row=0, column=0)
# singularities notebook
notebook_singular = ttk.Notebook(right_frame_frame)
notebook_singular.grid(row=1, column=0)
# plotting options notebook
notebook_small = ttk.Notebook(bot_frame_frame)
notebook_small.grid(row=0, column=1)

# singularities:
singular_frame = tk.LabelFrame(notebook_singular)
singular_frame.grid(row=0, column=1)
# main options:
right_frame = tk.LabelFrame(notebook)
right_frame.grid(row=0, column=0)
# field input
bot_frame = tk.LabelFrame(notebook_bottom)
bot_frame.grid(row=0, column=0)
# Line integrals
LI_frame = tk.LabelFrame(notebook)
LI_frame.grid(row=0, column=1)
# calculus
calculus_frame = tk.LabelFrame(notebook)
calculus_frame.grid(row=0, column=3)
# R3
r3_frame = tk.LabelFrame(notebook)
r3_frame.grid(row=0, column=4)
# dynamics
dynamics_frame = tk.LabelFrame(notebook)
dynamics_frame.grid(row=0, column=2)
# plotting options
small_frame = tk.LabelFrame(notebook_small)
small_frame.grid(row=0, column=0)

# finsalise them
notebook.add(right_frame, text='VF')
notebook.add(LI_frame, text='Line Integrals')
notebook.add(dynamics_frame, text='Dynamics')
notebook.add(calculus_frame, text='Ext. Alegebra')
notebook.add(r3_frame, text='R^3')
notebook_bottom.add(bot_frame, text='Fields')
notebook_singular.add(singular_frame, text='singularities')
notebook_small.add(small_frame, text='Plotting')

# bind the clicks on tabs to a function
notebook.bind_all('<<NotebookTabChanged>>', tab_selection)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set up initial parameters and plot the initial graph, put it in plot frame
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

'''  BASIC PARAMETERS  '''

# define scale of the graph
L = 5
pt_den = 21 # number of points on each axis

# Initialise auto scaling variable
ascale = tk.IntVar()
ascale.set(0)

# define x and y values
x = np.linspace(-L, L, pt_den)
y = np.linspace(-L, L, pt_den)
z = np.linspace(-L, L, pt_den)  # here only to begin slider in r3 tab

# create a grid on x-y plane
xg, yg = np.meshgrid(x, y)

# define an example vector field
u = yg*np.sin(xg)  # x component
v = -xg*np.cos(yg)  # y component
# for no dependance in any initial component, use : np.zeros(np.shape(xg)) or yg


'''LIST OF DEFAULT VECTOR FIELDS TO DISPLAY IN DROPDOWN MENU'''


# list of names of fields to display
field_name_list = ['Default: y*sin(x)dx - x*cos(y)dy',
                   'Simple pendulum: ydx  - sin(x)dy',
                   'Harmonic oscillator: ydx -xdy',
                   'Linear example 1: (x + 2*y)dx + (3*x - 4*y)dy',
                   'Linear example 2: xdx',
                   'Constant: 6dx + 3dy',
                   'Falling cat (Planar 3 link robot)',
                   'Electric Point Charge: -x/(x**2+y**2)dx + -y/(x**2+y**2)dy',
                   'H field of Current Carrying Wire: -y/(x**2+y**2)dx + x/(x**2+y**2)dy',
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
                'x + 2*y',
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
                '3*x - 4*y',
                '0',
                '3',
                '-(3*cos(x) + 4)/(15 + 6*cos(x) + 6*cos(y))',
                '-y/(x**2+y**2)',
                'x/(x**2+y**2)',
                'y/(sqrt(x**2 + y**2)*(1-2/(sqrt(x**2 + y**2)))) + x',
                '-2*y*((x^2+y^2)^(-1.5))*(1-(2/sqrt(x^2+y^2)))^(-2)'
                ]


# define 2-form names and equations for the pre-defined 2-forms dropdown
# its on R2, therefore all are only dx^dy
list_form_2_names = ['default: x**2*y',
                     'constant: 3',                     
                     ]

# define its equations:
list_form_2_equations = ['x**2*y',
                         '3',
                         ]


'''  INITIAL PLOT PARAMETERS and FIGURE '''

# set up quiver factors
arrows = False  # set up if arrows should be plotted on stacks or not.
orientation = 'mid'  # how the arrow rotates about its assigned grid point - options: tail, mid and tip as string
scale = 5  # the scale reduction factor, if None (as None-type, not str), automatically computed by average, if 1 = mag

# set up the delta_factor of additional axis space L/delta_factor gives extra space on axis
delta_factor = 10

# fraction of sheet length to graph length
fract = 0.05

# define the maximum number of stack to plot, dep. on magnitude (initialy)
s_max = 6

# set screen dpi
my_dpi = 100

# define denominator of fractional height and width of arrowhead based on stack size
w_head = 1/8
h_head = 1/4

# create a figure, use dpi to fit it more precisely to size of the frame
fig = plt.figure(figsize=(730/my_dpi, 573/my_dpi), dpi=my_dpi)

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

# define a variable to track showflux in LI tab
# showflux = tk.IntVar()
# showflux.set(0)

# showcirc = tk.IntVar()
# showcirc.set(0)

# define a variable to track shading areas in Calclulus area inetgrals
shadearea = tk.IntVar()
shadearea.set(1)

# set up line intergal enpty variables
LI_coord = []
LI_total = 0
flux = 0
shape_area = 0
ratio1 = 0
ratio2 = 0

# set up initial variables for 2_form integral
AI_coord = []
AI_area = 0
AI_result = 0
AI_verts = []


''' DEFINE CALCULUS PARAMETERS, 2-forms and R3 stuff '''

# set up initial strings for 2-forms window to display, for it to save properly after
to_wedge_x_1_str = ''
to_wedge_y_1_str = ''
to_wedge_x_2_str = ''
to_wedge_y_2_str = ''

# Initialise the click button selection
click_opt_int = 0

# variable to keep track of what is being plotted in the calculus code
# a 1 form or a 2 form
# starts set to 2 form, as that is the default when opening the tab
calculus_form_tracker = 2

# polar tracker
polar_tracker = False

# define initial pixes coordinates for the function to use
# when mouse is clicked, these are changed to the pressed position
x_m = float(0)
y_m = float(0)

# initial orientation for circle integration:
orient_int = 'ccw'

# define initial stack bool
stacks = True

# set the dimensionality
m = 2

# start the tab with a 2-form being supplied and plotted by blocks only
# define intiial variables for these:
# define the initial 2-form, in terms of a string, equation and numerically
form_2_str = 'x*y**2'  # dx^dy component
form_2_eq = format_eq(form_2_str)
form_2 = eval(form_2_eq)

''' initial variables for operations completed by the calculus GUI '''

# ##### GENERAL ONES FOR STACKS #######

# define colours to use for 2-form plots with components
# first string defines colour for positive (ccw), second for negative (cw)
# last one is an in case, for when the magnitude is exactly zero.
colour_str = ['red', 'blue', 'grey']

# set up a zero vector filed to plot x and y components as 2 separate fields:
zero_field = np.zeros(np.shape(xg))

# define a unit field so that code deals with constant fields
field_unit = np.ones(np.shape(xg))  # need it here cuz its used in eq_to_comps and in find_2_form

# ###### ONES NEEDED FOR 2 FORM FROM 1 FORMS #######

# define an example vector field, now - from string, initially
string_x = 'sin(x+y) - y'  # x component
string_y = 'sin(x+y)'  # y component
string_z = 'x*y*x'  # also define it for it to display correctly in R3 tab

# take the input strings and turn them into sympy expressions to be able to
# use sympy's partial differentiation
sympy_expr_x = parse_expr(string_x, evaluate=False)
sympy_expr_y = parse_expr(string_y, evaluate=False)
# for m > 2, need more components, and need these in 'expressions' too!

# combine the 2 into a list:
expressions = np.array([sympy_expr_x, sympy_expr_y])
# will use sympy partial derrivatives on these, as to get a 2-form on R2:
# need to differentiate each component w.r.t the coordinates that it's
# elementary 1-form does not contain.

# set up an array of coordinates that need to be used (in standard order)
coords = ['x', 'y']

# predefine vectors to store for the interior derivative
vector_ex_str = ''
vector_ey_str = ''

''' INITIALISE VARIABLES FOR R3 CODE '''

# to start with, set as viewing aling z axis onto x-y plane
axis_view = 'z'
# set the initial index of height of the viewed plane along viewing axis
h_index = 0
hvalue_string = '-5'


'''  Initialise variables for the Dynamics tab  '''

# store click coordinates
dyn_coord = []

# define array of points
dyn_point, = main_axis.plot([], [], 'ro', markersize=4)

# define initial variables needed by dynamics
tmax = 1
dyn_N = 100
dyn_time = np.linspace(0, tmax, dyn_N)

# bool to det if shapoes are to be drawn
dyn_join_shapes = tk.IntVar()
dyn_join_shapes.set(1)  # start by drawing shapes


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
        deriv_calc(x_m, y_m)
    # separately deal with needed click respponses in LI tab
    if tab_text == 'Line Integrals':
        if LI_shape_select.get() == 'Polygon':
            # Store the coordinates of the click in list
            LI_coord.append([x_m, y_m])
            line_int_poly(100000, string_x, string_y)
        # elif LI_shape_select.get() == 'square':
        #     print('not implemented')   # insert function later
        elif LI_shape_select.get() == 'Circle':
            # get the radius and call approperiate function
            Radius_LI_circ = eval(Radius_LI_circ_entry.get())
            # if its negative correct and put it into the entry box
            if Radius_LI_circ < 0:
                Radius_LI_circ *= -1
                Radius_LI_circ_entry.delete(0, 'end')
                Radius_LI_circ_entry.insert(0 , str(Radius_LI_circ))
            line_int_circ([x_m,y_m], Radius_LI_circ, 100000, string_x, string_y, orient_int)
    # and when want to use zoom or tools in calculus:
    if tab_text == 'Ext. Alegebra':
        if R2_tools_opt_int == 1:
            form_2_zoom(x_m, y_m)
        elif R2_tools_opt_int == 0:
            key_press_handler(event, canvas, toolbar)
        else:
            if R2_flux_shape.get() == 'Polygon':
                # only restart if shape has already been completed
                if shape_complete_tracker == 0:
                    # if the shape is not complete, just continue normally
                    pass
                else:
                    # if the shape has already been completed, restart first
                    # then continue
                    AI_restart()
                # call function to calculate
                area_finder_form_2_int(x_m, y_m)
            elif R2_flux_shape.get() == 'Circle':
                # restart
                AI_restart()
                # get the radius and call approperiate function
                Radius_R2_circ = eval(Radius_R2_circ_entry.get())
                # if its negative correct and put it into the entry box
                if Radius_R2_circ < 0:
                    Radius_R2_circ *= -1
                    Radius_R2_circ_entry.delete(0, 'end')
                    Radius_R2_circ_entry.insert(0 , str(Radius_LI_circ))
                integration_from_2_circle(x_m, y_m, Radius_R2_circ, AI_verts)
            else:
                pass
    if tab_text == 'Dynamics':
        # clicking should only register a new point and plot it:
        global dyn_coord
        if dyn_shape_select.get() == 'Polygon':
            dyn_coord.append([x_m, y_m])
            dot = patch.Circle(dyn_coord[-1], L*fract/6, color='red')
            main_axis.add_patch(dot)
            canvas.draw()
        elif dyn_shape_select.get() == 'Square':
            # clear
            clear_response()
            # get the input size
            dyn_pre_size = eval(dyn_pre_size_entry.get())
            if dyn_pre_size < 0:
                dyn_pre_size *= -1
                dyn_pre_size_entry.delete(0, 'end')
                dyn_pre_size_entry.insert(0 , str(Radius_LI_circ))
            # from that define a square
            # add these to the axis and draw it
            dyn_coord.append([x_m - dyn_pre_size, y_m - dyn_pre_size])
            dyn_coord.append([x_m - dyn_pre_size, y_m + dyn_pre_size])
            dyn_coord.append([x_m + dyn_pre_size, y_m + dyn_pre_size])
            dyn_coord.append([x_m + dyn_pre_size, y_m - dyn_pre_size])
            # draw these on:
            poly = mpl.patches.Polygon(dyn_coord, fill=True, color='blue')
            main_axis.add_artist(poly)
            canvas.draw()
        elif dyn_shape_select.get() == 'Human':
            # clear
            clear_response()
            # get the input size
            dyn_pre_size = eval(dyn_pre_size_entry.get())
            if dyn_pre_size < 0:
                dyn_pre_size *= -1
                dyn_pre_size_entry.delete(0, 'end')
                dyn_pre_size_entry.insert(0 , str(Radius_LI_circ))
            # from that define a square
            # add these to the axis and draw it
            dyn_coord.append([x_m + dyn_pre_size, y_m])
            dyn_coord.append([x_m, y_m + dyn_pre_size])
            dyn_coord.append([x_m - dyn_pre_size, y_m])
            dyn_coord.append([x_m, y_m - dyn_pre_size])
            dyn_coord.append([x_m, y_m - 2*dyn_pre_size])
            dyn_coord.append([x_m - 1.5*dyn_pre_size, y_m - dyn_pre_size])
            dyn_coord.append([x_m, y_m - 2*dyn_pre_size])
            dyn_coord.append([x_m + 1.5*dyn_pre_size, y_m - dyn_pre_size])
            dyn_coord.append([x_m, y_m - 2*dyn_pre_size])
            dyn_coord.append([x_m, y_m - 3*dyn_pre_size])
            dyn_coord.append([x_m - 1.5*dyn_pre_size, y_m - 4*dyn_pre_size])
            dyn_coord.append([x_m, y_m - 3*dyn_pre_size])
            dyn_coord.append([x_m + 1.5*dyn_pre_size, y_m - 4*dyn_pre_size])
            dyn_coord.append([x_m, y_m - 3*dyn_pre_size])
            dyn_coord.append([x_m, y_m - dyn_pre_size])
            # draw these on:
            poly = mpl.patches.Polygon(dyn_coord, fill=True, color='blue')
            main_axis.add_artist(poly)
            canvas.draw()
        elif dyn_shape_select.get() == 'Circle':
            # clear
            clear_response()
            # get the input size
            dyn_pre_size = eval(dyn_pre_size_entry.get())
            if dyn_pre_size < 0:
                dyn_pre_size *= -1
                dyn_pre_size_entry.delete(0, 'end')
                dyn_pre_size_entry.insert(0 , str(Radius_LI_circ))
            # from that define a Circle as a patch
            cirlce_dyn_patch = patch.CirclePolygon((x_m, y_m), dyn_pre_size, fill=False, color='red')
            # draw it on
            main_axis.add_artist(cirlce_dyn_patch)
            canvas.draw()
            # from that patch, extract points and add them to dyn_coord
            verts = cirlce_dyn_patch.get_path().vertices
            trans = cirlce_dyn_patch.get_patch_transform()
            points = trans.transform(verts)
            for i in range(len(points)):
                dyn_coord.append([points[i, 0], points[i, 1]])
        else:
            pass


# connect figure event to a function that responds to clicks, defined above
cid = fig.canvas.mpl_connect("button_press_event", on_key_press)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# define other needed functions, for input reponses
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


''' LINE INTEGRALS'''


# define a funciton to restart line integral calculation and lines
def LI_restart():
    global LI_total, LI_coord, shape_area, swaplistpoly, col_poly, flux
    # first, initialise variables again
    LI_coord = []
    LI_total = 0
    flux = 0
    shape_area = 0
    ratio1 = 0
    ratio2 = 0
    swaplistpoly = []
    col_poly = []
    
    # update the labels
    LI_total_label.configure(text=LI_total)
    flux_label.configure(text=flux)
    shape_area_label.configure(text=shape_area)
    ratio1_label.configure(text=ratio1)
    ratio2_label.configure(text=ratio2)
    # NOT IDEAL BUT HOPEFULLY TEMPORARY
    # call plot-respose to redraw (so lines are deleted)
    # ideally would want a way of deleting lines without restarting the plot
    PLOT_response()


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


# define a function that will change grid separation when in poly
def poly_grid_submit():
    # get the entry box input
    grid_separation_str = grid_sep_poly_entry.get()
    # evaluate it and make sure to avoid negatives and zero
    grid_separation = eval(grid_separation_str)
    if grid_separation < 0:
        grid_separation *= -1
    elif grid_separation == 0:
        tk.messagebox.showwarning('ERROR', 'Can\' t have zero separation')
        # set it to something managable
        grid_separation = 1
    else:
        pass
    # get an array of values for grid to plot at
    major_ticks = np.arange(-L, L + 0.01, grid_separation)
    # update the grid based on that
    main_axis.grid(False)  # remove initial
    main_axis.set_xticks(major_ticks)  # update ticks
    main_axis.set_yticks(major_ticks)
    main_axis.grid(True)
    canvas.draw()
    # in case original input was changed, update it
    grid_sep_poly_entry.delete(0, 'end')
    grid_sep_poly_entry.insert(0, str(round(grid_separation, 6)))


# define LI shape selection response to dropdown
def LI_shape_select_response(selected_shape):
    # deal with lines
    global LI_shape_selected, Radius_LI_circ_entry, Radius_LI_label, orient_int_btn
    global grid_sep_poly_entry, grid_LI_poly_label, submit_poly_sep_grid
    # get the chosen shape
    LI_shape_select.set(selected_shape)
    # depending on that selection, prepare rest of frame:
    if selected_shape == 'Circle':
        # restart the plot and the integral values
        LI_restart()
        # set up the grid
        poly_grid_submit()
        # if circle is selected, display an entry box for the radius of it
        Radius_LI_label = tk.Label(LI_frame, text='Circle Radius:')
        Radius_LI_label.grid(row=5, column=0)
        Radius_LI_circ_entry = tk.Entry(LI_frame, width=10)
        Radius_LI_circ_entry.grid(row=5, column=1)
        Radius_LI_circ_entry.insert(0, '1')
        # also define a button to change orientation
        # define a button that will change circle orientation
        orient_int_btn = tk.Button(LI_frame, text=orient_int, command=orient_int_response)
        orient_int_btn.grid(row=5, column=2)
        # enable the toggle switch
    else:
        # restart the plot and the integral values
        LI_restart()
        # set up the grid
        poly_grid_submit()
        # get rid of circle options
        try:
            Radius_LI_label.destroy()
            Radius_LI_circ_entry.destroy()
            orient_int_btn.destroy()
        except (UnboundLocalError, NameError):  
            pass
        # disable the toggle, as it does not work with polygons (at least yet)
    # update canvas
    canvas.draw()
        

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
    dx_norm = np.zeros(shape=(1, N))
    dy_norm = np.zeros(shape=(1, N))
    
    # Magnitude of increment equals that of the circumference subtended by dt for large N
    A = 2*np.pi*R/(N)
    
    # Loop through to assign coordinates to interval points and plot them
    # depending on chosen orientation
    intervals[0, :] = xc + R*np.cos(dt)
    intervals[1, :] = yc + R*np.sin(dt)

    # get the points along the circle and save
    uv_store[0, :] = eval(format_eq(u_str, LI=1))
    uv_store[1, :] = eval(format_eq(v_str, LI=1))
    
    # Increment vector components
    if orient_int == 'ccw':  
        dx = -A*np.sin(dt)
        dy = A*np.cos(dt)
    else:
        dx = A*np.sin(dt)
        dy = -A*np.cos(dt)
    
    dx_norm = A*np.cos(dt)
    dy_norm = A*np.sin(dt)
    
    # res = np.sum(dx[:-1]*uv_store[0, :-1] + dy[:-1]*uv_store[1, :-1])
    res = np.sum(dx[:-1]*uv_store[0, :-1] + dy[:-1]*uv_store[1, :-1])
    flux = np.sum(dx_norm[:-1]*uv_store[0, :-1] + dy_norm[:-1]*uv_store[1, :-1])

    
    # Colouring the circle based on flux
    if showcol.get() > 0:
        
        store = np.zeros(shape=(1,N))
        store_col = np.zeros(shape=(1,N))
        c = 0
        pstring = ''
        swaplist = []
        col_in = ['grey', 'blue', 'red']
        col_list = []
        
        if showcol.get() == 1:
            store[0,:] = dx[:]*uv_store[0, :] + dy[:]*uv_store[1, :]
        elif showcol.get() == 2:
            store[0,:] = dx_norm[:]*uv_store[0, :] + dy_norm[:]*uv_store[1, :]
        
        for i in range(N):
            if store[0,i] < 0:
                store_col[0,i] = 1
                
            elif store[0,i] > 0:
                store_col[0,i] = 2
                
            elif store[0,i] == 0 or store[0,i] == -0:
                store_col[0,i] = 0
                
            if i == 0:
                continue
            else:
                if store_col[0,i-1] != store_col[0,i]:
                    swaplist.append(dt[i])
                    col_list.append(col_in[int(store_col[0,i])])
                    c += 1
    
        if c > 0: 
            swap_ang = (180/np.pi)*np.array(swaplist)
            fe = np.array([swap_ang[0]])
            swap_ang = np.concatenate((swap_ang, fe))
            for a in range(c):
                exec('w' + str(a) + '= mpl.patches.Arc(cent, 2*R, 2*R, 0, swap_ang[a], swap_ang[a+1], fill=False, color=col_list[a], linewidth=3)')
                exec('main_axis.add_artist(w' + str(a) + ')')
        
        else:
            colour2 = col_in[int(store_col[0,0])]
            circle2 = mpl.patches.Circle(cent, R, fill=False, color=colour2, linewidth=3)
            main_axis.add_artist(circle2)
            
    else:
        # Plot the circle
        circle1 = mpl.patches.Circle(cent, R, fill=False, color='black', linewidth=2)
        main_axis.add_artist(circle1)
        
    fig.canvas.draw()
    
    if showcol.get() > 0:
        # Remove arcs from the plot
        if c > 0:
            for a in range(c):
                exec('w' + str(a) + '.remove()')
        else:
            circle2.remove()
    else:
        circle1.remove()
    
    # update the total
    LI_total = res
    
    # Update labels
    LI_total_label.configure(text=(str(round(LI_total, 2)) + ' (' + str(round((LI_total/np.pi), 2)) + ' pi)' ))
    
    flux_label.configure(text=(str(round(flux, 2)) + ' (' + str(round((flux/np.pi), 2)) + ' pi)' ))

    shape_area = np.pi*R**2
    shape_area_label.configure(text=str(round(shape_area, 2)))
    
    ratio1 = LI_total/abs(shape_area)
    ratio1_label.configure(text=str(round(ratio1, 2)))
    
    ratio2 = flux/shape_area
    ratio2_label.configure(text=str(round(ratio2, 2)))
    
    return res


# define a function that will complete the line integral
def line_int_poly(N, u_str, v_str):
    global LI_total, coord_array, coord_diff, LI_verts, flux, swaplistpoly, col_poly
    # set up a conuter to know how many time to complete the sum process
    c_count = len(LI_coord)
    
    # Tolerance for auto-joining lines together i.e. distance below which lines will join
    ctol = 0.02*L
    
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
        main_axis.add_line(Line2D((a[0], b[0]), (a[1], b[1]), linewidth=2, color='black'))
        
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
        
        flux_inc = 0
        
        # Evaluate line integral as sum of vector components multiplied by the small x and y displacements 
        res = dx*np.sum(uv_store[0, :]) + dy*np.sum(uv_store[1, :])
        flux_inc = dy*np.sum(uv_store[0, :]) - dx*np.sum(uv_store[1, :])
        
        if showcol.get() > 0:
        
            store = np.zeros(shape=(1,N))
            store_col = np.zeros(shape=(1,N))
            
            # Store flux at each point on line
            if showcol.get() == 1:
                store[0,:] = dx*uv_store[0, :] + dy*uv_store[1, :]
            elif showcol.get() == 2:
                store[0,:] = dy*uv_store[0, :] - dx*uv_store[1, :]

            c = 0
            col_in = ['grey', 'blue', 'red']
            
            for i in range(N):
                if store[0,i] < 0:
                    store_col[0,i] = 1
                    
                elif store[0,i] > 0:
                    store_col[0,i] = 2
                    
                elif store[0,i] == 0 or store[0,i] == -0:
                    store_col[0,i] = 0
                
                # Assign colour of first segment
                if i == 0:
                    swaplistpoly.append([intervals[0,i], intervals[1,i]])
                    col_poly.append(col_in[int(store_col[0,i])])
                    continue
                else:
                    if store_col[0,i-1] != store_col[0,i]:
                        swaplistpoly.append([intervals[0,i], intervals[1,i]])
                        col_poly.append(col_in[int(store_col[0,i])])
                        c += 1
        
        # update the total
        LI_total += res
        flux += flux_inc
        
        # update its label
        LI_total_label.configure(text=(str(round(LI_total, 2)) + ' (' + str(round((LI_total/np.pi), 2)) + ' pi)' ))    
        
        if len(LI_verts) > 3:
            shape_area = calc_area(LI_verts)
            shape_area_label.configure(text=str(round(abs(shape_area), 2)))
            
            if shape_area < 0:
                shape_area = (-1)*shape_area
                flux = (-1)*flux
                
                if showcol.get() == 2:
                    for b in range(len(col_poly)):
                        if col_poly[b] == 'blue':
                            col_poly[b] = col_poly[b].replace('blue', 'red')
                        elif col_poly[b] == 'red':
                            col_poly[b] = col_poly[b].replace('red', 'blue')
            else:
                pass
            
            if showcol.get() > 0:
                swaplistpoly.append(LI_verts[0])
                for a in range(len(col_poly)):
                    main_axis.add_line(Line2D((swaplistpoly[a][0], swaplistpoly[a+1][0]), (swaplistpoly[a][1], swaplistpoly[a+1][1]), linewidth=3, color=col_poly[a]))
            
            ratio1 = LI_total/shape_area
            ratio1_label.configure(text=str(round(ratio1, 2)))
            
            ratio2 = flux/shape_area
            ratio2_label.configure(text=str(round(ratio2, 2)))
            
            flux_label.configure(text=(str(round(flux, 2)) + ' (' + str(round((flux/np.pi), 2)) + ' pi)' ))
            
        # display the drawn on lines
        canvas.draw()
            
        return res
    
    else:
        # If only clicked once, plot small red circle
        circle = patch.Circle(LI_coord[0], L*fract/6, color='red')
        main_axis.add_patch(circle)
        canvas.draw()


# Calc polygon (signed) area when user completes a shape
# CCW is +ve, CW is -ve

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
    A = 0.5*(S1-S2)
    return A


''' RESPONSE FUNCTIONS TO PLOT '''


# define a function that will respond to radio buttons behind choosing vector types:
# Note, tensor is a local variable. So to get value in other functions, use tensor.get() to extract current radiobutton value
def vect_type_response(tensor):
    global click_option
    # clear the plot that is already there:
    main_axis.clear()
    # deal with grids if user is in the LI tab
    if tab_text == 'Line Integrals':
        LI_restart()
        # plot the grid
        poly_grid_submit()
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
    # now distinguish clearly between vector fields and 1 forms:
    # when stacks are selected, disable div and curl, these don't exist
    # also disable the derivative, it is not the same:
    # 1 form drrivative is the exterior derivative, not this geometrical one
    # change the labels of inputs too
    # when deselcted, return these
    if tensor == 0 or tensor == 2:
        # stacks included
        # unclick them first, by selecting the default, tools
        # or keeping tools
        if click_opt_int != 0 and click_opt_int != 1:
            click_option.set(0)
            click_option_handler(click_option.get())
        #disable spproperiately
        click_option_Deriv_btn.configure(state=tk.DISABLED)
        click_option_Div_btn.configure(state=tk.DISABLED)
        click_option_Curl_btn.configure(state=tk.DISABLED)
        # change labels to be inputs of 1 forms
        component_x_entry_label.configure(text='dx component')
        component_y_entry_label.configure(text='dy component')
        field_select_drop_label.configure(text='Select Pre-Defined 1-Form:')
        # change frame name too
        bot_frame_frame.configure(text='1-Form input frame')
    elif tensor == 1:
        # enable approperiately
        click_option_Deriv_btn.configure(state=tk.NORMAL)
        click_option_Div_btn.configure(state=tk.NORMAL)
        click_option_Curl_btn.configure(state=tk.NORMAL)
        # change labels to be inputs of 1 forms
        component_x_entry_label.configure(text='x component')
        component_y_entry_label.configure(text='y component')
        field_select_drop_label.configure(text='Select Pre-Defined Vector Field:')
        # change frame name too
        bot_frame_frame.configure(text='Vector Field input frame')


# define the PLOT button response function
def PLOT_response(test_for_clearing_dyn=1):
    # first, take from entry boxes, wanted parameters and make them global:
    # these must be changed globally for other functions to work with the new field.
    global u, v, arrows, stacks, polar_tracker, dyn_coord
    # take inputs and globally update them
    update_variables()
    # set radial tracker
    polar_tracker = False
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
    # clear the current axis
    main_axis.clear()
    # deal with grids if user is in the LI tab
    if tab_text == 'Line Integrals':
        global LI_total, LI_coord, shape_area
        # first, initialise variables again
        LI_coord = []
        LI_total = 0
        flux = 0
        shape_area = 0
        ratio1 = 0
        ratio2 = 0
        # update the labels
        LI_total_label.configure(text=LI_total)
        flux_label.configure(text=flux)
        shape_area_label.configure(text=shape_area)
        ratio1_label.configure(text=ratio1)
        ratio2_label.configure(text=ratio2)
        # plot the grid
        poly_grid_submit()
    if tab_text == 'Dynamics':
        if test_for_clearing_dyn == 1:
            try:
                for a in range(len(dyn_coord)):
                    exec('global ' + 'xy' + str(a) + '\n' + 'del ' + 'xy' + str(a))
            except NameError:
                pass
            # then clear coordinates
            dyn_coord = []
    # create a figure and display it
    stack_plot(xg, yg, main_axis, u, v, s_max, L, pt_den, fract, arrows, stacks, orientation, scale, w_head, h_head, 0)
    canvas.draw()
    # recolour pt_den to white, if it was red from polar plots
    pt_den_entry.configure(bg='white')
    # colour the x and y boxes green to show that these plot
    x_comp_entry.configure(bg='#C0F6BB')
    y_comp_entry.configure(bg='#C0F6BB')


# define a function that will respons to field selection in the drop down menu
def field_selection_response(event):
    global u, v, fract, calculus_form_tracker, polar_tracker, arrows, stacks
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
    # colour code to be able to distinguish what is being plotted
    x_comp_entry.configure(bg='#C0F6BB')
    y_comp_entry.configure(bg='#C0F6BB')
    # now call the plot function to finalise all these onto the plot
    # depending on tab ,use correct 1 form plotting fucntion
    if tab_text == 'Ext. Alegebra':
        # this plots 1 form always for these responses
        form_1_stacks_response()
    else:
        # check if the selected field is stricte a 1-form
        # if so, change representation.
        if selected_index == 7 or selected_index == 8 or selected_index == 9 or selected_index == 10:
            # select stacks to be plotted
            tensor.set(0)
            # respons to this by removing options unavaliable to 1-forms in
            # main tab:
            if click_opt_int != 0 and click_opt_int != 1:
                click_option.set(0)
                click_option_handler(click_option.get())
            click_option_Deriv_btn.configure(state=tk.DISABLED)
            click_option_Div_btn.configure(state=tk.DISABLED)
            click_option_Curl_btn.configure(state=tk.DISABLED)
            component_x_entry_label.configure(text='dx component')
            component_y_entry_label.configure(text='dy component')
            field_select_drop_label.configure(text='Select Pre-Defined 1-Form:')
        # then, with all these set, call the plot function.
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
    # DO not actually replot, just save these as globals
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
    submit_arr_btn = tk.Button(arrowH_opt_window, text='SAVE ALL', padx=20, pady=10, command=custom_submission)
    submit_arr_btn.grid(row=8, column=0, pady=10)


# define a response funcction to autoscale toggle button
def scale_toggle_response():
    global ascale
    if ascale.get() == 0:
        # the burron is off, and has been clicked therefore change the
        # variable to an and the image to on
        ascale.set(1)
        ascale_toggle.configure(image=toggle_image_on)
        ascale_toggle_LI.configure(image=toggle_image_on)
        # for it to update, reclick whatever radiobutton is selected
        # or, if stacks only is chosen, change it to both, to show some change
        vect_type_response(tensor.get())
    else:
        # the button is on and has been clicked
        # set it to off and change image
        ascale.set(0)
        ascale_toggle.configure(image=toggle_image_off)
        ascale_toggle_LI.configure(image=toggle_image_off)
        # for it to update, reclick whatever radiobutton is selected
        # or, if stacks only is chosen, change it to both, to show some change
        vect_type_response(tensor.get())


''' POLAR PLOTS '''

# define a function to repond to plotting apolar grid
# takes the same field, but plots it on a polar grid
def Polar_grid_plot_response(tensor):
    global xg, yg, u, v, s_max, pt_den_entry, polar_tracker
    # set the polar tracker
    polar_tracker = True
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
    # deal with grids if user is in the LI tab
    if tab_text == 'Line Integrals':
        global LI_total, LI_coord, shape_area
        # first, initialise variables again
        LI_coord = []
        LI_total = 0
        flux = 0
        shape_area = 0
        ratio1 = 0
        ratio2 = 0
        # update the labels
        LI_total_label.configure(text=LI_total)
        flux_label.configure(text=flux)
        shape_area_label.configure(text=shape_area)
        ratio1_label.configure(text=ratio1)
        ratio2_label.configure(text=ratio2)
        # plot the grid
        poly_grid_submit()
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
    # colour the x and y boxes green to show that these plot
    x_comp_entry.configure(bg='#C0F6BB')
    y_comp_entry.configure(bg='#C0F6BB')


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


'''

DIFFERENTIAL CALCULUS FUNCTIONS

'''


# define a function that will calculate the local, geometrical derivative
def deriv_calc(x_m, y_m):
    # Make u_s and v_s global for testing
    global u_s, v_s
    # Range and point density of the derivative plot
    d_range = 0.33*L/(zoom_slider.get())
    d_length = d_length_select.get()*L
    dpd = dpd_select.get()
    d_scale = scale*(zoom_slider.get())
    
    # Initialise arrays for storing components
    u_s = np.zeros(shape=(dpd, dpd))
    v_s = np.zeros(shape=(dpd, dpd))
    
    # define new axis in the derivative plot
    dx = np.linspace(-d_range+x_m, d_range+x_m, dpd)
    dy = np.linspace(-d_range+y_m, d_range+y_m, dpd)
    dxg, dyg = np.meshgrid(dx, dy)
    
    # define the vector field in these new axis
    u_zoom, v_zoom = eq_to_comps(string_x, string_y, dxg, dyg)

    # Define the components of the derivative field
    V = u_zoom - eval(format_eq_div(format_eq(string_x)))
    W = v_zoom - eval(format_eq_div(format_eq(string_y)))
    
    if click_opt_int == 1:
        u_s = u_zoom
        v_s = v_zoom
    
    elif click_opt_int == 2:
        u_s = V
        v_s = W
        
    if click_opt_int > 2:
        # get corrected dpd (for loops)
        N = dpd - 1
        
        # get number of points in quadrant
        if dpd % 2 == 1:
            quad_x = int(dpd/2)
            quad_y = int((dpd+1)/2)
        else:
            quad_x = int(dpd/2)
            quad_y = int(dpd/2)
        
        
        # loop over quadrants, making sure to preserve squares for pi/2 rotations
        # OUTER LOOP
        for i in range(quad_x):
            # get the l number, for projection of j on radial / i on tangent
            l = i - 0.5*(N)
            
            # INNER LOOP
            for j in range(quad_y):
                # get the k number of projection: i on radial / j on tangent
                k = j - 0.5*(N)
                
                # get the commuting parts of V and W for each square corner
                # (x and y components of the subtracted field)
                V_comm_1 = 0.25*(2*V[i, j] + W[j, N-i] - W[N-j, i])
                V_comm_2 = 0.25*(2*V[j, N-i] + W[N-i, N-j] - W[i, j])
                V_comm_3 = 0.25*(2*V[N-i, N-j] + W[N-j, i] - W[j, N-i])
                V_comm_4 = 0.25*(2*V[N-j, i] + W[i, j] - W[N-i, N-j])
                
                W_comm_1 = 0.25*(2*W[i, j] - V[j, N-i] + V[N-j, i])
                W_comm_2 = 0.25*(2*W[j, N-i] - V[N-i, N-j] + V[i, j])
                W_comm_3 = 0.25*(2*W[N-i, N-j] - V[N-j, i] + V[j, N-i])
                W_comm_4 = 0.25*(2*W[N-j, i] - V[i, j] + V[N-i, N-j])
                
                # gte a normalisation factor from l and k
                A = k**2 + l**2
                
                # for each corner of the square
                # find the div and curl field in u and in v
                # saving these into their arrays
                
                # Local divergence field
                if click_opt_int == 3:    
                    u_s[i, j] = (V_comm_1*k + W_comm_1*l)*k/A
                    v_s[i, j] = (V_comm_1*k + W_comm_1*l)*l/A
                    u_s[j, N-i] = (V_comm_2*l + W_comm_2*(-k))*l/A
                    v_s[j, N-i] = (V_comm_2*l + W_comm_2*(-k))*(-k)/A
                    u_s[N-i, N-j] = (V_comm_3*(-k) + W_comm_3*(-l))*(-k)/A
                    v_s[N-i, N-j] = (V_comm_3*(-k) + W_comm_3*(-l))*(-l)/A
                    u_s[N-j, i] = (V_comm_4*(-l) + W_comm_4*k)*(-l)/A
                    v_s[N-j, i] = (V_comm_4*(-l) + W_comm_4*k)*k/A
                    
                # Local curl field
                if click_opt_int == 4:        
                    u_s[i, j] = (V_comm_1*l + W_comm_1*(-k))*l/A
                    v_s[i, j] = (V_comm_1*l + W_comm_1*(-k))*(-k)/A
                    u_s[j, N-i] = (V_comm_2*(-k) + W_comm_2*(-l))*(-k)/A
                    v_s[j, N-i] = (V_comm_2*(-k) + W_comm_2*(-l))*(-l)/A
                    u_s[N-i, N-j] = (V_comm_3*(-l) + W_comm_3*k)*(-l)/A
                    v_s[N-i, N-j] = (V_comm_3*(-l) + W_comm_3*k)*k/A
                    u_s[N-j, i] = (V_comm_4*k + W_comm_4*l)*k/A
                    v_s[N-j, i] = (V_comm_4*k + W_comm_4*l)*l/A
        
        # correct for singular values
        for i in range(dpd):
            for j in range(dpd):
                if isnan(u_s[i, j]) is True or abs(u_s[i, j]) > 1e15 or  abs(u_s[i, j]) < 1e-10:
                    u_s[i, j] = 0
                if isnan(v_s[i, j]) is True or abs(v_s[i, j]) > 1e15 or abs(v_s[i, j]) < 1e-10:
                    v_s[i, j] = 0
    
    # Create axes at clicked position from supplied position and given axis sizes
    deriv_inset_ax = main_axis.inset_axes([(x_pix-116)/500 - (0.931*d_length/(2*L)), (y_pix-59)/500 - (0.931*d_length/(2*L)), 0.931*d_length/L, 0.931*d_length/L])

    if tensor.get() == 0:
        arrows = False
        stacks = True
    if tensor.get() == 1:
        arrows = True
        stacks = False
    if tensor.get() == 2:
        arrows = True
        stacks = True            
    
    stack_plot(dxg, dyg, deriv_inset_ax, u_s, v_s, 5, d_range, dpd, 0.1, arrows, stacks, orientation, d_scale, w_head, h_head, 1) 
    
    # Don't display the x and y axis values
    # if click_opt_int > 2:  
    #     deriv_inset_ax.set_xticks([])
    #     deriv_inset_ax.set_yticks([])
    
    # Redraw the figure canvas, showing the inset axis
    fig.canvas.draw()
    deriv_inset_ax.clear()
    deriv_inset_ax.remove()


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
    # loop over differentiating each, when differentiating w.r.t its coord, set to 0
    for coord_index in range(len(coords)):
        # loop over differentiating each component:
        for comp_index in range(len(expressions)):
            J[comp_index, coord_index] = str(diff(expressions[comp_index], coords[coord_index]))
    return J


# define a function that will update the variable that defines click action
# and deals with setting up responses to click option changes
def click_option_handler(click_option):
    global click_opt_int, toolbar, LI_coord, LI_total_label, LI_total, LI_instruction_label, LI_restart_btn, LI_shape_select, x_m, y_m
    global x_pix, y_pix, x_m, y_m
    # get the selected option
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
        # unless the canvas has not been clicked yet
        try:
            deriv_calc(x_m, y_m)
        except (TypeError):
            # make up coordinates
            x_m = 0
            y_m = 0
            x_pix = 427
            y_pix = 308
            deriv_calc(x_m, y_m)


# Additional formatting function used in divergence plots
def format_eq_div(string):
    string = string.replace('xg', 'x_m')
    string = string.replace('yg', 'y_m')
    return string


# Function for updating the inset plots when options are changed by the user.
def update_deriv(self):
    if 0 < click_opt_int < 5:
        deriv_calc(x_m, y_m)
    else:
        pass


# define a function to put in an inset axis at user specified position
def set_inset_target():
    global x_m, y_m, x_pix, y_pix
    # get inputs from the entry boxes
    x_m = eval(x_m_entry.get())
    y_m = eval(y_m_entry.get())
    # Could not figure out how to reliably get pixels from coordinates on plot
    # got approximate relations based on experimenting
    x_pix = (x_m/L * 229) + 427
    y_pix = (y_m/L * 229) + 308
    # call deriv calc
    if 0 < click_opt_int < 5:
        deriv_calc(x_m, y_m)
    else:
        pass


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
    global magnitude_check, L, pt_den
    # get L and pt_den
    L = float(L_entry.get())
    pt_den = int(pt_den_entry.get())
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
            # get the strings of fields
            # and format them to be evaluated on the plotted grids
            string_check_x = string_x.replace('x', 'x_sing')
            string_check_x = string_check_x.replace('y', 'y_vals_singular')
            string_check_y = string_y.replace('x', 'x_sing')
            string_check_y = string_check_y.replace('y', 'y_vals_singular')
            # also replace for ln exp and **
            string_check_x = string_check_x.replace('^', '**')
            string_check_x = string_check_x.replace('ln', 'log')
            string_check_x = string_check_x.replace('e**', 'exp')
            string_check_y = string_check_y.replace('^', '**')
            string_check_y = string_check_y.replace('ln', 'log')
            string_check_y = string_check_y.replace('e**', 'exp')
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
            # also replace for ln exp and **
            string_check_x = string_check_x.replace('^', '**')
            string_check_x = string_check_x.replace('ln', 'log')
            string_check_x = string_check_x.replace('e**', 'exp')
            string_check_y = string_check_y.replace('^', '**')
            string_check_y = string_check_y.replace('ln', 'log')
            string_check_y = string_check_y.replace('e**', 'exp')
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
            # also replace for ln exp and **
            string_check_x = string_check_x.replace('^', '**')
            string_check_x = string_check_x.replace('ln', 'log')
            string_check_x = string_check_x.replace('e**', 'exp')
            string_check_y = string_check_y.replace('^', '**')
            string_check_y = string_check_y.replace('ln', 'log')
            string_check_y = string_check_y.replace('e**', 'exp')
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


'''  2-forms on R2 function from the other code '''


# define a function that will plot stack components, coloured
# as per the orientation of the 2-form at that grid point
def form_2_components_plot(grid_x, grid_y, form_2_loc, angle, s_max, L, fract, colour_str, s_min=2, ax=main_axis, axis_check=1, form_2_glob=None):
    global ratio2
    # get axis lengths:
    x_len = len(grid_x[:, 0])
    y_len = len(grid_y[:, 0])
    
    if axis_check == 1:
        # Scaling of axes and setting equal proportions circles look like circles
        ax.set_aspect('equal')
        ax_L = L + L/delta_factor
        ax.set_xlim(-ax_L, ax_L)
        ax.set_ylim(-ax_L, ax_L)
    else:
        ax.set_xlim(-L + x_m - L/5, L + x_m + L/5)
        ax.set_ylim(-L + y_m - L/5, L + y_m + L/5)
    
    # get the signs of the input 2-form
    form_2_sgn = np.sign(form_2_loc)
    
    # define an empty array of magnitudes, to then fill with integer rel. mags
    R_int = np.zeros(shape=((x_len), (y_len)))
    
    # #########################################################################
    # get variables needed for the initial, simplified stack plot
    # #########################################################################
    
    # find direction of each arrow
    theta = angle*np.ones(np.shape(form_2_loc))
    
    # deal with sinularities in mag
    for i in range(x_len):
        for j in range(y_len):
            # set to zero points that are not defined or inf
            if isnan(form_2_loc[i, j]) is True or abs(form_2_loc[i, j]) == np.inf  or abs(form_2_loc[i, j]) > 1e15:
                # colour this region as a red dot, not square to
                # not confuse with nigh mag 2-forms in stacks. or worse, in
                # blocks
                circ = patch.Circle((grid_x[i, j], grid_y[i, j]), L*fract/3, color='red')
                ax.add_patch(circ)
                form_2_loc[i, j] = 0
            # ALso, since we got this lop anyway
            # correct for singularities in planar form 2:
            # set to zero points that are not defined or inf
            if isnan(form_2_loc[i, j]) is True:
                form_2_sgn[i, j] = 0
    
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
    max_size = abs(np.max(form_2_loc))   # careful with singularities, else ---> nan
    
    # find the relative magnitude of vectors to maximum, as an array
    R = abs(form_2_loc)/max_size
    
    # globally scale if plotting an iset axis
    if axis_check == 0:
        # not to change the global form 2 here, change it to local
        form_2_global = form_2_glob*1
        # clear the global form_2 from singularities first:
        for i in range(len(form_2_global[:, 0])):
            for j in range(len(form_2_global[0, :])):
                # set to zero points that are not defined or inf
                if isnan(form_2_global[i, j]) is True or abs(form_2_global[i, j]) == np.inf  or abs(form_2_global[i, j]) > 1e15:
                    form_2_global[i, j] = 0
        # from globals, get the maximum over the whole 2 form
        max_form_global2 = np.max(form_2_global)
        # find its ratio to the one plotted here:
        max_form_local2 = np.max(form_2_loc)  # use mag as always u or v is zero
        # find their ratio:
        ratio2 = max_form_local2/max_form_global2
        # multiply the relative scaling array approperiately
        R *= ratio2
    
    # define tigonometirc shifts
    I_sin = np.sin(theta)
    I_cos = np.cos(theta)
    
    # define the points that set out a line of the stack sheet (middle line)
    A_x = grid_x + (sheet_L/2)*np.sin(theta)
    A_y = grid_y - (sheet_L/2)*np.cos(theta)
    B_x = grid_x - (sheet_L/2)*np.sin(theta)
    B_y = grid_y + (sheet_L/2)*np.cos(theta)
    
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
                color_index = 2
            
            # linear scaling
            for t in range(s_min, s_max+2):
                if (t-2)/s_max <= R[i, j] <= (t-1)/s_max:
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


# define a function that will find the 2-form from given expressions
# in a given number of dimensions and in terms of given coordinate symbols
def find_2_form(expressions, coords, pt_den, m=2):
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
    merge the results into a 2-form (for 2-form on R^2, the result is a single component (dx^xy))
    do so by adding opposite elements along the diagonal ( / ) components of ext_ds
    this  includes taking elemets with switched i and j
    '''
    
    # set up a variable to count pairs (pairs because we are forming 2-forms):
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
    
    # set up a vector to store the 2-form numerically, from xg and yg
    # Note - need pt_den m times.
    if m == 2:
        form_2 = np.empty((1, pt_den, pt_den))
        form_2[0, :, :] = eval(result[0, 0])
    elif m == 3:
        form_2 = np.empty((3, pt_den, pt_den, pt_den))
        for d in range(3):
            form_2[d, :, :, :] = eval(result[d, 0])
    return form_2


'''

Define response functions to GUI interactions

'''

# define a function that will let the user zoom with a mouse, on R2
def form_2_zoom(x_m, y_m):
    global form_2_zoom_values
    # define axis to be used in the inset
    zoomR2_range = (1/3)*L/(zoom_slider_R2.get())
    zoomR2_length = zoomR2_length_select.get()*L
    zoomR2pd = zoomR2pd_select.get()
    fract_zoom = 2/(zoomR2pd - 1)
    R2x = np.linspace(-zoomR2_range + x_m, zoomR2_range + x_m, zoomR2pd)
    R2y = np.linspace(-zoomR2_range + y_m, zoomR2_range + y_m, zoomR2pd)
    R2xg, R2yg = np.meshgrid(R2x, R2y)
    
    # Create axes at clicked position from supplied position and axis sizes
    zoomR2_inset_ax = main_axis.inset_axes([(x_pix-116)/500 - (0.931*zoomR2_length/(2*L)), (y_pix-59)/500 - (0.931*zoomR2_length/(2*L)), 0.931*zoomR2_length/L, 0.931*zoomR2_length/L])
    
    # determine which one is being plotted
    # a 1-form or a 2-form
    # based on that, complete the plot in the inset axis
    if calculus_form_tracker == 1:
        # define the vector field in these new axis
        uzoomR2, vzoomR2 = eq_to_comps(string_x, string_y, R2xg, R2yg)
        # plot as stacks
        stack_plot(R2xg, R2yg, zoomR2_inset_ax, uzoomR2, vzoomR2, s_max, zoomR2_range, zoomR2pd, 0.1, False, True, 'mid', 1, w_head, h_head, 1)
    # complete the process for when 2 form is to be zoomed
    elif calculus_form_tracker == 2:
        # define the 2 form over this region, on these axis
        form_2_zoom_str = form_2_str.replace('x', 'R2xg')
        form_2_zoom_str = form_2_zoom_str.replace('y', 'R2yg')
        form_2_zoom_str = form_2_zoom_str.replace('^', '**')
        form_2_zoom_str = form_2_zoom_str.replace('ln', 'log')
        # evalue it numerically
        form_2_zoom_values = eval(form_2_zoom_str)
        # plot it
        form_2_components_plot(R2xg, R2yg, form_2_zoom_values/2, np.pi/2, s_max, zoomR2_range, fract_zoom, colour_str, ax=zoomR2_inset_ax, axis_check=0, form_2_glob=form_2)
        form_2_components_plot(R2xg, R2yg, form_2_zoom_values/2, 0, s_max, zoomR2_range, fract_zoom, colour_str, ax=zoomR2_inset_ax, axis_check=0, form_2_glob=form_2)
    # Don't display the x and y axis values
    zoomR2_inset_ax.set_xticks([])
    zoomR2_inset_ax.set_yticks([])
    
    # redraw
    fig.canvas.draw()
    zoomR2_inset_ax.clear()
    zoomR2_inset_ax.remove()


# define a funvtion to call to select R2 tools
def R2_tools_handler(R2_tools_opt_var):
    global R2_tools_opt_int, toolbar
    global label_AI_area_2, label_AI_result_2, restart_AI_btn
    global label_AI_area_1, label_AI_result_1
    global shadearea_toggle, label_shade_area
    global R2_flux_shape, R2_flux_shape_list
    global R2_flux_shape_instruction, R2_flux_shape_drop
    # get the variable as an integer, make it global not to have to repeat this
    R2_tools_opt_int = R2_tools_opt_var
    if R2_tools_opt_int == 0:
        # in tools option:
        fig.canvas.draw()
        # if the tools is selected again, add the zoom and pan buttons
        # get rid of the modified toolbar:
        toolbar.destroy()
        # put the default matplotlib toolbar, back on:
        toolbar = NavigationToolbar2Tk(canvas, plot_frame)
        toolbar.update()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        # tools, therefore disable the slider
        zoom_slider_R2.configure(state=tk.DISABLED)
        # return all options from AREa integral
        form_0_btn.configure(state=tk.NORMAL)
        form_1_btn.configure(state=tk.NORMAL)
        # get rid of restart integral button
        try:
            restart_AI_btn.destroy()
            label_AI_area_1.destroy()
            label_AI_area_2.destroy()
            label_AI_result_1.destroy()
            label_AI_result_2.destroy()
            label_shade_area.destroy()
            shadearea_toggle.destroy()
            R2_flux_shape.destroy()
            R2_flux_shape_list.destroy()
            R2_flux_shape_instruction.destroy()
            R2_flux_shape_drop.destroy()
        except (UnboundLocalError, NameError):  
            pass
    elif R2_tools_opt_int == 1:
        # in zoom option:
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
        # enable it again
        zoom_slider_R2.configure(state=tk.NORMAL)
        # return all options from AREa integral
        form_0_btn.configure(state=tk.NORMAL)
        form_1_btn.configure(state=tk.NORMAL)
        # get rid of restart integral button
        try:
            restart_AI_btn.destroy()
            label_AI_area_1.destroy()
            label_AI_area_2.destroy()
            label_AI_result_1.destroy()
            label_AI_result_2.destroy()
            label_shade_area.destroy()
            shadearea_toggle.destroy()
            R2_flux_shape.destroy()
            R2_flux_shape_list.destroy()
            R2_flux_shape_instruction.destroy()
            R2_flux_shape_drop.destroy()
        except (UnboundLocalError, NameError):  
            pass
    else:
        # Area Integrals
        # go back to home and disable tools:
        fig.canvas.draw()
        toolbar.home()
        state = fig.canvas.toolbar.mode
        if state == 'zoom rect':
            toolbar.zoom()
        if state == 'pan/zoom':
            toolbar.pan()
        toolbar.children['!button4'].pack_forget()
        toolbar.children['!button5'].pack_forget()
        # disable the slider too
        zoom_slider_R2.configure(state=tk.DISABLED)
        # block the other options for the time being:
        form_0_btn.configure(state=tk.DISABLED)
        form_1_btn.configure(state=tk.DISABLED)
        ''' NOT DONE YET '''  # !!!
        # set up a restart button
        restart_AI_btn = tk.Button(calculus_frame, text='Reset Integral', command=AI_restart)
        restart_AI_btn.grid(row=15, column=0)
        # set up labels for result
        label_AI_area_1 = tk.Label(calculus_frame, text='Area drawn:')
        label_AI_area_1.grid(row=13, column=0)
        label_AI_area_2 = tk.Label(calculus_frame, text='0')
        label_AI_area_2.grid(row=13, column=1)
        label_AI_result_1 = tk.Label(calculus_frame, text='Integration result:')
        label_AI_result_1.grid(row=14, column=0)
        label_AI_result_2 = tk.Label(calculus_frame, text='0')
        label_AI_result_2.grid(row=14, column=1)
        # Shade Area of Area integtral -  toggle button
        label_shade_area = tk.Label(calculus_frame, text='Shade Area:')
        label_shade_area.grid(row=15, column=1)
        shadearea_toggle = tk.Button(calculus_frame, image=toggle_image_on, bd=0, command=shadearea_response)
        shadearea_toggle.grid(row=15, column=2)
        # dropdown to pick different shapes
        R2_flux_shape = tk.StringVar()
        R2_flux_shape_list = ['Polygon', 'Circle']
        R2_flux_shape.set(R2_flux_shape_list[0])
        R2_flux_shape_instruction = tk.Label(calculus_frame, text='Shape:')
        R2_flux_shape_instruction.grid(row=16, column=0)
        R2_flux_shape_drop = tk.OptionMenu(calculus_frame, R2_flux_shape, *R2_flux_shape_list, command=R2_flux_shape_response)
        R2_flux_shape_drop.grid(row=16, column=1)
        # define a dropdown to select which 
        # make sure a 2-form is being plotted:
        form_2_response()
        


# upate the zooming tool
def update_2_form_zoom(self):
    if R2_tools_opt_int == 0:
        pass  # tools
    elif R2_tools_opt_int == 1:
        form_2_zoom(x_m, y_m)  # zooming
    else:
        pass  # area integrals


# gets 2-form from entry box and plots it as coloured blocks only
def form_2_response():
    # get globals
    global form_2_str, form_2_eq, form_2, comp_x, comp_y, u, v, fract
    global calculus_form_tracker, fract
    # set tracker
    calculus_form_tracker = 2
    update_variables()
    # restart the area calculation if chosen
    if R2_tools_opt_int == 2:
        AI_restart(1)
    else:
        pass
    # from these, establish the new fract, approperiate for 2-forms
    fract = 2/((pt_den-1))
    # get the input 2-form
    form_2_str = str(simplify(form_2_entry.get()))
    # format it to be python understood
    form_2_eq = format_eq(form_2_str)
    # check against constant and zero 2-forms being supplied
    form_2_eq = form_2_constant_correction(form_2_eq)
    # get the numerical evaluation of it
    form_2 = eval(form_2_eq)
    # clear the current plot
    main_axis.clear()
    # use plotting stacks to display these
    #ALWAYS HALF AND HALF SPLITTING NOW
    form_2_components_plot(xg, yg, form_2/2, np.pi/2, s_max, L, fract, colour_str)
    form_2_components_plot(xg, yg, form_2/2, 0, s_max, L, fract, colour_str)
    # display the new plot
    canvas.draw()
    # display a background green on the 2-form entry to show that
    # this entry is being displayed now.
    form_2_entry.configure(bg='#C0F6BB')
    # undo it for 1-forms
    x_comp_entry.configure(bg='#FFFFFF')
    y_comp_entry.configure(bg='#FFFFFF')
    # update the label to remove the zero form if int deriv of both was used
    form_0_entry.configure(bg='#FFFFFF')


# plots the vetor field with stacks only
def form_1_stacks_response():
    global u, v, fract, calculus_form_tracker
    # set the tracker
    calculus_form_tracker = 1
    # clear the current axis
    main_axis.clear()
    # get input values
    update_variables()
    # get the supplied 1-forms from entry boxes:
    x_comp_str = str(simplify(x_comp_entry.get()))
    y_comp_str = str(simplify(y_comp_entry.get()))
    # take all these values, and the input from field component bnoxes to set up the field:
    u, v = eq_to_comps(x_comp_str, y_comp_str, xg, yg)
    # redefine better stack size for 1 forms
    fract = 0.05
    # plot the new field
    stack_plot(xg, yg, main_axis, u, v, s_max, L, pt_den, fract)
    # put it onto the screen
    canvas.draw()
    # display a background green on 1-form components and get rid of the 2-form
    # colour, to show which one is being plotted
    x_comp_entry.configure(bg='#C0F6BB')
    y_comp_entry.configure(bg='#C0F6BB')
    # get rid of the 2-form colour:
    form_2_entry.configure(bg='#FFFFFF')
    # update the label to remove the zero form if int deriv of both was used
    form_0_entry.configure(bg='#FFFFFF')


# performs the interior derivative on supplied 2-form and plots it as stacks
# for 2-form only, not including the 1-form.
# if combined, use different fucntion.
def Int_deriv_2_form():
    global u, v, u_str, v_str, vector_ex_str, vector_ey_str, vector_ex_eq, vector_ey_eq, vector_ex, vector_ey
    global form_2, form_2_eq
    # get globals
    update_variables()
    # from these, establish the new fract, approperiate for 1-forms
    fract = 0.05
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
    # get the input 2-form
    form_2_str = str(simplify(form_2_entry.get()))
    # format it to be python understood
    form_2_eq = format_eq(form_2_str)
    # check against constant and zero 2-forms being supplied
    form_2_eq = form_2_constant_correction(form_2_eq)
    # get the numerical evaluation of it
    form_2 = eval(form_2_eq)
    # using interior product, get the u and v (dx and dy) components
    # of the resulting 1-form
    u = -form_2 * vector_ey
    v = form_2 * vector_ex
    # to be usable in ext_deriv, define strings of these variables
    # later to put them into their entry boxes.
    u_str = str(simplify('-(' + form_2_str + ')*(' + vector_ey_str + ')' ))
    v_str = str(simplify( '(' + form_2_str + ')*(' + vector_ex_str + ')' ))
    # use the stacks plotter to present this
    stack_plot(xg, yg, main_axis, u, v, s_max, L, pt_den, fract, s_min=1)


# define a function that will find the interior derivative of a given 1-form
def Int_deriv_1_form():
    global zero_form_str, zero_form, vector_ex_str, vector_ey_str, vector_ex_eq
    global vector_ey_eq, vector_ex, vector_ey, x_comp, y_comp, x_comp_eq, y_comp_eq
    # get globals
    update_variables()
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
    # get the input 1-forms
    x_comp_str = str(simplify(x_comp_entry.get()))
    y_comp_str = str(simplify(y_comp_entry.get()))
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
    # as a string
    zero_form_str = str(simplify('(' + x_comp_str + ')*(' + vector_ex_str + ')' + ' + (' + y_comp_str + ')*(' + vector_ey_str + ')'))
    #define a finer grid for these contours
    contour_x, contour_y = np.linspace(-L, L, pt_den*5), np.linspace(-L, L, pt_den*5)
    contour_x_grid, contour_y_grid = np.meshgrid(contour_x, contour_y)
    # format zero_form_str in those
    zero_form_str = zero_form_str.replace('x', 'contour_x_grid')
    zero_form_str = zero_form_str.replace('y', 'contour_y_grid')
    # get it numerically
    # and take care of evaluation when the result is zero
    if zero_form_str.find('contour_x_grid') & zero_form_str.find('contour_y_grid') == -1:
        zero_form = eval(zero_form_str)*np.ones(np.shape(contour_x_grid))
    else:
        zero_form = eval(zero_form_str)
    # plot the zero_form as contours with labeled levels
    CS = main_axis.contour(contour_x_grid, contour_y_grid, zero_form, 15)
    main_axis.clabel(CS, inline=True, fontsize=7)
    # unformat for the label
    zero_form_str = zero_form_str.replace('contour_x_grid', 'x')
    zero_form_str = zero_form_str.replace('contour_y_grid', 'y')


# define a function that will respond to plotting the 2-form only
def Int_deriv_22_form():
    global calculus_form_tracker
    # set the tracker
    calculus_form_tracker = 1
    # clear the axis
    main_axis.clear()
    # call the asked fucntion
    Int_deriv_2_form()
    # draw its result
    canvas.draw()
    # update the label to remove the zero form if int deriv of both was used
    form_0_entry.configure(bg='#FFFFFF')
    # change the entry box 1-form to the calculated ones
    x_comp_entry.delete(0, 'end')
    y_comp_entry.delete(0, 'end')
    x_comp_entry.insert(0, u_str)
    y_comp_entry.insert(0, v_str)
    # display that the 1-form is now being plotted, therefore get rid of
    # 2-form colour and show 1-form components in green:
    x_comp_entry.configure(bg='#C0F6BB')
    y_comp_entry.configure(bg='#C0F6BB')
    form_2_entry.configure(bg='#FFFFFF')
    # destroy the extra window
    int_vector_window.destroy()


# define a function that will find the interior derivative of both the 2-form
# and a 1-form, merged.
def Int_deriv_21_form():
    global calculus_form_tracker
    # set the tracker
    calculus_form_tracker = 1
    # clear the axis
    main_axis.clear()
    # first call the function to plot a 2-form and plot it
    Int_deriv_2_form()
    # then call the function that will do this for the 1-form
    # plot these together
    Int_deriv_1_form()
    # draw its result
    canvas.draw()
    # make a label for the found 0 form
    form_0_entry.delete(0, 'end')
    form_0_entry.insert(0, zero_form_str)
    form_0_entry.configure(bg='#C0F6BB')
    # change the entry box 1-form to the calculated ones
    x_comp_entry.delete(0, 'end')
    y_comp_entry.delete(0, 'end')
    x_comp_entry.insert(0, u_str)
    y_comp_entry.insert(0, v_str)
    # display that the 1-form is now being plotted, therefore get rid of
    # 2-form colour and show 1-form components in green:
    x_comp_entry.configure(bg='#C0F6BB')
    y_comp_entry.configure(bg='#C0F6BB')
    form_2_entry.configure(bg='#FFFFFF')
    # destroy the extra window
    int_vector_window.destroy()


# define a response function to plot int deriv of 1-form only
def Int_deriv_11_form():
    global calculus_form_tracker
    # set the tracker of forms
    calculus_form_tracker = 0
    # clear the axis
    main_axis.clear()
    # call the fucntion to complete the 1 form int deriv and plot
    Int_deriv_1_form()
    # draw its result
    canvas.draw()
    # make a label for the found 0 form
    form_0_entry.delete(0, 'end')
    form_0_entry.insert(0, zero_form_str)
    form_0_entry.configure(bg='#C0F6BB')
    # display that the 0-form is now being plotted, therefore get rid of
    # 2 and1-form colours and show 0-form components in green:
    x_comp_entry.configure(bg='#FFFFFF')
    y_comp_entry.configure(bg='#FFFFFF')
    form_2_entry.configure(bg='#FFFFFF')
    # destroy the extra window
    int_vector_window.destroy()


# Interior derivative response for 2-form only
# asks to give vectors w.r.t to which perform iota.
def Int_deriv_2_form_response(type_form):
    global int_vector_window, int_vect_ex_entry, int_vect_ey_entry
    # Just in case
    if type_form == 1:
        if calculus_form_tracker == 0:
            tk.messagebox.showinfo('', 'Interior derivative of a 0-form is None')
            main_axis.clear()
            canvas.draw()
        else:
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
            # for whatever is being plotted
        if calculus_form_tracker == 2:
            int_vector_load_btn = tk.Button(int_vector_window, text='PLOT', padx=20, pady=10, command=Int_deriv_22_form)
            int_vector_load_btn.grid(row=4, column=0, pady=10)
        if calculus_form_tracker == 1:
            int_vector_load_btn = tk.Button(int_vector_window, text='PLOT', padx=20, pady=10, command=Int_deriv_11_form)
            int_vector_load_btn.grid(row=4, column=0, pady=10)
        # for both 2-form and 1-form
    if type_form == 2:
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
        # for whatever is being plotted
        int_vector_load_btn = tk.Button(int_vector_window, text='PLOT', padx=20, pady=10, command=Int_deriv_21_form)
        int_vector_load_btn.grid(row=4, column=0, pady=10)


# define a function that will respond to the made choice reg. int deriv.
def int_deriv_choice(var):
    if var == 1:
        # deal with whatever is being plotted at the time
        Int_deriv_2_form_response(1)
    elif var == 2:
        # call the function, but make it deal with 2-form together with the
        # 1-form
        Int_deriv_2_form_response(2)
    # close the previous window
    int_option_window.destroy()


# define response to int deriv button
# shows window to select if should just use 2-form
# or a comnination of 2-form and 1-form
def Int_deriv_response():
    global int_option_window
    # show new window with label and two options.
    int_option_window = tk.Toplevel()
    int_option_window.title('input a vector for the interior derivative')
    int_option_window.geometry('460x200')
    # define Label
    tk.Label(int_option_window, text='Perform w.r.t 2-form only or combine given 2-form and 1-form ?').grid(row=0, column=0, columnspan=2)
    # define response buttons to the stated question
    int_plotted_form_btn = tk.Button(int_option_window, text='Plotted form', padx=50, pady=30, command=lambda: int_deriv_choice(1))
    int_form_2_and_1_btn = tk.Button(int_option_window, text='Both 2-form and 1-form', padx=2, pady=30, command=lambda: int_deriv_choice(2))
    int_plotted_form_btn.grid(row=1, column=0, pady=20)
    int_form_2_and_1_btn.grid(row=1, column=1, pady=20)


# perform ext deriv on the result of int_deriv and plots it as stacks
def Ext_deriv_response():
    global form_2, form_2_str
    global calculus_form_tracker
    # clear current axis
    main_axis.clear()
    if calculus_form_tracker == 1:
        # set the tracker
        calculus_form_tracker = 2
        # get globals
        update_variables()
        # from these, establish the new fract, approperiate for 2-forms
        fract = 2/((pt_den-1))
        # get the inpus from fields of x and u components
        x_comp_str = str(simplify(x_comp_entry.get()))
        y_comp_str = str(simplify(y_comp_entry.get()))
        # from found u and v in the interior derivative, set up sympy components
        sympy_expr_x = parse_expr(x_comp_str, evaluate=False)
        sympy_expr_y = parse_expr(y_comp_str, evaluate=False)
        # combine the 2 into a list:
        expressions = np.array([sympy_expr_x, sympy_expr_y])
        # set up an array of coordinates that need to be used (in standard order)
        coords = ['x', 'y']
        # from these, use the find_2_form function to get the 2-form
        form_2 = find_2_form(expressions, coords, pt_den, m)[0]
        # get the string of this new 2-form to use it in int deriv
        # also put it into the entry
        form_2_str = str(simplify(str(unformat(result[0][0]))))
        # format it
        form_2_eq = format_eq(form_2_str)
        # numerically evalue it, careful about constants etc.
        form_2_eq = form_2_constant_correction(form_2_eq)
        form_2 = eval(form_2_eq)
        # unformat it to display in the entry box, this way it does not
        # format twice if int deriv runs again
        form_2_str = unformat(form_2_str)
        form_2_entry.delete(0, 'end')
        form_2_entry.insert(0, form_2_str)
        # clear the current plot
        main_axis.clear()
        # use plotting stacks to display these
        form_2_components_plot(xg, yg, form_2/2, np.pi/2, s_max, L, fract, colour_str)
        form_2_components_plot(xg, yg, form_2/2, 0, s_max, L, fract, colour_str)
        # display the new plot
        canvas.draw()
        # display a background green on the 2-form entry to show that
        # this entry is being displayed now.
        form_2_entry.configure(bg='#C0F6BB')
        # undo it for 1-forms (even though they are hidden, they are not off)
        x_comp_entry.configure(bg='#FFFFFF')
        y_comp_entry.configure(bg='#FFFFFF')
        # update the label to remove the zero form if int deriv of both was used
        form_0_entry.configure(bg='#FFFFFF')
    elif calculus_form_tracker == 0:
        # change the tracker
        calculus_form_tracker = 1
        # get globals
        update_variables()
        # from these, establish the new fract, approperiate for a 1-form
        fract = 0.05
        # get the input from 0-form box
        form_0_str = str(simplify(form_0_entry.get()))
        # from this, need derivatives so set it as a SymPy object
        sympy_expr_form_0 = parse_expr(form_0_str, evaluate=False)
        # set up an array of coordinates that need to be used (in standard order)
        coords = ['x', 'y']
        # from these, find the derivatives
        form_0_x = str(diff(sympy_expr_form_0, coords[0]))
        form_0_y = str(diff(sympy_expr_form_0, coords[1]))
        # put these into their entry boxes
        x_comp_entry.delete(0, 'end')
        y_comp_entry.delete(0, 'end')
        x_comp_entry.insert(0, form_0_x)
        y_comp_entry.insert(0, form_0_y)
        # call a fucntion to plot these:
        form_1_stacks_response()
    else:
        tk.messagebox.showerror('', '3-form on R^2 is not defined')


# define a function that will wedge two 1-forms and plot them
def wedge_product_R2():
    global to_wedge_x_1_str, to_wedge_y_1_str, to_wedge_x_2_str, to_wedge_y_2_str
    global form_2_str, form_2_eq, form_2
    global calculus_form_tracker
    # set the tracker
    calculus_form_tracker = 2
    # get globals
    update_variables()
    # from these, establish the new fract, approperiate for 2-forms
    fract = 2/((pt_den-1))
    # first, get all entries out, save as string for these to display when
    # window is opened again
    to_wedge_x_1_str = str(simplify(str(to_wedge_x_1_entry.get())))
    to_wedge_y_1_str = str(simplify(str(to_wedge_y_1_entry.get())))
    to_wedge_x_2_str = str(simplify(str(to_wedge_x_2_entry.get())))
    to_wedge_y_2_str = str(simplify(str(to_wedge_y_2_entry.get())))
    # clear the axis:
    main_axis.clear()
    # first, find the result of the 2-form
    # this if, in terms of the above commented fields:
    # 2-form = f*m - g*h
    # get it mathematically, as a string
    form_2_str = str(simplify( '(' + to_wedge_x_1_str + ')*(' +  to_wedge_y_2_str + ')' + ' - (' + to_wedge_y_1_str + ')*(' +  to_wedge_x_2_str + ')' ))
    # put it into the entry box for 2-forms
    form_2_entry.delete(0, 'end')
    form_2_entry.insert(0, form_2_str)
    # plot these as stacks, with no arrowheads, on top of one another.
     # format it to be python understood
    form_2_eq = format_eq(form_2_str)
    # check against constant and zero 2-forms being supplied
    form_2_eq = form_2_constant_correction(form_2_eq)
    # get the numerical evaluation of it
    form_2 = eval(form_2_eq)
    # use plotting stacks to display these
    # ALWAYS HALF AND HALF SPLITTING NOW
    form_2_components_plot(xg, yg, form_2/2, np.pi/2, s_max, L, fract, colour_str)
    form_2_components_plot(xg, yg, form_2/2, 0, s_max, L, fract, colour_str)
    # display the new plot
    canvas.draw()
    # display a background green on the 2-form entry to show that
    # this entry is being displayed now.
    form_2_entry.configure(bg='#C0F6BB')
    # undo it for 1-forms
    x_comp_entry.configure(bg='#FFFFFF')
    y_comp_entry.configure(bg='#FFFFFF')
    # update the label to remove the zero form if int deriv of both was used
    form_0_entry.configure(bg='#FFFFFF')


# define a reponse function, opens new window where two 1-forms to be wedged can be entered
def wedge_2_response():
    global wedge_2_window, to_wedge_x_1_entry, to_wedge_y_1_entry, to_wedge_x_2_entry, to_wedge_y_2_entry
    # open a titled new window
    wedge_2_window = tk.Toplevel()
    wedge_2_window.title('input two 1-forms to wedge')
    # define all entry boxes
    tk.Label(wedge_2_window, text='first 1-form x component :').grid(row=0, column=0)
    to_wedge_x_1_entry = tk.Entry(wedge_2_window, width=30, borderwidth=1)
    to_wedge_x_1_entry.insert(0, to_wedge_x_1_str)
    to_wedge_x_1_entry.grid(row=1, column=0)
    tk.Label(wedge_2_window, text='first 1 from y component:').grid(row=2, column=0)
    to_wedge_y_1_entry = tk.Entry(wedge_2_window, width=30, borderwidth=1)
    to_wedge_y_1_entry.insert(0, to_wedge_y_1_str)
    to_wedge_y_1_entry.grid(row=3, column=0)
    tk.Label(wedge_2_window, text='second 1-form x component :').grid(row=4, column=0)
    to_wedge_x_2_entry = tk.Entry(wedge_2_window, width=30, borderwidth=1)
    to_wedge_x_2_entry.insert(0, to_wedge_x_2_str)
    to_wedge_x_2_entry.grid(row=5, column=0)
    tk.Label(wedge_2_window, text='second 1-form y component :').grid(row=6, column=0)
    to_wedge_y_2_entry = tk.Entry(wedge_2_window, width=30, borderwidth=1)
    to_wedge_y_2_entry.insert(0, to_wedge_y_2_str)
    to_wedge_y_2_entry.grid(row=7, column=0)
    # define a button that will plot these
    plot_wedge_btn = tk.Button(wedge_2_window, text='PLOT', padx=20, pady=10, command=wedge_product_R2)
    plot_wedge_btn.grid(row=8, column=0, pady=10)


# define a fucntion that will respond to finding the Hodge dual of the given froms
# 1-form or 2-form depending on chosen parameter
# f is zero form, A is dx component, B is dy component, D is (dx /\ dy) component. Return list of the components and print result
def Hodge2D(f='0',A='0',B='0',D='0'):
    # switch components as per rules of a Hidge on R2
    f_out = D
    A_out = '-(' + B + ')'
    B_out = A
    D_out = f
    # put them ito an array
    H = [f_out, A_out, B_out, D_out]
    # return to user
    return H


# define a fucntion that will respond to finding the Hodge dual of the given froms
# 1-form or 2-form depending on chosen parameter
def Hodge_1_form_response():
    global fract, calculus_form_tracker
    # update tracker for zoom to work
    calculus_form_tracker = 1
    # update fract to work best on 1 form (result)
    fract = 0.05
    # prepare variables
    main_axis.clear()
    update_variables()
    u_str = x_comp_entry.get()
    v_str = y_comp_entry.get()
    
    # use previously dfined function to take Hodge
    H = Hodge2D(A = u_str, B = v_str)
    
    # from that get components
    u, v = eq_to_comps(H[1], H[2], xg, yg)
    
    # simplify the result to input into result box
    H[1] = simplify((H[1]))
    H[2] = simplify((H[2]))
    
    # input into boxes and colour
    x_comp_entry.delete(0, 'end')
    y_comp_entry.delete(0, 'end')
    x_comp_entry.insert(0, H[1])
    y_comp_entry.insert(0, H[2])
    form_2_entry.configure(bg='#FFFFFF')
    x_comp_entry.configure(bg='#C0F6BB')
    y_comp_entry.configure(bg='#C0F6BB')                     
    # plot, always as a 1 form:
    arrows = False
    stacks = True
    stack_plot(xg, yg, main_axis, u, v, s_max, L, pt_den, fract, arrows, stacks, orientation, scale, w_head, h_head, 0) 
    # redraw
    canvas.draw()
    # update the label to remove the zero form if int deriv of both was used
    form_0_entry.configure(bg='#FFFFFF')


# Hodge the two form in the entry box, plotting the scalar function as a contour
def Hodge_2_form_response():
    # Note: can't approperiately set variable for zooming here
    # no zooming of contours is implemented, nor needed
    
    # prepare variables
    main_axis.clear()
    update_variables()
    
    # get the 0 form from the 2 form using prev. def. function
    form_2_to_hodge = form_2_entry.get()
    form_2_to_hodge = form_2_to_hodge.replace('ln', 'log')
    form_2_to_hodge = form_2_to_hodge.replace('^', '**')
    H = Hodge2D(D=form_2_to_hodge)
    zero_form_str = H[0]
    
    # set up contours, AS BEFORE, and plot:
    contour_x, contour_y = np.linspace(-L, L, pt_den*5), np.linspace(-L, L, pt_den*5)
    contour_x_grid, contour_y_grid = np.meshgrid(contour_x, contour_y)
    zero_form_str = zero_form_str.replace('x', 'contour_x_grid')
    zero_form_str = zero_form_str.replace('y', 'contour_y_grid')
    if zero_form_str.find('contour_x_grid') & zero_form_str.find('contour_y_grid') == -1:
        zero_form = eval(zero_form_str)*np.ones(np.shape(contour_x_grid))
    else:
        zero_form = eval(zero_form_str)
    CS = main_axis.contour(contour_x_grid, contour_y_grid, zero_form, 15)
    main_axis.clabel(CS, inline=True, fontsize=7)
    # colour approperiate boxes
    form_2_entry.configure(bg='#FFFFFF')
    x_comp_entry.configure(bg='#FFFFFF')
    y_comp_entry.configure(bg='#FFFFFF')
    # unformat the string of the zero_form to print to user:
    zero_form_str = zero_form_str.replace('contour_x_grid', 'x')
    zero_form_str = zero_form_str.replace('contour_y_grid', 'y')
    # input the 0 form into assigned (hidden) label
    form_0_entry.delete(0, 'end')
    form_0_entry.insert(0, zero_form_str)
    form_0_entry.configure(bg='#C0F6BB')
    canvas.draw()


# define a function to put in an inset axis at user specified position
# for the calculus tab
def set_inset_target_calc():
    global x_m, y_m, x_pix, y_pix
    # get inputs from the entry boxes
    x_m = eval(x_m_entry_calc.get())
    y_m = eval(y_m_entry_calc.get())
    # Could not figure out how to reliably get pixels from coordinates on plot
    # got approximate relations based on experimenting
    x_pix = (x_m/L * 229) + 427
    y_pix = (y_m/L * 229) + 308
    # call zooming tool
    if R2_tools_opt_int == 0:
        pass
    else:
        form_2_zoom(x_m, y_m)


# define a fucntion to plot a zero form when button is pressed.
def form_0_response():
    global zero_form_str, zero_form, calculus_form_tracker
    # clear the axis
    main_axis.clear()
    # set up the tracker
    calculus_form_tracker = 0
    # get the supplied form
    zero_form_str = str(simplify(form_0_entry.get()))
    # set up grids for contours
    contour_x, contour_y = np.linspace(-L, L, pt_den*5), np.linspace(-L, L, pt_den*5)
    contour_x_grid, contour_y_grid = np.meshgrid(contour_x, contour_y)
    # format the given ftring
    zero_form_str = zero_form_str.replace('x', 'contour_x_grid')
    zero_form_str = zero_form_str.replace('y', 'contour_y_grid')
    # evaluate bearing in mind zeros
    if zero_form_str.find('contour_x_grid') & zero_form_str.find('contour_y_grid') == -1:
        zero_form = eval(zero_form_str)*np.ones(np.shape(contour_x_grid))
    else:
        zero_form = eval(zero_form_str)
    # set up the contour plot
    CS = main_axis.contour(contour_x_grid, contour_y_grid, zero_form, 15)
    main_axis.clabel(CS, inline=True, fontsize=7)
    # colour approperiate boxes
    form_2_entry.configure(bg='#FFFFFF')
    x_comp_entry.configure(bg='#FFFFFF')
    y_comp_entry.configure(bg='#FFFFFF')
    # colour the 0-form box
    form_0_entry.configure(bg='#C0F6BB')
    # draw the plot on
    canvas.draw()


# define a fucntion that will respond to dropdown selection of 2-forms
def selection_form_2_response(event):
    # clear the 2-form box
    form_2_entry.delete(0, 'end')
    # get the index at which this entry is
    selected_index = list_form_2_names.index(str(select_form_2_drop.get()))
    # using that index, get the equation its list
    # and insert into 2-form box
    selected_form_2 = list_form_2_equations[selected_index]
    form_2_entry.insert(0, selected_form_2)
    # colour code to be able to distinguish what is being plotted
    x_comp_entry.configure(bg='#C0F6BB')
    y_comp_entry.configure(bg='#C0F6BB')
    # now call the plot function to finalise all these onto the plot
    form_2_response()


# define a function that will integrate 2-forms
# over given regions
def integration_form_2(AI_verts):
    global AI_area, AI_result, shape_complete_tracker, form_2_inside
    # Calculate the area of that shape:
    AI_area = calc_area(AI_verts)
    # check against negative areas:
    if AI_area < 0:
        AI_area *= -1
    # define a grid of points that covers the area the user has input:
    # change the list into an array
    coord_verts = np.array(AI_verts)
    # find the maxima and minima of shape in x and y
    AI_verts_x = coord_verts[:, 0]
    AI_verts_y = coord_verts[:, 1]
    shape_x_min = np.min(AI_verts_x)
    shape_x_max = np.max(AI_verts_x)
    shape_y_min = np.min(AI_verts_y)
    shape_y_max = np.max(AI_verts_y)
    # set up accuracy
    if shadearea.get() == 0:
        points_N = 120
    else:
        points_N = 100
    # based on these, define a superfine grid to integrate over
    x_shape_points = np.linspace(shape_x_min, shape_x_max, points_N)
    y_shape_points = np.linspace(shape_y_min, shape_y_max, points_N)
    # mesh these
    x_shape_g, y_shape_g = np.meshgrid(x_shape_points, y_shape_points)
    # CRUCIAL STEP:
    # need to test if points are inside the shape, if not, set them
    # to something to be excluded later
    # set up a list to store points coordinates, that are inside
    inside_list = []
    polygon_path = mplPath.Path(coord_verts)
    if shadearea.get() == 1:
        for i in range(points_N):
            for j in range(points_N):
                if polygon_path.contains_point((x_shape_g[i, j], y_shape_g[i, j])) is True:
                    # sppend to inside list
                    inside_list.append([x_shape_g[i, j], y_shape_g[i, j]])
                    # Was made for testing but I love it, colour the inside of the drawn shape
                    circle =  patch.Circle((x_shape_g[i, j], y_shape_g[i, j]), L*fract/20, color='#DADADA')
                    main_axis.add_patch(circle)
        # finisise the shading if it was chosen
        canvas.draw()
    else:
        for i in range(points_N):
            for j in range(points_N):
                if polygon_path.contains_point((x_shape_g[i, j], y_shape_g[i, j])) is True:
                    # sppend to inside list
                    inside_list.append([x_shape_g[i, j], y_shape_g[i, j]])
    # split the inside list to be now an array of x and y components
    inside_arr = np.array(inside_list)
    points_x = inside_arr[:, 0]
    points_y = inside_arr[:, 1]
    # grid these:
    points_xg, points_yg = np.meshgrid(points_x, points_y)
    # then need to evalue the 2_form at all of these
    form_2_str = str(simplify(form_2_entry.get()))
    form_2_inside_eq = form_2_str.replace('^', '**')
    form_2_inside_eq = form_2_inside_eq.replace('ln', 'log')
    form_2_inside_eq = form_2_inside_eq.replace('x', 'points_xg')
    form_2_inside_eq = form_2_inside_eq.replace('y', 'points_yg')
    # take care of constnts and zeros:
    if form_2_inside_eq.find('points_xg') & form_2_inside_eq.find('points_yg') == -1:
        form_2_inside_eq = '(' + str(form_2_inside_eq) + ')* np.ones((len(points_x), len(points_y)))'
    else:
        pass
    form_2_inside = eval(form_2_inside_eq)
    # then sum over all values in that array
    sum_inside = np.sum(form_2_inside)
    # get the elemental area
    dxdy_area = (x_shape_points[1] - x_shape_points[0])*(y_shape_points[1] - y_shape_points[0])
    # then evaluate the total of the integral
    AI_result = sum_inside * dxdy_area / len(inside_list)
    # show these results
    label_AI_result_2.configure(text=str(round(AI_result, 1)))
    label_AI_area_2.configure(text=str(round(AI_area, 1)))
    # set up a variable to establish that a shape has been ompleted now
    shape_complete_tracker = 1


# define a function to integrate 2-form flux over a circle
def integration_from_2_circle(x_m, y_m, radius_for_flux, AI_verts):
    global AI_area, AI_result, shape_complete_tracker, form_2_inside
    # globals for tests only
    global inside_arr, circle_for_flux, x_shape_points, y_shape_points
    global circle_as_polygon, points
    # define a circle patch
    circle_for_flux = patch.CirclePolygon((x_m, y_m), radius_for_flux, fill=False, color='red')
    # get its transformed path (found out this is needed from the 
    # circle patch object describtion and found a way of doing it
    # on stackoverflow):
    verts = circle_for_flux.get_path().vertices
    trans = circle_for_flux.get_patch_transform()
    points = trans.transform(verts)
    # set up an array these are to be sotred in
    coord_verts = np.zeros(shape=(len(points), 2))
    # get the points into an array
    for i in range(len(points)):
        row_to_append = np.array([[points[i, 0], points[i, 1]]])
        coord_verts[i, :] = row_to_append
    # get a list of these, put it into AI_verts
    AI_verts = coord_verts.tolist()
    # Calculate the area of that shape:
    AI_area = calc_area(AI_verts)
    # check against negative areas:
    if AI_area < 0:
        AI_area *= -1
    # add the array to a patch path
    polygon_path = mplPath.Path(coord_verts)
    # set up the grid and find points in the circle:
    shape_x_min = x_m - radius_for_flux
    shape_x_max = x_m + radius_for_flux
    shape_y_min = y_m - radius_for_flux
    shape_y_max = y_m + radius_for_flux
    # set up accuracy
    if shadearea.get() == 0:
        points_N = 100
    else:
        points_N = 80
    # based on these, define a superfine grid to integrate over
    x_shape_points = np.linspace(shape_x_min, shape_x_max, points_N)
    y_shape_points = np.linspace(shape_y_min, shape_y_max, points_N)
    # mesh these
    x_shape_g, y_shape_g = np.meshgrid(x_shape_points, y_shape_points)
    # CRUCIAL STEP:
    # need to test if points are inside the shape, if not, set them
    # to something to be excluded later
    # set up a list to store points coordinates, that are inside
    inside_list = []
    if shadearea.get() == 1:
        for i in range(points_N):
            for j in range(points_N):
                if polygon_path.contains_point((x_shape_g[i, j], y_shape_g[i, j])) is True:
                    # append to inside list
                    inside_list.append([x_shape_g[i, j], y_shape_g[i, j]])
                    # Was made for testing but I love it, colour the inside of the drawn shape
                    grey_patching =  patch.Circle((x_shape_g[i, j], y_shape_g[i, j]), L*fract/20, color='#DADADA')
                    main_axis.add_patch(grey_patching)
    else:
        for i in range(points_N):
            for j in range(points_N):
                if polygon_path.contains_point((x_shape_g[i, j], y_shape_g[i, j])) is True:
                    # sppend to inside list
                    inside_list.append([x_shape_g[i, j], y_shape_g[i, j]])
    # put the patch for the perimeter on the screen
    # after for this to not get covered
    main_axis.add_patch(circle_for_flux)
    # finilise drawing
    canvas.draw()
    # split the inside list to be now an array of x and y components
    inside_arr = np.array(inside_list)
    points_x = inside_arr[:, 0]
    points_y = inside_arr[:, 1]
    # grid these:
    points_xg, points_yg = np.meshgrid(points_x, points_y)
    # then need to evalue the 2_form at all of these
    form_2_str = str(simplify(form_2_entry.get()))
    form_2_inside_eq = form_2_str.replace('^', '**')
    form_2_inside_eq = form_2_inside_eq.replace('ln', 'log')
    form_2_inside_eq = form_2_inside_eq.replace('x', 'points_xg')
    form_2_inside_eq = form_2_inside_eq.replace('y', 'points_yg')
    # take care of constnts and zeros:
    if form_2_inside_eq.find('points_xg') & form_2_inside_eq.find('points_yg') == -1:
        form_2_inside_eq = '(' + str(form_2_inside_eq) + ')* np.ones((len(points_x), len(points_y)))'
    else:
        pass
    form_2_inside = eval(form_2_inside_eq)
    # then sum over all values in that array
    sum_inside = np.sum(form_2_inside)
    # get the elemental area
    dxdy_area = (x_shape_points[1] - x_shape_points[0])*(y_shape_points[1] - y_shape_points[0])
    # then evaluate the total of the integral
    AI_result = sum_inside * dxdy_area / len(inside_list)
    # show these results
    label_AI_result_2.configure(text=str(round(AI_result, 1)))
    label_AI_area_2.configure(text=str(round(AI_area, 1)))
    # set up a variable to establish that a shape has been ompleted now
    shape_complete_tracker = 1


# define a fucntion that will track user clicks in a list, and check if they
# form a closed area, if so, will call the above inetgration fucntion
def area_finder_form_2_int(x_click, y_click):
    global AI_verts
    # append the clikck coordinates  to the global AI_coords:
    AI_coord.append([x_m, y_m])
    # gte number of clicks
    N_clicks = len(AI_coord)
    # set a tolerance to when clicks are considered joined:
    tolerance = 0.02*L  # same as for LI
    # get an array from the AI list
    coord_array = np.array(AI_coord)
    if N_clicks == 1:
        # draw a small dot
        circle = patch.Circle(AI_coord[0], L*fract/18, color='red')
        main_axis.add_patch(circle)
        canvas.draw()
    else:
        # draw the line:
        # get coordinates from mouse clicks
        a = AI_coord[N_clicks - 2]
        b = AI_coord[N_clicks - 1]
        # Plot line between points a and b
        main_axis.add_line(Line2D((a[0], b[0]), (a[1], b[1]), linewidth=2, color='red'))
        # draw it on
        canvas.draw()
        # test if the user closed the area by click onto (very near)
        # a previously clicked point, if so, autojoin:
        # get disance between current and all prev. points
        coord_diff = coord_array - coord_array[N_clicks - 1, :]
        # cycle over them to check if click is near any
        for i in range(N_clicks):
            # check if distance is less than specified:
            if sqrt((coord_diff[i, 0])**2 + (coord_diff[i, 1])**2) < tolerance:
                # set the current cordiante to the previous one that was joined
                # onto:
                AI_coord[N_clicks-1] = AI_coord[i]
                # get verticies
                AI_verts = AI_coord[i:]
                break
        # now, if enough of these verticies were extracted, a shape was drawn
        # check for that:
        if len(AI_verts) > 3:  # as it included the return point
            # a closed shape was drawn, pass it onto the integral calculating
            # function
            integration_form_2(AI_verts)


# define a function that will restart the 2-form integration
def AI_restart(test=0):
    global AI_coord, AI_result, AI_area, AI_verts, shape_complete_tracker
    # first, initialise variables again
    AI_coord = []
    AI_result = 0
    AI_area = 0
    AI_verts = []
    # update labels
    label_AI_result_2.configure(text='0')
    label_AI_area_2.configure(text='0')
    # NOT IDEAL BUT HOPEFULLY TEMPORARY
    # call plot-respose to redraw (so lines are deleted)
    # ideally would want a way of deleting lines without restarting the plot
    if test == 0:
        form_2_response()
    else:
        pass
    # set up a tracker to account for the shape being undone
    shape_complete_tracker = 0


# define a function that will respond to toggle switch, to shade the area
def shadearea_response():
    global shadearea
    if shadearea.get() == 0:
        # the button is off, and has been clicked therefore change the
        # variable to an and the image to on
        shadearea.set(1)
        shadearea_toggle.configure(image=toggle_image_on)
    else:
        # the button is on and has been clicked
        # set it to off and change image
        shadearea.set(0)
        shadearea_toggle.configure(image=toggle_image_off)


def R2_flux_shape_response(selected_shape):
    global R2_flux_shape, Radius_R2_circ_entry, Radius_R2_label
    # get the chosen shape
    R2_flux_shape.set(selected_shape)
    # depending on that selection, prepare rest of frame:
    if selected_shape == 'Polygon':
        # restart the plot and the integral values
        AI_restart()
        # get rid of circle options
        try:
            Radius_R2_label.destroy()
            Radius_R2_circ_entry.destroy()
        except (UnboundLocalError, NameError):  
            pass
    if selected_shape == 'Circle':
        # restart the plot and the integral values
        AI_restart()
        # if circle is selected, display an entry box for the radius of it
        Radius_R2_label = tk.Label(calculus_frame, text='Circle Radius:')
        Radius_R2_label.grid(row=17, column=0)
        Radius_R2_circ_entry = tk.Entry(calculus_frame, width=10)
        Radius_R2_circ_entry.grid(row=17, column=1)
        Radius_R2_circ_entry.insert(0, '1')
    # update canvas
    canvas.draw()


''' DEFINE FUNCTIONS USED IN R3 CODE '''


# find the max value of each, excluding singularities
def find_global_max(string):
    global values
    # format
    equation = format_eq(string)
    # evaluate
    values = eval(equation)
    # get rid of infs and nans:
    try:
        for i in range(len(values[:, 0, 0])):
            for j in range(len(values[0, :, 0])):
                for k in range(len(values[0, 0, :])):
                    if isnan(values[i, j, k]) is True or abs(values[i, j, k]) == np.inf or values[i, j, k] > 1e15:
                        values[i, j, k] = 0
    except TypeError:
        pass
    # now find maximum
    maximum_in_it = np.max(values)
    return maximum_in_it


# define a function that will plot stack components, coloured
# as per the orientation of the 2-form at that grid point
def form_2_components_plot_3(grid_x, grid_y, h_index, axis_view, u, v, s_max, L, pt_den, fract, colour_str, s_min=0):
    global s_L
    
    # depending on axis_view and h_index, get the planar u and v from the given
    # 3D ones and the grids sorted
    # Also
    # get the 2-form signs and
    # change the 2_form signs to just be in terms of the selected plane
    if axis_view == 'z':
        grid_x = grid_x[:, :, h_index]
        grid_y = grid_y[:, :, h_index]
        u = u[:, :, h_index]
        v = v[:, :, h_index]
        form_2_sgn = np.sign(form_2[0])
        form_2_sgn_planar = form_2_sgn[:, :, h_index]
    elif axis_view == 'y':
        grid_x = grid_x[h_index, :, :]
        grid_y = grid_y[h_index, :, :]
        u = u[h_index, :, :]
        v = v[h_index, :, :]
        form_2_sgn = np.sign(form_2[1])
        form_2_sgn_planar = form_2_sgn[h_index, :, :]
    elif axis_view == 'x':
        grid_x = grid_x[:, h_index, :]
        grid_y = grid_y[:, h_index, :]
        u = u[:, h_index, :]
        v = v[:, h_index, :]
        form_2_sgn = np.sign(form_2[2])
        form_2_sgn_planar = form_2_sgn[:, h_index, :]
    else:
        print('Error can\'t find this axis')
    
    # get axis lengths:
    x_len = len(grid_x[:, 0])  # no need to change with axis_view
    y_len = len(grid_y[0, :])  # if grids are all same size
    
    # Scaling of axes and setting equal proportions circles look like circles
    main_axis.set_aspect('equal')
    ax_L = L + L/delta_factor
    main_axis.set_xlim(-ax_L, ax_L)
    main_axis.set_ylim(-ax_L, ax_L)
    
    # define an empty array of magnitudes, to then fill with integer rel. mags
    R_int = np.zeros(shape=((x_len), (y_len)))
    
    # #########################################################################
    # get variables needed for the initial, simplified stack plot
    # #########################################################################
    
    # find the arrow length corresponding to each point and store in mag array
    mag = np.sqrt(u**2 + v**2)
    
    # find direction of each arrow
    theta = np.arctan2(v, u)   # theta defined from positive x axis ccw
    
    # deal with sinularities in mag
    for i in range(x_len):
        for j in range(y_len):
            # set to zero points that are not defined or inf
            if isnan(mag[i, j]) is True or abs(mag[i, j]) == np.inf  or abs(mag[i, j]) > 1e15:
                # colour this region as a red dot, not square to
                # not confuse with nigh mag 2-forms in stacks. or worse, in
                # blocks
                circ = patch.Circle((grid_x[i, j], grid_y[i, j]), L*fract/3, color='red')
                main_axis.add_patch(circ)
                mag[i, j] = 0
            # ALso, since we got this lop anyway
            # correct for singularities in planar form 2:
            # set to zero points that are not defined or inf
            if isnan(mag[i, j]) is True:
                form_2_sgn_planar[i, j] = 0
            if mag[i, j] == np.inf  or mag[i, j] > 1e15:
                form_2_sgn_planar[i, j] = 1
            if mag[i, j] == -np.inf  or mag[i, j] < -1e15:
                form_2_sgn_planar[i, j] = -1
    
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
    max_size = np.max(mag)
    
    # find the relative magnitude of vectors to maximum, as an array
    R = mag/max_size
    
    # to get tubing right, scale with respect to
    # approperiate global max value
    # to scale, limit s_max based on its max mag.
    try:
        if axis_view == 'z':
            mag_loc = max_size/max_global_dxdy
            R *= mag_loc
        if axis_view == 'y':
            mag_loc = max_size/max_global_dxdz
            R *= mag_loc
        if axis_view == 'x':
            mag_loc = max_size/max_global_dydz
            R *= mag_loc
    except ValueError:
        print('NaN encountered in divide')
        pass
    # define tigonometirc shifts
    I_sin = np.sin(theta)
    I_cos = np.cos(theta)
    
    # define the points that set out a line of the stack sheet (middle line)
    A_x = grid_x + (sheet_L/2)*np.sin(theta)
    A_y = grid_y - (sheet_L/2)*np.cos(theta)
    B_x = grid_x - (sheet_L/2)*np.sin(theta)
    B_y = grid_y + (sheet_L/2)*np.cos(theta)
    
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
            
            for t in range(s_min, s_max+2):
                if (t-2)/s_max <= R[i, j] <= (t-1)/s_max:
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
                    main_axis.add_line(Line2D((Ax1, Bx1), (Ay1, By1), linewidth=0.5, color=colour_str[color_index]))
                    main_axis.add_line(Line2D((Ax2, Bx2), (Ay2, By2), linewidth=0.7, color=colour_str[color_index]))
                    
                    # update parameter to reapet and draw all needed arrows
                    s += 1
            # deal with the odd number of stacks:
            elif parity(n) is False:
                # Add the centre line for odd numbers of stacks
                main_axis.add_line(Line2D((A_x[i, j], B_x[i, j]), (A_y[i, j], B_y[i, j]), linewidth=0.7, color=colour_str[color_index]))
                
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
                 
                    main_axis.add_line(Line2D((Ax1,Bx1),(Ay1,By1), linewidth=0.7, color=colour_str[color_index]))
                    main_axis.add_line(Line2D((Ax2,Bx2),(Ay2,By2), linewidth=0.7, color=colour_str[color_index]))
                    
                    # change the parameter to loop over all changes in displacement for current magnitude
                    s += 1
    plt.close()


def form_1_to_2onR3():
    global string_x, string_y, string_z, form_2, F_xy_x, F_xy_y, F_xz_x, F_xz_z, F_yz_y, F_yz_z
    global max_global_dxdy, max_global_dxdz, max_global_dydz, fract
    global hvalue_string, h_index
    # get the updated variables
    update_variables(1)
    # update the h_index label based on the newly input arrays
    try:
        hvalue_string = str(z[h_index])
        axis_height_txt.configure(text=str(round(eval(hvalue_string), 2)))
    except IndexError:
        # indexes run out, set the index to the first
        h_index = 0
        # redo the label update
        hvalue_string = str(z[h_index])
        axis_height_txt.configure(text=str(round(eval(hvalue_string), 2)))
        # show warning about that
        tk.messagebox.showwarning('INDEX ERROR', 'you selected a value in slider that might no longer be avalaible, it has been changed avoid errors')
    # from these, establish the new fract, approperiate for 2-forms
    fract = 2/((pt_den-1))
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
    max_global_dxdy = find_global_max(form_2_str_dxdy)
    max_global_dxdz = find_global_max(form_2_str_dxdz)
    max_global_dydz = find_global_max(form_2_str_dydz)
    eq_1_xy = str(simplify('(' + form_2_str_dxdy + ')' + '*(' + str(1/2) + ')'))
    eq_2_xy = str(simplify('(' + form_2_str_dxdy + ')' + '*(' + str(1/2) + ')'))
    F_xy_x, F_xy_y = eq_to_comps(eq_1_xy, eq_2_xy, xg, yg)
    eq_1_xz = str(simplify('(' + form_2_str_dxdz + ')' + '*(' + str(1/2) + ')'))
    eq_2_xz = str(simplify('(' + form_2_str_dxdz + ')' + '*(' + str(1/2) + ')'))
    F_xz_x, F_xz_z = eq_to_comps(eq_1_xz, eq_2_xz, xg, yg)
    eq_1_yz = str(simplify('(' + form_2_str_dydz + ')' + '*(' + str(1/2) + ')'))
    eq_2_yz = str(simplify('(' + form_2_str_dydz + ')' + '*(' + str(1/2) + ')'))
    F_yz_y, F_yz_z = eq_to_comps(eq_1_yz, eq_2_yz, xg, yg)
    # call the slide function to easily plot them depending on h and on axis_view
    slide()
    # into the entry boxes for 2-forms, but the result down
    Entry_2_form_R3_dxdy.delete(0, 'end')
    Entry_2_form_R3_dxdz.delete(0, 'end')
    Entry_2_form_R3_dydz.delete(0, 'end')
    Entry_2_form_R3_dxdy.insert(0, str(form_2_str_dxdy))
    Entry_2_form_R3_dxdz.insert(0, str(form_2_str_dxdz))
    Entry_2_form_R3_dydz.insert(0, str(form_2_str_dydz))
    # colour the 2-forms green to show that these display now
    Entry_2_form_R3_dxdy.configure(bg='#C0F6BB')
    Entry_2_form_R3_dxdz.configure(bg='#C0F6BB')
    Entry_2_form_R3_dydz.configure(bg='#C0F6BB')


# define a function to update z Index and redraw the plot based on the slider
def slide():
    global hvalue_string, s_max
    # remove the currently displayed plot
    main_axis.clear()
    # CANT UPDATE VARIABLES HERE BECUASE THIS WILL NOT CHNAGE THE CALCULATED
    # RESULT FROM aBOVE TEHREFORE SHOW WARNING IF THIS CAUSES AN ERROR
    # check if that is needed, and if so, break
    if L != eval(L_entry.get()) or pt_den != int(pt_den_entry.get()):
        # show warning
        tk.messagebox.showwarning('PLOT', 'You changed the variables, replot first (on the right or below)')
        # end function
        return None
    # allow only s_max to be globally changed here
    s_max = int(s_max_entry.get())
    # replot the graph with that new h_index
    # and change the label under the slider to be the value of the chosen axis at that h_index
    if axis_view == 'z':
        form_2_components_plot_3(xg, yg, h_index, axis_view, F_xy_x, zero_field, s_max, L, pt_den, fract, colour_str, 2)
        form_2_components_plot_3(xg, yg, h_index, axis_view, zero_field, F_xy_y, s_max, L, pt_den, fract, colour_str, 2)
        main_axis.set_xlabel('$x$')
        main_axis.set_ylabel('$y$')
    elif axis_view == 'y':
        form_2_components_plot_3(xg, zg, h_index, axis_view, F_xz_x, zero_field, s_max, L, pt_den, fract, colour_str, 2)
        form_2_components_plot_3(xg, zg, h_index, axis_view, zero_field, F_xz_z, s_max, L, pt_den, fract, colour_str, 2)
        main_axis.set_xlabel('$x$')
        main_axis.set_ylabel('$z$')
    elif axis_view == 'x':
        form_2_components_plot_3(yg, zg, h_index, axis_view, F_yz_y, zero_field, s_max, L, pt_den, fract, colour_str, 2)
        form_2_components_plot_3(yg, zg, h_index, axis_view, zero_field, F_yz_z, s_max, L, pt_den, fract, colour_str, 2)
        main_axis.set_xlabel('$y$')
        main_axis.set_ylabel('$z$')
    # draw that onto the screen
    canvas.draw()


# deifne a function that will update the label
def label_update(var):
    global h_index, hvalue_string
    # update current height, remembering to account for ends of array
    if h_index == 0:
        if var == -1:
            tk.messagebox.showerror('WARNING', 'REACHED END INDEX')
        else:
            h_index += var
            hvalue_string = str(z[h_index])
    else:
        try:
            h_index += var
            hvalue_string = str(z[h_index])
        except IndexError:
            h_index -= var
            tk.messagebox.showerror('WARNING', 'REACHED END INDEX')
    try:
        # update the label
        axis_height_txt.configure(text=str(round(eval(hvalue_string), 2)))
    except UnboundLocalError:
        pass


# define a function that will repond to changing axis view with radiobuttons
def view_response(view_var):
    # get and make global the axis view variable
    global axis_view
    axis_view = view_tk.get()
    # call slide function to plot based on the wanted axis view and height.
    slide()
    # draw that onto the screen
    canvas.draw()


# define a fucntion that will plot 2-forms on R3 directly from user input
def form_2_R3_direct_plot():
    global max_global_dxdy, max_global_dxdz, max_global_dydz, F_xy_x, F_xy_y
    global F_xz_x, F_xz_z, F_yz_y, F_yz_z, form_2, form_2_sgn
    global form_2_str_dxdy, form_2_str_dxdz, form_2_str_dydz
    global hvalue_string, h_index, fract
    # get the updated variables
    update_variables(1)
    # update the h_index label based on the newly input arrays
    try:
        hvalue_string = str(z[h_index])
        axis_height_txt.configure(text=str(round(eval(hvalue_string), 2)))
    except IndexError:
        # indexes run out, set the index to the first
        h_index = 0
        # redo the label update
        hvalue_string = str(z[h_index])
        axis_height_txt.configure(text=str(round(eval(hvalue_string), 2)))
        # show warning about that
        tk.messagebox.showwarning('INDEX ERROR', 'you selected a value in slider that might no longer be avalaible, it has been changed avoid errors')
    # from these, establish the new fract, approperiate for 2-forms
    fract = 2/((pt_den-1))
    # get the strings from entry boxes
    form_2_str_dxdy = str(simplify(Entry_2_form_R3_dxdy.get()))
    form_2_str_dxdz = str(simplify(Entry_2_form_R3_dxdz.get()))
    form_2_str_dydz = str(simplify(Entry_2_form_R3_dydz.get()))
    # from these, get the signs of the new 2-form for colours
    form_2 = np.empty((3, pt_den, pt_den, pt_den))
    form_2[0, :, :, :] = eval(format_eq(form_2_str_dxdy))
    form_2[1, :, :, :] = eval(format_eq(form_2_str_dxdz))
    form_2[2, :, :, :] = eval(format_eq(form_2_str_dydz))
    # get global maxima
    max_global_dxdy = find_global_max(form_2_str_dxdy)
    max_global_dxdz = find_global_max(form_2_str_dxdz)
    max_global_dydz = find_global_max(form_2_str_dydz)
    eq_1_xy = str(simplify('(' + form_2_str_dxdy + ')' + '*(' + str(1/2) + ')'))
    eq_2_xy = str(simplify('(' + form_2_str_dxdy + ')' + '*(' + str(1/2) + ')'))
    F_xy_x, F_xy_y = eq_to_comps(eq_1_xy, eq_2_xy, xg, yg)
    eq_1_xz = str(simplify('(' + form_2_str_dxdz + ')' + '*(' + str(1/2) + ')'))
    eq_2_xz = str(simplify('(' + form_2_str_dxdz + ')' + '*(' + str(1/2) + ')'))
    F_xz_x, F_xz_z = eq_to_comps(eq_1_xz, eq_2_xz, xg, yg)
    eq_1_yz = str(simplify('(' + form_2_str_dydz + ')' + '*(' + str(1/2) + ')'))
    eq_2_yz = str(simplify('(' + form_2_str_dydz + ')' + '*(' + str(1/2) + ')'))
    F_yz_y, F_yz_z = eq_to_comps(eq_1_yz, eq_2_yz, xg, yg)
    # call the slide function to easily plot them depending on h and on axis_view
    slide()
    # colour the 2-forms green to show that these display now
    Entry_2_form_R3_dxdy.configure(bg='#C0F6BB')
    Entry_2_form_R3_dxdz.configure(bg='#C0F6BB')
    Entry_2_form_R3_dydz.configure(bg='#C0F6BB')
    # show other ones as white
    form_1_x_entry.configure(bg='#FFFFFF')
    form_1_y_entry.configure(bg='#FFFFFF')
    form_1_z_entry.configure(bg='#FFFFFF')


'''

Define function for the dynamics tab

'''


# function to get initial conditionsl ist
def ode1(xy, t):
    # Unpack initial conditions
    x, y = xy
    # List of expressions dx_dt and dy_dt
    string_x_dyn = str(x_comp_entry.get())
    string_x_dyn = string_x_dyn.replace('^', '**')
    string_x_dyn = string_x_dyn.replace('ln', 'log')
    string_y_dyn = str(y_comp_entry.get())
    string_y_dyn = string_y_dyn.replace('^', '**')
    string_y_dyn = string_y_dyn.replace('ln', 'log')
    L = [eval(string_x_dyn), eval(string_y_dyn)]
    return L


# response function to matplotlibs animate, supplied next points
def animate(i):
    global dyn_point, poly
    
    xplot = eval(x_dyn_str)
    yplot = eval(y_dyn_str)
    dyn_point.set_data(xplot, yplot)
    
    if dyn_join_shapes.get() == 1:
        poly_plot = eval(poly_str)
        poly = mpl.patches.Polygon(poly_plot, fill=True, color='blue')
        main_axis.add_artist(poly)
        return dyn_point, poly
    else:
        return dyn_point,


# function to carry out the animating
def animation_storing_function():
    global dyn_N
    ani = animation.FuncAnimation(fig, animate, dyn_N, interval=25, blit=True, repeat=False)
    return ani


# function to respond to button to begin the animation.
def animate_response():
    global dummy_variable_dyn, dyn_time
    global dyn_coord, x_dyn_str, y_dyn_str, poly_str, dyn_N, tmax
    # clear the axis and redraw
    PLOT_response(0)
    x_dyn_str = ''
    y_dyn_str = ''
    poly_str = ''
    dyn_N = int(round(dyn_N_slider.get(),0))
    tmax = tmax_slider.get()
    dyn_N = int(1000*tmax/50)  # int(round(dyn_N_slider.get(), 0))
    dyn_time = np.linspace(0, tmax, dyn_N)
    for a in range(len(dyn_coord)):
        exec('global ' +  'xy' + str(a) + '\n'
              'xy' + str(a) + ' = odeint(ode1, dyn_coord[a], dyn_time)')
        x_dyn_str += 'xy' + str(a) + '[i,0], '
        y_dyn_str += 'xy' + str(a) + '[i,1], '
        poly_str += '[xy'+ str(a) + '[i,0], xy' + str(a) + '[i,1]],'
    dummy_variable_dyn = animation_storing_function()


# pauses the animation
def pause_response():
    global dummy_variable_dyn
    dummy_variable_dyn.event_source.stop()


# plays it again if paused
def play_response():
    global dummy_variable_dyn
    dummy_variable_dyn.event_source.start()


# clears the grid and restes variables
def clear_response():
    global dummy_variable_dyn, dyn_coord
    try:
        dummy_variable_dyn.event_source.stop()
    except (NameError, AttributeError):
        pass
    PLOT_response(0)
    # delete the extra created variables:
    try:
        for a in range(len(dyn_coord)):
            exec('global ' + 'xy' + str(a) + '\n' + 'del ' + 'xy' + str(a))
    except NameError:
        pass
    # then clear coordinates
    dyn_coord = []


# Button to enable/ disable joining shapes in straight lines
def dyn_join_response():
    global dyn_join_shapes
    if dyn_join_shapes.get() == 0:
        # the button is off, and has been clicked therefore change the
        # variable to an and the image to on
        dyn_join_shapes.set(1)
        dyn_join_toggle.configure(image=toggle_image_on)

    else:
        # the button is on and has been clicked
        # set it to off and change image
        dyn_join_shapes.set(0)
        dyn_join_toggle.configure(image=toggle_image_off)


# funcction to respond to dropdown for dyn shapes
def dyn_shape_select_response(selected_shape):
    # deal with lines
    global dyn_shape_select, size_dyn_shape
    global dyn_pre_size_label, dyn_pre_size_entry
    # get the chosen shape
    dyn_shape_select.set(selected_shape)
    # depending on that selection, prepare rest of frame:
    if selected_shape == 'Polygon':
        # restart the plot and the integral values
        clear_response()
        # get rid of size options
        try:
            dyn_pre_size_label.destroy()
            dyn_pre_size_entry.destroy()
        except (UnboundLocalError, NameError):  
            pass
    else:
        # restart the plot and the integral values
        clear_response()
        # if circle is selected, display an entry box for the radius of it
        dyn_pre_size_label = tk.Label(dynamics_frame, text='Size')
        dyn_pre_size_label.grid(row=6, column=0)
        dyn_pre_size_entry = tk.Entry(dynamics_frame, width=10)
        dyn_pre_size_entry.grid(row=6, column=1)
        if selected_shape == 'Human':
            dyn_pre_size_entry.insert(0, '0.4')
        else:
            dyn_pre_size_entry.insert(0, '1')


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
tensor_label.grid(row=8, column=0)

# define each button and put them on the screen, in the right_frame
arrow_btn = tk.Radiobutton(right_frame, text='arrow', variable=tensor, value=1, command=lambda: vect_type_response(tensor.get())).grid(row=8, column=1)
arrow_stack_btn = tk.Radiobutton(right_frame, text='both', variable=tensor, value=2, command=lambda: vect_type_response(tensor.get())).grid(row=8, column=3)
stack_btn = tk.Radiobutton(right_frame, text='stack', variable=tensor, value=0, command=lambda: vect_type_response(tensor.get())).grid(row=8, column=2)

# DERIVATIVE FUNCTIONS

# Radiobuttons to select what happens when clicking the plot
click_option = tk.IntVar()
click_option.set(0)
click_option_Tools_btn = tk.Radiobutton(right_frame, text='Tools', variable=click_option, value=0, command=lambda: click_option_handler(click_option.get()))
click_option_Zoom_btn = tk.Radiobutton(right_frame, text='Zoom', variable=click_option, value=1, command=lambda: click_option_handler(click_option.get()))
click_option_Deriv_btn = tk.Radiobutton(right_frame, text='Deriv.', variable=click_option, value=2, command=lambda: click_option_handler(click_option.get()))
click_option_Div_btn = tk.Radiobutton(right_frame, text='Div.', variable=click_option, value=3, command=lambda: click_option_handler(click_option.get()))
click_option_Curl_btn = tk.Radiobutton(right_frame, text='Curl', variable=click_option, value=4, command=lambda: click_option_handler(click_option.get()))
click_option_Deriv_btn.configure(state=tk.DISABLED)
click_option_Div_btn.configure(state=tk.DISABLED)
click_option_Curl_btn.configure(state=tk.DISABLED)
click_option_Tools_btn.grid(row=1, column=0)
click_option_Zoom_btn.grid(row=1, column=1)
click_option_Deriv_btn.grid(row=1, column=2)
click_option_Div_btn.grid(row=2, column=0)
click_option_Curl_btn.grid(row=2, column=1)

# Zooming window zoom slider
tk.Label(right_frame, text='Zoom').grid(row=3, column=0)
zoom_slider = tk.Scale(right_frame, from_=1, to=100, orient=tk.HORIZONTAL, resolution=1)
zoom_slider.bind("<ButtonRelease-1>", update_deriv)
zoom_slider.grid(row=3, column=1)

# Drop down to select the derivative plot point density (dpd)
dpd_select = tk.IntVar()
dpd_select.set(5)
dpd_list = [5, 7, 9]

tk.Label(right_frame, text='Inset Plot Point Density:').grid(row=4, column=0)
dpd_drop = tk.OptionMenu(right_frame, dpd_select, *dpd_list, command = update_deriv)
dpd_drop.grid(row=4, column=1)

# Drop down to select inset axis size (d_length)
d_length_select = tk.DoubleVar()
d_length_list = [0.1, 0.2, 0.3, 0.4, 0.5]
d_length_select.set(d_length_list[2])
tk.Label(right_frame, text='Inset Fractional Size:').grid(row=5, column=0)
d_length_drop = tk.OptionMenu(right_frame, d_length_select, *d_length_list, command=update_deriv)
d_length_drop.grid(row=5, column=1)

# Autoscale Toggle
ascale_label = tk.Label(right_frame, text='Autoscale arrows:')
ascale_label.grid(row=7, column=0)
ascale_toggle = tk.Button(right_frame, image=toggle_image_off, bd=0, command=scale_toggle_response)
ascale_toggle.grid(row=7, column=1, pady=5)

# define entry boxes to allow user to input x_m and y_m
x_m_entry = tk.Entry(right_frame, width=12)
y_m_entry = tk.Entry(right_frame, width=12)
x_m_entry.grid(row=6, column=0)
y_m_entry.grid(row=6, column=1)
# and a button to submit these:
Set_target_btn = tk.Button(right_frame, text='Set Target', command=set_inset_target)
Set_target_btn.grid(row=6, column=2, padx=20)



'''

SINGULARITY NOTEBOOK

'''

# get a button to draw on singularities
singularity_button = tk.Button(singular_frame, text='search singularities', command=show_singularities)
singularity_button.grid(row=0, column=0)
# entry for N
tk.Label(singular_frame, text='<- sampling points').grid(row=0, column=2, columnspan=2)
fine_grid_N_entry = tk.Entry(singular_frame, width=5)
fine_grid_N_entry.grid(row=0, column=1)
fine_grid_N_entry.insert(0, 10)

# define an entry where the user can inpu known singularity equation
# this will be taken and plotted as a red, dotted line
tk.Label(singular_frame, text='singularity equation:').grid(row=1, column=0, columnspan=2)

# define a dropdown to select y= or x=
singular_var = tk.StringVar()
singular_list = ['y=', 'x=', 'point']
singular_var.set(singular_list[0])
dpd_drop = tk.OptionMenu(singular_frame, singular_var, *singular_list, command=singular_drop_response)
dpd_drop.grid(row=2, column=0)
# equation entry box
known_singularity_entry = tk.Entry(singular_frame, width=15)
known_singularity_entry.grid(row=2, column=1)
known_singularity_entry.insert(0, '')

# define asubmit button to that entry
submit_known_singularity_btn = tk.Button(singular_frame, text='show expression', command=known_singularity_response)
submit_known_singularity_btn.grid(row=3, column=0)


'''

set up all in BOTTOM FRAME

'''


# define entry boxes for the field equations in x and y
component_x_entry_label = tk.Label(bot_frame, text='dx component')
component_x_entry_label.grid(row=0, column=0)
x_comp_entry = tk.Entry(bot_frame, width=20, borderwidth=2)
x_comp_entry.grid(row=1, column=0)
x_comp_entry.insert(0, 'y*sin(x)')

component_y_entry_label = tk.Label(bot_frame, text='dy component')
component_y_entry_label.grid(row=0, column=1)
y_comp_entry = tk.Entry(bot_frame, width=20, borderwidth=2)
y_comp_entry.grid(row=1, column=1)
y_comp_entry.insert(0, '-x*cos(y)')


# define strings from initial components 
# these are needed by the derivative function, therefore for the derivative
# to work on the initial field, need to initially define them
# THESE ARE HERE TO APPEAR AFTER THE ENTRIES ARE ACTUALLY DEFINED
# NOT VITAL, BUT EASIER
string_x = str(x_comp_entry.get())
string_y = str(y_comp_entry.get())

# set up a dropdown box for 1-forms and VFs
field_select = tk.StringVar()
field_select.set(field_name_list[0])

field_select_drop_label = tk.Label(bot_frame, text='Select Pre-Defined 1-Form:')
field_select_drop_label.grid(row=2, column=0, columnspan=2)


field_select_drop = ttk.Combobox(bot_frame, value=field_name_list, width=40)
field_select_drop.current(0)
field_select_drop.grid(row=3, column=0, columnspan=2)
field_select_drop.bind("<<ComboboxSelected>>", field_selection_response)


'''

what was in small frame initally, now also in botton notebook

'''

# define the PLOT button
PLOT_btn = tk.Button(small_frame, text='PLOT', padx=40, pady=20, command=PLOT_response)
PLOT_btn.grid(row=0, column=2, rowspan=2)

# define a button that will just plot the given cartesian field
# on a polar grid
polar_grid_plot_btn = tk.Button(small_frame, text='Polar plot', padx=20, command= lambda: Polar_grid_plot_response(tensor.get()))
polar_grid_plot_btn.grid(row=1, column=1)

# define a button in small frame that will open new window to adjust arrowheads
custom_btn = tk.Button(small_frame, text='Visuals Customise', padx=1, pady=1, command=custom_btn_reponse)
custom_btn.grid(row=0, column=0)

# define a button to customise the polar grids
polar_grid_custom_btn = tk.Button(small_frame, text='Polar customise', padx=7, pady=1, command=polar_grid_custom_reponse)
polar_grid_custom_btn.grid(row=1, column=0)

# define entry boxes for each (in order): L, pt_den, s_max and a ; and info txt
# Also input into them the initial values
tk.Label(small_frame, text='Size').grid(row=2, column=0)
L_entry = tk.Entry(small_frame, width=5, borderwidth=1)
L_entry.grid(row=3, column=0, padx=2)
L_entry.insert(0, L)

tk.Label(small_frame, text='grid').grid(row=2, column=1)
pt_den_entry = tk.Entry(small_frame, width=5, borderwidth=1)
pt_den_entry.grid(row=3, column=1, padx=2)
pt_den_entry.insert(0, pt_den)

tk.Label(small_frame, text='max sheets').grid(row=2, column=2)
s_max_entry = tk.Entry(small_frame, width=5, borderwidth=1)
s_max_entry.grid(row=3, column=2, padx=2)
s_max_entry.insert(0, s_max)



'''

set up all in LI tab

'''


# define a label that will display it
tk.Label(LI_frame, text='Shape Area:').grid(row=0, column=0)
shape_area_label = tk.Label(LI_frame, text=shape_area)
shape_area_label.grid(row=0, column=1)

tk.Label(LI_frame, text='LI Total:').grid(row=1, column=0, padx=10)
LI_total_label = tk.Label(LI_frame, text=LI_total)
LI_total_label.grid(row=1, column=1)

tk.Label(LI_frame, text='Ratio:').grid(row=1, column=2)
ratio1_label = tk.Label(LI_frame, text=ratio1)
ratio1_label.grid(row=1, column=3)

tk.Label(LI_frame, text='Flux:').grid(row=2, column=0)
flux_label = tk.Label(LI_frame, text=flux)
flux_label.grid(row=2, column=1)

tk.Label(LI_frame, text='Ratio:').grid(row=2, column=2)
ratio2_label = tk.Label(LI_frame, text=ratio2)
ratio2_label.grid(row=2, column=3)

# display a restart button that will clear the lines
# and restart the variables.
LI_restart_btn = tk.Button(LI_frame, text='LI Restart', padx=20, command=LI_restart)
LI_restart_btn.grid(row=3, column=0, columnspan=2)

# define a drop down to draw: connected lines, square or circle.
LI_shape_select = tk.StringVar()
LI_shape_list = ['Polygon', 'Circle']
LI_shape_select.set(LI_shape_list[0])
LI_shape_instruction = tk.Label(LI_frame, text='Shape:')
LI_shape_instruction.grid(row=4, column=0)
LI_shape_drop = tk.OptionMenu(LI_frame, LI_shape_select, *LI_shape_list, command=LI_shape_select_response)
LI_shape_drop.grid(row=4, column=1)


# input the radiobuttons for arrows, stacks or both again here
arrow_btn = tk.Radiobutton(LI_frame, text='arrow', variable=tensor, value=1, command=lambda: vect_type_response(tensor.get())).grid(row=7, column=1)
stack_btn = tk.Radiobutton(LI_frame, text='stack', variable=tensor, value=0, command=lambda: vect_type_response(tensor.get())).grid(row=7, column=2)

# Radiobutton for showing flux/circulation
showcol = tk.IntVar()
showcol.set(0)
tk.Label(LI_frame, text='Colour Curve:').grid(row=8, column=0)
shownone_btn = tk.Radiobutton(LI_frame, text='None', variable=showcol, value=0).grid(row=8, column=1)
showcirc_btn = tk.Radiobutton(LI_frame, text='Show Circ.', variable=showcol, value=1).grid(row=8, column=2)
showflux_btn = tk.Radiobutton(LI_frame, text='Show Flux', variable=showcol, value=2).grid(row=8, column=3)

# add an entry so the user can choose what grid they want
# by choosing grid lines separation
grid_LI_poly_label = tk.Label(LI_frame, text='Grid Separation:').grid(row=9, column=0)
grid_sep_poly_entry = tk.Entry(LI_frame, width=10)
grid_sep_poly_entry.grid(row=9, column=1)
grid_sep_poly_entry.insert(0, '2')
# define a button to submit these changes
submit_poly_sep_grid = tk.Button(LI_frame, text='Submit Grid', command=poly_grid_submit)
submit_poly_sep_grid.grid(row=9, column=2)

# Autoscale Arrows Toggle
tk.Label(LI_frame, text='Autoscale arrows:').grid(row=10, column=0)
ascale_toggle_LI = tk.Button(LI_frame, image=toggle_image_off, bd=0, command=scale_toggle_response)
ascale_toggle_LI.grid(row=10, column=1, pady=5)


'''

DEFINE ALL WINDGETS IN CALCULUS TAB

'''

# define a window to supply the 2-form
tk.Label(calculus_frame, text='2-form on R2').grid(row=0, column=1)
form_2_entry = tk.Entry(calculus_frame, width=15, borderwidth=2)
form_2_entry.grid(row=0, column=0)
form_2_entry.insert(0, form_2_str)
# begin  displaying it on green colour to show that this is ebing displayed to
# beign with
form_2_entry.configure(bg='#C0F6BB')

# extra label for 0 form.
# for now not an entry box because it doesn't get used anywhere
# I make it red to remember to come back to it later and use it.
Label_zero_form = tk.Label(calculus_frame, text='Zero form:')
Label_zero_form.grid(row=4, column=0)

# set up an entry for zero forms:
form_0_entry = tk.Entry(calculus_frame, width=15, borderwidth=2)
form_0_entry.grid(row=4, column=1)
form_0_entry.insert(0, '')

# set up a button to plot the 0-form
form_0_btn = tk.Button(calculus_frame, text='0-form plot', padx=3, pady=5, command=form_0_response)
form_0_btn.grid(row=4, column=2)

# define a button to submit the supplied 2-form and plot it as blocks
form_2_btn = tk.Button(calculus_frame, text='2-form plot', padx=3, pady=5, command=form_2_response)
form_2_btn.grid(row=3, column=1)

# define a button that will just plot the 1-form
# this will not be needed when its it merged with main GUI
# as there will already be a plot button there
form_1_btn = tk.Button(calculus_frame, text='1-form plot', padx=3, pady=5, command=form_1_stacks_response)
form_1_btn.grid(row=3, column=0)

# add a button to plot the interior derivative as superposing stack fields
INT_btn = tk.Button(calculus_frame, text='Int Deriv', padx=0, pady=2, command=Int_deriv_response)
INT_btn.grid(row=5, column=0)

# define a button to plot the exterior derivative from given u and v
# Note, it will get the 2-form first, then return back down
# to a one form to avoid concellations
# therefore, it will also just be one possible representation
# not the only possible one
EXT_int_btn = tk.Button(calculus_frame, text='Ext Deriv', padx=0, pady=2, command=Ext_deriv_response)
EXT_int_btn.grid(row=5, column=1)

# define a wedge product button that will let the user input TWO 1-forms
# in a new window to be wedged to gice a 2-form
wedge_btn = tk.Button(calculus_frame, text='Wedge', padx=0, pady=2, command=wedge_2_response)
wedge_btn.grid(row=6, column=0)

# define ab utton that will Find the Hodge dual
Hodge_1_form_btn = tk.Button(calculus_frame, text='Hodge 1-Form', padx=0, pady=2, command=Hodge_1_form_response)
Hodge_1_form_btn.grid(row=7, column=0)

Hodge_2_form_btn = tk.Button(calculus_frame, text='Hodge 2-Form', padx=0, pady=2, command=Hodge_2_form_response)
Hodge_2_form_btn.grid(row=7, column=1)

# define radiobuttons button to choose zooming with the mouse on 2-forms on R2
# and as opposed to tools
R2_tools_opt = tk.IntVar()
R2_tools_opt.set(0)
R2_tools_Tools_btn = tk.Radiobutton(calculus_frame, text='Tools', variable=R2_tools_opt, value=0, command=lambda: R2_tools_handler(R2_tools_opt.get()))
R2_tools_Zoom_btn = tk.Radiobutton(calculus_frame, text='Zoom', variable=R2_tools_opt, value=1, command=lambda: R2_tools_handler(R2_tools_opt.get()))
R2_tools_int_btn = tk.Radiobutton(calculus_frame, text='Area Int', variable=R2_tools_opt, value=2, command=lambda: R2_tools_handler(R2_tools_opt.get()))
R2_tools_Tools_btn.grid(row=8, column=0)
R2_tools_Zoom_btn.grid(row=8, column=1)
R2_tools_int_btn.grid(row=8, column=2)

# set up a zooming tool for that too
tk.Label(calculus_frame, text='Zoom').grid(row=9, column=0)
zoom_slider_R2 = tk.Scale(calculus_frame, from_=1, to=20, orient=tk.HORIZONTAL)
zoom_slider_R2.bind("<ButtonRelease-1>", update_2_form_zoom)
zoom_slider_R2.grid(row=9, column=1)
zoom_slider_R2.configure(state=tk.DISABLED)


# Drop down to select the R2 2 form zoom plot point density
zoomR2pd_select = tk.IntVar()
zoomR2pd_select.set(11)
zoomR2pd_list = [5, 6, 10, 11, 15, 16, 20, 21]

tk.Label(calculus_frame, text='Inset Plot Point Density:').grid(row=10, column=0)
zoomR2pd_drop = tk.OptionMenu(calculus_frame, zoomR2pd_select, *zoomR2pd_list, command=update_2_form_zoom)
zoomR2pd_drop.grid(row=10, column=1)

# Drop down to select inset axis size for R2 2 forms
zoomR2_length_select = tk.DoubleVar()
zoomR2_length_list = [0.1, 0.2, 0.3, 0.4, 0.5]
zoomR2_length_select.set(zoomR2_length_list[2])
tk.Label(calculus_frame, text='Inset Fractional Size:').grid(row=11, column=0)
zoomR2_length_drop = tk.OptionMenu(calculus_frame, zoomR2_length_select, *zoomR2_length_list, command=update_2_form_zoom)
zoomR2_length_drop.grid(row=11, column=1)

# define entry boxes to allow user to input x_m and y_m
x_m_entry_calc = tk.Entry(calculus_frame, width=12)
y_m_entry_calc = tk.Entry(calculus_frame, width=12)
x_m_entry_calc.grid(row=12, column=0)
y_m_entry_calc.grid(row=12, column=1)
# and a button to submit these:
Set_target_btn_calc = tk.Button(calculus_frame, text='Set Target', command=set_inset_target_calc)
Set_target_btn_calc.grid(row=12, column=2, padx=20)


# add a dropdown menu for 2-forms, for now at the end of this tab
select_form_2 = tk.StringVar()
select_form_2.set(list_form_2_names[0])

select_form_2_drop_label = tk.Label(calculus_frame, text='Select Pre-Defined 2-Form:')
select_form_2_drop_label.grid(row=1, column=0, columnspan=3)

select_form_2_drop = ttk.Combobox(calculus_frame, value=list_form_2_names, width=40)
select_form_2_drop.current(0)
select_form_2_drop.grid(row=2, column=0, columnspan=3)
select_form_2_drop.bind("<<ComboboxSelected>>", selection_form_2_response)


'''

DEFINE WIDGETS USED IN R3 CODE

'''

height_frame = tk.LabelFrame(r3_frame, text='viewing frame', padx=2, pady=2)
height_frame.grid(row=0, column=0)

# Label to show current axis value
axis_height_txt = tk.Label(height_frame, text=str(z[0]))
axis_height_txt.grid(row=1, column=0)

# on the left, make a 'move down' button
down_height = tk.Button(height_frame, text=' \/ ', command=lambda: label_update(-1))
down_height.grid(row=2, column=0)

# on the right, make a 'move up' button
up_height = tk.Button(height_frame, text=' /\ ', command=lambda: label_update(1))
up_height.grid(row=0, column=0)

# define a button to submit the currently chosen value:
Submit_h_btn = tk.Button(height_frame, text='SUBMIT', padx=2, pady=50, command=slide)
Submit_h_btn.grid(row=0, column=1, rowspan=3, padx=5)


# define rediobuttons to chose from which axis the user is looking:
view_tk = tk.StringVar()
view_tk.set('z')
view_z_btn = tk.Radiobutton(height_frame, text='z', variable=view_tk, value='z', command=lambda: view_response(view_tk.get())).grid(row=0, column=2)
view_z_btn = tk.Radiobutton(height_frame, text='y', variable=view_tk, value='y', command=lambda: view_response(view_tk.get())).grid(row=1, column=2)
view_z_btn = tk.Radiobutton(height_frame, text='x', variable=view_tk, value='x', command=lambda: view_response(view_tk.get())).grid(row=2, column=2)


# NOTE  NOT GREAT I KNOW BUT TEMPORARY:
# define a new frame for the fields to be input
field_input_frame = tk.LabelFrame(r3_frame, text='Fields frame', padx=5, pady=5)
field_input_frame.grid(row=0, column=1)


# define a button that will let the user chose the splitting option
# for 2-forms plotted as stacks.
# define entry boxes for the three 1-forms that are being plotted
# define entries for a 1-form
tk.Label(field_input_frame, text='dx comp. 1-form').grid(row=1, column=0)
form_1_x_entry = tk.Entry(field_input_frame, width=20, borderwidth=2)
form_1_x_entry.grid(row=2, column=0, columnspan=2)
form_1_x_entry.insert(0, string_x)

tk.Label(field_input_frame, text='dy comp. 1-form').grid(row=3, column=0)
form_1_y_entry = tk.Entry(field_input_frame, width=20, borderwidth=2)
form_1_y_entry.grid(row=4, column=0, columnspan=2)
form_1_y_entry.insert(0, string_y)

tk.Label(field_input_frame, text='dz comp. 1-form').grid(row=5, column=0)
form_1_z_entry = tk.Entry(field_input_frame, width=20, borderwidth=2)
form_1_z_entry.grid(row=6, column=0, columnspan=2)
form_1_z_entry.insert(0, string_z)

# deifne a button to plot from these:
form_2_from_1_R3_btn = tk.Button(field_input_frame, text='Ext. Deriv. R3', padx=3, pady=5, command=form_1_to_2onR3)
form_2_from_1_R3_btn.grid(row=7, column=0, columnspan=2)

# define a frame for R3 2-form results and input by the user
# It will be filled when R3 is opened (in tab_selection)
# as it needs 'result' to fill.
form_2_frame = tk.LabelFrame(r3_frame, text='2-form frame', padx=5, pady=5)
form_2_frame.grid(row=2, column=0, columnspan=2)
 

# set up the elemental 2-form labels there
Label_2_form_xy = tk.Label(form_2_frame, text='  dx^dy').grid(row=0, column=1)
Label_2_form_xz = tk.Label(form_2_frame, text='  dx^dz').grid(row=1, column=1)
Label_2_form_yz = tk.Label(form_2_frame, text='  dy^dz').grid(row=2, column=1)

# set up entry boxes for 2-forms, keep them empty for now
# they will fill when R3 is opened (look in tab_selection())
Entry_2_form_R3_dxdy = tk.Entry(form_2_frame, width=20)
Entry_2_form_R3_dxdy.grid(row=0, column=0)
Entry_2_form_R3_dxdz = tk.Entry(form_2_frame, width=20)
Entry_2_form_R3_dxdz.grid(row=1, column=0)
Entry_2_form_R3_dydz = tk.Entry(form_2_frame, width=20)
Entry_2_form_R3_dydz.grid(row=2, column=0)


# define a button that will plot the supplied 2-forms directly
form_2_R3_direct_plot_btn = tk.Button(form_2_frame, text='PLOT 2-form', padx=10, command=form_2_R3_direct_plot)
form_2_R3_direct_plot_btn.grid(row=3, column=0, columnspan=2)


'''

Dynamics Frame

'''

arrow_btn = tk.Radiobutton(dynamics_frame, text='arrow', variable=tensor, value=1, command=lambda: vect_type_response(tensor.get())).grid(row=0, column=0)
stack_btn = tk.Radiobutton(dynamics_frame, text='stack', variable=tensor, value=0, command=lambda: vect_type_response(tensor.get())).grid(row=0, column=1)

# add a button that will respond with animating the clikced points
animate_btn = tk.Button(dynamics_frame, text='Animate', padx=10, pady=10, command=animate_response)
animate_btn.grid(row=1, column=0, padx=5, pady=5)

pause_btn = tk.Button(dynamics_frame, text='Stop', command=pause_response)
pause_btn.grid(row=1, column=1)

play_btn = tk.Button(dynamics_frame, text='Play', command=play_response)
play_btn.grid(row=1, column=2)

clear_btn = tk.Button(dynamics_frame, text='Clear', command=clear_response)
clear_btn.grid(row=1, column=3)

tk.Label(dynamics_frame, text='Animation Frames:').grid(row=2, column=0)
dyn_N_slider = tk.Scale(dynamics_frame, from_=dyn_N, to=500, orient=tk.HORIZONTAL)
dyn_N_slider.grid(row=2, column=1)

tk.Label(dynamics_frame, text='Tmax:').grid(row=3, column=0)
tmax_slider = tk.Scale(dynamics_frame, from_=tmax, to=30, orient=tk.HORIZONTAL)
tmax_slider.grid(row=3, column=1)

# Button to enable to disable straight line shape joining
tk.Label(dynamics_frame, text='Join to shapes').grid(row=4, column=0)
dyn_join_toggle = tk.Button(dynamics_frame, image=toggle_image_on, bd=0, command=dyn_join_response)
dyn_join_toggle.grid(row=4, column=1)

# define a dropdown menu for predefined shape options
dyn_shape_select = tk.StringVar()
dyn_shape_list = ['Polygon', 'Square', 'Human', 'Circle']
dyn_shape_select.set(dyn_shape_list[0])
dyn_shape_instruction = tk.Label(dynamics_frame, text='Shape:').grid(row=5, column=0)
dyn_shape_drop = tk.OptionMenu(dynamics_frame, dyn_shape_select, *dyn_shape_list, command=dyn_shape_select_response)
dyn_shape_drop.grid(row=5, column=1)


# return time to run
stop = timeit.default_timer()
print('Time: ', stop - start)

# keep the program running until otherwise specified by user
tk.mainloop()
