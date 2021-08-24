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
    global ax_L, main_axis, form_2, form_2_sgn
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
        # and get the signs of the 2-form
        form_2_sgn = np.sign(form_2)
    else:
        x = np.linspace(-L, L, pt_den)
        y = np.linspace(-L, L, pt_den)
        z = np.linspace(-L, L, pt_den)
        xg, yg, zg = np.meshgrid(x, y, z)
        # define zero field
        zero_field = np.zeros(np.shape(xg))
        # define the unit field
        field_unit = np.ones(np.shape(xg))
        # get 2-form to update too
        #form_2 = eval(form_2_eq)
        # and get the signs of the 2-form
        #form_2_sgn = np.sign(form_2)


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
    
    # Scaling of axes and setting equal proportions circles look like circles
    # axis.set_aspect('equal')
    # ax_L = L + L/delta_factor
    # axis.set_xlim(-ax_L, ax_L)
    # axis.set_ylim(-ax_L, ax_L)
    
    # Account for change to grid centre for divergence plot
    if axis_check == 1:
        # if click_opt_int > 2:       
        #     axis.set_xlim(-L-L/5, L+L/5)
        #     axis.set_ylim(-L-L/5, L+L/5)
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


# define a function to take care of tab changes
def tab_selection(event):
    global x_m, y_m, R3_use_track, click_opt_int
    global fract, form_2, notebook_bottom, m
    global expressions, coords, form_2_str, form_2_eq, form_2_sgn
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
        form_2_sgn = np.sign(form_2)
        # 2-form was returned to defualt on R2 so put that in the entry box
        form_2_entry.delete(0, 'end')
        form_2_entry.insert(0, form_2_str)
        # restore the fields input frame
        notebook_bottom.add(bot_frame, text='fields')
        notebook_bottom.select(0)
        # enable the plot and polar buttons again
        PLOT_btn['state'] = tk.NORMAL
        polar_grid_plot_btn['state'] = tk.NORMAL
    if tab_text == 'Line Integrals':
        global toolbar, LI_coord, LI_total, LI_shape_select, tensor
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
        click_opt_int = 5
        # chenage fract to a better one for 1-forms
        fract = 0.05
        # deal with consistancy if both was initially selected
        if tensor.get() == 2:
            tensor.set(0)
            vect_type_response(tensor.get())
        # return the option to default
        LI_shape_select.set(LI_shape_list[0])
        LI_shape_select_response('Polygon')
        # restrat the LI calculations and replot to start from default
        LI_restart()
        # set up grid initally as polygon is selected
        main_axis.grid(True)
        # draw it on
        canvas.draw()
        # enable plot buttons again
        PLOT_btn['state'] = tk.NORMAL
        polar_grid_plot_btn['state'] = tk.NORMAL#
    elif tab_text == 'Main':
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
    elif tab_text == 'Calculus':
        main_axis.clear()
        # get globals form entry boxes
        update_variables()
        # from these, establish the new fract, approperiate for 2-forms
        fract = 2/((pt_den-1))
        # put the initial plot onto the canvas
        form_2_components_plot(xg, yg, form_2/2, zero_field, form_2_sgn, s_max, L, fract, colour_str, 2)
        form_2_components_plot(xg, yg, zero_field, form_2/2, form_2_sgn, s_max, L, fract, colour_str, 2)
        canvas.draw()
        # disable the PLOT button and the ploar grid button
        PLOT_btn['state'] = tk.DISABLED
        polar_grid_plot_btn['state'] = tk.DISABLED
        # clear the zero form label as plot comes back to default
        Label_zero_form.configure(text='')
        # this is now never Vector fields and never arrows therefore
        # set the labels as such
        component_x_entry_label.configure(text='dx component')
        component_y_entry_label.configure(text='dy component')
        field_select_drop_label.configure(text='Select Pre-Defined 1-Form:')
        # change frame name too
        bot_frame_frame.configure(text='1-Form input frame')
    elif tab_text == '\mathbb{R}^{3}':
        global form_2_frame
        global F_xy_x, F_xy_y, F_xz_x, F_xz_z, F_yz_y, F_yz_z
        global max_global_dxdy, max_global_dxdz, max_global_dydz
        global field_select, h_index, hvalue_string
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
        # hide the VFA input frame for now
        notebook_bottom.hide(0)
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
    # if anything but the main window is selected, change to tools
    if tab_text != 'Main':
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

# right frame:
right_frame_frame = tk.LabelFrame(root, text='Options Frame', padx=5, pady=5)
right_frame_frame.grid(row=1, column=1)

# bot frame:
bot_frame_frame = tk.LabelFrame(root, text='1-Form Input Frame', padx=5, pady=5)
bot_frame_frame.grid(row=2, column=0)

# plot frame:
plot_frame = tk.LabelFrame(root, text='Plot Frame', padx=5, pady=5)
plot_frame.grid(row=1, column=0)

# plot characteristics frame and plot button
small_frame = tk.LabelFrame(root, text='Plot Customisation Frame', padx=5, pady=5)
small_frame.grid(row=2, column=1)

# define notebook for tabs
notebook = ttk.Notebook(right_frame_frame)
notebook.grid(row=0, column=0)
# notebook for bottom field input, when needed to disappear.
notebook_bottom = ttk.Notebook(bot_frame_frame)
notebook_bottom.grid(row=0, column=0)

# main options:
right_frame = tk.LabelFrame(notebook)
right_frame.grid(row=0, column=1)
# field input
bot_frame = tk.LabelFrame(notebook_bottom)
bot_frame.grid(row=0, column=0)
# Line integrals
LI_frame = tk.Frame(notebook)
LI_frame.grid(row=0, column=2)
# calculus
calculus_frame = tk.Frame(notebook)
calculus_frame.grid(row=0, column=2)
# R3
r3_frame = tk.Frame(notebook)
r3_frame.grid(row=0, column=3)

# finsalise them
notebook.add(right_frame, text='Main')
notebook.add(LI_frame, text='Line Integrals')
notebook.add(calculus_frame, text='Calculus')
notebook.add(r3_frame, text='\mathbb{R}^{3}')
notebook_bottom.add(bot_frame, text='fields')

# bind the clicks on tabs to a function
notebook.bind_all('<<NotebookTabChanged>>', tab_selection)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# set up initial parameters and plot the initial graph, put it in plot frame
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

'''  BASIC PARAMETERS  '''

# define scale of the graph
L = 5
pt_den = 11 # number of points on each axis

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
                   'Gravitational/Electric Point Charge: -x/(x**2+y**2)dx + -y/(x**2+y**2)dy',
                   'Magnetism of Current Carrying Wire: -y/(x**2+y**2)dx + x/(x**2+y**2)dy',
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

# set up line intergal enpty variables
LI_coord = []
LI_total = 0
flux = 0
shape_area = 0
ratio1 = 0
ratio2 = 0


''' DEFINE CALCULUS PARAMETERS, 2-forms and R3 stuff '''

# set up initial strings for 2-forms window to display, for it to save properly after
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

# get the signs of the 2-form
form_2_sgn = np.sign(form_2)

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
            # if its negative correct and put it into the entry box
            if Radius_LI_circ < 0:
                Radius_LI_circ *= -1
                Radius_LI_circ_entry.delete(0, 'end')
                Radius_LI_circ_entry.insert(0 , str(Radius_LI_circ))
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
    LI_coord = []
    LI_total = 0
    flux = 0
    shape_area = 0
    ratio1 = 0
    ratio2 = 0
    # update the label
    LI_total_label.configure(text=LI_total)
    flux_label.configure(text=flux)
    shape_area_label.configure(text=shape_area)
    ratio1_label.configure(text=ratio1)
    ratio2_label.configure(text=ratio2)
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
        # get rid of grid specifications also
        # no need for try and except, frame opens with these, these will
        # always be there
        grid_LI_poly_label.destroy()
        grid_sep_poly_entry.destroy()
        submit_poly_sep_grid.destroy()
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
        # get rid of grid on plot
        main_axis.grid(False)
        # update canvas
        canvas.draw()
    else:
        # restart the plot and the integral values
        LI_restart()
        # set up a grid
        main_axis.grid(True)
        # add an entry so the user can choose what grid they want
        # by choosing grid lines separation
        grid_LI_poly_label = tk.Label(LI_frame, text='grid separation:')
        grid_LI_poly_label.grid(row=5, column=0)
        grid_sep_poly_entry = tk.Entry(LI_frame, width=10)
        grid_sep_poly_entry.grid(row=5, column=1)
        grid_sep_poly_entry.insert(0, '')
        # define a button to submit these changes
        submit_poly_sep_grid = tk.Button(LI_frame, text='Submit grid', command=poly_grid_submit)
        submit_poly_sep_grid.grid(row=5, column=2)
        # update canvas
        canvas.draw()
        # get rid of circle options
        try:
            Radius_LI_label.destroy()
            Radius_LI_circ_entry.destroy()
            orient_int_btn.destroy()
        except (UnboundLocalError, NameError):  
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
    
    # Plot the circle
    circle1 = mpl.patches.Circle(cent, R, fill=False, color='red')
    main_axis.add_artist(circle1)
    fig.canvas.draw()
    circle1.remove()
    
    # update the total
    LI_total = res
    
    # Update labels
    LI_total_label.configure(text=str(round(LI_total, 4)))
    
    flux_label.configure(text=str(round(flux, 4)))

    shape_area = np.pi*R**2
    shape_area_label.configure(text=str(round(shape_area, 4)))
    
    ratio1 = LI_total/abs(shape_area)
    ratio1_label.configure(text=str(round(ratio1, 4)))
    
    ratio2 = flux/shape_area
    ratio2_label.configure(text=str(round(ratio2, 4)))
    
    return res


# define a function that will complete the line integral
def line_int_poly(N, u_str, v_str):
    global LI_total, coord_array, coord_diff, LI_verts, flux
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
        
        flux_inc = 0
        
        # Evaluate line integral as sum of vector components multiplied by the small x and y displacements 
        res = dx*np.sum(uv_store[0, :]) + dy*np.sum(uv_store[1, :])
        flux_inc = dy*np.sum(uv_store[0, :]) - dx*np.sum(uv_store[1, :])
        
        # display the drawn on lines
        canvas.draw()
        
        # update the total
        LI_total += res
        flux += flux_inc
        
        # update its label
        LI_total_label.configure(text=str(round(LI_total, 4)))
        
        if len(LI_verts) > 3:
            shape_area = calc_area(LI_verts)
            shape_area_label.configure(text=str(round(abs(shape_area), 4)))
            
            if shape_area < 0:
                shape_area = (-1)*shape_area
                flux = (-1)*flux
            else:
                pass
            
            ratio1 = LI_total/shape_area
            ratio1_label.configure(text=str(round(ratio1, 4)))
            
            ratio2 = flux/shape_area
            ratio2_label.configure(text=str(round(ratio2, 4)))
            
            flux_label.configure(text=str(round(flux, 4)))
            
            flux = 0
            
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
def PLOT_response():
    # first, take from entry boxes, wanted parameters and make them global:
    # these must be changed globally for other functions to work with the new field.
    global u, v, arrows, stacks
    # take inputs and globally update them
    update_variables()
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
    # Range and point density of the derivative plot
    d_range = 0.33*L/(zoom_slider.get())
    d_length = d_length_select.get()
    dpd = dpd_select.get()
    d_scale = scale*(zoom_slider.get())
    
    # define new axis in the derivative plot
    dx = np.linspace(-d_range+x_m, d_range+x_m, dpd)
    dy = np.linspace(-d_range+y_m, d_range+y_m, dpd)
    dxg, dyg = np.meshgrid(dx, dy)
    # define the vector field in these new axis
    uzoom, vzoom = eq_to_comps(string_x, string_y, dxg, dyg)

    # Define the components of the derivative field
    V = uzoom - eval(format_eq_div(format_eq(string_x)))
    W = vzoom - eval(format_eq_div(format_eq(string_y)))
    
    if click_opt_int > 2:
        # prepare grids to store the components
        u_div = np.zeros(shape=(dpd, dpd))
        v_div = np.zeros(shape=(dpd, dpd))
        u_curl = np.zeros(shape=(dpd, dpd))
        v_curl = np.zeros(shape=(dpd, dpd))
        
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
                
                if click_opt_int == 3:    
                    u_div[i, j] = (V_comm_1*k + W_comm_1*l)*k/A
                    v_div[i, j] = (V_comm_1*k + W_comm_1*l)*l/A
                    u_div[j, N-i] = (V_comm_2*l + W_comm_2*(-k))*l/A
                    v_div[j, N-i] = (V_comm_2*l + W_comm_2*(-k))*(-k)/A
                    u_div[N-i, N-j] = (V_comm_3*(-k) + W_comm_3*(-l))*(-k)/A
                    v_div[N-i, N-j] = (V_comm_3*(-k) + W_comm_3*(-l))*(-l)/A
                    u_div[N-j, i] = (V_comm_4*(-l) + W_comm_4*k)*(-l)/A
                    v_div[N-j, i] = (V_comm_4*(-l) + W_comm_4*k)*k/A
                    
                if click_opt_int == 4:        
                    u_curl[i, j] = (V_comm_1*l + W_comm_1*(-k))*l/A
                    v_curl[i, j] = (V_comm_1*l + W_comm_1*(-k))*(-k)/A
                    u_curl[j, N-i] = (V_comm_2*(-k) + W_comm_2*(-l))*(-k)/A
                    v_curl[j, N-i] = (V_comm_2*(-k) + W_comm_2*(-l))*(-l)/A
                    u_curl[N-i, N-j] = (V_comm_3*(-l) + W_comm_3*k)*(-l)/A
                    v_curl[N-i, N-j] = (V_comm_3*(-l) + W_comm_3*k)*k/A
                    u_curl[N-j, i] = (V_comm_4*k + W_comm_4*l)*k/A
                    v_curl[N-j, i] = (V_comm_4*k + W_comm_4*l)*l/A
        
        # correct for singular values
        for i in range(dpd):
            for j in range(dpd):
                if isnan(u_div[i, j]) is True or abs(u_div[i, j]) > 1e15 or  abs(u_div[i, j]) < 1e-10:
                    u_div[i, j] = 0
                if isnan(v_div[i, j]) is True or abs(v_div[i, j]) > 1e15 or abs(v_div[i, j]) < 1e-10:
                    v_div[i, j] = 0
    
    # Create axes at clicked position from supplied position and given axis sizes
    deriv_inset_ax = main_axis.inset_axes([(x_pix-178)/500 - (0.931*d_length/(2*L)), (y_pix-59)/500 - (0.931*d_length/(2*L)), 0.931*d_length/L, 0.931*d_length/L])
    
    # Check radiobutton selection
    if click_opt_int == 1:
        u_s = uzoom
        v_s = vzoom
        scale_s = d_scale
        
    elif click_opt_int == 2:
        u_s = V
        v_s = W
        scale_s = d_scale
        
    elif click_opt_int == 3:
        u_s = u_div
        v_s = v_div
        scale_s = d_scale
        
    elif click_opt_int == 4:
        u_s = u_curl
        v_s = v_curl
        scale_s = d_scale

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
    elif 0 < click_opt_int < 5 :
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
        # unless the canvas has not been clicked yet, in that case, make up
        # values
        if str(x_m) == 'None':
            x_m = 3.62101167
            y_m = 1.60383546
            x_pix = 592
            y_pix = 391
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
def form_2_components_plot(grid_x, grid_y, u, v, form_2_sgn, s_max, L, fract, colour_str, arrowheads=False, w_head=1/8, h_head=1/4, s_min=2):
    global s_L
    # get axis lengths:
    x_len = len(grid_x[:, 0])
    y_len = len(grid_y[:, 0])
    
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
                form_2_sgn[i, j] = 0
            if mag[i, j] == np.inf  or mag[i, j] > 1e15:
                form_2_sgn[i, j] = 1
            if mag[i, j] == -np.inf  or mag[i, j] < -1e15:
                form_2_sgn[i, j] = -1
    
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
    
    # define points of stack arrowheads as arrays for all stacks
    p_sh1x = grid_x + (s_L/2)*I_cos + (sheet_L*w_head)*I_sin
    p_sh1y = grid_y + (s_L/2)*I_sin - (sheet_L*w_head)*I_cos
    p_sh2x = grid_x + (s_L/2)*I_cos - (sheet_L*w_head)*I_sin
    p_sh2y = grid_y + (s_L/2)*I_sin + (sheet_L*w_head)*I_cos
    p_sh3x = grid_x + (s_L*0.5 + s_L*h_head)*I_cos
    p_sh3y = grid_y + (s_L*0.5 + s_L*h_head)*I_sin
    
    # define these for when there is only 1 line in the stack plot:
    P_sh1x = grid_x + (sheet_L*w_head)*I_sin
    P_sh1y = grid_y - (sheet_L*w_head)*I_cos
    P_sh2x = grid_x - (sheet_L*w_head)*I_sin
    P_sh2y = grid_y + (sheet_L*w_head)*I_cos
    P_sh3x = grid_x + (s_L*h_head)*I_cos
    P_sh3y = grid_y + (s_L*h_head)*I_sin
    
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
                    
                    main_axis.add_line(Line2D((Ax1, Bx1), (Ay1, By1), linewidth=0.7, color=colour_str[color_index]))
                    main_axis.add_line(Line2D((Ax2, Bx2), (Ay2, By2), linewidth=0.7, color=colour_str[color_index]))
                    
                    # change the parameter to loop over all changes in displacement for current magnitude
                    s += 1
            if arrowheads is True:
                # plot lines of arrowheads from central sheet for n = 1 or on top sheet for n>1 
                if n > 1:   # for all lines ubt the single sheet one
                    main_axis.add_line(Line2D((p_sh1x[i, j], p_sh3x[i, j]), (p_sh1y[i, j], p_sh3y[i, j]), linewidth=1, color=colour_str[color_index]))
                    main_axis.add_line(Line2D((p_sh2x[i, j], p_sh3x[i, j]), ((p_sh2y[i, j], p_sh3y[i, j])), linewidth=1, color=colour_str[color_index]))
                # then define it for the stacks with only 1 sheet:
                else:
                    main_axis.add_line(Line2D((P_sh1x[i, j], P_sh3x[i, j]), (P_sh1y[i, j], P_sh3y[i, j]), linewidth=1, color=colour_str[color_index]))
                    main_axis.add_line(Line2D((P_sh2x[i, j], P_sh3x[i, j]), ((P_sh2y[i, j], P_sh3y[i, j])), linewidth=1, color=colour_str[color_index]))
            else:
                pass


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

# gets 2-form from entry box and plots it as coloured blocks only
def form_2_response():
    global form_2_str, form_2_eq, form_2_sgn, form_2, comp_x, comp_y, u, v, fract
    # get globals
    update_variables()
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
    # get the signs of thsi new 2-form
    form_2_sgn = np.sign(form_2)
    # clear the current plot
    main_axis.clear()
    # use plotting stacks to display these
    #ALWAYS HALF AND HALF SPLITTING NOW
    form_2_components_plot(xg, yg, form_2/2, zero_field, form_2_sgn, s_max, L, fract, colour_str)
    form_2_components_plot(xg, yg, zero_field, form_2/2, form_2_sgn, s_max, L, fract, colour_str)
    # display the new plot
    canvas.draw()
    # display a background green on the 2-form entry to show that
    # this entry is being displayed now.
    form_2_entry.configure(bg='#C0F6BB')
    # undo it for 1-forms
    x_comp_entry.configure(bg='#FFFFFF')
    y_comp_entry.configure(bg='#FFFFFF')
    # update the label to remove the zero form if int deriv of both was used
    Label_zero_form.configure(text='')


# plots the vetor field with stacks only
def form_1_stacks_response():
    global u, v, fract
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
    Label_zero_form.configure(text='')


# performs the interior derivative on supplied 2-form and plots it as stacks
# for 2-form only, not including the 1-form.
# if combined, use different fucntion.
def Int_deriv_2_form():
    global u, v, u_str, v_str, vector_ex_str, vector_ey_str, vector_ex_eq, vector_ey_eq, vector_ex, vector_ey
    global form_2, form_2_eq, form_2_sgn
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
    # get the signs of this new 2-form
    form_2_sgn = np.sign(form_2)
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
    # clear the axis
    main_axis.clear()
    # call the asked fucntion
    Int_deriv_2_form()
    # draw its result
    canvas.draw()
    # update the label to remove the zero form if int deriv of both was used
    Label_zero_form.configure(text='')
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


# define a function that will find the interior derivative of both the 2-form
# and a 1-form, merged.
def Int_deriv_21_form():
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
    Label_zero_form.configure(text=zero_form_str)
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
    return


# Interior derivative response for 2-form only
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
    # for a 2-form only
    if type_form == 2:
        int_vector_load_btn = tk.Button(int_vector_window, text='PLOT', padx=20, pady=10, command=Int_deriv_22_form)
        int_vector_load_btn.grid(row=4, column=0, pady=10)
    # for 2-form with 1-form included:
    elif type_form == 1:
        int_vector_load_btn = tk.Button(int_vector_window, text='PLOT', padx=20, pady=10, command=Int_deriv_21_form)
        int_vector_load_btn.grid(row=4, column=0, pady=10)


# define a function that will respond to the made choice reg. int deriv.
def int_deriv_choice(var):
    if var == 0:
        # only use 2-form, therefore call previous functions for this
        Int_deriv_2_form_response(2)
    elif var == 1:
        # call the function, but make it deal with 2-form together with the
        # 1-form
        Int_deriv_2_form_response(1)
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
    int_option_window.geometry('425x200')
    # define Label
    tk.Label(int_option_window, text='Perform w.r.t 2-form only or combine given 2-form and 1-form ?').grid(row=0, column=0, columnspan=2)
    # define response buttons to the stated question
    form_2_only_btn = tk.Button(int_option_window, text='2-form', padx=30, pady=30, command=lambda: int_deriv_choice(0))
    from_2_and_1_btn = tk.Button(int_option_window, text='Both', padx=30, pady=30, command=lambda: int_deriv_choice(1))
    form_2_only_btn.grid(row=1, column=0, pady=20)
    from_2_and_1_btn.grid(row=1, column=1, pady=20)


# perform ext deriv on the result of int_deriv and plots it as stacks
def Ext_deriv_response():
    global form_2, form_2_str, form_2_sgn
    # get globals
    update_variables()
    # from these, establish the new fract, approperiate for 2-forms
    fract = 2/((pt_den-1))
    # celar current axis
    main_axis.clear()
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
    # get the signs of the 2-form
    form_2_sgn = np.sign(form_2)
    # get the string of this new 2-form to use it in int deriv
    # also put it into the entry
    form_2_str = str(simplify(str(unformat(result[0][0]))))
    # unformat it to display in the entry box, this way it does not
    # format twice if int deriv runs again
    form_2_str = unformat(form_2_str)
    form_2_entry.delete(0, 'end')
    form_2_entry.insert(0, form_2_str)
    # clear the current plot
    main_axis.clear()
    # use plotting stacks to display these
    form_2_components_plot(xg, yg, form_2/2, zero_field, form_2_sgn, s_max, L, fract, colour_str)
    form_2_components_plot(xg, yg, zero_field, form_2/2, form_2_sgn, s_max, L, fract, colour_str)
    # display the new plot
    canvas.draw()
    # display a background green on the 2-form entry to show that
    # this entry is being displayed now.
    form_2_entry.configure(bg='#C0F6BB')
    # undo it for 1-forms
    x_comp_entry.configure(bg='#FFFFFF')
    y_comp_entry.configure(bg='#FFFFFF')
    # update the label to remove the zero form if int deriv of both was used
    Label_zero_form.configure(text='')


# define a function that will wedge two 1-forms and plot them
def wedge_product_R2():
    global to_wedge_x_1_str, to_wedge_y_1_str, to_wedge_x_2_str, to_wedge_y_2_str
    global form_2_str, form_2_eq, form_2, form_2_sgn
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
    # get the signs of thsi new 2-form
    form_2_sgn = np.sign(form_2)
    # use plotting stacks to display these
    # ALWAYS HALF AND HALF SPLITTING NOW
    form_2_components_plot(xg, yg, form_2/2, zero_field, form_2_sgn, s_max, L, fract, colour_str)
    form_2_components_plot(xg, yg, zero_field, form_2/2, form_2_sgn, s_max, L, fract, colour_str)
    # display the new plot
    canvas.draw()
    # display a background green on the 2-form entry to show that
    # this entry is being displayed now.
    form_2_entry.configure(bg='#C0F6BB')
    # undo it for 1-forms
    x_comp_entry.configure(bg='#FFFFFF')
    y_comp_entry.configure(bg='#FFFFFF')
    # update the label to remove the zero form if int deriv of both was used
    Label_zero_form.configure(text='')


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
# 1-form or 2-form depending on chosen parameter  - to be implemented later
def Hodge_response():
    tk.messagebox.showwarning('Hodge error', 'This has not yet been implemented, await further updates')


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
zoom_slider = tk.Scale(right_frame, from_=1, to=200, orient=tk.HORIZONTAL)
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
component_x_entry_label = tk.Label(bot_frame, text='dx component')
component_x_entry_label.grid(row=1, column=0)
x_comp_entry = tk.Entry(bot_frame, width=20, borderwidth=2)
x_comp_entry.grid(row=2, column=0)
x_comp_entry.insert(0, 'y*sin(x)')

component_y_entry_label = tk.Label(bot_frame, text='dy component')
component_y_entry_label.grid(row=1, column=1)
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

field_select_drop_label = tk.Label(bot_frame, text='Select Pre-Defined 1-Form:')
field_select_drop_label.grid(row=3, column=0)

field_select_drop = ttk.Combobox(bot_frame, value=field_name_list, width=40)
field_select_drop.current(0)
field_select_drop.grid(row=4, column=0)
field_select_drop.bind("<<ComboboxSelected>>", field_selection_response)


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
LI_shape_instruction = tk.Label(LI_frame, text='Select what to draw:')
LI_shape_instruction.grid(row=4, column=0)
LI_shape_drop = tk.OptionMenu(LI_frame, LI_shape_select, *LI_shape_list, command=LI_shape_select_response)
LI_shape_drop.grid(row=4, column=1)


# input the radiobuttons for arrows, stacks or both again here
arrow_btn = tk.Radiobutton(LI_frame, text='arrow', variable=tensor, value=1, command=lambda: vect_type_response(tensor.get())).grid(row=7, column=1)
stack_btn = tk.Radiobutton(LI_frame, text='stack', variable=tensor, value=0, command=lambda: vect_type_response(tensor.get())).grid(row=7, column=2)



'''

DEFINE ALL WINDGETS IN CALCULUS TAB

'''

# define a window to supply the 2-form
tk.Label(calculus_frame, text='2-form on R2').grid(row=0, column=1)
form_2_entry = tk.Entry(calculus_frame, width=20, borderwidth=2)
form_2_entry.grid(row=0, column=0)
form_2_entry.insert(0, form_2_str)
# begin  displaying it on green colour to show that this is ebing displayed to
# beign with
form_2_entry.configure(bg='#C0F6BB')

# extra label for 0 form.
# for now not an entry box because it doesn't get used anywhere
# I make it red to remember to come back to it later and use it.
Label_zero_form = tk.Label(calculus_frame, text='', fg='red')
Label_zero_form.grid(row=3, column=0)

# define a button to submit the supplied 2-form and plot it as blocks
form_2_btn = tk.Button(calculus_frame, text='2-form plot', padx=3, pady=5, command=form_2_response)
form_2_btn.grid(row=1, column=0)

# define a button that will just plot the 1-form
# this will not be needed when its it merged with main GUI
# as there will already be a plot button there
form_1_stacks_btn = tk.Button(calculus_frame, text='1-form plot', padx=3, pady=5, command=form_1_stacks_response)
form_1_stacks_btn.grid(row=2, column=0)

# add a button to plot the interior derivative as superposing stack fields
INT_btn = tk.Button(calculus_frame, text='Int Deriv', padx=63, pady=10, command=Int_deriv_response)
INT_btn.grid(row=8, column=0)

# define a button to plot the exterior derivative from given u and v
# Note, it will get the 2-form first, then return back down
# to a one form to avoid concellations
# therefore, it will also just be one possible representation
# not the only possible one
EXT_int_btn = tk.Button(calculus_frame, text='Ext Deriv', padx=62, pady=10, command=Ext_deriv_response)
EXT_int_btn.grid(row=9, column=0)

# define a wedge product button that will let the user input TWO 1-forms
# in a new window to be wedged to gice a 2-form
wedge_btn = tk.Button(calculus_frame, text='wedge two 1-forms', padx=27, pady=10, command=wedge_2_response)
wedge_btn.grid(row=10, column=0)

# define ab utton that will Find the Hodge dual
Hodge_btn = tk.Button(calculus_frame, text='Hodge', padx=67, pady=10, command=Hodge_response)
Hodge_btn.grid(row=11, column=0)

'''

DEFINE WIDGETS USED IN R3 CODE

'''

height_frame = tk.LabelFrame(r3_frame, text='viewing frame', padx=32, pady=5)
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
field_input_frame = tk.LabelFrame(r3_frame, text='Fields frame', padx=32, pady=5)
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
form_2_frame = tk.LabelFrame(r3_frame, text='2-form frame', padx=32, pady=5)
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

# return time to run
stop = timeit.default_timer()
print('Time: ', stop - start)

# keep the program running until otherwise specified by user
tk.mainloop()
