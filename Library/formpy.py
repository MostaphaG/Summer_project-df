# Differential form python module attempt - 1

# import needed modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
import matplotlib as mpl
from sympy import diff, simplify
from sympy.parsing.sympy_parser import parse_expr
from math import isnan
from matplotlib import patches as patch
import matplotlib.path as mplPath
from matplotlib import animation
from scipy.integrate import odeint
import matplotlib.path as mPath

# input many numpy functions to deal with user input
from numpy import sin, cos, tan, sqrt, log, arctan, arcsin, arccos, tanh
from numpy import sinh, cosh, arcsinh, arccosh, arctanh, exp, pi, e

# define some functions that are not very important to user but are useful
# for other functions we write:

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
    '''
    Takes in 3 arguments:
    2 numbers: det number of sheets to draw, and which one is to be drawn now
    1: int bool, as 0 or 1, defines parity.
    defines coefficints needed to displace stack sheets along direction perp.
    to form, depending  on how many are to be plotted
    '''
    if c == 0:
        return ((2*s + 1)/(2*(n-1)))
    else:
        return (s/(n-1))


# %%

'''

function to create a 1-form object and define methods for it

'''


# define a function taht will set up a 1-form object that can be customised and
# plotted
def form_1(xg, yg, F_x, F_y, F_x_eqn=None, F_y_eqn=None):
    '''
    defines a 1-form object and returns it to user
    Takes 4 arguments, these are the 2 grids in 2D, which muse be square
    and of equal sizes. Then 2 arguments for the dx component and dy component
    bbased on the same grids
    '''
    # define the 1-form object and all its methods
    class form_set_up():
        # set up all variables
        def __init__(self, xg, yg, F_x, F_y, s_max=6, s_min=2, fract=0.05, deltafactor=10, fig_size=[7, 7]):
            self.xg = xg
            self.yg = yg
            self.F_x = F_x
            self.F_y = F_y
            self.figure = plt.figure(figsize=(fig_size[0], fig_size[1]))
            self.axis = self.figure.gca()
            self.s_max = s_max
            self.s_min = s_min
            self.pt_den = len(xg[:, 0])# + 1  # assume square grids
            self.fract = fract
            self.arrow_bool = False
            self.stack_bool = True
            self.orientation = 'mid'
            self.scale = 1
            self.w_head = 1/8
            self.h_head = 1/4
            self.arrowheads = True
            self.color = 'green'
            self.logarithmic_scale_bool = 0
            self.scale_bool = True
            self.delta_factor = deltafactor
            self.form_1_str_x = F_x_eqn  # to start with, use rmust change to access some methods
            # Note, the string must be given with x and y as variables
            self.form_1_str_y = F_y_eqn
        
        # #####################################################################
        # write some methods that will allow the user to chenge some of the
        # above variables
        # #####################################################################
        
        # define a mehtod to allow user to supply the string equation
        # of the 0-form
        def give_eqn(self, equation_str_x, equation_str_y):
            '''
            Takes in 1-argument, string
            This must be the equation of the supplied numerical 0-form
            It must be in terms of x and y.
            Has to be given, for some methods to be calculatable.
            '''
            self.form_1_str_x = equation_str_x
            self.form_1_str_y = equation_str_y
        
        # deifne a function to return the string equations to the user
        def return_string(self):
            '''
            Takes in no arguments, returns the unformatted strings back to user
            This is done in case user wants to access strings
            that got here not by input but by ext. alg.
            '''
            return self.form_1_str_x, self.form_1_str_x
        
        # define a method to change figure size
        def fig_size(self, n, m):
            ''' Takes two inputs, float or int numbers, sets the figure
            size to these dimensions in inches. Uses set_size_inches from
            matploitlib so can just use that on
            the atribute figure, this function is here just for
            easier nameing'''
            self.figure.set_size_inches(n, m)
        
        # change colour
        def colour(self, color):
            '''
            Takes input of a single string.String must be formatted
            as to be accepted by maplotlib colors
            changes the colour of plotted stacks.
            '''
            self.color = str(color)
        
        # change arrowsheads
        def arrow_heads(self):
            '''
            Takes no arguments.
            Changes the boolean that determines if arrowheads are plotted
            on stacks. Whenever it is called, it changes that boolean to opposite
            The form object is initialised with this as True
            '''
            self.arrowheads = not self.arrowheads
        
        # change w_head
        def head_width(self, wide):
            '''
            Takes one argument, needs to be a float or int.
            Sets the width of the arrowhead on a stacks to the desired float
            as a fraction of the stack length in the direction perp. to form
            '''
            self.w_head = float(wide)
        
        # change h_head
        def head_height(self, high):
            '''
            Takes one argument, needs to be a float or int.
            Sets the height of the arrowhead on a stacks to the desired float
            as a fraction of the stack length in the direction parall. to form
            '''
            self.h_head = float(high)
        
        # change orientation:
        def orient(self, string):
            '''
            Takes one input, needs to be a string understood by matplotlib
            quiver to orient arrows
            eg. 'tip', 'tail', 'mid' etc.
            Orients arrows on quiver plot depending on this
            '''
            self.orientation = str(string)
        
        # plot stacks boolean
        def stacks(self):
            '''
            Takes no arguments
            Changes the boolean that determines if stcaks are plotted
            Whenever it is called, it changes that boolean to opposite
            The form object is initialised with this as True
            '''
            self.stack_bool = not self.stack_bool
        
        # plot arrows boolean
        def arrows(self):
            '''
            Takes no arguments
            Changes the boolean that determines if arrows are plotted
            Whenever it is called, it changes that boolean to opposite
            The form object is initialised with this as False
            '''
            self.arrow_bool = not self.arrow_bool
        
        # change boolean that det. if to sclae logarithmically
        def log_scaling(self):
            '''
            Takes no arguments
            Changes the boolean that determines if scaling is logarithmic
            Whenever it is called, it changes that boolean to opposite
            The form object is initialised with this as False
            '''
            self.logarithmic_scale_bool = not self.logarithmic_scale_bool
        
        # define a method to be able to change bool that det. if arrows autoscale
        def autoscale(self):
            '''
            Takes no arguments
            Changes the boolean that determines if arrows are autoscaled
            Whenever it is called, it changes that boolean to opposite
            The form object is initialised with this as False
            '''
            self.scale_bool = not self.scale_bool
        
        # define methods to change s_max
        def max_sheets(self, maximum):
            '''
            Takes one argument, must be int
            Changes maximum number of sheets to draw on a stack.
            These still scale relative to max magnitude.
            '''
            self.s_max = maximum
        
        # define method to change fraction of sheetsize w.r.t graoh size:
        def sheet_size(self, fraction):
            '''
            Takes a single argument, float
            Changes the size of stack in direction perp. to form
            it is done in in terms of the fraction of graph size
            the graph size is extracted form maximum values of grid
            Grids are assumed to be square and origin centered.
            '''
            self.fract = fraction
        
        # define a method to change spare spacing around figure
        def surround_space(self, delta_denominator):
            '''
            Takes in one argument, float or int
            Sets the extra blank space around the domain of grids in axis
            The input number defines the denominator or fraction to use
            eg. supplying 3 will make the white space 1/3 of the width
            of the domain of the grid.
            '''
            self.delta_factor = delta_denominator
        
        # define a method to change the density of grids in same range
        # requires string input of 1-form:
        def same_range_density(self, points_number):
            '''
            takes in one argument, requires the string equation to be
            supplied
            Changes the desnity of points in the same range to the input value
            '''
            if self.form_1_str_x == None or self.form_1_str_y == None:
                # Error
                print('Error: You need to supply the 1-form equation to do this, look at \'give_eqn\' method')
            else:
                # redefine the grids
                v = np.linspace(-self.xg[0, -1], self.xg[0, -1], points_number)
                self.xg, self.yg = np.meshgrid(v, v)
                # based on these change other, dependant variables
                self.pt_den = len(self.xg[:, 0])
                # substitute these into the equation:
                # but keep it local
                str_x = self.form_1_str_x + ''
                str_y = self.form_1_str_y + ''
                str_x = str_x.replace('x', '(self.xg)')
                str_x = str_x.replace('y', '(self.yg)')
                str_y = str_y.replace('x', '(self.xg)')
                str_y = str_y.replace('y', '(self.yg)')
                # re-evaluate the 2-form numerically
                self.F_x = eval(str_x)
                self.F_y = eval(str_y)
        
        # #####################################################################
        # Write more complicated methods. That will use this form object
        # eg. plot, exterior derivative, Hodge etc.
        # #####################################################################
        
        # define a fucntion that will use the set up 1-form and plot it
        # stcakplot: but it takes the above defined variables:
        def plot(self, keep=True):
            '''
            Finilises the plotting
            Uses the attribues of the object as set originally and as customised
            with methods to create a plot of the 1-form
            '''
            # from self, get axis
            axis = self.axis
            
            # depending on input, clear the axis or don't
            if keep is True:
                pass
            else:
                axis.clear()
            
            # get the lengths of x and y from their grids
            x_len = len(self.xg[:, 0])
            y_len = len(self.yg[0, :])
            
            # get L from largest entry in the array, assume they are square:
            L = self.xg[0, -1]
            
            # define axis limits based on supplied arrays
            ax_L = L + L/self.delta_factor
            axis.set_xlim(-ax_L, ax_L)
            axis.set_ylim(-ax_L, ax_L)
            
            # find the distance between neightbouring points on the grid
            dist_points = self.xg[0, 1] - self.xg[0, 0]
            
            # define an empty array of magnitudes, to then fill with integer rel. mags
            R_int = np.zeros(shape=((x_len), (y_len)))
        
            # #########################################################################
            # get variables needed for the initial, simplified stack plot
            # #########################################################################
            
            # find the arrow length corresponding to each point and store in mag array
            mag = np.sqrt(self.F_x**2 + self.F_y**2)
            
            # find direction of each arrow
            angles = np.arctan2(self.F_y, self.F_x)   # theta defined from positive x axis ccw
            
            # find regions ON GRID that are nan or inf as a bool array
            #bool_array = undef_region(mag)
            
            # deal with infs and nans in mag
            for i in range(x_len):
                for j in range(y_len):
                    # set to zero points that are not defined or inf
                    if isnan(mag[i, j]) is True:
                        #colour this region as a shaded square
                        rect = patch.Rectangle((self.xg[i, j] - dist_points/2, yg[i, j]  - dist_points/2), dist_points, dist_points, color='#B5B5B5')
                        axis.add_patch(rect)
                        mag[i, j] = 0
                    if abs(mag[i, j]) == np.inf  or abs(mag[i, j]) > 1e15:
                        # colour this point as a big red dot
                        circ = patch.Circle((self.xg[i, j], self.yg[i, j]), L*self.fract/3, color='red')
                        axis.add_patch(circ)
                        mag[i, j] = 0
            
            # #########################################################################
            # use the the direction of arrows to define stack properties
            # #########################################################################
            
            # s_L and sheet_L are defined the same? 
            
            # define length of sheet as a fraction of total graph scale
            # sheet_L = L * self.fract
            # set up the max, total height of stack (along arrow)
            
            s_L = self.fract * L
            
            # #########################################################################
            # define the stacks based on geometrical arguments
            # to be perp. to arrow. shifted parallel to it, their density porp to mag
            # of the arrow and with an arrowhead on top.
            # #########################################################################
            
            # find the maximum magnitude for scaling
            max_size = np.max(mag)   # careful with singularities, else ---> nan
            
            # Define scaling factor
            #ScaleFactor = max_size/(0.9*(2*L/self.pt_den))
            
            # find the relative magnitude of vectors to maximum, as an array
            R = mag/max_size
            
            # logarithmic attempt
            if self.logarithmic_scale_bool == 1:
                log_a = 1000000
                R = np.where(R<=1/log_a, 1/log_a, R)  # Remove the values less than critical for log
                R = log(log_a*R)/log(log_a)
            else:
                pass
            
            if self.scale_bool is False:
                ScaleFactor = self.scale
            elif self.scale_bool is True:
                ScaleFactor = max_size/(0.9*(2*L/self.pt_den))
            
            # for arrows to work, with nan and infs
            # make a local variable of F_x and F_y
            # so that thye don't alter globally
            F_x_local = self.F_x * 1
            F_y_local = self.F_y * 1
            
            # plot the quiver plot on grid points if chosen in original function
            if self.arrow_bool is True:
                # prevent any magnitudes from being inf or nan
                # only here, need to do it to u and v not just mag
                for i in range(x_len):
                    for j in range(y_len):
                        if isnan(F_x_local[i,j]) == True or isnan(F_y_local[i,j]) == True or abs(F_x_local[i, j]) == np.inf or abs(F_y_local[i, j]) == np.inf or abs(F_y_local[i, j]) > 1e15 or abs(F_x_local[i, j]) > 1e15:
                            F_x_local[i,j] = F_y_local[i,j] = 0
                axis.quiver(self.xg, self.yg, F_x_local, F_y_local, pivot=self.orientation, scale=ScaleFactor, scale_units='xy') 
            else:
                pass
            
            if self.stack_bool is True:
                # define tigonometirc shifts
                I_sin = np.sin(angles)
                I_cos = np.cos(angles)
                
                # define the points that set out a line of the stack sheet (middle line)
                A_x = self.xg + (s_L/2)*I_sin
                A_y = self.yg - (s_L/2)*I_cos
                B_x = self.xg - (s_L/2)*I_sin
                B_y = self.yg + (s_L/2)*I_cos
                
                # define points of stack arrowheads as arrays for all stacks
                p_sh1x = self.xg + (s_L/2)*I_cos + (s_L*self.w_head)*I_sin
                p_sh1y = self.yg + (s_L/2)*I_sin - (s_L*self.w_head)*I_cos
                p_sh2x = self.xg + (s_L/2)*I_cos - (s_L*self.w_head)*I_sin
                p_sh2y = self.yg + (s_L/2)*I_sin + (s_L*self.w_head)*I_cos
                p_sh3x = self.xg + (s_L*0.5 + s_L*self.h_head)*I_cos
                p_sh3y = self.yg + (s_L*0.5 + s_L*self.h_head)*I_sin
                
                # define these for when there is only 1 line in the stack plot:
                P_sh1x = self.xg + (s_L*self.w_head)*I_sin
                P_sh1y = self.yg - (s_L*self.w_head)*I_cos
                P_sh2x = self.xg - (s_L*self.w_head)*I_sin
                P_sh2y = self.yg + (s_L*self.w_head)*I_cos
                P_sh3x = self.xg + (s_L*self.h_head)*I_cos
                P_sh3y = self.yg + (s_L*self.h_head)*I_sin
                
                # loop over each arrow coordinate in x and y
                for i in range(x_len):
                    for j in range(y_len):
                        
                        # Label each element with the number of stacks required: linear scaling
                        for t in range(1, self.s_max+1):
                            if (t-1)/self.s_max <= R[i, j] <= t/self.s_max:
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
                                axis.add_line(Line2D((Ax1, Bx1), (Ay1, By1), linewidth=1, color=self.color))
                                axis.add_line(Line2D((Ax2, Bx2), (Ay2, By2), linewidth=1, color=self.color))
                                
                                # update parameter to reapet and draw all needed arrows
                                s += 1
                        # deal with the odd number of stacks:
                        elif parity(n) is False:
                            # Add the centre line for odd numbers of stacks
                            axis.add_line(Line2D((A_x[i, j], B_x[i, j]), (A_y[i, j], B_y[i, j]), linewidth=1, color=self.color))
                            
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
                                axis.add_line(Line2D((Ax1,Bx1),(Ay1,By1), linewidth=1, color=self.color))
                                axis.add_line(Line2D((Ax2,Bx2),(Ay2,By2), linewidth=1, color=self.color))
                                
                                # change the parameter to loop over all changes in displacement for current magnitude
                                s += 1
                        if self.arrowheads == True:
                            # plot lines of arrowheads from central sheet for n = 1 or on top sheet for n>1 
                            if n > 1:   # for all lines ubt the single sheet one
                                axis.add_line(Line2D((p_sh1x[i, j],p_sh3x[i, j]),(p_sh1y[i, j],p_sh3y[i, j]), linewidth=1, color = self.color))
                                axis.add_line(Line2D((p_sh2x[i, j],p_sh3x[i, j]),((p_sh2y[i, j],p_sh3y[i, j])), linewidth=1, color = self.color))
                            # then define it for the stacks with only 1 sheet:
                            else:
                                axis.add_line(Line2D((P_sh1x[i, j], P_sh3x[i, j]), (P_sh1y[i, j], P_sh3y[i, j]), linewidth=1, color = self.color))
                                axis.add_line(Line2D((P_sh2x[i, j], P_sh3x[i, j]), ((P_sh2y[i, j], P_sh3y[i, j])), linewidth=1, color = self.color))
                        else:
                            pass
    
    # now call that object to create it:
    form_1_object = form_set_up(xg, yg, F_x, F_y)
    # return it to user to store
    return form_1_object


# %%

'''

function to create a 2-form object and define methods for it

'''

# define a function that will set up a 2-form object that can be customised and
# plotted
def form_2(xg, yg, form_2, form_2_eq=None):
    '''
    defines a 2-form object and returns it to user
    Takes 3 arguments basic, these are the 2 grids in 2D, which muse be square
    and of equal sizes. Then 1 argument for the dx^dy component
    based on the same grids.
    '''
    # define the 1-form object and all its methods
    class form_set_up():
        # set up all variables
        def __init__(self, xg, yg, form_2, s_max=6, s_min=2, fig_size=[7, 7], deltafactor=10):
            self.xg = xg
            self.yg = yg
            self.form_2 = form_2
            self.figure = plt.figure(figsize=(fig_size[0], fig_size[1]))
            self.axis = self.figure.gca()
            self.s_max = s_max
            self.s_min = s_min
            self.pt_den = len(xg[:, 0])  # + 1  # assume square grids
            self.fract = 2/((self.pt_den - 1))
            self.colour_list = ['red', 'blue', 'grey']
            self.logarithmic_scale_bool = 0
            self.delta_factor = deltafactor
            self.form_2_str = form_2_eq  # to start with, use rmust change to access some methods
            # Note, the string must be given with x and y as variables
        
        # #####################################################################
        # Define basic methods to customise this object
        # #####################################################################
        
        # define a mehtod to allow user to supply the string equation
        # of the 2-form
        def give_eqn(self, equation_str):
            '''
            Takes in 1-argument, string
            This must be the equation of the supplied numerical 0-form
            It must be in terms of x and y.
            Has to be given, for some methods to be calculatable.
            '''
            self.form_2_str = equation_str
        
        # define a method to change figure size
        def fig_size(self, n, m):
            '''
            Takes two inputs, float or int numbers, sets the figure
            size to these dimensions in inches. Uses set_size_inches from
            matploitlib so can just use that on
            the atribute figure, this function is here just for
            easier nameing
            '''
            self.figure.set_size_inches(n, m)
        
        # change colour list
        def colour(self, colours):
            '''
            Takes input of a single string. String must be formatted
            as to be accepted by maplotlib colors
            changes the colour of plotted stacks.
            '''
            self.colour_list = str(colours)
        
        # change boolean that det. if to sclae logarithmically
        def log_scaling(self):
            '''
            Takes no arguments
            Changes the boolean that determines if scaling is logarithmic
            Whenever it is called, it changes that boolean to opposite
            The form object is initialised with this as False (as 0)
            '''
            self.logarithmic_scale_bool = not self.logarithmic_scale_bool
        
        # define methods to change s_max
        def max_sheets(self, maximum):
            '''
            Takes one argument, must be int
            Changes maximum number of sheets to draw on a stack.
            These still scale relative to max magnitude.
            '''
            self.s_max = maximum
        
        # define method to change fraction of sheetsize w.r.t graoh size:
        def sheet_size(self, fraction):
            '''
            Takes a single argument, float
            Changes the size of stack in direction perp. to form
            it is done in in terms of the fraction of graph size
            the graph size is extracted form maximum values of grid
            Grids are assumed to be square and origin centered.
            Note, for 2-forms, the size in directions parall. and perp.
            are always set to be the same.
            '''
            self.fract = fraction
        
        #define a method to change spare spacing around figure
        def surround_space(self, delta_denominator):
            '''
            Takes in one argument, float or int
            Sets the extra blank space around the domain of grids in axis
            The input number defines the denominator or fraction to use
            eg. supplying 3 will make the white space 1/3 of the width
            of the domain of the grid.
            '''
            self.delta_factor = delta_denominator
        
        # define a method to change the density of grids in same range
        # requires string input of 1-form:
        def same_range_density(self, points_number):
            '''
            takes in one argument, requires the string equation to be
            supplied
            Changes the desnity of points in the same range to the input value
            '''
            if self.form_2_str == None:
                # Error
                print('Error: You need to supply the 2-form equation to do this, look at \'give_eqn\' method')
            else:
                # redefine the grids
                v = np.linspace(-self.xg[0, -1], self.xg[0, -1], points_number)
                self.xg, self.yg = np.meshgrid(v, v)
                # based on these change other, dependant variables
                self.pt_den = len(self.xg[:, 0])
                self.fract = 2/(self.pt_den - 1)
                # substitute these into the equation:
                # but keep it local
                str_2 = self.form_2_str + ''
                str_2 = str_2.replace('x', 'self.xg')
                str_2 = str_2.replace('y', 'self.yg')
                # re-evaluate the 2-form numerically
                self.form_2 = eval(str_2)
        
        # #####################################################################
        # Write more complicated methods. That will use this form object
        # eg. plot, exterior derivative, Hodge etc.
        # #####################################################################
        
        # define a function to plot the set up 2-form
        # originally form_2_components_plot
        def plot(self, keep=True):
            '''
            Finilises the plotting
            Takes in 1 input, that determines if axis should be cleared before
            default is True
            Uses the attribues of the object as set originally and as customised
            with methods to create a plot of the 2-form
            '''
            
            # for ease of later writting:
            axis = self.axis  # from self, get axis
            form_2 = self.form_2  # from self, get 2-form
            
            # depending on input, clear the axis or don't
            if keep is True:
                pass
            else:
                axis.clear()
            
            # get the lengths of x and y from their grids
            x_len = len(self.xg[:, 0])
            y_len = len(self.yg[0, :])
            
            # get L from largest entry in the array, assume they are square:
            L = self.xg[0, -1]
            
            # define axis limits based on supplied arrays
            ax_L = L + L/self.delta_factor
            axis.set_xlim(-ax_L, ax_L)
            axis.set_ylim(-ax_L, ax_L)
            
            # get the signs of the input 2-form
            form_2_sgn = np.sign(form_2)
            
            # define an empty array of magnitudes, to then fill with integer rel. mags
            R_int = np.zeros(shape=((x_len), (y_len)))
            
            # #########################################################################
            # get variables needed for the initial, simplified stack plot
            # #########################################################################
            
            # set up directions
            angles =[0*np.ones(np.shape(form_2)), (np.pi/2)*np.ones(np.shape(form_2))]
            
            # deal with sinularities that appear on evaluated points
            for i in range(x_len):
                for j in range(y_len):
                    # set to zero points that are not defined or inf
                    if isnan(form_2[i, j]) is True or abs(form_2[i, j]) == np.inf  or abs(form_2[i, j]) > 1e15:
                        # colour this region as a red dot, not square to
                        # not confuse with nigh mag 2-forms in stacks. or worse, in
                        # blocks
                        circ = patch.Circle((self.xg[i, j], self.yg[i, j]), L*self.fract/3, color='red')
                        axis.add_patch(circ)
                        form_2[i, j] = 0
                    # ALso, since we got this lop anyway
                    # correct for singularities in planar form 2:
                    # set to zero points that are not defined or inf
                    if isnan(form_2[i, j]) is True:
                        form_2_sgn[i, j] = 0
            
            # #########################################################################
            # use the the direction of arrows to define stack properties
            # #########################################################################
            
            # set up the max, total height of stack (along arrow)
            s_L = self.fract*L
            
            # #########################################################################
            # define the stacks based on geometrical arguments
            # to be perp. to arrow. shifted parallel to it, their density porp to mag
            # of the arrow and with an arrowhead on top.
            # #########################################################################
            # find the maximum magnitude for scaling
            max_size = abs(np.max(form_2))   # careful with singularities, else ---> nan
            
            # find the relative magnitude of vectors to maximum, as an array
            R = abs(form_2)/max_size
            
            # logarithmic attempt on 2-forms:
            if self.logarithmic_scale_bool == 1:
                log_a = 1000000
                R = np.where(R<=1/log_a, 1/log_a, R)  # Remove the values less than critical for log
                R = log(log_a*R)/log(log_a)
            else:
                pass
            
            
            # Now, for both values of theta, complete plotting:
            for theta in angles:
                # define tigonometirc shifts
                I_sin = np.sin(theta)
                I_cos = np.cos(theta)
                
                # define the points that set out a line of the stack sheet (middle line)
                A_x = self.xg + (s_L/2)*I_sin
                A_y = self.yg - (s_L/2)*I_cos
                B_x = self.xg - (s_L/2)*I_sin
                B_y = self.yg + (s_L/2)*I_cos
                
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
                        for t in range(self.s_min, self.s_max+2):
                            if (t-2)/self.s_max <= R[i, j] <= (t-1)/self.s_max:
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
                                axis.add_line(Line2D((Ax1, Bx1), (Ay1, By1), linewidth=0.5, color=self.colour_list[color_index]))
                                axis.add_line(Line2D((Ax2, Bx2), (Ay2, By2), linewidth=0.7, color=self.colour_list[color_index]))
                                
                                # update parameter to reapet and draw all needed arrows
                                s += 1
                        # deal with the odd number of stacks:
                        elif parity(n) is False:
                            # Add the centre line for odd numbers of stacks
                            axis.add_line(Line2D((A_x[i, j], B_x[i, j]), (A_y[i, j], B_y[i, j]), linewidth=0.7, color=self.colour_list[color_index]))
                            
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
                                
                                axis.add_line(Line2D((Ax1, Bx1), (Ay1, By1), linewidth=0.7, color=self.colour_list[color_index]))
                                axis.add_line(Line2D((Ax2, Bx2), (Ay2, By2), linewidth=0.7, color=self.colour_list[color_index]))
                                
                                # change the parameter to loop over all changes in displacement for current magnitude
                                s += 1
    
    # now call that object to create it:
    form_2_object = form_set_up(xg, yg, form_2)
    # return it to user to store
    return form_2_object


# %%

'''

function to create a 0-form object and define methods for it

'''


# define a function that will set up a 2-form object that can be customised and
# plotted
def form_0(xg, yg, form_0):
    '''
    defines a 0-form object and returns it to user
    Takes 3 arguments basic, these are the 2 grids in 2D, which muse be square
    and of equal sizes. Then 1 argument 0-form based on the same grids.
    '''
    # define the 1-form object and all its methods
    class form_set_up():
        # set up all initial, defualt variables
        def __init__(self, xg, yg, form_0, fig_size=[7, 7], deltafactor=10):
            self.xg = xg
            self.yg = yg
            self.form_0 = form_0
            self.figure = plt.figure(figsize=(fig_size[0], fig_size[1]))
            self.axis = self.figure.gca()
            self.pt_den = len(xg[:, 0])  # + 1  # assume square grids
            self.logarithmic_scale_bool = 0
            self.delta_factor = deltafactor
            self.denser = 1
            self.lines = 15
            self.fontsize = 7
            self.inline_bool = True
            self.form_0_str = None  # to start with, use rmust change to access some methods
            # Note, the string must be given with x and y as variables
            self.form_0_contour = None  # Initialise with that, will be changed, if user
            # gets contour plot with new density.
        
        # #####################################################################
        # Define basic methods to customise this object
        # #####################################################################
        
        # define a mehtod to allow user to supply the string equation
        # of the 0-form
        def give_eqn(self, equation_str):
            '''
            Takes in 1-argument, string
            This must be the equation of the supplied numerical 0-form
            It must be in terms of x and y.
            Has to be given, for some methods to be calculatable.
            '''
            self.form_0_str = equation_str
        
        # define a method to change figure size
        def fig_size(self, n, m):
            '''
            Takes two inputs, float or int numbers, sets the figure
            size to these dimensions in inches. Uses set_size_inches from
            matploitlib so can just use that on
            the atribute figure, this function is here just for
            easier nameing
            '''
            self.figure.set_size_inches(n, m)
        
        #define a method to change spare spacing around figure
        def surround_space(self, delta_denominator):
            '''
            Takes in one argument, float or int
            Sets the extra blank space around the domain of grids in axis
            The input number defines the denominator or fraction to use
            eg. supplying 3 will make the white space 1/3 of the width
            of the domain of the grid.
            '''
            self.delta_factor = delta_denominator
        
        
        def density_increase(self, factor):
            '''
            Takes 1 float/int argument
            sets increase in density between form grids and contour grids
            needed if this was accessed by other forms via ext.alg methods
            Note, This cannot be set to anything but 1, if the 0-form
            equation as string is not also supplied correctly.
            '''
            self.denser = factor
        
        def lines_number(self, number):
            '''
            Takes 1 int argument
            changes number of contour lines that get drawn
            supplied to contour plot from matplotlib via levels
            '''
            self.lines = number
        
        def fonts_size(self, size):
            '''
            Takes 1 float/int argument
            Changes fontsize for contour labels
            '''
            self.fontsize = size
        
        def labels(self):
            '''
            Takes no arguments
            determines if hight labels are put on contours
            Starts off as True, calling changes it each time
            '''
            self.inline_bool = not self.inline_bool
        
        # define a method to change the density of grids in same range
        # requires string input of 1-form:
        # technically not so needed here as all plotting is done with denser
        # which is similar enough. But it might be useful to change
        # the grids as stored and not just locally for a plot
        def same_range_density(self, points_number):
            '''
            takes in one argument, requires the string equation to be
            supplied
            Changes the desnity of points in the same range to the input value
            '''
            if self.form_0_str == None:
                # Error
                print('Error: You need to supply the 0-form equation to do this, look at \'give_eqn\' method')
            else:
                # redefine the grids
                v = np.linspace(-self.xg[0, -1], self.xg[0, -1], points_number)
                self.xg, self.yg = np.meshgrid(v, v)
                # based on these change other, dependant variables
                self.pt_den = len(self.xg[:, 0])
                # substitute these into the equation:
                # but keep it local
                str_0 = self.form_0_str + ''
                str_0 = str_0.replace('x', 'self.xg')
                str_0 = str_0.replace('y', 'self.yg')
                # re-evaluate the 2-form numerically
                self.form_0 = eval(str_0)
        
        # #####################################################################
        # Write more complicated methods. That will use this form object
        # eg. plot, exterior derivative, Hodge etc.
        # #####################################################################
        
        
        # define a fucntion to plot a zero form when button is pressed.
        def plot(self, keep=True):
            '''
            Finilises the plotting
            Uses the attribues of the object as set originally and as customised
            with methods to create a plot of the 2-form.
            Can take one parameter: bool. Default is True
            determines if axis should be cleared before plotting.
            '''
            
            # check if user wants to clear first:
            if keep is True:
                pass
            else:
                self.axis.clear()
            
            # for ease of later writting:
            axis = self.axis  # from self, get axis
            form_0 = self.form_0  # from self, get 2-form
            
            # get L from largest entry in the array, assume they are square:
            L = self.xg[0, -1]
            
            # define axis limits based on supplied arrays
            ax_L = L + L/self.delta_factor
            axis.set_xlim(-ax_L, ax_L)
            axis.set_ylim(-ax_L, ax_L)
            
            if self.denser != 1:
                if self.form_0_str == None:
                    # This cannot be done if a string has not been supplied
                    # ERROR
                    print('Error: You need to supply the 0-form equation to do this, look at \'give_eqn\' method')
                else:
                    # get the supplied form as a string
                    zero_form_str = str(simplify(self.form_0_str))
                    # set up grids for contours
                    contour_x, contour_y = np.linspace(-L, L, self.pt_den*self.denser), np.linspace(-L, L, self.pt_den*self.denser)
                    contour_x_grid, contour_y_grid = np.meshgrid(contour_x, contour_y)
                    # format the given ftring
                    zero_form_str = zero_form_str.replace('x', 'contour_x_grid')
                    zero_form_str = zero_form_str.replace('y', 'contour_y_grid')
                    # evaluate bearing in mind zeros
                    if zero_form_str.find('contour_x_grid') & zero_form_str.find('contour_y_grid') == -1:
                        form_0_contour = eval(zero_form_str)*np.ones(np.shape(contour_x_grid))
                    else:
                        form_0_contour = eval(zero_form_str)
                    # set up the contour plot
                    CS = axis.contour(contour_x_grid, contour_y_grid, form_0_contour, levels=self.lines)
                    axis.clabel(CS, inline=self.inline_bool, fontsize=self.fontsize)
            else:
                # set up the contour plot with given grids
                CS = axis.contour(self.xg, self.yg, form_0, levels=self.lines)
                axis.clabel(CS, inline=self.inline_bool, fontsize=self.fontsize)
        
        # define a method to compute the exterior derivative
        def ext_d(self):
            '''
            Takes in no argument
            Returns 1 form object and the strings of its x comp. and y-comp.
            computes the exterior derivative and returns it as the 1-form object
            '''
            
            # first make sure that the string has been supplied
            if self.form_0_str == None:
                    # ERROR
                    print('Error: You need to supply the 0-form equation to do this, look at \'give_eqn\' method')
            else:
                # can compute the exterior derivative:
                form_0_str = str(simplify(self.form_0_str))
                # from this, need derivatives so set it as a SymPy object
                sympy_expr_form_0 = parse_expr(form_0_str, evaluate=False)
                # set up an array of coordinates that need to be used (in standard order)
                coords = ['x', 'y']
                # from these, find the derivatives
                form_1_x_str = str(diff(sympy_expr_form_0, coords[0]))
                form_1_y_str = str(diff(sympy_expr_form_0, coords[1]))
                # need to uspply these unformatted, so save those:
                form_1_x_unformated, form_1_y_unformated = form_1_x_str*1, form_1_y_str*1
                # from these strings, get the numerical 1-form:
                form_1_x_str = form_1_x_str.replace('x', 'self.xg')
                form_1_x_str = form_1_x_str.replace('y', 'self.yg')
                form_1_y_str = form_1_y_str.replace('x', 'self.xg')
                form_1_y_str = form_1_y_str.replace('y', 'self.yg')
                form_0_x = eval(form_1_x_str)
                form_0_y = eval(form_1_y_str)
                # supply these to the 1-form object function
                result_1_form = form_1(self.xg, self.yg, form_0_x, form_0_y, form_1_x_unformated, form_1_y_unformated)
                return result_1_form
    
    # now call that object to create it:
    form_0_object = form_set_up(xg, yg, form_0)
    # return it to user to store
    return form_0_object
