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


# define a function that will find the 2-form from given expressions
# in a given number of dimensions and in terms of given coordinate symbols
def find_2_form(expressions, coords, xg, yg, zg=None, m=2):
    
    # from the grids, find pt_den
    # again, assume that they are square
    pt_den = len(xg[:, 0])
    
    # define a sympy expression for string 0
    sympy_expr_zero = parse_expr('0*x', evaluate=False)
    
    # set up an array to store derrivatives.
    ext_ds = np.empty((m, m), dtype='object')
    
    # set up an array to store the results
    # in 2D only dx^dy, in 3D (m=3) (in order): dx^dy, dx^dz, dy^dz
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
    
    
    # merge the results into a 2-form (for 2-form on R^2, the result is a single component (dx^xy))
    # do so by adding opposite elements along the diagonal ( / ) components of ext_ds
    # this  includes taking elemets with switched i and j
    
    
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
    
    # create a local result, that will be used to evaluate the resulting string
    loc_res = result + ''
    
    # format string in each result row
    for d in range(pair):
        # format the result to be 'python understood' to be able to use the eval()
        loc_res[d, 0] = loc_res[d, 0].replace('x', 'xg')
        loc_res[d, 0] = loc_res[d, 0].replace('y', 'yg')
        loc_res[d, 0] = loc_res[d, 0].replace('z', 'zg')
        
        # check against constant result, to be of correct shape before eval
        if loc_res[d, 0].find('x') & loc_res[d, 0].find('y') == -1:
            loc_res[d, 0] = '(' + str(loc_res[d, 0]) + ')* np.ones(np.shape(xg))'
        if loc_res[d, 0].find('x') & loc_res[d, 0].find('y') == -1:
            loc_res[d, 0] = '(' + str(loc_res[d, 0]) + ')* np.ones(np.shape(yg))'

    # set up a vector to store the 2-form numerically, from xg and yg
    # Note - need pt_den m times.
    if m == 2:
        form_2 = np.empty((1, pt_den, pt_den))
        form_2[0, :, :] = eval(loc_res[0, 0])
    elif m == 3:
        form_2 = np.empty((3, pt_den, pt_den, pt_den))
        for d in range(3):
            form_2[d, :, :, :] = eval(loc_res[d, 0])
    
    # return useful findings to the user
    return form_2, result, ext_ds


# %%

'''

function to create a 1-form object and define methods for it

'''


# define a function taht will set up a 1-form object that can be customised and
# plotted
def form_1(xg, yg, F_x, F_y, F_x_eqn=None, F_y_eqn=None, fig=None, subplots=False, sub_axis_list=[]):
    '''
    defines a 1-form object and returns it to user
    Takes 9 arguments, these are the 2 grids in 2D, which muse be square
    and of equal sizes. Then 2 arguments for the dx component and dy component
    based on the same grids
    Then 2 equations, for x and y (not always needed)
    Then, can supply a figure if user doesn't wat this object
    to create a new one for itself
    Can also supply a bool input to define if subplots are to be allowed
    these can be added using a method (add_subplot)
    also, user can supply sub_axis_list, to provide the axis they have set
    these only work if a figure has been supplied and if subplots is True
    '''
    # define the 1-form object and all its methods
    class form_set_up():
        # set up all variables
        def __init__(self, xg, yg, F_x, F_y, s_max=6, s_min=2, fract=0.05, deltafactor=10, fig_size=[7, 7]):
            self.xg = xg
            self.yg = yg
            self.F_x = F_x
            self.F_y = F_y
            # set up a figure if one was not given:
            if fig == None:
                self.figure = plt.figure(figsize=(fig_size[0], fig_size[1]))
                self.axis = self.figure.gca()
            else:
                if subplots is False:
                    self.figure = fig
                    self.axis  = self.figure.gca()
                elif subplots is True:
                    if len(sub_axis_list) == 0:
                        self.figure = fig
                        self.axis = []
                    else:
                        self.figure = fig
                        self.axis = sub_axis_list
                else:
                    raise TypeError('Error, incorrect input for \'subplots\'')
            self.subplots = subplots
            self.s_max = s_max
            self.s_min = s_min
            self.pt_den = len(xg[:, 0])# + 1  # assume square grids
            self.fract = fract
            self.orientation = 'mid'
            self.scale = 1
            self.w_head = 1/8
            self.h_head = 1/4
            self.arrowheads = True
            self.color = 'green'
            self.logarithmic_scale_bool = 0
            self.scale_bool = True
            self.delta_factor = deltafactor
            if F_x_eqn is not None:
                self.form_1_str_x = str(simplify(F_x_eqn))  # to start with, use rmust change to access some methods
                # Note, the string must be given with x and y as variables
            else:
                self.form_1_str_x = None
            
            if F_y_eqn is not None:
                self.form_1_str_y = str(simplify(F_y_eqn))
            else:
                self.form_1_str_x = None
            
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
            # set equation parameters to simplified inputs
            self.form_1_str_x = str(simplify(equation_str_x))
            self.form_1_str_y = str(simplify(equation_str_y))
            # make the values match automatically to limit how often mismatch occurs
            # substitute these into the equation:
            # but keep it local
            str_x = self.form_1_str_x + ''
            str_y = self.form_1_str_y + ''
            str_x = str_x.replace('x', '(self.xg)')
            str_x = str_x.replace('y', '(self.yg)')
            str_y = str_y.replace('x', '(self.xg)')
            str_y = str_y.replace('y', '(self.yg)')
            
            # check against constant forms, to have correct shape
            if str_x.find('x') & str_x.find('y') == -1:
                str_x = '(' + str(str_x) + ')* np.ones(np.shape(self.xg))'
            if str_y.find('x') & str_y.find('y') == -1:
                str_y = '(' + str(str_y) + ')* np.ones(np.shape(self.yg))'
            # re-evaluate the 2-form numerically
            self.F_x = eval(str_x)
            self.F_y = eval(str_y)
        
        # deifne a function to return the string equations to the user
        def return_string(self):
            '''
            Takes in no arguments, returns the unformatted strings back to user
            This is done in case user wants to access strings
            that got here not by input but by ext. alg.
            '''
            return self.form_1_str_x, self.form_1_str_y
        
        # define a method to add a subplot
        def add_subplot(self, order):
            '''
            Takes in one argument, as as matplotlib add_subplot
            It determines the shape of the subplot structure and which axis is
            being set (eg. 231 gives a set of plots 2 rows by 3 columns
            and this sets the first axis in that configuration)
            Adds a subplot to the figure that this form occupies
            Returns nothing, saves the subplot axis to the axis self param.
            '''
            # chck if the correct option was chosen to allow for this:
            if type(self.axis) != list:
                raise ValueError('Error, set up the object allowing for subplots')
            else:
                sub_axis = self.figure.add_subplot(order)
                self.axis.append(sub_axis)
        
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
                raise ValueError('Error: You need to supply the 1-form equation to do this, look at \'give_eqn\' method')
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
                
                # cehck against constant forms, to have correct shape
                if str_x.find('x') & str_x.find('y') == -1:
                    str_x = '(' + str(str_x) + ')* np.ones(np.shape(self.xg))'
                if str_y.find('x') & str_y.find('y') == -1:
                    str_y = '(' + str(str_y) + ')* np.ones(np.shape(self.yg))'
            
                # re-evaluate the 2-form numerically
                self.F_x = eval(str_x)
                self.F_y = eval(str_y)
        
        # #####################################################################
        # Write more complicated methods. That will use this form object
        # eg. plot, exterior derivative, Hodge etc.
        # #####################################################################
        
        # define a fucntion that will use the set up 1-form and plot it
        # stcakplot: but it takes the above defined variables:
        def plot(self, keep=True, subplot_index=None):
            '''
            Finilises the plotting
            Uses the attribues of the object as set originally and as customised
            with methods to create a plot of the 1-form
            Takes in 2 arguments:
            --- \'keep\', which allows the user to plot
            on top of the previous pltos they created without clearing axis
            When its set to False, the axis are cleared first
            Default is True.
            --- \' subplot_index \', default set to None, can be input if
            the user has selected subplots to be allowed when creating the
            object. Determines which aixs to draw on, indecies are in order
            that they were added to the object
            '''
            # from self, get axis
            # depending on if subplots are wanted:
            if type(self.axis) != list:
                axis = self.axis
            else:
                axis = self.axis[subplot_index]
            
            # depending on input, clear the axis or don't
            if keep is True:
                pass
            else:
                axis.clear()
            
            # get the lengths of x and y from their grids
            x_len = len(self.xg[:, 0])
            y_len = len(self.yg[0, :])
            
            # Extract L from the x and y grids. Assumes they are square.
            L = 0.5*(self.xg[0, -1] - self.xg[0, 0])
            x0 = self.xg[0,0] + L
            y0 = self.yg[0,0] + L
            
            ax_L = L + L/self.delta_factor
            axis.set_xlim(-ax_L + x0, ax_L + x0)
            axis.set_ylim(-ax_L + y0, ax_L + y0)
            
            # find the distance between neightbouring points on the grid
            dist_points = self.xg[0, 1] - self.xg[0, 0]
            
            # define an empty array of magnitudes, to then fill with integer rel. mags
            R_int = np.zeros(shape=((x_len), (y_len)))
        
            # #########################################################################
            # get variables needed for the initial, simplified stack plot
            # #########################################################################
            
            # set all insignificant values to zero:
            self.F_x[np.abs(self.F_x) < 1e-15] = 0
            self.F_x[np.abs(self.F_x) < 1e-15] = 0
            
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
        
        # define a method to find its exterior derivative
        def ext_d(self, pass_on_figure=False):
            '''
            Takes in 1 argument:
            determines if the figure should be passed onto the 2-form
               object that is returned, or if it should create one of its own
            Returns 2 form object.
            Computes the exterior derivative and returns it
            as the 2-form object
            '''
            if self.form_1_str_x == None or self.form_1_str_y == None:
                    # ERROR
                    raise ValueError('Error: You need to supply the 1-form equations to do this, look at \'give_eqn\' method')
            else:
                # the strings have been correctly given, compute the
                # exterior derivative
                # get the inpus from fields of x and u components
                x_comp_str = self.form_1_str_x
                y_comp_str = self.form_1_str_y
                # from found u and v in the interior derivative, set up sympy components
                sympy_expr_x = parse_expr(x_comp_str, evaluate=False)
                sympy_expr_y = parse_expr(y_comp_str, evaluate=False)
                # combine the 2 into a list:
                expressions = np.array([sympy_expr_x, sympy_expr_y])
                # set up an array of coordinates that need to be used (in standard order)
                coords = ['x', 'y']
                # set up dimensionality
                m = 2
                
                
                # ################## from these get the 2-form ################
                
                
                # define a sympy expression for string 0
                sympy_expr_zero = parse_expr('0*x', evaluate=False)
                
                # set up an array to store derrivatives.
                ext_ds = np.empty((m, m), dtype='object')
                
                # set up an array to store the results
                # in 2D only dx^dy, in 3D (m=3) (in order): dx^dy, dx^dz, dy^dz
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
                
                
                # merge the results into a 2-form (for 2-form on R^2, the result is a single component (dx^xy))
                # do so by adding opposite elements along the diagonal ( / ) components of ext_ds
                # this  includes taking elemets with switched i and j
                
                
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
                
                
                # ################### done dinfding 2-form ####################
                
                # get the string of this new 2-form
                form_2_str = str(simplify(result[0][0]))
                
                # keep a local, unformatted version of this
                # to supply to form_2
                form_2_str_loc = form_2_str*1
                
                # numerically evaluate it, careful about constants
                # to evaluate it, make sure to use grids
                form_2_str = form_2_str.replace('x', '(self.xg)')
                form_2_str = form_2_str.replace('y', '(self.yg)')
                if form_2_str.find('x') & form_2_str.find('y') == -1:
                    form_2_str = '(' + str(form_2_str) + ')* np.ones(np.shape(self.xg))'
                else:
                    pass
                
                form_2_result = eval(form_2_str)
                
                if pass_on_figure is False:
                    # supply these to the 2-form object creator
                    result_form = form_2(self.xg, self.yg, form_2_result, form_2_str_loc)
                elif pass_on_figure is True:
                    if self.subplots is False:
                        result_form = form_2(self.xg, self.yg, form_2_result, form_2_str_loc, fig=self.figure, subplots=False)
                    elif self.subplots is True:
                        result_form = form_2(self.xg, self.yg, form_2_result, form_2_str_loc, fig=self.figure, subplots=True, sub_axis_list=self.axis)
                else:
                    raise ValueError('Error, Incorrect input for \' pass_on_figure \'')
                
                # return it to the user
                return result_form
        
        # define a funciton to complete numerical only curl
        def num_ext_d(self, pass_on_figure=False):
            '''
            Takes in 1 argument:
            --- pass on figure: Determies if figure should be passed onto the
               new object if it is to be created
             
            returns 2-form object
            computes the exterior derivative numerically only
            The equations do not need to be given
            If given, they do not get passed onto the 2-form object anyway
            NUMERICAL ONLY, they will be lost!
            '''
            
            # get steps in dx and dy:
            dx = self.xg[0, :]
            dy = self.yg[:, 0]
            
            # copy F_x and F_y, locally
            fx = self.F_x + np.zeros(np.shape(self.xg))
            fy = self.F_y + np.zeros(np.shape(self.xg))
            
            #np.set_printoptions(True)
            
            # clean up F_x and F_y from nan
            # keep inf and large values, for gradient to be found still
            for i in range(len(self.xg[:, 0])):
                for j in range(len(self.yg[0, :])):
                    # correct for ill defined values
                    if isnan(fx[i, j]):
                        fx[i, j] = 0
                    if isnan(fy[i, j]):
                        fy[i, j] = 0
                    if abs(fx[i, j]) == np.inf  or abs(fx[i, j]) > 1e15:
                        fx[i, j] = 1e10
                    if abs(fy[i, j]) == np.inf  or abs(fy[i, j]) > 1e15:
                        fy[i, j] = 1e10
            
            
            # Calculate deirvatvies as needed, using numpy gradient.
            dy_F_x, _ = np.gradient(fx, dx, dy)
            _, dx_F_y = np.gradient(fy, dx, dy)
            
#            dy_F_x, _ = np.gradient(0.5*(fx + fx.transpose()), dx, dy)
#            _, dx_F_y = np.gradient(0.5*(fy + fy.transpose()), dx, dy)
            
#            dy_F_x, _ = np.gradient(fx - fy, dx, dy)
#            _, dx_F_y = np.gradient(fy - fx, dx, dy)
            
            # from these, get the 2-form
            form_2_result = dx_F_y - dy_F_x
            #form_2_result = 0.5*(form_2_result + form_2_result.transpose())
            
            # return 2-form object to user
            if pass_on_figure is False:
                # supply these to the 2-form object creator
                result_form = form_2(self.xg, self.yg, form_2_result)
            elif pass_on_figure is True:
                if self.subplots is False:
                    result_form = form_2(self.xg, self.yg, form_2_result, fig=self.figure, subplots=False)
                elif self.subplots is True:
                    result_form = form_2(self.xg, self.yg, form_2_result, fig=self.figure, subplots=True, sub_axis_list=self.axis)
            else:
                raise ValueError('Error, Incorrect input for \' pass_on_figure \'')
                
            # return it to the user
            return result_form
        
        
        # define a method to Hodge it
        def Hodge(self, numerical_only=True, keep_object=True, pass_on_figure=False):
            '''
            Takes in three bool arguments:
            1) determines if the calculation should be numerical or analytic
            True if numerical, False if analytic and numerical. Default is True
            For the analytic one, the equations as string must be supplied to
            the object
            Note, choosing analytical, changes the equations AND the numerical answers
            2) determines if the result should be returned as a new 1-form or
                if the current one need to be changed. Default is True
            3) Determies if figure should be passed onto the new object
                if it is to be created
            It calulates the Hodge on R^2 by the standard definition:
            dx -> dy and dy -> -dx
            return nothing 1-form if keep_object is False, else returns nothing
            '''
            # distinguish between doing it numerically and alaytically
            if numerical_only is True:
                # check if equations have been given:
                # if they have, doing it only numerically would create
                # a mismatch, avoid that
                if self.form_1_str_x == None or self.form_1_str_y == None:
                    pass
                else:
                    # equations have been given, a mismatch may occur
                    # warn the user
                    print('Warning: You supplied equations, doing it numerically only will result in a mismacth between numerical values and equations')
                # now complete the process numerically save as instructed
                # check keep_object:
                if keep_object is True:
                    # change the object self properties accoringly
                    new_x = -self.F_y
                    new_y = self.F_x
                    self.F_x = new_x
                    self.F_y = new_y
                elif keep_object is False:
                    # pass these in to the object to create a new one:
                    # DEPENDING ON FIGURES AND SUBPLOTS:
                    # N.B no equations to supply
                    if pass_on_figure is False:
                        # supply these to the 2-form object creator
                        new_object = form_1(self.xg, self.yg, -self.F_y, self.F_x)
                    elif pass_on_figure is True:
                        if self.subplots is False:
                            new_object = form_1(self.xg, self.yg, -self.F_y, self.F_x, fig=self.figure, subplots=False)
                        elif self.subplots is True:
                            new_object = form_1(self.xg, self.yg, -self.F_y, self.F_x, fig=self.figure, subplots=True, sub_axis_list=self.axis)
                    else:
                        raise ValueError('Error, Incorrect input for \' pass_on_figure \'')
                
                    # return the new one to the user:
                    return new_object
                else:
                    raise ValueError('Error, Invalid input for \'keep_object\'')
            
            elif numerical_only is False:
                # can only be done if equations have been given, check:
                if self.form_1_str_x == None or self.form_1_str_y == None:
                    # ERROR
                    raise TypeError('Error: You need to supply the 1-form equation to do this, look at \'give_eqn\' method')
                else:
                    # some equations are there, compute the Hodge on these:
                    # Note: Upto user to make sure their equations match their
                    # numerical input
                    new_str_x = '-(' + self.form_1_str_y + ')'
                    new_str_y = self.form_1_str_x
                    # from these, get numerical solutions, evaulated on local
                    # strings changed to relate to the self grids
                    # need to uspply these unformatted, so save those:
                    form_1_x_unformated, form_1_y_unformated = new_str_x*1, new_str_y*1
                    # from these strings, get the numerical 1-form:
                    new_str_x = new_str_x.replace('x', '(self.xg)')
                    new_str_x = new_str_x.replace('y', '(self.yg)')
                    new_str_y = new_str_y.replace('x', '(self.xg)')
                    new_str_y = new_str_y.replace('y', '(self.yg)')
                    
                    if new_str_x.find('x') & new_str_x.find('y') == -1:
                        new_str_x = '(' + str(new_str_x) + ')* np.ones(np.shape(self.xg))'
                    if new_str_y.find('x') & new_str_y.find('y') == -1:
                        new_str_y = '(' + str(new_str_y) + ')* np.ones(np.shape(self.yg))'
                    
                    form_1_x = eval(new_str_x)
                    form_1_y = eval(new_str_y)
                    
                    # depending on keep_object, return:
                    if keep_object is True:
                        self.F_x = form_1_x
                        self.F_y = form_1_y
                        self.form_1_str_x = form_1_x_unformated
                        self.form_1_str_y = form_1_y_unformated
                    elif keep_object is False:
                        if pass_on_figure is True:
                            # pass these in to the object to create a new one:
                            if self.subplots is False:
                                new_object = form_1(self.xg, self.yg, form_1_x, form_1_y, F_x_eqn=form_1_x_unformated, F_y_eqn=form_1_y_unformated, fig=self.figure, subplots=False)
                            elif self.subplots is True:
                                new_object = form_1(self.xg, self.yg, form_1_x, form_1_y, F_x_eqn=form_1_x_unformated, F_y_eqn=form_1_y_unformated, fig=self.figure, subplots=True, sub_axis_list=self.axis)
                        elif pass_on_figure is False:
                            new_object = form_1(self.xg, self.yg, -self.F_y, self.F_x)
                        else:
                            raise ValueError('Error, Incorrect input for \' pass_on_figure \'')
                        # return the new one to the user:
                        return new_object
                    else:
                        raise ValueError('Error, Invalid input for \'keep_object\'')
            else:
                # Error
                raise ValueError('ERROR: Invalid input for \'numerical_only\'')
            
        
        # define a fucntion to compute a wedge product of two 1 forms
        def wedge_analytical(self, equation_2_x, equation_2_y, pass_on_figure=False):
            '''
            Takes in 3 arguments, 2 are for the strings of the other form to
            wedge with this one.
            Third determines if figure from this object is to be passed on
            to the 2-form object
            To do so here, strings for the form must be supplied.
            Computes the Wedge product using strings, ANALYTICALLY
            Returns a 2-form object
            '''
            # first, get all entries out, save as string for these to display when
            # window is opened again
            to_wedge_x_1_str = self.form_1_str_x
            to_wedge_y_1_str = self.form_1_str_y
            to_wedge_x_2_str = equation_2_x
            to_wedge_y_2_str = equation_2_y
            # first, find the result of the 2-form
            # this if, in terms of the above commented fields:
            # 2-form = f*m - g*h
            # get it mathematically, as a string
            form_2_str = str(simplify( '(' + to_wedge_x_1_str + ')*(' +  to_wedge_y_2_str + ')' + ' - (' + to_wedge_y_1_str + ')*(' +  to_wedge_x_2_str + ')' ))
            # keep it as it is locally to supply it to object maker later
            form_2_str_loc = form_2_str + ''
            # format it to be in terms of grids and:
            # check against constant and zero 2-forms being supplied
            # get the numerical evaluation of it
            form_2_str = form_2_str.replace('x', 'self.xg')
            form_2_str = form_2_str.replace('y', 'self.yg')
            if form_2_str.find('x') & form_2_str.find('y') == -1:
                form_2_str = '(' + str(form_2_str) + ')* np.ones(np.shape(self.xg))'
            else:
                pass
            
            # evaluate it numerically on the grid supplied
            form_2_result = eval(form_2_str)
            # create a 2-form object from this; to return.
            # do this depending on what return figures was chosen
            if pass_on_figure is False:
                ret_object = form_2(self.xg, self.yg, form_2_result, form_2_str_loc)
            elif pass_on_figure is True:
                if self.subplots is False:
                    ret_object = form_2(self.xg, self.yg, form_2_result, form_2_str_loc, fig=self.figure, subplots=False)
                elif self.subplots is True:
                    ret_object = form_2(self.xg, self.yg, form_2_result, form_2_str_loc, fig=self.figure, subplots=True, sub_axis_list=self.axis)
            else:
                raise ValueError('Error, Incorrect input for \' pass_on_figure \'')
            
            # return the resulting object
            return ret_object
        
        
        # define a method for numerical wedge product
        def wedge_num(self, form_1_second=None, pass_on_figure=False):
            '''
            Takes in 2 arguments:
            The first one is the form_1_second, either as:
                --- 1-form object
                --- numerical components
            If none supplied, assumed to be a (1, 1) 1-form
            Third determines if figure from this object is to be passed on
            to the 2-form object
            Computes the Wedge product numerically
            Returns a 2-form object
            '''
            
            # test if equations were given first:
            if self.form_1_str_x == None or self.form_1_str_y == None:
                pass
            else:
                # Warn user that these will be lost
                print('The first 1-form you are completing the wedge with has equations supplied, these will be lost')
            
            # if the vector field was supplied, extract its equations, if possible
            if form_1_second is None:
                # if none was given, do it with respect to uniform 1, 1
                f12_x = np.ones(np.shape(self.xg))
                f12_y = np.ones(np.shape(self.xg))
            elif type(form_1_second) == tuple:
                # if numerical grids were given, take these, is equations were given here, break!
                if type(form_1_second[0]) == str:
                    raise ValueError('for numerical calulation, supply 1-form arrays, not equations')
                else:
                    f12_x = form_1_second[0]
                    f12_y = form_1_second[1]
            else:
                # object supplied, get numericals
                f12_x = form_1_second.F_x
                f12_y = form_1_second.F_y
                
                # warn user if equations were given too:
                if form_1_second.form_1_str_x == None or form_1_second.form_1_str_y == None:
                    pass
                else:
                    # Warn user that these will be lost
                    print('The second 1-form in wedge has equations supplied, these will be lost')
            
            # from these get the numerical 2-form
            result = self.F_x * f12_y - self.F_y * f12_x
            
            # return it to user:
            if pass_on_figure is False:
                # supply these to the 0-form object creator
                result_form = form_2(self.xg, self.yg, result)
            elif pass_on_figure is True:
                if self.subplots is False:
                    result_form = form_2(self.xg, self.yg, result, fig=self.figure, subplots=False)
                elif self.subplots is True:
                    result_form = form_2(self.xg, self.yg, result, fig=self.figure, subplots=True, sub_axis_list=self.axis)
            else:
                raise ValueError('Error, Incorrect input for \' pass_on_figure \'')
            
            # return it to the user
            return result_form
            
        
        def zoom(self, target=[0,0], zoom=2, dpd=9):
            '''
            Create a new window which displays the field zoomed at a certain point
            User gives arguments
            Target: Determines the zoom location, coordinates
            Zoom: +ve float, determines zooming amount
            dpd: +int, determines how many points on each axis
            '''
            
            # Requires user to provide eqn of the 1-form they are zooming on.
            
            if self.form_1_str_x == None or self.form_1_str_y == None:
                # ERROR
                raise TypeError('Error: No equation provided')
            else:
                
                # Target coordinates
                x_m = target[0]
                y_m = target[1]
                
                d_range = self.xg[0, -1]/zoom
                #d_length = 
                
                # Set up zoom window grids
                dx = np.linspace(-d_range + x_m, d_range + x_m, dpd)
                dy = np.linspace(-d_range + y_m, d_range + y_m, dpd)
                dxg, dyg = np.meshgrid(dx, dy)
                
                # Create variables for the user provided equation strings
                u_str = self.form_1_str_x
                v_str = self.form_1_str_y
                
                # Check if the equations provided contain x and y terms
                if u_str.find('x') & u_str.find('y') == -1:
                    u_str = '(' + str(u_str) + ')* np.ones(np.shape(dxg))'
                else:
                    u_str = u_str.replace('x', 'dxg')
                    u_str = u_str.replace('y', 'dyg')
          
                if v_str.find('x') & v_str.find('y') == -1:
                    v_str = '(' + str(v_str) + ')* np.ones(np.shape(dyg))'
                else:
                    v_str = v_str.replace('x', 'dxg')
                    v_str = v_str.replace('y', 'dyg')
                    
                # Generate arrays for the components of the zoom field
                u_zoom = eval(u_str)
                v_zoom = eval(v_str)
                
                zoom_form = form_1(dxg, dyg, u_zoom, v_zoom)
                
                return zoom_form
        
        # define a mehtod to evaluate the interior derivative of the 1-form
        # with respect to a given vector field object or without.
        def interior_d(self, vector_field=None, pass_on_figure=False, numerical_only=False):
            '''
            Computes the interior derivative of the 1-form
            Takes in:
            -- Vector_field = vector field object of formpy library to do the
            derivative with respect to, needs equations to work with
            nuymerical_only being False. Can also supply equations in a tuple:
            (eqn_x, eqn_y). If using numerical only, can supply object or
            tuple of numpy arrays (array_x, atrray_y). If nothing is supplied
            for it, it assumes F_x = 1 and F_y = 1, with correct form and shape
            
            -- pass_on_figure = determines if figure from this object is to be
            passed on to the 0-form object.
            
            --- numerical_only = bool, if true, it calculates only numerically
            otherwise, calculates it based on given equations, evaluates
            it numerically and supplies all to 0-form obejct creator
            
            Returns 0-form object
            '''
            
            # split up the code depending if numerical only or analytical too:
            if numerical_only is False:
                # test if equations were given first:
                if self.form_1_str_x == None or self.form_1_str_y == None:
                    # ERROR
                    raise ValueError('Error: You need to supply the 1-form equations to do this, look at \'give_eqn\' method')
                # if the vector field was supplied, extract its equations, if possible
                if vector_field is None:
                    # if none was given, do it with respect to uniform 1, 1
                    vf_x_str = '1'
                    vf_y_str = '1'
                elif type(vector_field) == tuple:
                    # if equations were given, take these, is numericals were given here, break!
                    if type(vector_field[0]) == str:
                        vf_x_str = vector_field[0]
                        vf_y_str = vector_field[1]
                    else:
                        raise ValueError('for analytical result, supply VF equations')
                else:
                    if vector_field.str_x == None or vector_field.str_y == None:
                        # ERROR
                        raise ValueError('Error: You need to supply the VF equations to do this, look at \'give_eqn\' method')
                    else:
                        vf_x_str = str(simplify(vector_field.str_x))
                        vf_y_str = str(simplify(vector_field.str_y))
                
                # combine them correctly with the 1-form strings:
                zero_form_str = str(simplify('(' + self.form_1_str_x + ')*(' + vf_x_str + ')' + ' + (' + self.form_1_str_y + ')*(' + vf_y_str + ')'))
                
                # keep an unformatted version to supply to the 0-form
                zero_form_str_unformatted = zero_form_str + ''
                
                # format the expression to be evluated
                zero_form_str = zero_form_str.replace('x', 'self.xg')
                zero_form_str = zero_form_str.replace('y', 'self.yg')
                
                # check against constants in the expression to be evaluated
                if zero_form_str.find('x') & zero_form_str.find('y') == -1:
                    zero_form_str = '(' + str(zero_form_str) + ')* np.ones(np.shape(self.xg))'
                else:
                    pass
                
                # evaulate the numerical zero form:
                zero_form_result = eval(zero_form_str)
                
                # return it, with equations, to user, depending on their figure
                # preferances
                
                if pass_on_figure is False:
                    # supply these to the 0-form object creator
                    result_form = form_0(self.xg, self.yg, zero_form_result, zero_form_str_unformatted)
                elif pass_on_figure is True:
                    if self.subplots is False:
                        result_form = form_0(self.xg, self.yg, zero_form_result, zero_form_str_unformatted, fig=self.figure, subplots=False)
                    elif self.subplots is True:
                        result_form = form_0(self.xg, self.yg, zero_form_result, zero_form_str_unformatted, fig=self.figure, subplots=True, sub_axis_list=self.axis)
                else:
                    raise ValueError('Error, Incorrect input for \' pass_on_figure \'')
                
                # return it to the user
                return result_form
            
            # deal with it if user wants to only do it numerically
            elif numerical_only is True:
                # check if equations have been given:
                # if they have, doing it only numerically would create
                # a mismatch, avoid that
                if self.form_1_str_x == None or self.form_1_str_y == None:
                    pass
                else:
                    # equations have been given, a mismatch may occur
                    # warn the user
                    print('Warning: You supplied equations, doing it numerically only will not pass equations to the 0-form and these will be lost')
                # now complete the process numerically save as instructed
                
                # Take the vector field components, checking what was input!
                if vector_field is None:
                    # if none was given, do it with respect to uniform 1, 1
                    vf_x = np.ones(np.shape(xg))
                    vf_y = np.ones(np.shape(xg))
                elif type(vector_field) == tuple:
                    # if equations were given, take these, is numericals were given here, break!
                    if type(vector_field[0]) == str:
                        raise ValueError('for numerical calulation, supply VF arrays, not equations')
                    else:
                        vf_x = vector_field[0]
                        vf_y = vector_field[1]
                else:
                    # extract needed properties from the object supplied
                    vf_x = vector_field.F_x
                    vf_y = vector_field.F_y
                
                # Complete the interior derivative 1-form --> 0-form:
                zero_form_result = self.F_x * vf_x + self.F_y * vf_y
                
                # return it to user:
                if pass_on_figure is False:
                    # supply these to the 0-form object creator
                    result_form = form_0(self.xg, self.yg, zero_form_result)
                elif pass_on_figure is True:
                    if self.subplots is False:
                        result_form = form_0(self.xg, self.yg, zero_form_result, fig=self.figure, subplots=False)
                    elif self.subplots is True:
                        result_form = form_0(self.xg, self.yg, zero_form_result, fig=self.figure, subplots=True, sub_axis_list=self.axis)
                else:
                    raise ValueError('Error, Incorrect input for \' pass_on_figure \'')
                
                # return it to the user
                return result_form
        
        # define a method to change a supplied Vector filed to the 1-form
        def vectorise(self, pass_on_figure=False, g=[['1', '0'], ['0', '1']]):
            '''
            Passes in everything it can (all it has been supplied)
            to the VF object.
            Works via the ('inverse') metric on R2
            Can supply the metric in as equations or as evaluated arrays
            Format of the metric is a list of numpy arrays
            0th array is the top row, its 0th component is 11, 1st is 12
            1st array is the botton row, its 0th comp is 21 and 1st is 22.
            Note, if it is supplied as arrays, they must come from numpy grids
            via meshgrid, if it is supplied as strings, needs to be in terms of
            x and y, and contain no special funtions, apart from ones imported
            here automatically and listed in the documentation #!!!
            Apart form the metric, takes in pass_on_figure:
                ----determines if figure from this object is to be
                passed on to the 1-form object.
            
            Returns a single object (VF object)
            '''
            
            # extract what is needed form the metric depending on what the user
            # supplied
            # check if its has string components
            if type(g[0][0]) == str and type(g[0][1]) == str and type(g[1][0]) == str and type(g[1][1]) == str:
                # deal with supplied string metric
                # need to format it, correct it for constants and evaluate it's numerical equivalent
                str_comp_00 = g[0][0] + ''
                str_comp_01 = g[0][1] + ''
                str_comp_10 = g[1][0] + ''
                str_comp_11 = g[1][1] + ''
                str_comp_00 = str_comp_00.replace('x', '(self.xg)')
                str_comp_00 = str_comp_00.replace('y', '(self.yg)')
                str_comp_01 = str_comp_01.replace('x', '(self.xg)')
                str_comp_01 = str_comp_01.replace('y', '(self.yg)')
                str_comp_10 = str_comp_10.replace('x', '(self.xg)')
                str_comp_10 = str_comp_10.replace('y', '(self.yg)')
                str_comp_11 = str_comp_11.replace('x', '(self.xg)')
                str_comp_11 = str_comp_11.replace('y', '(self.yg)')
                # check against constant form components:
                if str_comp_00.find('x') & str_comp_00.find('y') == -1:
                    str_comp_00 = '(' + str(str_comp_00) + ')* np.ones(np.shape(self.xg))'
                if str_comp_01.find('x') & str_comp_01.find('y') == -1:
                    str_comp_01 = '(' + str(str_comp_01) + ')* np.ones(np.shape(self.yg))'
                if str_comp_10.find('x') & str_comp_10.find('y') == -1:
                    str_comp_10 = '(' + str(str_comp_10) + ')* np.ones(np.shape(self.yg))'
                if str_comp_11.find('x') & str_comp_11.find('y') == -1:
                    str_comp_11 = '(' + str(str_comp_11) + ')* np.ones(np.shape(self.yg))'
                
                # evaluate the components numerically, inputting them into a
                # store numerical metric
                comp_00 = eval(str_comp_00)
                comp_01 = eval(str_comp_01)
                comp_10 = eval(str_comp_10)
                comp_11 = eval(str_comp_11)
                g_num = [[comp_00, comp_01], [comp_10, comp_11]]
                
                # set up a dummy variable to store the fact that numericals were given
                # not to check again later
                analytics = True
                
            elif type(g[0][0]) == np.ndarray and type(g[0][1]) == np.ndarray and type(g[1][0]) == np.ndarray and type(g[1][1]) == np.ndarray:
                # deal with the metric being supplied as components
                # if the user has 1-form equations, warn that these can't
                # be passed anymore, because we don't have equations for this
                # metric
                if self.form_1_str_x == None and self.form_1_str_y == None:
                    pass
                else:
                    print('The 1-form has equations, but the metric does not, these will be lost and the resulting VF will only have numerical values, not equations supplied')
                # No need to do anythng more to the metric, upto the user to make sure its
                # correctly sized, as with other code in this library
                # just rename the metric here
                g_num = g
                
                # set up a dummy variable to store the fact that numericals were
                # not given, not to check again later
                analytics = False
                
            else:
                # Inconsistant metric components
                raise TypeError('Metric components are inconcisstant')
            
            # from 1-form components, get VF components by the metric
            # first, do so numerically, as this must always happen
            form_x = self.F_x * g_num[0][0] + self.F_y * g_num[0][1]
            form_y = self.F_y * g_num[1][1] + self.F_x * g_num[1][0]
            
            # if the equations were given, evaluate these analytically too:
            # only if vector file doriginally has equations
            if analytics:
                if self.form_1_str_x == None and self.form_1_str_y == None:
                    print('You supplied the metric as equations (or it was default), but did not give 1-form equations, therefore only numericals will be completed')
                    analytics = False
                else:
                    x_str_form = '(' + self.form_1_str_x + ')*(' + g[0][0] + ') + (' + self.form_1_str_y + ')*(' + g[0][1] + ')'
                    y_str_form = '(' + self.form_1_str_y + ')*(' + g[1][1] + ') + (' + self.form_1_str_x + ')*(' + g[1][0] + ')'
                    # simplify them
                    x_str_form = str(simplify(x_str_form))
                    y_str_form = str(simplify(y_str_form))
            else:
                pass

            # based on what was given into the Vector field, return a 1-form object with these parameters
            if analytics:
                if pass_on_figure is False:
                    # supply these to the 1-form object creator
                    result_field = vector_field(self.xg, self.yg, form_x, form_y, x_str_form, y_str_form)
                elif pass_on_figure is True:
                    if self.subplots is False:
                        result_field = vector_field(self.xg, self.yg, form_x, form_y, x_str_form, y_str_form, fig=self.figure, subplots=False)
                    elif self.subplots is True:
                        result_field = vector_field(self.xg, self.yg, form_x, form_y, x_str_form, y_str_form, fig=self.figure, subplots=True, sub_axis_list=self.axis)
                else:
                    raise ValueError('Error, Incorrect input for \' pass_on_figure \'')
            elif not analytics:
                if pass_on_figure is False:
                    # supply these to the 1-form object creator
                    result_field = vector_field(self.xg, self.yg, form_x, form_y)
                elif pass_on_figure is True:
                    if self.subplots is False:
                        result_field = vector_field(self.xg, self.yg, form_x, form_y, fig=self.figure, subplots=False)
                    elif self.subplots is True:
                        result_field = vector_field(self.xg, self.yg, form_x, form_y, fig=self.figure, subplots=True, sub_axis_list=self.axis)
                else:
                    raise ValueError('Error, Incorrect input for \' pass_on_figure \'')
        
            # return the found object
            return result_field

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
def form_2(xg, yg, form2, form_2_eq=None, fig=None, subplots=False, sub_axis_list=[]):
    '''
    defines a 2-form object and returns it to user
    Takes 3 arguments basic, these are the 2 grids in 2D, which muse be square
    and of equal sizes. Then 1 argument for the dx^dy component
    based on the same grids. Also takes in an equation which is needed for some
    operaions
    Takes in a figure if one is to be supplied. Can take axis for subplots in
    The subplots only occur if subplots input is set to True, default is False
    '''
    # define the 1-form object and all its methods
    class form_set_up():
        # set up all variables
        def __init__(self, xg, yg, form2, s_max=6, s_min=2, fig_size=[7, 7], deltafactor=10):
            self.xg = xg
            self.yg = yg
            self.form_2 = form2
            # set up a figure if one was not given:
            if fig == None:
                self.figure = plt.figure(figsize=(fig_size[0], fig_size[1]))
                self.axis = self.figure.gca()
            else:
                if subplots is False:
                    self.figure = fig
                    self.axis  = self.figure.gca()
                elif subplots is True:
                    if len(sub_axis_list) == 0:
                        self.figure = fig
                        self.axis = []
                    else:
                        self.figure = fig
                        self.axis = sub_axis_list
                else:
                    raise ValueError('Error, incorrect input for \'subplots\'')
            self.subplots = subplots
            self.s_max = s_max
            self.s_min = s_min
            self.pt_den = len(xg[:, 0])  # + 1  # assume square grids
            self.fract = 2/((self.pt_den - 1))
            self.colour_list = ['red', 'blue', 'grey']
            self.logarithmic_scale_bool = 0
            self.delta_factor = deltafactor
            if form_2_eq is not None:
                self.form_2_str = str(simplify(form_2_eq))  # to start with, user must change to access some methods
                # Note, the string must be given with x and y as variables
            else:
                self.form_2_str = None
        
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
            
            # update the numerical values to always match
            string = self.form_2_str + ''
            string = string.replace('x', '(self.xg)')
            string = string.replace('y', '(self.yg)')
            
            # correct for consatnt form before evaluating
            if string.find('x') & string.find('y') == -1:
                string = '(' + str(string) + ')* np.ones(np.shape(self.xg))'
            else:
                pass
            
            # re-evaluate the 2-form numerically
            self.form_2 = eval(string)
        
        # deifne a function to return the string equation to the user
        def return_string(self):
            '''
            Takes in no arguments, returns the unformatted string back to user
            This is done in case user wants to access strings
            that got here not by input but by ext. alg.
            '''
            return self.form_2_str
        
        # define a method to add a subplot
        def add_subplot(self, order):
            '''
            Takes in one argument, as as matplotlib add_subplot
            It determines the shape of the subplot structure and which axis is
            being set (eg. 231 gives a set of plots 2 rows by 3 columns
            and this sets the first axis in that configuration)
            Adds a subplot to the figure that this form occupies
            Returns nothing, saves the subplot axis to the axis self param.
            '''
            # chck if the correct option was chosen to allow for this:
            if type(self.axis) != list:
                raise TypeError('Error, set up the object allowing for subplots')
            else:
                sub_axis = self.figure.add_subplot(order)
                self.axis.append(sub_axis)
        
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
                raise TypeError('Error: You need to supply the 2-form equation to do this, look at \'give_eqn\' method')
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
                str_2 = str_2.replace('x', '(self.xg)')
                str_2 = str_2.replace('y', '(self.yg)')
                
                # correct for consatnt form before evaluating
                if str_2.find('x') & str_2.find('y') == -1:
                    str_2 = '(' + str(str_2) + ')* np.ones(np.shape(self.xg))'
                else:
                    pass
                # re-evaluate the 2-form numerically
                self.form_2 = eval(str_2)
        
        # #####################################################################
        # Write more complicated methods. That will use this form object
        # eg. plot, exterior derivative, Hodge etc.
        # #####################################################################
        
        # define a function to plot the set up 2-form
        # originally form_2_components_plot
        def plot(self, keep=True, subplot_index=None):
            '''
            Finilises the plotting
            Takes in 2 inputs:
            1) \'keep\'determines if axis should be cleared before.
                Default is True
            2) \' subplot_index \', default set to None, can be input if
            the user has selected subplots to be allowed when creating the
            object. Determines which aixs to draw on, indecies are in order
            that they were added to the object.
            
            Uses the attribues of the object as set originally and as customised
            with methods to create a plot of the 2-form
            '''
            # for ease of later writting:
            # from self, get axis
            # depending on if subplots are wanted:
            if type(self.axis) != list:
                axis = self.axis
            else:
                axis = self.axis[subplot_index]
            
            
            form2 = self.form_2  # from self, get 2-form too
            
            # set all insignificant values to zero:
            form2[np.abs(form2) < 1e-15] = 0
            
            # depending on input, clear the axis or don't
            if keep is True:
                pass
            else:
                axis.clear()
            
            # get the lengths of x and y from their grids
            x_len = len(self.xg[:, 0])
            y_len = len(self.yg[0, :])
            
            # find L based on the origin of given grid is
            L = 0.5*(self.xg[0, -1] - self.xg[0, 0])
            x0 = self.xg[0,0] + L
            y0 = self.yg[0,0] + L
            
            # rescale axis
            ax_L = L + L/self.delta_factor
            axis.set_xlim(-ax_L + x0, ax_L + x0)
            axis.set_ylim(-ax_L + y0, ax_L + y0)
            
            # get the signs of the input 2-form
            form_2_sgn = np.sign(form2)
            
            # define an empty array of magnitudes, to then fill with integer rel. mags
            R_int = np.zeros(shape=((x_len), (y_len)))
            
            # #########################################################################
            # get variables needed for the initial, simplified stack plot
            # #########################################################################
            
            # set up directions
            angles =[0*np.ones(np.shape(form2)), (np.pi/2)*np.ones(np.shape(form2))]
            
            # deal with sinularities that appear on evaluated points
            for i in range(x_len):
                for j in range(y_len):
                    # set to zero points that are not defined or inf
                    if isnan(form2[i, j]) is True or abs(form2[i, j]) == np.inf  or abs(form2[i, j]) > 1e15:
                        # colour this region as a red dot, not square to
                        # not confuse with nigh mag 2-forms in stacks. or worse, in
                        # blocks
                        circ = patch.Circle((self.xg[i, j], self.yg[i, j]), L*self.fract/3, color='red')
                        axis.add_patch(circ)
                        form2[i, j] = 0
                    # ALso, since we got this lop anyway
                    # correct for singularities in planar form 2:
                    # set to zero points that are not defined or inf
                    if isnan(form2[i, j]) is True:
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
            max_size = np.max(abs(form2))   # careful with singularities, else ---> nan
            
            # find the relative magnitude of vectors to maximum, as an array
            R = abs(form2)/max_size
            
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
        
        # define a fucntion to Hodge the 2-form (into a 0-form)
        def Hodge(self, numerical_only=True, pass_on_figure=False):
            '''
            Takes in two bool arguments:
            
            numerical_only
            Determines if the calculation should be numerical or analytic
            True if numerical, False if analytic and numerical. Default is True
            For the analytic one, the equations as string must be supplied to
            the object
            Note, choosing analytical, changes the equations AND the numerical answers
            It calulates the Hodge on R^2 by the standard definition:
            *(dx^dy) = 1
            
            pass_on_figure
            Determies if figure should be passed onto the new object
            if it is to be created
            
            returns a 0-form
            '''
            # distinguish between doing it numerically and alaytically
            if numerical_only is True:
                # check if equations have been given:
                # if they have, doing it only numerically would create
                # a mismatch, avoid that
                if self.form_2_str == None:
                    pass
                else:
                    # equations have been given, a mismatch may occur
                    # warn the user
                    print('Warning: You supplied equations, doing it numerically only will result in a mismacth between numerical values and equations')
                # now complete the process numerically
                # pass these in to the object to create a new one:
                if pass_on_figure is False:
                    new_object = form_0(self.xg, self.yg, self.form_2, self.form_2_str)  # N.B no equations to supply
                elif pass_on_figure is True:
                     new_object = form_0(self.xg, self.yg, self.form_2, self.form_2_str, fig=self.figure, subplots=self.subplots, sub_axis_list=self.axis)
                else:
                    raise ValueError('Error, Incorrect input for \' pass_on_figure \'')
                # return the new one to the user:
                return new_object
            
            elif numerical_only is False:
                # can only be done if equations have been given, check:
                if self.form_2_str == None:
                    # ERROR
                    raise TypeError('Error: You need to supply the 2-form equation to do this, look at \'give_eqn\' method')
                else:
                    # some equations are there, compute the Hodge on these:
                    # Note: Upto user to make sure their equations match their
                    # numerical input, unless using give eqn, then its updates
                    # numerical values to match
                    
                    # get numerical solutions, evaulated on local
                    # strings changed to relate to the self grids
                    # need to uspply these unformatted, so save those:
                    form_0_str_unformated = self.form_2_str + '' 
                    string_0_form = self.form_2_str  # formated
                    # from these strings, get the numerical 0-form:
                    string_0_form = string_0_form.replace('x', '(self.xg)')
                    string_0_form = string_0_form.replace('y', '(self.yg)')
                    
                    # correct for constant forms
                    if string_0_form.find('x') & string_0_form.find('y') == -1:
                        string_0_form = '(' + str(string_0_form) + ')* np.ones(np.shape(self.xg))'
                    
                    # evaulated numerically
                    form_0_result = eval(string_0_form)
                    
                    # return object, depending on option for figure passage:
                    # pass these in to the object to create a new one:
                    if pass_on_figure is True:
                        if self.subplots is False:
                            new_object = form_0(self.xg, self.yg, form_0_result, form_0_eqn=form_0_str_unformated, fig=self.figure, subplots=False)
                        if self.subplots is True:
                            new_object = form_0(self.xg, self.yg, form_0_result, form_0_eqn=form_0_str_unformated, fig=self.figure, subplots=True, sub_axis_list=self.axis)
                    elif pass_on_figure is False:
                        new_object = form_0(self.xg, self.yg, form_0_result, form_0_eqn=form_0_str_unformated)
                    else:
                        raise ValueError('Error, Incorrect input for \' pass_on_figure \'')
                    
                    # return the new one to the user:
                    return new_object
            else:
                # Error
                raise ValueError('ERROR: Invalid input for \'numerical_only\'')
    
        # define a method to create a zoomed in 2-form
        def zooming(self, target=[0, 0], zoom=2, dpd=9):
            '''
            Creates a new window which displays the 2-form zoomed at a certain point
            User gives arguments:
            Target: Determines the zoom location, coordinates
            Zoom: +ve float, determines zooming amount
            dpd: +int, determines how many points on each axis
            '''
            
            # Requires user to provide eqn of the 1-form they are zooming on.
            if self.form_2_str == None:
                # ERROR
                raise TypeError('Error: No equation provided')
            else:
                
                # Target coordinates
                x_m = target[0]
                y_m = target[1]
                
                d_range = self.xg[0, -1]/zoom
                
                # Set up zoom window grids
                dx = np.linspace(-d_range + x_m, d_range + x_m, dpd)
                dy = np.linspace(-d_range + y_m, d_range + y_m, dpd)
                dxg, dyg = np.meshgrid(dx, dy)
                
                # Create variables for the user provided equation strings
                zoom_str = self.form_2_str + ''
                
                # Check if the equations provided contain x and y terms
                if zoom_str.find('x') & zoom_str.find('y') == -1:
                    zoom_str = '(' + str(zoom_str) + ')* np.ones(np.shape(dxg))'
                else:
                    zoom_str = zoom_str.replace('x', '(dxg)')
                    zoom_str = zoom_str.replace('y', '(dyg)')
                
                # Generate arrays for the components of the zoom field
                zoom_2form = eval(zoom_str)
                
                # set up a new 2-form, that is the form, after zooming in
                # that will be returned to the user
                zoom_form = form_2(dxg, dyg, zoom_2form, form_2_eq=self.form_2_str)
                
                return zoom_form
        
        # define a mehtod to evaluate the interior derivative of the 2-form
        # with respect to a given vector field object or without.
        def interior_d(self, vector_field=None, pass_on_figure=False, numerical_only=False):
            '''
            Computes the interior derivative of the 2-form
            Takes in:
            -- Vector_field = vector field object of formpy library to do the
            derivative with respect to, needs equations to work with
            nuymerical_only being False. Can also supply equations in a tuple:
            (eqn_x, eqn_y). If using numerical only, can supply object or
            tuple of numpy arrays (array_x, atrray_y). If nothing is supplied
            for it, it assumes F_x = 1 and F_y = 1, with correct form and shape
            
            -- pass_on_figure = determines if figure from this object is to be
            passed on to the 1-form object.
            
            --- numerical_only = bool, if true, it calculates only numerically
            otherwise, calculates it based on given equations, evaluates
            it numerically and supplies all to 1-form obejct creator
            
            '''
            
            # split up the code depending if numerical only or analytical too:
            if numerical_only is False:
                # test if the equation was given first:
                if self.form_2_str == None:
                    # ERROR
                    raise ValueError('Error: You need to supply the 2-form equations to do this, look at \'give_eqn\' method')
                
                # if the vector field was supplied, extract its equations, if possible
                if vector_field is None:
                    # if none was given, do it with respect to uniform 1, 1
                    vf_x_str = '1'
                    vf_y_str = '1'
                elif type(vector_field) == tuple:
                    # if equations were given, take these, is numericals were given here, break!
                    if type(vector_field[0]) == str:
                        vf_x_str = vector_field[0]
                        vf_y_str = vector_field[1]
                    else:
                        raise ValueError('for analytical result, supply VF equations')
                else:
                    if vector_field.str_x == None or vector_field.str_y == None:
                        # ERROR
                        raise ValueError('Error: You need to supply the VF equations to do this, look at \'give_eqn\' method')
                    else:
                        vf_x_str = str(simplify(vector_field.str_x))
                        vf_y_str = str(simplify(vector_field.str_y))
                
                
                # define strings of the resulting 1-form components
                u_str = str(simplify('-(' + self.form_2_str + ')*(' + vf_y_str + ')' ))
                v_str = str(simplify( '(' + self.form_2_str + ')*(' + vf_x_str + ')' ))
                
                # keep an unformatted version to supply to the 1-form
                u_str_unformatted = u_str + ''
                v_str_unformatted = v_str + ''
                
                u_str = u_str.replace('x', '(self.xg)')
                u_str = u_str.replace('y', '(self.yg)')
                v_str = v_str.replace('x', '(self.xg)')
                v_str = v_str.replace('y', '(self.yg)')
                if u_str.find('x') & u_str.find('y') == -1:
                    u_str = '(' + str(u_str) + ')* np.ones(np.shape(self.xg))'
                if v_str.find('x') & v_str.find('y') == -1:
                    v_str = '(' + str(v_str) + ')* np.ones(np.shape(self.yg))'
                
                # evaulate the numerical 1-form components form:
                form_x = eval(u_str)
                form_y = eval(v_str)
                
                # return it, with equations, to user, depending on their figure
                # preferances
                
                if pass_on_figure is False:
                    # supply these to the 1-form object creator
                    result_form = form_1(self.xg, self.yg, form_x, form_y, u_str_unformatted, v_str_unformatted)
                elif pass_on_figure is True:
                    if self.subplots is False:
                        result_form = form_1(self.xg, self.yg, form_x, form_y, u_str_unformatted, v_str_unformatted, fig=self.figure, subplots=False)
                    elif self.subplots is True:
                        result_form = form_1(self.xg, self.yg, form_x, form_y, u_str_unformatted, v_str_unformatted, fig=self.figure, subplots=True, sub_axis_list=self.axis)
                else:
                    raise ValueError('Error, Incorrect input for \' pass_on_figure \'')
                
                # return it to the user
                return result_form
            
            # deal with it if user wants to only do it numerically
            elif numerical_only is True:
                # check if equations have been given:
                # if they have, doing it only numerically would create
                # a mismatch, avoid that
                if self.form_2_str == None:
                    pass
                else:
                    # equations have been given, a mismatch may occur
                    # warn the user
                    print('Warning: You supplied equations, doing it numerically only will not pass equations to the 1-form and these will be lost')
                # now complete the process numerically save as instructed
                
                # Take the vector field components, checking what was input!
                if vector_field is None:
                    # if none was given, do it with respect to uniform 1, 1
                    vf_x = np.ones(np.shape(xg))
                    vf_y = np.ones(np.shape(xg))
                elif type(vector_field) == tuple:
                    # if equations were given, take these, is numericals were given here, break!
                    if type(vector_field[0]) == str:
                        raise ValueError('for numerical calulation, supply VF arrays, not equations')
                    else:
                        vf_x = vector_field[0]
                        vf_y = vector_field[1]
                else:
                    # extract needed properties from the object supplied
                    vf_x = vector_field.F_x
                    vf_y = vector_field.F_y
                
                # Complete the interior derivative 2-form --> 1-form:
                form_x = -self.form_2 * vf_y
                form_y = self.form_2 * vf_x
                
                if pass_on_figure is False:
                    # supply these to the 1-form object creator
                    result_form = form_1(self.xg, self.yg, form_x, form_y)
                elif pass_on_figure is True:
                    if self.subplots is False:
                        result_form = form_1(self.xg, self.yg, form_x, form_y, fig=self.figure, subplots=False)
                    elif self.subplots is True:
                        result_form = form_1(self.xg, self.yg, form_x, form_y, fig=self.figure, subplots=True, sub_axis_list=self.axis)
                else:
                    raise ValueError('Error, Incorrect input for \' pass_on_figure \'')
                
                # return it to the user
                return result_form
        
    # now call that object to create it:
    form_2_object = form_set_up(xg, yg, form2)
    # return it to user to store
    return form_2_object


# %%

'''

function to create a 0-form object and define methods for it

'''


# define a function that will set up a 0-form object that can be customised and
# plotted
def form_0(xg, yg, form_0, form_0_eqn=None, fig=None, subplots=False, sub_axis_list=[]):
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
            # set up a figure if one was not given:
            if fig == None:
                self.figure = plt.figure(figsize=(fig_size[0], fig_size[1]))
                self.axis = self.figure.gca()
            else:
                if subplots is False:
                    self.figure = fig
                    self.axis  = self.figure.gca()
                elif subplots is True:
                    if len(sub_axis_list) == 0:
                        self.figure = fig
                        self.axis = []
                    else:
                        self.figure = fig
                        self.axis = sub_axis_list
                else:
                    raise ValueError('Error, incorrect input for \'subplots\'')
            
            self.subplots = subplots
            self.pt_den = len(xg[:, 0])  # + 1  # assume square grids
            self.logarithmic_scale_bool = 0
            self.delta_factor = deltafactor
            self.denser = 1
            self.lines = 15
            self.fontsize = 7
            self.inline_bool = True
            if form_0_eqn is not None:
                self.form_0_str = str(simplify(form_0_eqn))  # to start with, use rmust change to access some methods
            else:
                self.form_0_str = None
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
            
            # update the numerical values to always match
            string = self.form_0_str + ''
            
            # Check if the equations provided contain x and y terms
            # and format them to be evaluated
            if string.find('x') & string.find('y') == -1:
                string = '(' + str(string) + ')* np.ones(np.shape(xg))'
            else:
                string = string.replace('x', '(self.xg)')
                string = string.replace('y', '(self.yg)')
            
            # re-evaluate the 2-form numerically, preventing mismatch
            self.form_0 = eval(string)
            
        # deifne a function to return the string equation to the user
        def return_string(self):
            '''
            Takes in no arguments, returns the unformatted string back to user
            This is done in case user wants to access strings
            that got here not by input but by ext. alg.
            '''
            return self.form_0_str
        
        # define a method to add a subplot
        def add_subplot(self, order):
            '''
            Takes in one argument, as as matplotlib add_subplot
            It determines the shape of the subplot structure and which axis is
            being set (eg. 231 gives a set of plots 2 rows by 3 columns
            and this sets the first axis in that configuration)
            Adds a subplot to the figure that this form occupies
            Returns nothing, saves the subplot axis to the axis self param.
            '''
            # chck if the correct option was chosen to allow for this:
            if type(self.axis) != list:
                raise ValueError('Error, set up the object allowing for subplots')
            else:
                sub_axis = self.figure.add_subplot(order)
                self.axis.append(sub_axis)
        
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
                raise TypeError('Error: You need to supply the 0-form equation to do this, look at \'give_eqn\' method')
            else:
                # redefine the grids
                v = np.linspace(-self.xg[0, -1], self.xg[0, -1], points_number)
                self.xg, self.yg = np.meshgrid(v, v)
                # based on these change other, dependant variables
                self.pt_den = len(self.xg[:, 0])
                # substitute these into the equation:
                # but keep it local
                str_0 = self.form_0_str + ''
                str_0 = str_0.replace('x', '(self.xg)')
                str_0 = str_0.replace('y', '(self.yg)')
                # correct for constant forms
                if str_0.find('x') & str_0.find('y') == -1:
                    str_0 = '(' + str(str_0) + ')* np.ones(np.shape(self.xg))'
                # re-evaluate the 2-form numerically
                self.form_0 = eval(str_0)
        
        # #####################################################################
        # Write more complicated methods. That will use this form object
        # eg. plot, exterior derivative, Hodge etc.
        # #####################################################################
        
        
        # define a fucntion to plot a zero form when button is pressed.
        def plot(self, keep=True, subplot_index=None):
            '''
            Finilises the plotting
            Uses the attribues of the object as set originally and as customised
            with methods to create a plot of the 2-form.
            Can take one parameter: bool. Default is True
            determines if axis should be cleared before plotting.
            Another parameter it takes:
                index of subplot on which to plot it, if object was set up with
                axis
            '''
            
            # for ease of later writting:
            # from self, get axis
            # depending on if subplots are wanted:
            if type(self.axis) != list:
                axis = self.axis
            else:
                axis = self.axis[subplot_index]
            
            # check if user wants to clear first:
            if keep is True:
                pass
            else:
                axis.clear()
            
            form_0 = self.form_0  # from self, get 2-form
            
            # set all insignificant values to zero:
            form_0[np.abs(form_0) < 1e-15] = 0
            
            # find L based on the origin of given grid is
            L = 0.5*(self.xg[0, -1] - self.xg[0, 0])
            x0 = self.xg[0,0] + L
            y0 = self.yg[0,0] + L
            
            # rescale axis
            ax_L = L + L/self.delta_factor
            axis.set_xlim(-ax_L + x0, ax_L + x0)
            axis.set_ylim(-ax_L + y0, ax_L + y0)
            
            if self.denser != 1:
                if self.form_0_str == None:
                    # This cannot be done if a string has not been supplied
                    # ERROR
                    raise TypeError('Error: You need to supply the 0-form equation to do this, look at \'give_eqn\' method')
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
        def ext_d(self, pass_on_figure=False):
            '''
            Takes in 1 argument:
            -- pass on figure: Determies if figure should be passed onto the
               new object if it is to be created
            Returns 1 form object
            computes the exterior derivative and returns it as the 1-form object
            '''
            
            # first make sure that the string has been supplied
            if self.form_0_str == None:
                    # ERROR
                    raise TypeError('Error: You need to supply the 0-form equation to do this, look at \'give_eqn\' method')
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
                form_1_x_str = form_1_x_str.replace('x', '(self.xg)')
                form_1_x_str = form_1_x_str.replace('y', '(self.yg)')
                form_1_y_str = form_1_y_str.replace('x', '(self.xg)')
                form_1_y_str = form_1_y_str.replace('y', '(self.yg)')
                if form_1_x_str.find('x') & form_1_x_str.find('y') == -1:
                    form_1_x_str = '(' + str(form_1_x_str) + ')* np.ones(np.shape(self.xg))'
                if form_1_y_str.find('x') & form_1_y_str.find('y') == -1:
                    form_1_y_str = '(' + str(form_1_y_str) + ')* np.ones(np.shape(self.yg))'
                form_1_x = eval(form_1_x_str)
                form_1_y = eval(form_1_y_str)
                # supply these to the 1-form object function
                if pass_on_figure is True:
                    if self.subplots is False:
                        result_1_form = form_1(self.xg, self.yg, form_1_x, form_1_y, form_1_x_unformated, form_1_y_unformated, fig=self.figure, subplots=False)
                    elif self.subplots is True:
                        result_1_form = form_1(self.xg, self.yg, form_1_x, form_1_y, form_1_x_unformated, form_1_y_unformated, fig=self.figure, subplots=True, sub_axis_list=self.axis)
                    return result_1_form
                elif pass_on_figure is False:
                    result_1_form = form_1(self.xg, self.yg, form_1_x, form_1_y, form_1_x_unformated, form_1_y_unformated)
                else:
                     raise ValueError('Error, Incorrect input for \' pass_on_figure \'')
                return result_1_form
        
        # deifne a method to complete the exterior derivative numerically
        def num_ext_d(self, edge_order=1, pass_on_figure=False):
            '''
            Takes in 2 arguments
            -- edge_order: determines order same as in numpy gradient {1 or 2}
            -- pass on figure: Determies if figure should be passed onto the
               new object if it is to be created
            
            Return 1 object - 1-form
            computes the exterior derivative numerically only
            The equations do not need to be given
            If given, they do not get passed onto the 1-form object anyway
            NUMERICAL ONLY
            '''
            
            # from numpy gradient, get the gradient array
            fy, fx = np.gradient(form_0, edge_order=edge_order)
            # supply these to the 1-form object function
            if pass_on_figure is True:
                if self.subplots is False:
                    result_1_form = form_1(self.xg, self.yg, fx, fy, fig=self.figure, subplots=False)
                elif self.subplots is True:
                    result_1_form = form_1(self.xg, self.yg, fx, fy, fig=self.figure, subplots=True, sub_axis_list=self.axis)
            elif pass_on_figure is False:
                result_1_form = form_1(self.xg, self.yg, fx, fy)
            else:
                raise ValueError('Error, Incorrect input for \' pass_on_figure \'')
            
            # return the new object to user
            return result_1_form
        
        # deinfe a method for Hodge of a 0-form
        def Hodge(self, numerical_only=True, pass_on_figure=False):
            '''
            Takes in two bool arguments:
            
            numerical_only
            Determines if the calculation should be numerical or analytic
            True if numerical, False if analytic and numerical. Default is True
            For the analytic one, the equations as string must be supplied to
            the object
            Note, choosing analytical, changes the equations AND the numerical answers
            It calulates the Hodge on R^2 by the standard definition:
            1* = (dx^dy)
            
            pass_on_figure
            Determies if figure should be passed onto the new object
            if it is to be created
            
            returns a 2-form
            '''
            # distinguish between doing it numerically and alaytically
            if numerical_only is True:
                # check if equations have been given:
                # if they have, doing it only numerically would create
                # a mismatch, avoid that
                if self.form_0_str == None:
                    pass
                else:
                    # equations have been given, a mismatch may occur
                    # warn the user
                    print('Warning: You supplied equations, doing it numerically only will result in a mismacth between numerical values and equations')
                # now complete the process numerically
                # pass these in to the object to create a new one:
                if pass_on_figure is False:
                    new_object = form_2(self.xg, self.yg, self.form_0, self.form_0_str)  # N.B no equations to supply
                elif pass_on_figure is True:
                     new_object = form_2(self.xg, self.yg, self.form_0, self.form_0_str, fig=self.figure, subplots=self.subplots, sub_axis_list=self.axis)
                else:
                    raise ValueError('Error, Incorrect input for \' pass_on_figure \'')
                # return the new one to the user:
                return new_object
            
            elif numerical_only is False:
                # can only be done if equations have been given, check:
                if self.form_0_str == None:
                    # ERROR
                    raise TypeError('Error: You need to supply the 2-form equation to do this, look at \'give_eqn\' method')
                else:
                    # some equations are there, compute the Hodge on these:
                    # Note: Upto user to make sure their equations match their
                    # numerical input, unless using give eqn, then its updates
                    # numerical values to match
                    
                    # get numerical solutions, evaulated on local
                    # strings changed to relate to the self grids
                    # need to uspply these unformatted, so save those:
                    form_2_str_unformated = self.form_0_str + '' 
                    string_2_form = self.form_0_str  # to be formated
                    # from these strings, get the numerical 0-form:
                    string_2_form = string_2_form.replace('x', '(self.xg)')
                    string_2_form = string_2_form.replace('y', '(self.yg)')
                    
                    if string_2_form.find('x') & string_2_form.find('y') == -1:
                        string_2_form = '(' + str(string_2_form) + ')* np.ones(np.shape(self.xg))'
                    
                    # evaulated numerically
                    form_2_result = eval(string_2_form)
                    
                    # return object, depending on option for figure passage:
                    # pass these in to the object to create a new one:
                    if pass_on_figure is True:
                        if self.subplots is False:
                            new_object = form_2(self.xg, self.yg, form_2_result, form_2_eq=form_2_str_unformated, fig=self.figure, subplots=False)
                        if self.subplots is True:
                            new_object = form_2(self.xg, self.yg, form_2_result, form_2_eq=form_2_str_unformated, fig=self.figure, subplots=True, sub_axis_list=self.axis)
                    elif pass_on_figure is False:
                        new_object = form_2(self.xg, self.yg, form_2_result, form_2_eq=form_2_str_unformated)
                    else:
                        raise ValueError('Error, Incorrect input for \' pass_on_figure \'')
                    
                    # return the new one to the user:
                    return new_object
            else:
                # Error
                raise ValueError('ERROR: Invalid input for \'numerical_only\'')
    
    # now call that object to create it:
    form_0_object = form_set_up(xg, yg, form_0)
    # return it to user to store
    return form_0_object


# %%

'''

function to create a vector field object and define methods for it

'''

# define a function that will set up a vector field object that can be customised and
# plotted
def vector_field(xg, yg, F_x, F_y, F_x_eqn=None, F_y_eqn=None, fig=None, subplots=False, sub_axis_list=[]):
    '''
    defines a vector field object and returns it to user
    Takes 9 arguments, these are the 2 grids in 2D, which muse be square
    and of equal sizes. Then 2 arguments for the i component and j component
    based on the same grids
    Then 2 equations, for x and y (not always needed)
    Then, can supply a figure if user doesn't wat this object
    to create a new one for itself
    Can also supply a bool input to define if subplots are to be allowed
    these can be added using a method (add_subplot)
    also, user can supply sub_axis_list, to provide the axis they have set
    these only work if a figure has been supplied and if subplots is True
    '''
    # define the 1-form object and all its methods
    class field_set_up():
        # set up all variables
        def __init__(self, xg, yg, F_x, F_y, deltafactor=10, fig_size=[7, 7]):
            self.xg = xg
            self.yg = yg
            self.F_x = F_x
            self.F_y = F_y
            # set up a figure if one was not given:
            if fig == None:
                self.figure = plt.figure(figsize=(fig_size[0], fig_size[1]))
                self.axis = self.figure.gca()
            else:
                if subplots is False:
                    self.figure = fig
                    self.axis  = self.figure.gca()
                elif subplots is True:
                    if len(sub_axis_list) == 0:
                        self.figure = fig
                        self.axis = []
                    else:
                        self.figure = fig
                        self.axis = sub_axis_list
                else:
                    raise TypeError('Error, incorrect input for \'subplots\'')
            self.subplots = subplots
            self.pt_den = len(xg[:, 0])# + 1  # assume square grids
            self.orientation = 'mid'
            self.scale = 1
            self.color = 'black'
            self.logarithmic_scale_bool = 0
            self.scale_bool = True
            self.delta_factor = deltafactor
            if F_x_eqn is not None:
                self.str_x = str(simplify(F_x_eqn))  # to start with, use rmust change to access some methods
                # Note, the string must be given with x and y as variables
            else:
                self.str_x = None
            
            if F_y_eqn is not None:
                self.str_y = str(simplify(F_y_eqn))
            else:
                self.str_y = None
        
        # #####################################################################
        # write some methods that will allow the user to chenge some of the
        # above variables
        # #####################################################################
        # deifne a function to return the string equations to the user
        def give_eqn(self, equation_str_x, equation_str_y):
            '''
            Takes in 1-argument, string
            This must be the equation of the supplied numerical 0-form
            It must be in terms of x and y.
            Has to be given, for some methods to be calculatable.
            '''
            # set equation parameters to simplified inputs
            self.str_x = str(simplify(equation_str_x))
            self.str_y = str(simplify(equation_str_y))
            # make the values match automatically to limit how often mismatch occurs
            # substitute these into the equation:
            # but keep it local
            str_x = self.str_x + ''
            str_y = self.str_y + ''
            str_x = str_x.replace('x', '(self.xg)')
            str_x = str_x.replace('y', '(self.yg)')
            str_y = str_y.replace('x', '(self.xg)')
            str_y = str_y.replace('y', '(self.yg)')
            # check kagainst constant form components:
            if str_x.find('x') & str_x.find('y') == -1:
                str_x = '(' + str(str_x) + ')* np.ones(np.shape(self.xg))'
            if str_y.find('x') & str_y.find('y') == -1:
                str_y = '(' + str(str_y) + ')* np.ones(np.shape(self.yg))'
            
            # re-evaluate the 2-form numerically
            self.F_x = eval(str_x)
            self.F_y = eval(str_y)
            
        def return_string(self):
            '''
            Takes in no arguments, returns the unformatted strings back to user
            This is done in case user wants to access strings
            that got here not by input but by ext. alg.
            '''
            return self.str_x, self.str_y
        
        # define a method to add a subplot
        def add_subplot(self, order):
            '''
            Takes in one argument, as matplotlib add_subplot
            It determines the shape of the subplot structure and which axis is
            being set (eg. 231 gives a set of plots 2 rows by 3 columns
            and this sets the first axis in that configuration)
            Adds a subplot to the figure that this form occupies
            Returns nothing, saves the subplot axis to the axis self param.
            '''
            # chck if the correct option was chosen to allow for this:
            if type(self.axis) != list:
                raise ValueError('Error, set up the object allowing for subplots')
            else:
                sub_axis = self.figure.add_subplot(order)
                self.axis.append(sub_axis)
        
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
        
        # change orientation:
        def orient(self, string):
            '''
            Takes one input, needs to be a string understood by matplotlib
            quiver to orient arrows
            eg. 'tip', 'tail', 'mid' etc.
            Orients arrows on quiver plot depending on this
            '''
            self.orientation = str(string)
        
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
            if self.str_x == None or self.str_y == None:
                # Error
                raise ValueError('Error: You need to supply the 1-form equation to do this, look at \'give_eqn\' method')
            else:
                # redefine the grids
                v = np.linspace(-self.xg[0, -1], self.xg[0, -1], points_number)
                self.xg, self.yg = np.meshgrid(v, v)
                # based on these change other, dependant variables
                self.pt_den = len(self.xg[:, 0])
                # substitute these into the equation:
                # but keep it local
                str_x_l = self.str_x + ''
                str_y_l = self.str_y + ''
                str_x_l = str_x_l.replace('x', '(self.xg)')
                str_x_l = str_x_l.replace('y', '(self.yg)')
                str_y_l = str_y_l.replace('x', '(self.xg)')
                str_y_l = str_y_l.replace('y', '(self.yg)')
                # check kagainst constant form components:
                if str_x_l.find('x') & str_x_l.find('y') == -1:
                    str_x_l = '(' + str(str_x_l) + ')* np.ones(np.shape(self.xg))'
                if str_y_l.find('x') & str_y_l.find('y') == -1:
                    str_y_l = '(' + str(str_y_l) + ')* np.ones(np.shape(self.yg))'
                
                # re-evaluate the 2-form numerically
                self.F_x = eval(str_x_l)
                self.F_y = eval(str_y_l)
        
        # define a method to plot the vector field using quiver
        def plot(self, keep=True, subplot_index=None):
            '''
            Finilises the plotting
            Uses the attribues of the object as set originally and as customised
            with methods to create a plot of the 1-form
            Takes in 2 arguments:
            --- \'keep\', which allows the user to plot
            on top of the previous pltos they created without clearing axis
            When its set to False, the axis are cleared first
            Default is True.
            --- \' subplot_index \', default set to None, can be input if
            the user has selected subplots to be allowed when creating the
            object. Determines which aixs to draw on, indecies are in order
            that they were added to the object
            '''
            
            # from self, get axis
            # depending on if subplots are wanted:
            if type(self.axis) != list:
                axis = self.axis
            else:
                axis = self.axis[subplot_index]
            
            # depending on input, clear the axis or don't
            if keep is True:
                pass
            else:
                axis.clear()
            
            # get the lengths of x and y from their grids
            x_len = len(self.xg[:, 0])
            y_len = len(self.yg[0, :])
            
            # Extract L from the x and y grids. Assumes they are square.
            L = 0.5*(self.xg[0, -1] - self.xg[0, 0])
            x0 = self.xg[0,0] + L
            y0 = self.yg[0,0] + L
            
            # adjust axis limits based on that.
            ax_L = L + L/self.delta_factor
            axis.set_xlim(-ax_L + x0, ax_L + x0)
            axis.set_ylim(-ax_L + y0, ax_L + y0)
            
            # for arrows to work, with nan and infs
            # make a local variable of F_x and F_y
            # so that thye don't alter globally
            F_x_local = self.F_x * 1
            F_y_local = self.F_y * 1
            
            # prevent any magnitudes from being inf or nan
            # only here, need to do it to u and v not just mag
            for i in range(x_len):
                for j in range(y_len):
                    if isnan(F_x_local[i,j]) == True or isnan(F_y_local[i,j]) == True or abs(F_x_local[i, j]) == np.inf or abs(F_y_local[i, j]) == np.inf or abs(F_y_local[i, j]) > 1e15 or abs(F_x_local[i, j]) > 1e15:
                        F_x_local[i,j] = F_y_local[i,j] = 0
            
            
            # set all insignificant values to zero:
            F_x_local[np.abs(F_x_local) < 1e-15] = 0
            F_y_local[np.abs(F_y_local) < 1e-15] = 0
            
            
            # find the magnitude corresponding to each point and store in mag array
            mag = np.sqrt(F_x_local**2 + F_y_local**2)
            
            # find the maximum magnitude for scaling
            max_size = np.max(mag)   # careful with singularities, else ---> nan
            
            # deal with requested autoscaling
            if self.scale_bool is False:
                ScaleFactor = self.scale
            elif self.scale_bool is True:
                ScaleFactor = max_size/(0.9*(2*L/self.pt_den))
            
            
            axis.quiver(self.xg, self.yg, F_x_local, F_y_local, pivot=self.orientation, scale=ScaleFactor, scale_units='xy', color=self.color) 
        
        def zoom(self, target=[0,0], zoom=2, dpd=9, pass_on_figure=False):
            '''
            Create a new window which displays the field zoomed at a certain point
            User gives arguments
            Target: Determines the zoom location, coordinates
            Zoom: +ve float, determines zooming amount
            dpd: +int, determines how many points on each axis
            '''
            
            # Requires user to provide eqn of the 1-form they are zooming on.
            
            if self.str_x == None or self.str_y == None:
                # ERROR
                raise TypeError('Error: No equation provided')
            else:
                
                # Target coordinates
                x_m = target[0]
                y_m = target[1]
                
                d_range = self.xg[0, -1]/zoom
                #d_length = 
                
                # Set up zoom window grids
                dx = np.linspace(-d_range + x_m, d_range + x_m, dpd)
                dy = np.linspace(-d_range + y_m, d_range + y_m, dpd)
                dxg, dyg = np.meshgrid(dx, dy)
                
                # Create variables for the user provided equation strings
                u_str = self.str_x
                v_str = self.str_y
                
                # Check if the equations provided contain x and y terms
                if u_str.find('x') & u_str.find('y') == -1:
                    u_str = '(' + str(u_str) + ')* np.ones(np.shape(dxg))'
                else:
                    u_str = u_str.replace('x', 'dxg')
                    u_str = u_str.replace('y', 'dyg')
          
                if v_str.find('x') & v_str.find('y') == -1:
                    v_str = '(' + str(v_str) + ')* np.ones(np.shape(dyg))'
                else:
                    v_str = v_str.replace('x', 'dxg')
                    v_str = v_str.replace('y', 'dyg')
                    
                # Generate arrays for the components of the zoom field
                u_zoom = eval(u_str)
                v_zoom = eval(v_str)
                
                if pass_on_figure is False:
                    # New figure
                    zoom_vf = vector_field(dxg, dyg, u_zoom, v_zoom)
                elif pass_on_figure is True:
                    if self.subplots is False:
                        zoom_vf = vector_field(dxg, dyg, u_zoom, v_zoom, fig=self.figure, subplots=False)
                    elif self.subplots is True:
                        zoom_vf = vector_field(dxg, dyg, u_zoom, v_zoom, fig=self.figure, subplots=True, sub_axis_list=self.axis)
                else:
                    raise ValueError('Error, Incorrect input for \' pass_on_figure \'')
                    
                return zoom_vf
            
        def zoom_inset(self, target=[0,0], zoom=2, dpd=9):
            '''
            Zoom, but display in an inset plot
            '''
            
            # Requires user to provide eqn of the 1-form they are zooming on.
            
            if self.str_x == None or self.str_y == None:
                # ERROR
                raise TypeError('Error: No equation provided')
            else:
                
                # Target coordinates
                x_m = target[0]
                y_m = target[1]
                
                # Get the size of the original VF
                L = 0.5*(self.xg[0, -1] - self.xg[0, 0])
                
                # Range of the points in the new zoomed grid
                # Assumes the original VF is centred on (0,0)
                d_range = self.xg[0, -1]/zoom
                
                # Size of the inset plot (default as 0.3)
                d_length = 0.3

                # Set up zoom window grids
                dx = np.linspace(-d_range + x_m, d_range + x_m, dpd)
                dy = np.linspace(-d_range + y_m, d_range + y_m, dpd)
                dxg, dyg = np.meshgrid(dx, dy)
                
                # Create variables for the user provided equation strings
                u_str = self.str_x
                v_str = self.str_y
                
                # Check if the equations provided contain x and y terms
                if u_str.find('x') & u_str.find('y') == -1:
                    u_str = '(' + str(u_str) + ')* np.ones(np.shape(dxg))'
                else:
                    u_str = u_str.replace('x', 'dxg')
                    u_str = u_str.replace('y', 'dyg')
          
                if v_str.find('x') & v_str.find('y') == -1:
                    v_str = '(' + str(v_str) + ')* np.ones(np.shape(dyg))'
                else:
                    v_str = v_str.replace('x', 'dxg')
                    v_str = v_str.replace('y', 'dyg')
                    
                # Generate arrays for the components of the zoom field
                u_zoom = eval(u_str)
                v_zoom = eval(v_str)
                
                # Specify the position and the size of the inset plot
                # [x0, y0, width, height], (x0,y0) are coords of the lower left corner
                # Given as fraction of the main plot size. (e.g. 0,0,0.5,0.5 covers the bottom left quadrant)
                q = 0.92
                zoom_inset_ax = self.axis.inset_axes([0.5*(1 + q*x_m/L - d_length), 0.5*(1 + q*y_m/L - d_length), d_length, d_length])
                
                # if pass_on_figure is False:
                #     # New figure
                #     zoom_vf = vector_field(dxg, dyg, u_zoom, v_zoom)
                # elif pass_on_figure is True:
                #     if self.subplots is False:
                #         zoom_vf = vector_field(dxg, dyg, u_zoom, v_zoom, fig=self.figure, subplots=False)
                #     elif self.subplots is True:
                #         zoom_vf = vector_field(dxg, dyg, u_zoom, v_zoom, fig=self.figure, subplots=True, sub_axis_list=self.axis)
                # else:
                #     raise ValueError('Error, Incorrect input for \' pass_on_figure \'')
                
                
                
                # return zoom_vf
            
        def DF(self, target=[0,0], zoom=2, dpd=9, pass_on_figure=False):
            '''
            Creates new vector field object at a target location, showing the derivative field at this point.
            User gives arguments:
            Target - derivative plot location
            Zoom - Magnification level
            dpd - New plot point density
            '''
            if self.str_x == None or self.str_y == None:
                # ERROR
                raise TypeError('Error: No equation provided')
            else:
                
                # Target coordinates
                x_m = target[0]
                y_m = target[1]
                
                d_range = self.xg[0, -1]/zoom
                #d_length = 
                
                # Set up zoom window grids
                dx = np.linspace(-d_range + x_m, d_range + x_m, dpd)
                dy = np.linspace(-d_range + y_m, d_range + y_m, dpd)
                dxg, dyg = np.meshgrid(dx, dy)
                
                # Create variables for the user provided equation strings
                u_str = self.str_x
                v_str = self.str_y

                # Create string to evaluate the field at the target location
                u_str_point = u_str.replace('x', 'x_m')
                u_str_point = u_str_point.replace('y', 'y_m')
                
                v_str_point = v_str.replace('x', 'x_m')
                v_str_point = v_str_point.replace('y', 'y_m')
                
                # Check if the equations provided contain x and y terms
                if u_str.find('x') & u_str.find('y') == -1:
                    u_str_grid = '(' + str(u_str) + ')* np.ones(np.shape(dxg))'
                else:
                    u_str_grid = u_str.replace('x', 'dxg')
                    u_str_grid = u_str_grid.replace('y', 'dyg')
          
                if v_str.find('x') & v_str.find('y') == -1:
                    v_str_grid = '(' + str(v_str) + ')* np.ones(np.shape(dyg))'
                else:
                    v_str_grid = v_str.replace('x', 'dxg')
                    v_str_grid = v_str_grid.replace('y', 'dyg')
                    
                # Generate arrays for the components of the derivative field          
                U = eval(u_str_grid) - eval(u_str_point)
                V = eval(v_str_grid) - eval(v_str_point)
                
                # Return object to the user
                if pass_on_figure is False:
                    # New figure
                    DF_vf = vector_field(dxg, dyg, U, V)
                elif pass_on_figure is True:
                    if self.subplots is False:
                        DF_vf = vector_field(dxg, dyg, U, V, fig=self.figure, subplots=False)
                    elif self.subplots is True:
                        DF_vf = vector_field(dxg, dyg, U, V, fig=self.figure, subplots=True, sub_axis_list=self.axis)
                else:
                    raise ValueError('Error, Incorrect input for \' pass_on_figure \'')
                
                return DF_vf
            
        def Div(self, target=[0,0], zoom=2, dpd=9, pass_on_figure=False):
            '''
            Creates new vector field object at a target location, showing the Divergence of the field at this point.
            User gives arguments:
            Target - derivative plot location
            Zoom - Magnification level
            dpd - New plot point density
            '''
            if self.str_x == None or self.str_y == None:
                # ERROR
                raise TypeError('Error: No equation provided')
            else:
                
                # Target coordinates
                x_m = target[0]
                y_m = target[1]
                
                d_range = self.xg[0, -1]/zoom
                #d_length = 
                
                # Set up zoom window grids
                dx = np.linspace(-d_range + x_m, d_range + x_m, dpd)
                dy = np.linspace(-d_range + y_m, d_range + y_m, dpd)
                dxg, dyg = np.meshgrid(dx, dy)
                
                # Create variables for the user provided equation strings
                u_str = self.str_x
                v_str = self.str_y

                # Create string to evaluate the field at the target location
                u_str_point = u_str.replace('x', 'x_m')
                u_str_point = u_str_point.replace('y', 'y_m')
                
                v_str_point = v_str.replace('x', 'x_m')
                v_str_point = v_str_point.replace('y', 'y_m')
                
                # Check if the equations provided contain x and y terms
                if u_str.find('x') & u_str.find('y') == -1:
                    u_str_grid = '(' + str(u_str) + ')* np.ones(np.shape(dxg))'
                else:
                    u_str_grid = u_str.replace('x', 'dxg')
                    u_str_grid = u_str_grid.replace('y', 'dyg')
          
                if v_str.find('x') & v_str.find('y') == -1:
                    v_str_grid = '(' + str(v_str) + ')* np.ones(np.shape(dyg))'
                else:
                    v_str_grid = v_str.replace('x', 'dxg')
                    v_str_grid = v_str_grid.replace('y', 'dyg')
                    
                # Generate arrays for the components of the derivative field          
                U = eval(u_str_grid) - eval(u_str_point)
                V = eval(v_str_grid) - eval(v_str_point)
                
                # =============================================================================
                # Geometric Divergence Method - See Documentation                
                # =============================================================================
                
                U_div = np.zeros(shape=(dpd, dpd))
                V_div = np.zeros(shape=(dpd, dpd))
                
                # Looping Constant
                N = dpd - 1
        
                # get number of points in quadrant
                if dpd % 2 == 1:
                    quad_x = int(dpd/2)
                    quad_y = int((dpd+1)/2)
                else:
                    quad_x = int(dpd/2)
                    quad_y = int(dpd/2)
                    
                for i in range(quad_x):
                    # get the l number, for projection of j on radial / i on tangent
                    l = i - 0.5*N
                    
                    # INNER LOOP
                    for j in range(quad_y):
                        # get the k number of projection: i on radial / j on tangent
                        k = j - 0.5*N
                        
                        # get the commuting parts of V and W for each square corner
                        # (x and y components of the subtracted field)
                        U_comm_1 = 0.25*(2*U[i, j] + V[j, N-i] - V[N-j, i])
                        U_comm_2 = 0.25*(2*U[j, N-i] + V[N-i, N-j] - V[i, j])
                        U_comm_3 = 0.25*(2*U[N-i, N-j] + V[N-j, i] - V[j, N-i])
                        U_comm_4 = 0.25*(2*U[N-j, i] + V[i, j] - V[N-i, N-j])
                        
                        V_comm_1 = 0.25*(2*V[i, j] - U[j, N-i] + U[N-j, i])
                        V_comm_2 = 0.25*(2*V[j, N-i] - U[N-i, N-j] + U[i, j])
                        V_comm_3 = 0.25*(2*V[N-i, N-j] - U[N-j, i] + U[j, N-i])
                        V_comm_4 = 0.25*(2*V[N-j, i] - U[i, j] + U[N-i, N-j])
                        
                        # gte a normalisation factor from l and k
                        A = k**2 + l**2
                        
                        U_div[i, j] = (U_comm_1*k + V_comm_1*l)*k/A
                        V_div[i, j] = (U_comm_1*k + V_comm_1*l)*l/A
                        U_div[j, N-i] = (U_comm_2*l + V_comm_2*(-k))*l/A
                        V_div[j, N-i] = (U_comm_2*l + V_comm_2*(-k))*(-k)/A
                        U_div[N-i, N-j] = (U_comm_3*(-k) + V_comm_3*(-l))*(-k)/A
                        V_div[N-i, N-j] = (U_comm_3*(-k) + V_comm_3*(-l))*(-l)/A
                        U_div[N-j, i] = (U_comm_4*(-l) + V_comm_4*k)*(-l)/A
                        V_div[N-j, i] = (U_comm_4*(-l) + V_comm_4*k)*k/A
                
                # Return object to the user
                if pass_on_figure is False:
                    # New figure
                    Div_vf = vector_field(dxg, dyg, U_div, V_div)
                elif pass_on_figure is True:
                    if self.subplots is False:
                        Div_vf = vector_field(dxg, dyg, U_div, V_div, fig=self.figure, subplots=False)
                    elif self.subplots is True:
                        Div_vf = vector_field(dxg, dyg, U_div, V_div, fig=self.figure, subplots=True, sub_axis_list=self.axis)
                else:
                    raise ValueError('Error, Incorrect input for \' pass_on_figure \'')
                
                return Div_vf
            
        def Curl(self, target=[0,0], zoom=2, dpd=9, pass_on_figure=False):
            '''
            Creates new vector field object at a target location, showing local rotation (Curl)
            User gives arguments:
            Target - derivative plot location
            Zoom - Magnification level
            dpd - New plot point density
            '''
            if self.str_x == None or self.str_y == None:
                # ERROR
                raise TypeError('Error: No equation provided')
            else:
                
                # Target coordinates
                x_m = target[0]
                y_m = target[1]
                
                d_range = self.xg[0, -1]/zoom
                #d_length = 
                
                # Set up zoom window grids
                dx = np.linspace(-d_range + x_m, d_range + x_m, dpd)
                dy = np.linspace(-d_range + y_m, d_range + y_m, dpd)
                dxg, dyg = np.meshgrid(dx, dy)
                
                # Create variables for the user provided equation strings
                u_str = self.str_x
                v_str = self.str_y

                # Create string to evaluate the field at the target location
                u_str_point = u_str.replace('x', 'x_m')
                u_str_point = u_str_point.replace('y', 'y_m')
                
                v_str_point = v_str.replace('x', 'x_m')
                v_str_point = v_str_point.replace('y', 'y_m')
                
                # Check if the equations provided contain x and y terms
                if u_str.find('x') & u_str.find('y') == -1:
                    u_str_grid = '(' + str(u_str) + ')* np.ones(np.shape(dxg))'
                else:
                    u_str_grid = u_str.replace('x', 'dxg')
                    u_str_grid = u_str_grid.replace('y', 'dyg')
          
                if v_str.find('x') & v_str.find('y') == -1:
                    v_str_grid = '(' + str(v_str) + ')* np.ones(np.shape(dyg))'
                else:
                    v_str_grid = v_str.replace('x', 'dxg')
                    v_str_grid = v_str_grid.replace('y', 'dyg')
                    
                # Generate arrays for the components of the derivative field          
                U = eval(u_str_grid) - eval(u_str_point)
                V = eval(v_str_grid) - eval(v_str_point)
                
                # =============================================================================
                # Geometric Curl Method - See Documentation                
                # =============================================================================
                
                U_curl = np.zeros(shape=(dpd, dpd))
                V_curl = np.zeros(shape=(dpd, dpd))
                
                # Looping Constant
                N = dpd - 1
        
                # Quadrant Points
                if dpd % 2 == 1:
                    quad_x = int(dpd/2)
                    quad_y = int((dpd+1)/2)
                else:
                    quad_x = int(dpd/2)
                    quad_y = int(dpd/2)
                    
                for i in range(quad_x):
                    # get the l number, for projection of j on radial / i on tangent
                    l = i - 0.5*N
                    
                    # INNER LOOP
                    for j in range(quad_y):
                        # get the k number of projection: i on radial / j on tangent
                        k = j - 0.5*N
                        
                        # get the commuting parts of V and W for each square corner
                        # (x and y components of the subtracted field)
                        U_comm_1 = 0.25*(2*U[i, j] + V[j, N-i] - V[N-j, i])
                        U_comm_2 = 0.25*(2*U[j, N-i] + V[N-i, N-j] - V[i, j])
                        U_comm_3 = 0.25*(2*U[N-i, N-j] + V[N-j, i] - V[j, N-i])
                        U_comm_4 = 0.25*(2*U[N-j, i] + V[i, j] - V[N-i, N-j])
                        
                        V_comm_1 = 0.25*(2*V[i, j] - U[j, N-i] + U[N-j, i])
                        V_comm_2 = 0.25*(2*V[j, N-i] - U[N-i, N-j] + U[i, j])
                        V_comm_3 = 0.25*(2*V[N-i, N-j] - U[N-j, i] + U[j, N-i])
                        V_comm_4 = 0.25*(2*V[N-j, i] - U[i, j] + U[N-i, N-j])
                        
                        # gte a normalisation factor from l and k
                        A = k**2 + l**2
                        
                        U_curl[i, j] = (U_comm_1*l + V_comm_1*(-k))*l/A
                        V_curl[i, j] = (U_comm_1*l + V_comm_1*(-k))*(-k)/A
                        U_curl[j, N-i] = (U_comm_2*(-k) + V_comm_2*(-l))*(-k)/A
                        V_curl[j, N-i] = (U_comm_2*(-k) + V_comm_2*(-l))*(-l)/A
                        U_curl[N-i, N-j] = (U_comm_3*(-l) + V_comm_3*k)*(-l)/A
                        V_curl[N-i, N-j] = (U_comm_3*(-l) + V_comm_3*k)*k/A
                        U_curl[N-j, i] = (U_comm_4*k + V_comm_4*l)*k/A
                        V_curl[N-j, i] = (U_comm_4*k + V_comm_4*l)*l/A
                    
                # Return object to the user
                if pass_on_figure is False:
                    # New figure
                    Curl_vf = vector_field(dxg, dyg, U_curl, V_curl)
                elif pass_on_figure is True:
                    if self.subplots is False:
                        Curl_vf = vector_field(dxg, dyg, U_curl, V_curl, fig=self.figure, subplots=False)
                    elif self.subplots is True:
                        Curl_vf = vector_field(dxg, dyg, U_curl, V_curl, fig=self.figure, subplots=True, sub_axis_list=self.axis)
                else:
                    raise ValueError('Error, Incorrect input for \' pass_on_figure \'')
                
                return Curl_vf
        
        # define a method to change a supplied Vector filed to the 1-form
        def formalise(self, pass_on_figure=False, g=[['1', '0'], ['0', '1']]):
            '''
            Passes in everything it can (all it has been supplied)
            to the 1-form object.
            Works via the metric on R2
            Can supply the metric in as equations or as evaluated arrays
            Format of the metric is a list of numpy arrays
            0th array is the top row, its 0th component is 11, 1st is 12
            1st array is the botton row, its 0th comp is 21 and 1st is 22.
            Note, if it is supplied as arrays, they must come from numpy grids
            via meshgrid, if it is supplied as strings, needs to be in terms of
            x and y, and contain no special funtions, apart from ones imported
            here automatically and listed in the documentation #!!!
            Apart form the metric, takes in pass_on_figure:
                ----determines if figure from this object is to be
                passed on to the 1-form object.
            
            Returns a single object (1-form object)
            '''
            
            # extract what is needed form the metric depending on what the user
            # supplied
            # check if its has string components
            if type(g[0][0]) == str and type(g[0][1]) == str and type(g[1][0]) == str and type(g[1][1]) == str:
                # deal with supplied string metric
                # need to format it, correct it for constants and evaluate it's numerical equivalent
                str_comp_00 = g[0][0] + ''
                str_comp_01 = g[0][1] + ''
                str_comp_10 = g[1][0] + ''
                str_comp_11 = g[1][1] + ''
                str_comp_00 = str_comp_00.replace('x', '(self.xg)')
                str_comp_00 = str_comp_00.replace('y', '(self.yg)')
                str_comp_01 = str_comp_01.replace('x', '(self.xg)')
                str_comp_01 = str_comp_01.replace('y', '(self.yg)')
                str_comp_10 = str_comp_10.replace('x', '(self.xg)')
                str_comp_10 = str_comp_10.replace('y', '(self.yg)')
                str_comp_11 = str_comp_11.replace('x', '(self.xg)')
                str_comp_11 = str_comp_11.replace('y', '(self.yg)')
                # check against constant form components:
                if str_comp_00.find('x') & str_comp_00.find('y') == -1:
                    str_comp_00 = '(' + str(str_comp_00) + ')* np.ones(np.shape(self.xg))'
                if str_comp_01.find('x') & str_comp_01.find('y') == -1:
                    str_comp_01 = '(' + str(str_comp_01) + ')* np.ones(np.shape(self.yg))'
                if str_comp_10.find('x') & str_comp_10.find('y') == -1:
                    str_comp_10 = '(' + str(str_comp_10) + ')* np.ones(np.shape(self.yg))'
                if str_comp_11.find('x') & str_comp_11.find('y') == -1:
                    str_comp_11 = '(' + str(str_comp_11) + ')* np.ones(np.shape(self.yg))'
                
                # evaluate the components numerically, inputting them into a
                # store numerical metric
                comp_00 = eval(str_comp_00)
                comp_01 = eval(str_comp_01)
                comp_10 = eval(str_comp_10)
                comp_11 = eval(str_comp_11)
                g_num = [[comp_00, comp_01], [comp_10, comp_11]]
                
                # set up a dummy variable to store the fact that numericals were given
                # not to check again later
                analytics = True
                
            elif type(g[0][0]) == np.ndarray and type(g[0][1]) == np.ndarray and type(g[1][0]) == np.ndarray and type(g[1][1]) == np.ndarray:
                # deal with the metric being supplied as components
                # if the user has vector field equations, warn that these can't
                # be passed anymore, because we don't have equations for this
                # metric
                if self.str_x == None and self.str_y == None:
                    pass
                else:
                    print('The Vector field has equations, but the metric does not, these will be lost and the resulting 1-form will only have numerical values, not equations supplied')
                # No need to do anythng more to the metric, upto the user to make sure its
                # correctly sized, as with other code in this library
                # just rename the metric here
                g_num = g
                
                # set up a dummy variable to store the fact that numericals were
                # not given, not to check again later
                analytics = False
                
            else:
                # Inconsistant metric components
                raise TypeError('Metric components are inconcisstant')
            
            # from vector field components, get 1-form components by the metric
            # first, do so numerically, as this must always happen
            form_x = self.F_x * g_num[0][0] + self.F_y * g_num[0][1]
            form_y = self.F_y * g_num[1][1] + self.F_x * g_num[1][0]
            
            # if the equations were given, evaluate these analytically too:
            # only if vector file doriginally has equations
            if analytics:
                if self.str_x == None and self.str_y == None:
                    print('You supplied the metric as equations (or it was default), but did not give VF equations, therefore only numericals will be completed')
                    analytics = False
                else:
                    x_str_form = '(' + self.str_x + ')*(' + g[0][0] + ') + (' + self.str_y + ')*(' + g[0][1] + ')'
                    y_str_form = '(' + self.str_y + ')*(' + g[1][1] + ') + (' + self.str_x + ')*(' + g[1][0] + ')'
                    # simplify them
                    x_str_form = str(simplify(x_str_form))
                    y_str_form = str(simplify(y_str_form))
            else:
                pass

            # based on what was given into the Vector field, return a 1-form object with these parameters
            if analytics:
                if pass_on_figure is False:
                    # supply these to the 1-form object creator
                    result_form = form_1(self.xg, self.yg, form_x, form_y, x_str_form, y_str_form)
                elif pass_on_figure is True:
                    if self.subplots is False:
                        result_form = form_1(self.xg, self.yg, form_x, form_y, x_str_form, y_str_form, fig=self.figure, subplots=False)
                    elif self.subplots is True:
                        result_form = form_1(self.xg, self.yg, form_x, form_y, x_str_form, y_str_form, fig=self.figure, subplots=True, sub_axis_list=self.axis)
                else:
                    raise ValueError('Error, Incorrect input for \' pass_on_figure \'')
            elif not analytics:
                if pass_on_figure is False:
                    # supply these to the 1-form object creator
                    result_form = form_1(self.xg, self.yg, form_x, form_y)
                elif pass_on_figure is True:
                    if self.subplots is False:
                        result_form = form_1(self.xg, self.yg, form_x, form_y, fig=self.figure, subplots=False)
                    elif self.subplots is True:
                        result_form = form_1(self.xg, self.yg, form_x, form_y, fig=self.figure, subplots=True, sub_axis_list=self.axis)
                else:
                    raise ValueError('Error, Incorrect input for \' pass_on_figure \'')
        
            # return the found object
            return result_form

    # now call that object to create it:
    v_object = field_set_up(xg, yg, F_x, F_y)
    # return it to user to store
    return v_object

