# Differential form python module attempt - 1

# import needed modules
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sympy import diff, simplify
from sympy.parsing.sympy_parser import parse_expr
from math import isnan
from matplotlib import patches as patch
from matplotlib import cm

# input many numpy functions to deal with user input
from numpy import sin, cos, tan, sqrt, log, arctan, arcsin, arccos, tanh
from numpy import sinh, cosh, arcsinh, arccosh, arctanh, exp, pi, e

# define function that sets the recursion constant for the loop to plot stacks
# pre-define the displacements from mid point needed
# c is the Truth value from parity (odd or even number n)
def G(s, n, c):
    '''
    G(s, n, c)
    
    Defines coefficints needed to displace stack sheets along direction perp.
    to form, depending  on how many are to be plotted.
    
    Parameters:
    --------
    s - det. number of sheets to draw
    n - which sheet is sequence one is to be drawn now
    c - int bool, as 0 or 1, defines parity of n
    
    Returns:
    --------
    Coefficient to fractional sheet displacement.
    
    '''
    if c == 0:
        return ((2*s + 1)/(2*(n-1)))
    else:
        return (s/(n-1))


# define a function that will analytically find the 2-form from given expressions
# in a given number of dimensions and in terms of given coordinate symbols
def find_2_form(expressions, coords, xg, yg, zg=None, m=2):
    '''
    find_2_form(expressions, coords, xg, yg, zg=None, m=2)
    
    Finds the analytical 2 form using sympy experssion handling.
    
    Parameters:
    ---------------
    expressions - list of sympy experssions for the 1 form scaling fucntions
    coords - list of coordinate names as strings, that were used in experssions
    xg, yg - grids
    zg - possible grid
    m - number of dimensions
    
    Returns:
    ---------------
    result  - analytical, unformatted 2-form equation
    '''
    
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
    
    return result
    
#    # create a local result, that will be used to evaluate the resulting string
#    # and format it
#    loc_res = result + ''
#    
#    # format string in each result row
#    for d in range(pair):
#        # format the result to be 'python understood' to be able to use the eval()
#        loc_res[d, 0] = loc_res[d, 0].replace('x', 'xg')
#        loc_res[d, 0] = loc_res[d, 0].replace('y', 'yg')
#        loc_res[d, 0] = loc_res[d, 0].replace('z', 'zg')
#        
#        # check against constant result, to be of correct shape before eval is used
#        if loc_res[d, 0].find('x') & loc_res[d, 0].find('y') == -1:
#            loc_res[d, 0] = '(' + str(loc_res[d, 0]) + ')* np.ones(np.shape(xg))'
#        if loc_res[d, 0].find('x') & loc_res[d, 0].find('y') == -1:
#            loc_res[d, 0] = '(' + str(loc_res[d, 0]) + ')* np.ones(np.shape(yg))'
#
#    # set up a vector to store the 2-form numerically, from xg and yg and possibly further
#    # Note - need pt_den being supplied m times.
#    # not overall generalised, as not needed past m=3.
#    if m == 2:
#        form_2 = np.empty((1, pt_den, pt_den))
#        form_2[0, :, :] = eval(loc_res[0, 0])
#    elif m == 3:
#        form_2 = np.empty((3, pt_den, pt_den, pt_den))
#        for d in range(3):
#            form_2[d, :, :, :] = eval(loc_res[d, 0])
#    
#    # return useful findings to the user
#    return form_2, result, ext_ds


# %%

'''

function to create a 1-form object and define methods for it

'''

# define the 1-form object and all its methods
class form_1():
    '''
    form_1(xg, yg, F_x, F_y, F_x_eqn=None, F_y_eqn=None)
    
    Defines a 1-form object and returns it to user. 
    
    Parameters:
    ---------------
    xg - grid of x values (2D numpy.ndarray)
    yg - grid of y values (2D numpy.ndarray)
    F_x - grid of dx form components (2D numpy.ndarray)
    F_y - grid of dy form components (2D numpy.ndarray)
    
    Optional:
    F_x_eqn - expression for dx form component f(x,y) (string)
    F_y_eqn - expression for dy form component f(x,y) (string)
    
    
    Instance variables:
    ---------------
    xg, yg, F_x, F_y
    s_max - int - maximum number of sheets per stack
    s_min - int - minimum number of sheets per stack
    pt_den - int - number of points on grids, extracted from grids, assumes square grid
    fract - float/int - length of sheet in stack as fraction of whole plot size
    scale - float/int - constant multpilier to change scaling
    w_head - float/int - width of arrowghead on stack as size of sheet
    h_head - float/int - height of arrowghead on stack as size of sheet
    arrowheads - bool - determines of arrowheads showld be drawn on stacks
    color - str - colour to draw stacks with, can be Hex when using '#FFFFFF'
    logarithmic_scale_bool - bool - determines if log scaling is used
    delta_factor - float/int - determined size of blank boarder in figure
                                as fraction of whole plot size
    
    Methods:
    ---------------
    give_eqn
    return_string
    colour
    arrow_heads
    head_width
    head_height
    log_scaling
    max_sheets
    sheet_size
    surround_space
    set_density
    plot
    ext_d
    num_ext_d
    hodge
    wedge_analytical
    wedge_num
    zoom 
    interior_d
    contravariant
    
    '''
    
    def __init__(self, xg, yg, F_x, F_y, F_x_eqn=None, F_y_eqn=None):
        self.xg = xg
        self.yg = yg
        self.F_x = F_x
        self.F_y = F_y
        self.s_max = 6
        self.s_min = 1
        self.fract = 0.05
        self.scale = 1
        self.w_head = 1/8
        self.h_head = 1/4
        self.arrowheads = True
        self.color = '#8B14F3'
        self.logarithmic_scale_bool = 0
        self.delta_factor = 10
        # define equations if given:
        # user must change to access some methods, will indicate when needed
        # Note, the string must be given with x and y as variables
        if F_x_eqn is not None:
            self.form_1_str_x = str(simplify(F_x_eqn))
        else:
            self.form_1_str_x = None
        
        if F_y_eqn is not None:
            self.form_1_str_y = str(simplify(F_y_eqn))
        else:
            self.form_1_str_y = None
    
    # #####################################################################
    # write customising methods
    # #####################################################################
    
    # define a mehtod to allow user to supply the string equation
    # of the 1-form
    def give_eqn(self, equation_str_x, equation_str_y):
        '''
        give_eqn(equation_str_x, equation_str_y)
        
        This must be the equation of the supplied numerical 1-form
        in terms of variables x and y.
        All formatting is as in numpy, but no library calling in string
        exception - exponential - call it as e**(expression)
        
        it re-evaluates the numerical values to match the new equations
        A warning is shown if any differences are detected, not rigorous though
        Will often show for most minor changes
        
        Has to be given, for some methods to be computable
        Methods will indicate when needed
        
        Parameters:
        ---------------
        equation_str_x - string of the dx component, with x and y as variables
        equation_str_y - string of the dy component, with x and y as variables
        
        Returns: None
        
        '''
        # set equation parameters to simplified inputs
        self.form_1_str_x = str(simplify(equation_str_x))
        self.form_1_str_y = str(simplify(equation_str_y))
        # make the values match automatically to limit how often mismatch occurs
        # substitute these into the equation, but keep it local: 
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
        
        # re-evaluate the 2-form numerically, warn user if changed
        if not ((self.F_x is eval(str_x)) and (self.F_y is eval(str_y))):
            print('Warning: Equations did not exactly match numerical values, and these were change to agree with equations')
        
        # evaluate formatted equations and save
        self.F_x = eval(str_x)
        self.F_y = eval(str_y)
    
    # deifne a function to return the string equations to the user
    def return_string(self):
        '''
        Returns unformatted strings for component equations back to user
        Done in case user wants to access strings that got here by ext. alg.
        
        Parmateres: None
        Returns: None
        
        '''
        return self.form_1_str_x, self.form_1_str_y
    
    # change colour
    def colour(self, color):
        '''
        Changes the colour that stacks plot in
        Note, not strictly needed, can change it by instance.color(colour)
        
        Parmaeters:
        -------------
        color - string - string to define a colour of stacks
                    can be any matplotlib understood colour or Hex in #FFFFFF
        
        Returns: None
        '''
        self.color = str(color)
    
    # change arrowsheads
    def arrow_heads(self):
        '''
        Changes the boolean that determines if arrowheads are plotted
        on stacks. Whenever it is called, it changes that boolean to opposite
        The form object is initialised with this as True
        Note, not strictly needed, can change it by instance.arrowheads(bool)
        
        Parmaeters: None
        Returns: None
        '''
        self.arrowheads = not self.arrowheads
    
    # change w_head
    def head_width(self, wide):
        '''
        Sets the width of the arrowhead on a stacks to the desired float
        as a fraction of the stack length in the direction perp. to form
        Note, not strictly needed, can change it by instance.w_head(width)
        
        Parmaeters:
        ---------------
        wide - float/int - Sets the width as a fraction of the stack length
        
        Returns: None
        
        '''
        self.w_head = float(wide)
    
    # change h_head
    def head_height(self, high):
        '''
        Sets the height of the arrowhead on a stacks to the desired float
        as a fraction of the stack length in the direction parall. to form
        Note, not strictly needed, can change it by instance.h_head(height)
        
        Parmaeters:
        ---------------
        high - float/int - Sets the height as a fraction of the stack length
        
        Returns: None
        
        '''
        self.h_head = float(high)
    
    # change boolean that det. if to sclae logarithmically
    def log_scaling(self):
        '''
        Changes the boolean that determines if scaling is logarithmic
        Whenever it is called, it changes that boolean to opposite
        The form object is initialised with this as False
        Note, not strictly needed, can change it by instance.logarithmic_scale_bool(bool)
        
        Parmaeters: None
        Returns: None
        '''
        self.logarithmic_scale_bool = not self.logarithmic_scale_bool
        # self.base = base
    
    # define methods to change s_max
    def max_sheets(self, maximum):
        '''
        Changes maximum number of sheets to draw on a stack.
        These still scale relative to max magnitude.
        Note, not strictly needed, can change it by instance.s_max(maximum)
        
        Parmaeters:
        ---------------
        maximum - int - Max number of sheets to plot per stack
        
        Returns: None
        
        '''
        self.s_max = maximum
    
    # define method to change fraction of sheetsize w.r.t graph size:
    def sheet_size(self, fraction):
        '''
        Changes the size of stack in direction perp. to form.
        It is done in in terms of the fraction of plot size
        Note, not strictly needed, can change it by instance.fract(fraction)
        
        Parmaeters:
        ---------------
        fraction - float/int - size of stack in terms of the fraction of plot size
        
        Returns: None
        '''
        self.fract = fraction
    
    # define a method to change spare spacing around figure
    def surround_space(self, delta_denominator):
        '''
        Sets the extra blank space around the domain of grids in axis
        Note, not strictly needed, can change it by instance.delta_factor(delta_denominator)
        
        Parmaeters:
        ---------------
        delta_denominator - float/int - denominator or fraction to use
            eg. supplying 3 will make the white space 1/3 of the width
                of the domain of the grid.
        
        Returns: None
        
        '''
        self.delta_factor = delta_denominator
    
    # define a method to change the density of grids in same range
    # requires string input of 1-form:
    def set_density(self, points_number):
        '''
        Changes the desnity of points in the same range to the input value
        Requires the string equation to be supplied to not 'extrapolate'
        
        Only creates 2 axis with same number of points each
        cannot be used for any custom grids
        
        Parameters:
        --------------
        points_number - new number of points to use per axis
        
        Returns: None
        '''
        if self.form_1_str_x == None or self.form_1_str_y == None:
            # Error
            raise ValueError('Error: You need to supply the 1-form equation to do this, see \'give_eqn\' method')
        else:
            # redefine the grids
            x = np.linspace(self.xg[0,0], self.xg[0, -1], points_number)
            y = np.linspace(self.yg[0,0], self.yg[-1, 0], points_number)
            self.xg, self.yg = np.meshgrid(x, y)
            # substitute these into the equation, but keep it local:
            str_x = self.form_1_str_x + ''
            str_y = self.form_1_str_y + ''
            str_x = str_x.replace('x', '(self.xg)')
            str_x = str_x.replace('y', '(self.yg)')
            str_y = str_y.replace('x', '(self.xg)')
            str_y = str_y.replace('y', '(self.yg)')
            
            # check against constant forms, to have correct array shape
            if str_x.find('x') & str_x.find('y') == -1:
                str_x = '(' + str(str_x) + ')* np.ones(np.shape(self.xg))'
            if str_y.find('x') & str_y.find('y') == -1:
                str_y = '(' + str(str_y) + ')* np.ones(np.shape(self.yg))'
        
            # re-evaluate the 1-form numerically
            self.F_x = eval(str_x)
            self.F_y = eval(str_y)
    
    # #####################################################################
    # More useful methods (plotting, zooming and ext. alg.)
    # #####################################################################
    
    # define a fucntion that will use the set up 1-form and plot it
    def plot(self, axis):
        '''
        
        plot(axis)
        
        Uses the attribues of the object as set originally and as customised
        with methods to create a plot of the 1-form
        
        Parameters:
        -------------
        axis: matplotlib axes that the plot it to be put on
        
        Returns: None
        '''
        
        # get the lengths of x and y from their grids
        x_len = len(self.xg[:, 0])
        y_len = len(self.yg[0, :])
        
        # Extract L from the x and y grids
        Lx = 0.5*(self.xg[0, -1] - self.xg[0, 0])
        Ly = 0.5*(self.yg[-1, 0] - self.yg[0, 0])
        L = 0.5*(Lx + Ly)  # average, needed for stack sizes only
        x0 = self.xg[0, 0] + Lx
        y0 = self.yg[0, 0] + Ly
        
        # reset axis limits
        ax_Lx = Lx + Lx/self.delta_factor
        ax_Ly = Ly + Ly/self.delta_factor
        axis.set_xlim(-ax_Lx + x0, ax_Lx + x0)
        axis.set_ylim(-ax_Ly + y0, ax_Ly + y0)
        
        # find the distance between neightbouring points on the grid
        # for drawing extra arefacts
        dist_points = self.xg[0, 1] - self.xg[0, 0]
        
        # define an empty array of magnitudes, to then fill with integer rel. mags
        R_int = np.zeros(shape=((x_len), (y_len)))
    
        # #########################################################################
        # get variables needed for the initial stack plot
        # #########################################################################
        
        # set all insignificant values to zero:
        self.F_x[np.abs(self.F_x) < 1e-15] = 0
        self.F_x[np.abs(self.F_x) < 1e-15] = 0
        
        # find the arrow length corresponding to each point and store in mag array
        mag = np.sqrt(self.F_x**2 + self.F_y**2)
        
        # find direction of each arrow
        angles = np.arctan2(self.F_y, self.F_x)   # theta defined from positive x axis ccw
        
        # find regions ON GRID that are nan or inf as a bool array
        
        # deal with infs and nans in mag
        # set to zero points that are not defined or inf
        # and mark them on axis
        isnan_arr = np.isnan(mag)
        for i in range(x_len):
            for j in range(y_len):
                if isnan_arr[i, j]:
                    # colour this region as a shaded square
                    rect = patch.Rectangle((self.xg[i, j] - dist_points/2, self.yg[i, j]  - dist_points/2), dist_points, dist_points, color='#B5B5B5')
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
        
        # define length of sheet as a fraction of total graph scale
        # this also sets max, total height of stack (along its direction)
        s_L = self.fract * L
        
        # #########################################################################
        # define stack based on geometrical arguments
        # sheets perp. to hypothetical arrow, shifted along it
        # their density porp to mag, + arrowhead on top
        # #########################################################################
        
        # find the maximum magnitude for scaling
        max_size = np.max(mag)
        
        # setrelative scaling, linear or logarithmic
        if self.logarithmic_scale_bool:
            mag1 = mag + 1
            logmag1 = np.log(mag1)
            R = logmag1/np.max(logmag1)  # Re-assign R
        else:
            R = mag/max_size

        # define tigonometirc shifts
        I_sin = np.sin(angles)
        I_cos = np.cos(angles)
        
        # precalculate heavy operations
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
        
        #  special case, when there is only 1 line in the stack plot:
        P_sh1x = self.xg + (s_L*self.w_head)*I_sin
        P_sh1y = self.yg - (s_L*self.w_head)*I_cos
        P_sh2x = self.xg - (s_L*self.w_head)*I_sin
        P_sh2y = self.yg + (s_L*self.w_head)*I_cos
        P_sh3x = self.xg + (s_L*self.h_head)*I_cos
        P_sh3y = self.yg + (s_L*self.h_head)*I_sin
        
        # array of number of sheets for each stack
        for i in range(self.s_max - self.s_min + 1):
            t = self.s_max - i
            R_int[R <= t/self.s_max] = t
        
        # loop over each coordinate plotting
        for i in range(x_len):
            for j in range(y_len):
                # varible for current considered magnitude as it is reused
                # avoids extracting from R many times.
                n = R_int[i, j]
                
                # do not plot anything if magnitude is exactly zero
                if mag[i,j] == 0:
                    continue
                
                # deal with even number of sheets from magnitudes:
                if n % 2 == 0:
                    # parameter to loop over in the recursion equation
                    s = 0
                    
                    # points for sheets required for the given magnitude
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
                else:
                    # Add the centre line for odd numbers of stacks
                    axis.add_line(Line2D((A_x[i, j], B_x[i, j]), (A_y[i, j], B_y[i, j]), linewidth=1, color=self.color))
                    
                    # then loop over the remaining lines as per the recursion formula:
                    s = 1  # exclude already completed 0
                    
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
                        
                        # update parameter
                        s += 1
                
                # dela with arrowheads
                if self.arrowheads:
                    # from central sheet for n=1 or on top sheet for n>1 
                    if n > 1:   # for all lines but the single sheet one
                        axis.add_line(Line2D((p_sh1x[i, j],p_sh3x[i, j]),(p_sh1y[i, j],p_sh3y[i, j]), linewidth=1, color = self.color))
                        axis.add_line(Line2D((p_sh2x[i, j],p_sh3x[i, j]),((p_sh2y[i, j],p_sh3y[i, j])), linewidth=1, color = self.color))
                    else:
                        # when only 1-sheet is drawn
                        axis.add_line(Line2D((P_sh1x[i, j], P_sh3x[i, j]), (P_sh1y[i, j], P_sh3y[i, j]), linewidth=1, color = self.color))
                        axis.add_line(Line2D((P_sh2x[i, j], P_sh3x[i, j]), ((P_sh2y[i, j], P_sh3y[i, j])), linewidth=1, color = self.color))
                else:
                    pass
    
    # method to find its exterior derivative
    def ext_d(self):
        '''
        
        ext_d()
        
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
            
            # from these get the 2-form
            result = find_2_form(expressions, coords, self.xg, self.yg, zg=None, m=m)
            
            # format, and evaluate
            
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
            
            # evaluate, set up new object and return
            form_2_result = eval(form_2_str)
            result_form = form_2(self.xg, self.yg, form_2_result, form_2_str_loc)
            
            # return it to the user
            return result_form
    
    # define a funciton to complete numerical only curl
    def num_ext_d(self):
        '''
        Takes in no arguments
        
        computes the exterior derivative numerically only
        The equations do not need to be given
        If given, they do not get passed onto the 2-form object anyway
        NUMERICAL ONLY, they will be lost!
         
        returns 2-form object
        '''
        
        # get steps in dx and dy:
        dx = self.xg[0, :]
        dy = self.yg[:, 0]
        
        # copy F_x and F_y, locally
        fx = self.F_x + np.zeros(np.shape(self.xg))
        fy = self.F_y + np.zeros(np.shape(self.xg))
        
        # clean up F_x and F_y from nan etc
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
        
        # from these, get the 2-form
        form_2_result = dx_F_y - dy_F_x
        
        # return 2-form object to user
        result_form = form_2(self.xg, self.yg, form_2_result)
        
        # return it to the user
        return result_form
    
    # define a method to Hodge it
    def hodge(self, keep_object=False):
        '''
        
        hodge(keep_object=False)
        
        Parameters:
        -------------
        keep_object - determines if the result should be returned as a new
                      1-form or if current one need to be changed.
                      Default is False. When False, a new object is created
                      When true, the acted on is modified.
        
        It calulates the Hodge on R^2 by the standard definition:
        dx -> dy and dy -> -dx
        Does no analytically using the equations provided in the instance
        
        returns: 1-form if keep_object is False, else returns nothing
        '''
        
        # check for equations:
        if self.form_1_str_x == None or self.form_1_str_y == None:
            # ERROR
            raise TypeError('Error: You need to supply the 1-form equation to do this, look at \'give_eqn\' method')
        else:
            # some equations are there, compute the Hodge on these:
            new_str_x = '-(' + self.form_1_str_y + ')'
            new_str_y = self.form_1_str_x
            # from these, get numerical solutions, evaulated on local
            # strings changed to relate to the self grids
            # need to supply these unformatted, so save those:
            form_1_x_unformated, form_1_y_unformated = new_str_x*1, new_str_y*1
            # from these strings, get the numerical 1-form:
            new_str_x = new_str_x.replace('x', '(self.xg)')
            new_str_x = new_str_x.replace('y', '(self.yg)')
            new_str_y = new_str_y.replace('x', '(self.xg)')
            new_str_y = new_str_y.replace('y', '(self.yg)')
            # correct for constants
            if new_str_x.find('x') & new_str_x.find('y') == -1:
                new_str_x = '(' + str(new_str_x) + ')* np.ones(np.shape(self.xg))'
            if new_str_y.find('x') & new_str_y.find('y') == -1:
                new_str_y = '(' + str(new_str_y) + ')* np.ones(np.shape(self.yg))'
            
            # evaluate
            form_1_x = eval(new_str_x)
            form_1_y = eval(new_str_y)
            
            # depending on keep_object, return:
            if keep_object:
                self.F_x = form_1_x
                self.F_y = form_1_y
                self.form_1_str_x = form_1_x_unformated
                self.form_1_str_y = form_1_y_unformated
            elif not keep_object:
                new_object = form_1(self.xg, self.yg, form_1_x, form_1_y, F_x_eqn=form_1_x_unformated, F_y_eqn=form_1_y_unformated)
                return new_object
            else:
                raise ValueError('Error, Invalid input for \'keep_object\'')
    
    def num_hodge(self, keep_object=False):
        '''
        
        num_hodge(keep_object=False)
        
        Parameters:
        -------------
        keep_object - determines if the result should be returned as a new
                      1-form or if current one need to be changed.
                      Default is False. When False, a new object is created
                      When true, the acted on is modified.
        
        It calulates the Hodge on R^2 by the standard definition:
        dx -> dy and dy -> -dx
        
        Does no numerically using only component arrays.
        If equations have been previously provided, this method will
        loose them
        
        returns: 1-form if keep_object is False, else returns nothing
        '''
        # check if equations have been given:
        # if they have, doing it only numerically would create
        # a mismatch, warn user
        if self.form_1_str_x != None or self.form_1_str_y != None:
            print('Warning: You supplied equations, doing it numerically only will result in a mismacth between numerical values and equations')
        
        # now complete the process numerically save as instructed
        # check keep_object:
        if keep_object:
            # change the object self properties accoringly
            new_x = -self.F_y
            new_y = self.F_x
            self.F_x = new_x
            self.F_y = new_y
        elif not keep_object:
            # pass these in to the object to create a new one:
            # N.B no equations to supply
            new_object = form_1(self.xg, self.yg, -self.F_y, self.F_x)
            # return the new one to the user:
            return new_object
        else:
            raise ValueError('Error, Invalid input for \'keep_object\'')
    
    # define a fucntion to compute a wedge product of two 1 forms
    def wedge(self, form_second, degree=1, keep_object=False):
        '''
        
        wedge(form_second, degree=1, keep_object=False)
        
        Parameters:
        ----------------
        form_second - the form to wedge the 1-form with.
                    Can be supplied as a DFormPy instance, a tuple of equations,
                    or a single string equation depending on what form is to be
                    wedged.
                    To wedge with 1-form, supply 1-form instance, or tuple of
                    component equations as strings in terms of x and y.
                    To wedge with 0-form or 2-form, supply corresponding
                    instances or a single equation. When using equations,
                    to distinguish between them, provide parmater 'degree'.
        degree - default is 1. Only used when a single string is supplied
                    as form_second, to distinguish betwen 0-form and 2-form
                    for 0-form, degree=0, for 2-form, degree=2.
                    Determines what form is to be wegded with the
                    given 1-form.
        keep_object - bool -default=False - only used when 1-form is wedged
                    with a 0-form. If False, a new object is created as 
                    a result of the wedge. If True, the 1-form acted on
                    is modified to be the result of the wedge. 
        
        To do so here, strings for the form must be supplied.
        Computes the Wedge product using strings, ANALYTICALLY
        
        Returns:
        --------------
        Wedged with 0-form returns a 1-form object if keep_object is False
                    (default), and returns nothing when it is True
        Wedged with a 1-form, returns a 2-form instance
        Wedged with a 2-form, operation makes a 3-form, which on R^2 is
                    always = zero, only message displays.
        
        '''
        
        # test if equations were given first:
        if self.form_1_str_x == None or self.form_1_str_y == None:
            raise ValueError('Error: You need to supply the 1-form equation to do this, look at \'give_eqn\' method')
        
        # set up variable to store order of supplied form, initially assume 1-form
        order = 1
        
        # get needed second obejct strings dep. on input
        if isinstance(form_second, tuple):
            # if equations were given here take these, if numerical grids were given - error!
            # check size , should be a 1-form
            if len(form_second) == 2:
                # 1-form/\1-form, check if strings supplied
                if isinstance(form_second[0], str) and isinstance(form_second[1], str):
                    to_wedge_x_2_str = form_second[0]
                    to_wedge_y_2_str = form_second[1]
                    order = 1
                else:
                    raise ValueError('for analytical calulation, supply 1-form equations as strings')
            else:
                raise ValueError('too many or too little equations given in tuple')
        elif isinstance(form_second, str):
            # single string, could be 0-form or 2-form, check given degree:
            if degree == 0:
                to_wedge_0_form_str = form_second
                order = 0
            elif degree == 2:
                # Error, gives 3 form = 0 on R2
                order = None
                print('This operation makes a 3-form, which on R^2 is always = zero')
            else:
                raise ValueError('not possible digree given or supplied one string for a 1-form')
        else:
            # object supplied, get numericals checking which object is given:
            if isinstance(form_second, form_1):
                if form_second.form_1_str_x is None or form_second.form_1_str_y is None:
                     raise ValueError('supplied 1-form instance must contain equations for analytical calculation')
                else:
                    to_wedge_x_2_str = form_second.form_1_str_x
                    to_wedge_y_2_str = form_second.form_1_str_y
                    order = 1
            elif isinstance(form_second, form_0):
                if form_second.form_0_str is None:
                    raise ValueError('supplied 0-form instance must contain equations for analytical calculation')
                else:
                    to_wedge_0_form_str = form_second.form_0_str
                    order = 0       
            elif isinstance(form_second, form_2):
                order = None
                print('This operation makes a 3-form, which on R^2 is always = zero')
            else:
                raise TypeError('Supplied form to wedge with is not recognised')
        
        # Deal with 1-form/\1-form:
        if order == 1:
            # first, mathematically:  2-form = f*m - g*h
            form_2_str = str(simplify( '(' + self.form_1_str_x + ')*(' +  to_wedge_y_2_str + ')' + ' - (' + self.form_1_str_y + ')*(' +  to_wedge_x_2_str + ')' ))
            # keep it as it is locally to supply it to object maker later
            form_2_str_loc = form_2_str + ''
            # format it to be in terms of grids and:
            # check against constant and zero 2-forms being supplied
            # get the numerical evaluation of it
            form_2_str = form_2_str.replace('x', 'self.xg')
            form_2_str = form_2_str.replace('y', 'self.yg')
            if form_2_str.find('x') & form_2_str.find('y') == -1:
                form_2_str = '(' + str(form_2_str) + ')* np.ones(np.shape(self.xg))'
            
            # evaluate it numerically on the grid supplied
            form_2_result = eval(form_2_str)
            
            # create a 2-form object from this; to return and do so
            ret_object = form_2(self.xg, self.yg, form_2_result, form_2_str_loc)
            return ret_object
        
        elif order == 0:
            # first, find the result of the 1-form:
            new_str_x = str(simplify('(' + self.form_1_str_x + ')*(' +  to_wedge_0_form_str + ')'))
            new_str_y = str(simplify('(' + self.form_1_str_y + ')*(' +  to_wedge_0_form_str + ')'))
            # keep it as it is locally to supply it to object maker later
            form_1_str_x_loc = new_str_x + ''
            form_1_str_y_loc = new_str_y + ''
            # format it to be in terms of grids and:
            # check against constant and zero 1-forms being supplied
            # get the numerical evaluation of it
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
            if keep_object:
                self.F_x = form_1_x
                self.F_y = form_1_y
                self.form_1_str_x = form_1_str_x_loc
                self.form_1_str_y = form_1_str_y_loc
            elif not keep_object:
                new_object = form_1(self.xg, self.yg, form_1_x, form_1_y, F_x_eqn=form_1_str_x_loc, F_y_eqn=form_1_str_y_loc)
                # return the new one to the user:
                return new_object
            else:
                raise ValueError('Error, Invalid input for \'keep_object\'')
        elif order is None:
            # made a form that is always zero on R2, no need to make it
            # Warning already shown, when degree was set
            pass
        else:
            # should never happen, but in case
            raise ValueError('Variable change during code running, look at \'order\' parameter')
    
    
    # define a method for numerical wedge product
    def num_wedge(self, form_second, degree=1, keep_object=False):
        '''
        
        num_wedge(form_second, degree=1, keep_object=False)
        
        Parameters:
        ----------------
        form_second - the form to wedge the 1-form with.
                    Can be supplied as a DFormPy instance, a tuple of grids of
                    same size and dimensions as this 1-form,
                    or a single grid of scaling function values depending on
                    what form is to be wedged.
                    To wedge with 1-form, supply 1-form instance, or tuple of
                    component grids of same size as 1-form acted on.
                    To wedge with 0-form or 2-form, supply corresponding
                    instances or a single grid. When using grids,
                    to distinguish between them, provide parmater 'degree'.
        degree - default is 1. Only used when a single grid is supplied
                    as form_second, to distinguish betwen 0-form and 2-form
                    for 0-form, degree=0, for 2-form, degree=2.
                    Determines what form is to be wegded with the
                    given 1-form.
        keep_object - bool -default=False - only used when 1-form is wedged
                    with a 0-form. If False, a new object is created as 
                    a result of the wedge. If True, the 1-form acted on
                    is modified to be the result of the wedge. 
        
        Computes the Wedge product numerically
        
        Returns:
        --------------
        Wedged with 0-form returns a 1-form object if keep_object is False
                    (default), and returns nothing when it is True
        Wedged with a 1-form, returns a 2-form instance
        Wedged with a 2-form, operation makes a 3-form, which on R^2 is
                    always = zero, only message displays.
        
        '''
        
        # test if equations were given first:
        if isinstance(self.form_1_str_x, str) or isinstance(self.form_1_str_y, str):
            print('The first 1-form you are completing the wedge with has equations supplied, these will be lost')
        
        # set up variable to store order of supplied form, initially assume 1-form
        order = 1
        
        # get needed second obejct grids dep. on input
        if isinstance(form_second, tuple):
            # check size to see what it is to be wedged with.
            # tuple should only be length 2 --> 1-form/\1-form
            if len(form_second) == 2:
                # 1-form/\1-form, extract components
                # if numerical grids were given, take these, if equations, change to values on grids:
                if isinstance(form_second[0], str) and isinstance(form_second[1], str):
                    new_str_x = form_second[0].replace('x', '(self.xg)')
                    new_str_x = new_str_x.replace('y', '(self.yg)')
                    new_str_y = form_second[1].replace('x', '(self.xg)')
                    new_str_y = new_str_y.replace('y', '(self.yg)')
                    if new_str_x.find('x') & new_str_x.find('y') == -1:
                        new_str_x = '(' + str(new_str_x) + ')* np.ones(np.shape(self.xg))'
                    if new_str_y.find('x') & new_str_y.find('y') == -1:
                        new_str_y = '(' + str(new_str_y) + ')* np.ones(np.shape(self.yg))'
                    f12_x = eval(new_str_x)
                    f12_y = eval(new_str_y)
                    order = 1
                elif isinstance(form_second[0], np.ndarray) and isinstance(form_second[1], np.ndarray):
                    f12_x = form_second[0]
                    f12_y = form_second[1]
                    order = 1
                else:
                    raise ValueError('Not recognised input tuple')
            else:
                raise ValueError('too many or too little equations given in tuple')
        
        elif isinstance(form_second, np.ndarray):
            # check degree:
            if degree == 0:
                to_wedge_0_form = form_second
                order = 0
            elif degree == 1:
                raise ValueError('for degree 1, supply a 1-form, not a single grid')
            elif degree == 2:
                # Error, gives 3 form = 0 on R2
                order = None
                print('This operation makes a 3-form, which on R^2 is always = zero')
        
        elif isinstance(form_second, str):
            # single string, could be 0-form or 2-form, check given degree:
            if degree == 0:
                    str_0_form = form_second.replace('x', '(self.xg)')
                    str_0_form = str_0_form.replace('y', '(self.yg)')
                    if str_0_form.find('x') & str_0_form.find('y') == -1:
                        str_0_form = '(' + str(str_0_form) + ')* np.ones(np.shape(self.xg))'
                    
                    to_wedge_0_form = eval(str_0_form)
                    order = 0
            elif degree == 2:
                # Error, gives 3 form = 0 on R2
                order = None
                print('This operation makes a 3-form, which on R^2 is always = zero')
            else:
                raise ValueError('not possible digree given or supplied one string for a 1-form')
        
        # object supplied, get grids checking which object is given:
        
        elif isinstance(form_second, form_1):
            f12_x = form_second.F_x
            f12_y = form_second.F_y
            order = 1
        elif isinstance(form_second, form_0):
            to_wedge_0_form = form_second.form_0
            order = 0
        elif isinstance(form_second, form_2):
            order = None
            print('This operation makes a 3-form, which on R^2 is always = zero')
        else:
            raise TypeError('Supplied form to wedge with is not recognised')
        
        # USe given inputs to evaluate the result:
        
        # Deal with 1-form/\1-form:
        if order == 1:
            # from these get the numerical 2-form
            result = self.F_x * f12_y - self.F_y * f12_x
            
            # return it to user:
            ret_object = form_2(self.xg, self.yg, result)
            return ret_object
        
        elif order == 0:
            # first, find the result of the 1-form
            new_form_1_x = to_wedge_0_form * self.F_x
            new_form_1_y = to_wedge_0_form * self.F_y
            
            # depending on keep_object, return:
            if keep_object:
                self.F_x = new_form_1_x
                self.F_y = new_form_1_y
            elif not keep_object:
                new_object = form_1(self.xg, self.yg, new_form_1_x, new_form_1_y)
                # return the new one to the user:
                return new_object
            else:
                raise ValueError('Error, Invalid input for \'keep_object\'')
        elif order is None:
            # made a form that is always zero on R2, no need to make it
            # Warning already shown, when degree was set
            pass
        else:
            # should never happen, but in case
            raise ValueError('Variable change during code running, look at \'order\' parameter')
        
    
    def zoom(self, target=[0, 0], mag=2, dpd=9, inset=False, axis=None, insize=0.3):
        '''
        
        zoom(target=[0, 0], mag=2, dpd=9, inset=False, axis=None, insize=0.3)
        
        Parameters:
        -------------- 
        Create a new window which displays the field zoomed at a certain point
        User gives arguments
        Target: Determines the zoom location, coordinates
        mag: +ve float, determines zooming amount
        dpd: +int, determines how many points on each axis
        inset - bool - determines if zoom is to plotted as an inset
                if True, need to also give axis on which to plot
        axis - matplotlib axes instance - on it, the instance will plot.
        insize - float - size of inset as fraction of total figure
        
        returns:
        --------------
            if inset is False, returns the zoomed in insatnce as a 0-form
            object
            if inset if True, returns the inset axis, with the plot on them
            on top of the given axis and the 0-form instance
        '''
        
        # Requires user to provide eqn of the 1-form they are zooming on.
        
        if self.form_1_str_x == None or self.form_1_str_y == None:
            # ERROR
            raise TypeError('No equation provided, see \'give_eqn\' method')
        else:
             # Zoom must be one or greater
            if mag < 1:
                raise ValueError('Mag must be greater than one')
            else:
                
                if insize > 1 or insize < 0:
                    raise ValueError('Insize must be +ve and less than one')
                else:
                    
                    # If no inset, set the size of the zoom axis to allow normal plotting
                    if inset == False:
                        insize = 1
            
                    # Target coordinates
                    x_m = target[0]
                    y_m = target[1]
                    
                    # Get the size of the original VF
                    Lx = 0.5*(self.xg[0, -1] - self.xg[0, 0])
                    Ly = 0.5*(self.yg[-1, 0] - self.yg[0, 0])
                    
                    # Zoom axis range
                    d_range_x = insize*Lx/mag
                    d_range_y = insize*Ly/mag
                    
                    # Set up zoom window grids
                    dx = np.linspace(-d_range_x + x_m, d_range_x + x_m, dpd)
                    dy = np.linspace(-d_range_y + y_m, d_range_y + y_m, dpd)
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
                    
                    # crate the zoomed in form
                    zoom_form = form_1(dxg, dyg, u_zoom, v_zoom, self.form_1_str_x, self.form_1_str_y)
                    zoom_form.sheet_size(1/dpd)
                    
                    q = 1
                    
                    xi = (x_m - self.xg[0,0])/(2*Lx)
                    yi = (y_m - self.yg[0,0])/(2*Ly)
                    
                    if inset == True:
                        if axis != None:
                            # Create inset axis in the current axis.
                            
                            zoom_inset_ax = axis.inset_axes([(xi - 0.5*insize), (yi - 0.5*insize), insize, insize])
                            zoom_form.plot(zoom_inset_ax)
                            
                            # return the zoomed on axis
                            # also return zoomed in form in case user wants that.
                            return zoom_inset_ax, zoom_form
                        else:
                            raise ValueError('Cannot inset without supplied axis')
                    else:
                        # inset is false, just return the new zoomed in instance
                        return zoom_form
    
    
    # define a mehtod to evaluate the interior derivative of the 1-form
    # with respect to a given vector field object or without.
    def interior_d(self, vector_field=None):
        '''
        
        interior_d(vector_field=None)
        
        Computes the interior derivative of the 1-form
        
        Parameters:
        ------------------
        Vector_field = vector field object of DFormPy library to do the
            derivative with respect to, needs equations to work with
            nuymerical_only being False. Can also supply equations in a tuple:
            (eqn_x, eqn_y). If using numerical only, can supply object or
            tuple of numpy arrays (array_x, atrray_y). If nothing is supplied
            for it, it assumes F_x = 1 and F_y = 1, with correct form and shape
        
        Does no analytically using equations provided in instance
        
        Returns 0-form object
        '''
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
        result_form = form_0(self.xg, self.yg, zero_form_result, zero_form_str_unformatted)
        
        # return it to the user
        return result_form
        
    # numerical interior derivaitve
    def num_interior_d(self, vector_field=None):
        '''
        
        num_interior_d(self, vector_field=None)
        
        Computes the interior derivative of the 1-form
        
        Parameters:
        --------------
        Vector_field = vector field object of DFormPy library to do the
            derivative with respect to, needs equations to work with
            nuymerical_only being False. Can also supply equations in a tuple:
            (eqn_x, eqn_y). If using numerical only, can supply object or
            tuple of numpy arrays (array_x, atrray_y). If nothing is supplied
            for it, it assumes F_x = 1 and F_y = 1, with correct form and shape
        
        Does no numerically using arrays provided in instance
        If equations were proivided, this method will lose them
        
        Returns 0-form object
        '''
        # check if equations have been given:
        # if they have, doing it only numerically would create
        # a mismatch, Warn user
        if self.form_1_str_x == None or self.form_1_str_y == None:
            pass
        else:
            # equations have been given, a mismatch may occur
            # warn the user
            print('Warning: You supplied equations, doing it numerically only will not pass equations to the 0-form and these will be lost')
        
        # Take the vector field components, checking what was input
        if vector_field is None:
            # if none was given, do it with respect to uniform 1, 1
            vf_x = np.ones(np.shape(self.xg))
            vf_y = np.ones(np.shape(self.yg))
        elif type(vector_field) == tuple:
            # if numerical grids were given, take these
            # if equations were given here, evaulate them to grids
            if type(vector_field[0]) == str:
                new_str_x = vector_field[0].replace('x', '(self.xg)')
                new_str_x = new_str_x.replace('y', '(self.yg)')
                new_str_y = vector_field[1].replace('x', '(self.xg)')
                new_str_y = new_str_y.replace('y', '(self.yg)')
                if new_str_x.find('x') & new_str_x.find('y') == -1:
                    new_str_x = '(' + str(new_str_x) + ')* np.ones(np.shape(self.xg))'
                if new_str_y.find('x') & new_str_y.find('y') == -1:
                    new_str_y = '(' + str(new_str_y) + ')* np.ones(np.shape(self.yg))'
                vf_x = eval(new_str_x)
                vf_y = eval(new_str_y)
            else:
                vf_x = vector_field[0]
                vf_y = vector_field[1]
        else:
            # extract needed properties from the object supplied
            vf_x = vector_field.F_x
            vf_y = vector_field.F_y
        
        # Complete the interior derivative 1-form --> 0-form:
        zero_form_result = self.F_x * vf_x + self.F_y * vf_y
        
        # supply these to the 0-form object creator
        result_form = form_0(self.xg, self.yg, zero_form_result)
        
        # return it to the user
        return result_form
    
    # define a method to change a supplied Vector filed to the 1-form
    def contravariant(self, g=[['1', '0'], ['0', '1']]):
        '''
        
        contravariant(g=[['1', '0'], ['0', '1']])
        
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
        here automatically and listed in the documentation
        
        Returns a single object (VF object)
        '''
        
        # extract what is needed form the metric depending on what the user
        # supplied
        # check if it has string components
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
            # stored numerical metric
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
            # if the user has 1-form equations, warn that these 
            # will be lost, due to numerical calculations
            if self.form_1_str_x == None and self.form_1_str_y == None:
                pass
            else:
                print('The 1-form has equations, but the metric does not, these will be lost and the resulting VF will only have numerical values, not equations supplied')
            # No need to do anythng more to the metric
            # just rename the metric here
            g_num = g
            
            # set up the dummy variable
            analytics = False
        else:
            # Inconsistant metric components
            raise TypeError('Metric components are inconcisstant')
        
        # from 1-form components, get VF components by the metric
        # first, do so numerically, as this must always happen
        form_x = self.F_x * g_num[0][0] + self.F_y * g_num[0][1]
        form_y = self.F_y * g_num[1][1] + self.F_x * g_num[1][0]
        
        # if the equations were given, evaluate these analytically too:
        # only if vector filed originally has equations
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
            result_field = vector_field(self.xg, self.yg, form_x, form_y, x_str_form, y_str_form)
        elif not analytics:
            result_field = vector_field(self.xg, self.yg, form_x, form_y)
        
        # return the found object
        return result_field

# %%

'''

function to create a 2-form object and define methods for it

'''

# define 2-form object that can be customised and plotted
class form_2():
    '''
    defines a 2-form object and returns it to user
    Takes 3 arguments basic, these are the 2 grids in 2D, which muse be square
    and of equal sizes. Then 1 argument for the dx^dy component
    based on the same grids. Also takes in an equation which is needed for some
    operaions
    Takes in a figure if one is to be supplied. Can take axis for subplots in
    The subplots only occur if subplots input is set to True, default is False
    '''
    # set up all variables
    def __init__(self, xg, yg, form2, form_2_eq=None):
        self.xg = xg
        self.yg = yg
        self.form_2 = form2
        self.s_max = 6
        self.s_min = 2
        self.pt_den_x = len(xg[0, :])
        self.pt_den_y = len(yg[:, 0])
        self.fract_x = 2/((self.pt_den_x - 1))
        self.fract_y = 2/((self.pt_den_y - 1))
        self.colour_list = ['red', 'blue', 'grey']
        self.logarithmic_scale_bool = 0
        # self.base = 10
        self.delta_factor = 10
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
    
    # change colour list
    def colours(self, colour_list):
        '''
        Takes input of a list of three string. String must be formatted
        as to be accepted by maplotlib colors
        changes the colours for 2-form orientation.
        Order: [clockwise, counterclosckwise, zero]
        '''
        
        # make sure input was a list of strings:
        if not(isinstance(colour_list[0], str) and isinstance(colour_list[1], str) and isinstance(colour_list[2], str)):
            raise TypeError('Wrongly formatted string list, chech required inputs')
        
        # change stored colour list
        self.colour_list = colour_list
    
    # change boolean that det. if to sclae logarithmically
    def log_scaling(self):
        '''
        Takes no arguments
        Changes the boolean that determines if scaling is logarithmic
        Whenever it is called, it changes that boolean to opposite
        The form object is initialised with this as False (as 0)
        '''
        self.logarithmic_scale_bool = not self.logarithmic_scale_bool
        # self.base = base
    
    # define methods to change s_max
    def max_sheets(self, maximum):
        '''
        Takes one argument, must be int
        Changes maximum number of sheets to draw on a stack.
        These still scale relative to max magnitude.
        '''
        self.s_max = maximum
    
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
    def set_density2(self, points_number_x, points_number_y):
        '''
        
        Changes number of points on grids to given, if equations have been given
        
        Parameters:
        -------------
        points_number_x - int - number of points to put along the x axis
        points_number_y - int - number of points to put along the y axis
        
        Returns: None
        '''
        if self.form_2_str == None:
            # Error
            raise TypeError('Error: You need to supply the 2-form equation to do this, look at \'give_eqn\' method')
        else:
            # redefine the grids
            x = np.linspace(self.xg[0,0], self.xg[0,-1], points_number_x)
            y = np.linspace(self.yg[0,0], self.yg[-1,0], points_number_y)
            self.xg, self.yg = np.meshgrid(x, y)
            # based on these change other, dependant variables
            self.pt_den_x = len(self.xg[0, :])
            self.pt_den_y = len(self.yg[:, 0])
            self.fract_x = 2/(self.pt_den_x - 1)
            self.fract_y = 2/(self.pt_den_y - 1)
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
    def plot(self, axis):
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
        
        
        form2 = self.form_2 * 1  # from self, get 2-form too
        
        # set all insignificant values to zero:
        form2[np.abs(form2) < 1e-12] = 0
        
        # get the lengths of x and y from their grids
        x_len = len(self.xg[0, :])
        y_len = len(self.yg[:, 0])
        
        # Extract L from the x and y grids
        Lx = 0.5*(self.xg[0, -1] - self.xg[0, 0])
        Ly = 0.5*(self.yg[-1, 0] - self.yg[0, 0])
        L = 0.5*(Lx + Ly)
        x0 = self.xg[0, 0] + Lx
        y0 = self.yg[0, 0] + Ly
        
        # reset axis limits
        ax_Lx = Lx + Lx/self.delta_factor
        ax_Ly = Ly + Ly/self.delta_factor
        axis.set_xlim(-ax_Lx + x0, ax_Lx + x0)
        axis.set_ylim(-ax_Ly + y0, ax_Ly + y0)
        
        # get the signs of the input 2-form
        form_2_sgn = np.sign(form2)
        
        # define an empty array of magnitudes, to then fill with integer rel. mags
        R_int = np.zeros(shape=((y_len), (x_len)))
        
        # #########################################################################
        # get variables needed for the initial, simplified stack plot
        # #########################################################################
        
        # set up directions
        angles =[0*np.ones(np.shape(form2)), (np.pi/2)*np.ones(np.shape(form2))]
        
        # deal with sinularities that appear on evaluated points
        isnan_arr = np.isnan(form2)
        for i in range(y_len):
            for j in range(x_len):
                # set to zero points that are not defined or inf
                if isnan_arr[i, j] or abs(form2[i, j]) == np.inf  or abs(form2[i, j]) > 1e15:
                    # colour this region as a red dot, not square to
                    # not confuse with nigh mag 2-forms in stacks. or worse, in
                    # blocks
                    circ = patch.Circle((self.xg[i, j], self.yg[i, j]), L*(self.fract_x + self.fract_y)/6, color='red')
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
        s_L_x = self.fract_x*Lx
        s_L_y = self.fract_y*Ly
        
        # #########################################################################
        # define the stacks based on geometrical arguments
        # to be perp. to arrow. shifted parallel to it, their density porp to mag
        # of the arrow and with an arrowhead on top.
        # #########################################################################
        # find the maximum magnitude for scaling
        
        mag = abs(form2)
        
        max_size = np.max(mag)   # careful with singularities, else ---> nan

        if self.logarithmic_scale_bool:
            # Add 1 to each magnitude
            mag1 = mag + 1
            # Calculate the appropriate scaling factor
            # a = max_size**(1/self.s_max)
            # a = self.base
            # Take log(base=a) of mag1
            logmag1 = np.log(mag1)
            # Re-assign R
            R = logmag1/np.max(logmag1)    
        else:
            # find the relative magnitude of vectors to maximum, as an array
            R = mag/max_size
        
        # if self.logarithmic_scale_bool:
        #     mag1 = mag + 1 
        #     form_2_norm = form2/mag1
        #     logmag = np.log10(mag1)
        #     form2 = form2_norm*logmag
        #     mag = np.abs(form2)
        #     max_size = np.max(mag)

        # Now, for both values of theta, complete plotting:
        for theta in angles:
            # define tigonometirc shifts
            I_sin = np.sin(theta)
            I_cos = np.cos(theta)
            
            # define the points that set out a line of the stack sheet (middle line)
            A_x = self.xg + (s_L_x/2)*I_sin
            A_y = self.yg - (s_L_y/2)*I_cos
            B_x = self.xg - (s_L_x/2)*I_sin
            B_y = self.yg + (s_L_y/2)*I_cos
            
            
            for i in range(self.s_max - self.s_min + 1):
                t = self.s_max - i
                R_int[R <= t/self.s_max] = t

            # loop over each arrow coordinate in x and y
            for i in range(y_len):
                for j in range(x_len):
                    # define it for all magnitudes. Separately for odd and even corr. number of sheets:
                    
                    # Label each element with the number of stacks required: linear scaling
                    if form_2_sgn[i, j] == +1:
                        color_index = 0
                    elif form_2_sgn[i, j] == -1:
                        color_index = 1
                    else:
                        color_index = 2
                    
                    # # linear scaling
                    # for t in range(self.s_min, self.s_max+2):
                    #     if (t-2)/self.s_max <= R[i, j] <= (t-1)/self.s_max:
                    #         R_int[i, j] = t
                    
                    # set a varible for current considered magnitude as it is reused
                    # avoids extracting from R many times.
                    n = R_int[i, j]
                    
                    # deal with even number of sheets from magnitudes:
                    if n % 2 == 0:
                        # define a parameter to loop over in the recursion equation
                        s = 0
                        
                        # Define the points for sheets required for the given magnitude
                        # from these define all the needed lines and plot them
                        while s <= 0.5*(n-2):  # maximum set by equations (documentation)
                            # define all the points for the 2 currently looped +- sheets in while loop
                            Ax1 = A_x[i, j] + G(s, n, 0)*s_L_x*I_cos[i, j]
                            Ay1 = A_y[i, j] + G(s, n, 0)*s_L_y*I_sin[i, j]
                            Bx1 = B_x[i, j] + G(s, n, 0)*s_L_x*I_cos[i, j]
                            By1 = B_y[i, j] + G(s, n, 0)*s_L_y*I_sin[i, j]
                            Ax2 = A_x[i, j] - G(s, n, 0)*s_L_x*I_cos[i, j]
                            Ay2 = A_y[i, j] - G(s, n, 0)*s_L_y*I_sin[i, j]
                            Bx2 = B_x[i, j] - G(s, n, 0)*s_L_x*I_cos[i, j]
                            By2 = B_y[i, j] - G(s, n, 0)*s_L_y*I_sin[i, j]
                            
                            # from these, define the 2 lines, for this run
                            axis.add_line(Line2D((Ax1, Bx1), (Ay1, By1), linewidth=0.5, color=self.colour_list[color_index]))
                            axis.add_line(Line2D((Ax2, Bx2), (Ay2, By2), linewidth=0.7, color=self.colour_list[color_index]))
                            
                            # update parameter to reapet and draw all needed arrows
                            s += 1
                    # deal with the odd number of stacks:
                    else:
                        # Add the centre line for odd numbers of stacks
                        axis.add_line(Line2D((A_x[i, j], B_x[i, j]), (A_y[i, j], B_y[i, j]), linewidth=0.7, color=self.colour_list[color_index]))
                        
                        # then loop over the remaining lines as per the recursion formula:
                        s = 1  # change the looping parametr to exclude already completed 0 (corr. to middle sheet here)
                        
                        # define all remaining sheets for the magnitude:
                        while s <= 0.5*(n-1):  # maximum set by equations (documentation)
                            # define all the points for the current +- displacement in while loop
                            Ax1 = A_x[i, j] + G(s, n, 1)*s_L_x*I_cos[i, j]
                            Ay1 = A_y[i, j] + G(s, n, 1)*s_L_y*I_sin[i, j]
                            Bx1 = B_x[i, j] + G(s, n, 1)*s_L_x*I_cos[i, j]
                            By1 = B_y[i, j] + G(s, n, 1)*s_L_y*I_sin[i, j]
                            Ax2 = A_x[i, j] - G(s, n, 1)*s_L_x*I_cos[i, j]
                            Ay2 = A_y[i, j] - G(s, n, 1)*s_L_y*I_sin[i, j]
                            Bx2 = B_x[i, j] - G(s, n, 1)*s_L_x*I_cos[i, j]
                            By2 = B_y[i, j] - G(s, n, 1)*s_L_y*I_sin[i, j]
                            
                            # from these, define the 2 displaced lines
                            
                            axis.add_line(Line2D((Ax1, Bx1), (Ay1, By1), linewidth=0.7, color=self.colour_list[color_index]))
                            axis.add_line(Line2D((Ax2, Bx2), (Ay2, By2), linewidth=0.7, color=self.colour_list[color_index]))
                            
                            # change the parameter to loop over all changes in displacement for current magnitude
                            s += 1
    
    def ext_d(self):
        '''
        No inputs, no outputs, exterior derivative of a 2-form gives
        a 3-form, which on R2 is always =0
        '''
        print('This operation makes a 3-form, which on R^2 is always = zero')
    
    # define a fucntion to Hodge the 2-form (into a 0-form)
    def num_hodge(self):
        '''
        Takes in no arguments
        Does the hodge numerically based on instance provieded arrays
        If equations were provided, it will lose them.
        
        It calulates the Hodge on R^2 by the standard definition:
        *(dx^dy) = 1
        
        returns a 0-form
        '''
        # check if equations have been given:
        # if they have, doing it only numerically would create
        # a mismatch, avoid that
        if self.form_2_str != None:
            # equations have been given, a mismatch may occur
            # warn the user
            print('Warning: You supplied equations, doing it numerically only will lose these')
        
        # now complete the process numerically
        # pass these in to the object to create a new one:
        new_object = form_0(self.xg, self.yg, self.form_2)  # N.B no equations to supply
        
        # return the new one to the user:
        return new_object
    
    
    def hodge(self):
        '''
        Takes in no arguments
        Does the hodge analuically based on instance provieded equations
        changes the equations AND the numerical answers
        
        It calulates the Hodge on R^2 by the standard definition:
        *(dx^dy) = 1
        
        returns a 0-form
        '''
        # can only be done if equations have been given, check:
        if self.form_2_str != None:
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
            new_object = form_0(self.xg, self.yg, form_0_result, form_0_eqn=form_0_str_unformated)
            
            # return the new one to the user:
            return new_object
        else:
            # ERROR
            raise TypeError('You need to supply the 2-form equation to do this, look at \'give_eqn\' method')

    # define a method to create a zoomed in 2-form
    def zoom(self, target=[0, 0], mag=2, dpd=9, inset=False, axis=None, insize=0.3):
        
        '''
        Creates a new window which displays the 2-form zoomed at a certain point
        User gives arguments:
        Target: Determines the zoom location, coordinates
        mag: +ve float, determines zooming amount
        dpd: +int, determines how many points on each axis
        
        inset - bool - determies if the zoom is plotted on the given axis
        as an inset
        axis - matplotlib axis, only supply if inset is True, plots intset on these
        insize - float - size of inset as fraction of total figure
        
        returns:
        --------------
            if inset is False, returns the zoomed in insatnce as a 0-form
            object
            if inset if True, returns the inset axis, with the plot on them
            on top of the given axis and the 0-form instance
        '''
        
        # Requires user to provide eqn of the 1-form they are zooming on.
        if self.form_2_str == None:
            # ERROR
            raise TypeError('No equation provided')
        else:
             # Zoom must be one or greater
            if mag < 1:
                raise ValueError('Mag must be greater than one')
            else:
                
                if insize > 1 or insize < 0:
                    raise ValueError('Insize must be +ve and less than one')
                else:
                    
                    # If no inset, set the size of the zoom axis to allow normal plotting
                    if inset == False:
                        insize = 1
                        
                    # Target coordinates
                    x_m = target[0]
                    y_m = target[1]
                    
                    # Get the size of the original VF
                    Lx = 0.5*(self.xg[0, -1] - self.xg[0, 0])
                    Ly = 0.5*(self.yg[-1, 0] - self.yg[0, 0])
                    
                    # Zoom axis range
                    d_range_x = insize*Lx/mag
                    d_range_y = insize*Ly/mag
                    
                    # Set up zoom window grids
                    dx = np.linspace(-d_range_x + x_m, d_range_x + x_m, dpd)
                    dy = np.linspace(-d_range_y + y_m, d_range_y + y_m, dpd)
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
                    
                    # from that create 2-form instance
                    zoomform2 = form_2(dxg, dyg, zoom_2form, self.form_2_str)
                    
                    q = 1
                    
                    xi = (x_m - self.xg[0,0])/(2*Lx)
                    yi = (y_m - self.yg[0,0])/(2*Ly)
                    
                    if inset == True:
                        if axis != None:
                            # Create inset axis in the current axis.
                            zoom_inset_ax = axis.inset_axes([(xi - 0.5*insize), (yi - 0.5*insize), insize, insize])
                            zoomform2.plot(zoom_inset_ax)
                            # return the zoomed on axis
                            # also return zoomed in form in case user wants that.
                            return zoom_inset_ax, zoomform2
                        else:
                            raise ValueError('Cannot inset without supplied axis')
                    else:
                        # inset is false, just return the new zoomed in instance
                        return zoomform2
            
    
    # define a mehtod to evaluate the interior derivative of the 2-form
    # with respect to a given vector field object or without.
    def interior_d(self, vector_field=None):
        '''
        Computes the interior derivative of the 2-form
        Takes in:
        -- Vector_field = vector field object of DFormPy library to do the
        derivative with respect to, needs equations to work with
        nuymerical_only being False. Can also supply equations in a tuple:
        (eqn_x, eqn_y). If using numerical only, can supply object or
        tuple of numpy arrays (array_x, atrray_y). If nothing is supplied
        for it, it assumes F_x = 1 and F_y = 1, with correct form and shape
            
        Does no analytically via equations in instance
        
        Returns: 0-form
        
        '''
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
        
        # create the object to return
        result_form = form_1(self.xg, self.yg, form_x, form_y, u_str_unformatted, v_str_unformatted)
        
        # return it to the user
        return result_form
        
    
    def num_interior_d(self, vector_field=None):
        '''
        Computes the interior derivative of the 2-form
        Takes in:
        -- Vector_field = vector field object of DFormPy library to do the
        derivative with respect to, needs equations to work with
        nuymerical_only being False. Can also supply equations in a tuple:
        (eqn_x, eqn_y). If using numerical only, can supply object or
        tuple of numpy arrays (array_x, atrray_y). If nothing is supplied
        for it, it assumes F_x = 1 and F_y = 1, with correct form and shape
            
        Does no numerically via arrays in instance
        If equations were provided, these will be lost
        
        Returns: 0-form
        
        '''
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
            vf_x = np.ones(np.shape(self.xg))
            vf_y = np.ones(np.shape(self.xg))
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
        
        # supply these to the 1-form object creator
        result_form = form_1(self.xg, self.yg, form_x, form_y)
        
        # return it to the user
        return result_form
    
    # define a fucntion to compute a wedge product
    def wedge(self, form_second, degree=0, keep_object=False):
        '''
        Parameters:
        ----------------
        form_second - the form to wedge the 2-form with.
                    Can be supplied as a DFormPy instance, a tuple of equations,
                    or a single string equation depending on what form is to be
                    wedged.
                    To wedge with 1-form, supply 1-form instance, or tuple of
                    component equations as strings in terms of x and y.
                    To wedge with 0-form or 2-form, supply corresponding
                    instances or a single equation. When using equations,
                    to distinguish between them, provide parmater 'degree'.
        degree - default is 0. Only used when a single string is supplied
                    as form_second, to distinguish betwen 0-form and 2-form
                    for 0-form, degree=0, for 2-form, degree=2.
                    Determines what form is to be wegded with the
                    given 2-form.
        keep_object - bool - default=False - only used when 2-form is wedged
                    with a 0-form. If False, a new object is created as 
                    a result of the wedge. If True, the 1-form acted on
                    is modified to be the result of the wedge. 
        
        To do so here, strings for the form must be supplied.
        Computes the Wedge product using strings, ANALYTICALLY
        
        Returns:
        --------------
        Wedged with 0-form returns a 2-form object if keep_object is False
                    (default), and returns nothing when it is True
        Wedged with a 1-form, operation makes a 3-form, which on R^2 is
                    always = zero, only message displays.
        Wedged with a 2-form, operation makes a 4-form, which on R^2 is
                    always = zero, only message displays.
        
        '''
        
        # test if equations were given first:
        if self.form_2_str == None:
            raise ValueError('Error: You need to supply the 2-form equation to do this, look at \'give_eqn\' method')
        
        # set up variable to store order of supplied form, initially assume 1-form
        order = 0
        
        # get needed second obejct strings dep. on input
        if isinstance(form_second, tuple):
            # if equations were given here take these, if numerical grids were given - error!
            # check size , should be a 1-form
            if len(form_second) == 2:
                # 2-form/\1-form attempt, error
                if isinstance(form_second[0], str) and isinstance(form_second[1], str):
                    order = None
                    print('This operation makes a 3-form, which on R^2 is always = zero')
                else:
                    raise ValueError('for analytical calulation, supply 1-form equations as strings')
            else:
                raise ValueError('too many or too little equations given in tuple')
        elif isinstance(form_second, str):
            # single string, could be 0-form or 2-form, check given degree:
            if degree == 0:
                to_wedge_0_form_str = form_second
                order = 0
            elif degree == 2:
                # Error, gives 4 form = 0 on R2
                order = None
                print('This operation makes a 4-form, which on R^2 is always = zero')
            else:
                raise ValueError('not possible digree given or supplied one string for a 1-form')
        else:
            # object supplied, get numericals checking which object is given:
            if isinstance(form_second, form_1):
                print('This operation makes a 3-form, which on R^2 is always = zero')
                order = None
            elif isinstance(form_second, form_0):
                if form_second.form_0_str is None:
                    raise ValueError('supplied 0-form instance must contain equations for analytical calculation')
                else:
                    to_wedge_0_form_str = form_second.form_0_str
                    order = 0
            elif isinstance(form_second, form_2):
                order = None
                print('This operation makes a 4-form, which on R^2 is always = zero')
            else:
                raise TypeError('Supplied form to wedge with is not recognised')
        
        # Deal with 2-form/\0-form:
        
        
        if order == 0:
            # first, mathematically:  2-form = f*m - g*h
            form_2_str = str(simplify('(' + self.form_2_str + ')*(' +  to_wedge_0_form_str + ')'))
            # keep it as it is locally to supply it to object maker later
            form_2_str_loc = form_2_str + ''
            # format it to be in terms of grids and:
            # check against constant and zero 2-forms being supplied
            # get the numerical evaluation of it
            form_2_str = form_2_str.replace('x', 'self.xg')
            form_2_str = form_2_str.replace('y', 'self.yg')
            if form_2_str.find('x') & form_2_str.find('y') == -1:
                form_2_str = '(' + str(form_2_str) + ')* np.ones(np.shape(self.xg))'
            
            # evaluate it numerically on the grid supplied
            form_2_result = eval(form_2_str)
            
            # depending on keep_object, return:
            if keep_object:
                self.form_2 = form_2_result
                self.form_2_str = form_2_str_loc
            elif not keep_object:
                new_object = form_2(self.xg, self.yg, form_2_result, form_2_eq=form_2_str_loc)
                # return the new one to the user:
                return new_object
            else:
                raise ValueError('Error, Invalid input for \'keep_object\'')
        
        elif order is None:
            # made a form that is always zero on R2, no need to make it
            # Warning already shown, when degree was set
            pass
        else:
            # should never happen, but in case
            raise ValueError('Variable change during code running, look at \'order\' parameter')
    
    
    # define a method for numerical wedge product
    def num_wedge(self, form_second, degree=0, keep_object=False):
        '''
        Parameters:
        ----------------
        form_second - the form to wedge the 2-form with.
                    Can be supplied as a DFormPy instance, a tuple of component
                    grids, or a single string equation depending on what form
                    is to be wedged.
                    To wedge with 1-form, supply 1-form instance, or tuple of
                    component grids.
                    To wedge with 0-form or 2-form, supply corresponding
                    instances or a single grid. When using grids,
                    to distinguish between them, provide parmater 'degree'.
        degree - default is 0. Only used when a grid string is supplied
                    as form_second, to distinguish betwen 0-form and 2-form.
                    For 0-form, degree=0 and for 2-form, degree=2.
                    Determines what form is to be wegded with the
                    given 2-form.
        keep_object - bool - default=False - only used when 2-form is wedged
                    with a 0-form. If False, a new object is created as 
                    a result of the wedge. If True, the 1-form acted on
                    is modified to be the result of the wedge. 
        
        Computes the Wedge product using strings, numerically
        
        Returns:
        --------------
        Wedged with 0-form returns a 2-form object if keep_object is False
                    (default), and returns nothing when it is True
        Wedged with a 1-form, operation makes a 3-form, which on R^2 is
                    always = zero, only message displays.
        Wedged with a 2-form, operation makes a 4-form, which on R^2 is
                    always = zero, only message displays.
        '''
        
        # test if equations were given first, warn user of losses:
        if self.form_2_str == None:
            print('The first 1-form you are completing the wedge with has equations supplied, these will be lost')
        
        # set up variable to store order of supplied form, initially assume 1-form
        order = 0
        
        # get needed second obejct grids dep. on input
        if isinstance(form_second, tuple):
            # check size to see what it is to be wedged with.
            # tuple should only be length 2 --> 1-form/\1-form
            if len(form_second) == 2:
                # 2-form/\1-form attempt, error
                order = None
                print('This operation makes a 3-form, which on R^2 is always = zero')
            else:
                raise ValueError('too many or too little equations given in tuple')
        
        elif isinstance(form_second, np.ndarray):
            # check degree:
            if degree == 0:
                to_wedge_0_form = form_second
                order = 0
            elif degree == 1:
                raise ValueError('for degree 1, supply a 1-form, not a single grid')
            elif degree == 2:
                # Error, gives 3 form = 0 on R2
                order = None
                print('This operation makes a 4-form, which on R^2 is always = zero')
        
        elif isinstance(form_second, str):
            # single string, could be 0-form or 2-form, check given degree:
            if degree == 0:
                    str_0_form = form_second.replace('x', '(self.xg)')
                    str_0_form = str_0_form.replace('y', '(self.yg)')
                    if str_0_form.find('x') & str_0_form.find('y') == -1:
                        str_0_form = '(' + str(str_0_form) + ')* np.ones(np.shape(self.xg))'
                    
                    to_wedge_0_form = eval(str_0_form)
                    order = 0
            elif degree == 1:
                raise ValueError('for degree 1, supply a 1-form, not a single equation')
            elif degree == 2:
                # Error, gives 4 form = 0 on R2
                order = None
                print('This operation makes a 4-form, which on R^2 is always = zero')
            else:
                raise ValueError('not possible digree given')
        
        # object supplied, get grids checking which object is given:
        
        elif isinstance(form_second, form_1):
            # Error, gives 3 form = 0 on R2
            order = None
            print('This operation makes a 3-form, which on R^2 is always = zero')
        elif isinstance(form_second, form_0):
            to_wedge_0_form = form_second.form_0
            order = 0
        elif isinstance(form_second, form_2):
            order = None
            print('This operation makes a 4-form, which on R^2 is always = zero')
        else:
            raise TypeError('Supplied form to wedge with is not recognised')
        
        # Use given inputs to evaluate the result:
        
        if order == 0:
            # depending on keep_object, return:
            if keep_object:
                self.form_2 = to_wedge_0_form * self.form_2
            elif not keep_object:
                new_object = form_2(self.xg, self.yg, to_wedge_0_form * self.form_2)
                # return the new one to the user:
                return new_object
            else:
                raise ValueError('Error, Invalid input for \'keep_object\'')
        elif order is None:
            # made a form that is always zero on R2, no need to make it
            # Warning already shown, when degree was set
            pass
        else:
            # should never happen, but in case
            raise ValueError('Variable change during code running, look at \'order\' parameter')
    

# %%

'''

function to create a 0-form object and define methods for it

'''

# define a function that will set up a 0-form object that can be customised and
# plotted
class form_0():
    '''
    
    form_0(xg, yg, form_0, form_0_eqn=None)
    
    
    Defines a 0-form object and returns it to user. 
    
    Parameters:
    ---------------
    xg - grid of x values (2D numpy.ndarray)
    yg - grid of y values (2D numpy.ndarray)
    form_0 - sclar form grid (2D numpy.ndarray)
    
    Optional:
    form_0_eqn - expression for scalar form f(x,y) (string)
    
    
    Instance variables:
    ---------------
    xg, yg, form_0, form_0_eqn
    pt_den - int - number of points on grids, extracted from grids, assumes square grid
    color - str - colour to draw stacks with, can be Hex when using '#FFFFFF'
    logarithmic_scale_bool - bool - determines if log scaling is used
    N - int - base for log scaling
    delta_factor - float/int - determined size of blank boarder in figure
                                as fraction of whole plot size
    inline_bool - bool - if labels on contours are put on contour lines
    denser - int, default is 1 - if equations are given, increases density
                                of contours
    lines - int - number of contour lines to draw
    cmap - matplotlib colourmap - colour mapping to use
    
    
    Methods:
    ---------------
    give_eqn
    return_string
    colour
    log_scaling
    surround_space
    set_density
    plot
    ext_d
    num_ext_d
    hodge
    wedge_analytical
    wedge_num
    
    '''
    # set up all initial, defualt variables
    def __init__(self, xg, yg, form_0, form_0_eqn=None):
        self.xg = xg
        self.yg = yg
        self.form_0 = form_0
        self.pt_den_x = len(xg[0, :])
        self.pt_den_y = len(xg[:, 0])
        self.delta_factor = 10
        self.denser = 1
        self.lines = 15
        self.fontsize = 7
        self.inline_bool = True
        # Log scaling parameters
        self.logarithmic_scale_bool = 0
        self.N = 30
        
        if form_0_eqn is not None:
            self.form_0_str = str(simplify(form_0_eqn))  # user must change to access some methods
        else:
            self.form_0_str = None
        # Note, the string must be given with x and y as variables
        # gets contour plot with new density.
        self.cmap = cm.viridis
    
    # #####################################################################
    # Define basic methods to customise this object
    # #####################################################################
    
    # define a mehtod to allow user to supply the string equation
    # of the 0-form
    def give_eqn(self, equation_str):
        '''
        
        Allows user to supply equation to instance, if not initially done so
        
        Parameters:
        ------------
        equation_str - str - equation of the supplied numerical 0-form
                        It must be in terms of x and y.
                        Has to be given, for some methods to be calculatable.
        
        Returns: None
        
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
    
    def levels(self, values):
        '''
        Takes 1 argument: values - int or list
        if int: changes number of contour lines that get drawn
                the values are set automatically by matplotlib
        if list: sets values to draw level lines (ascending order)
        
        supplied to contour plot from matplotlib via levels
        '''
        
        if isinstance(values, int) or isinstance(values, list):
            self.lines = values
        else:
            raise TypeError('Require input to be integer or list, if you used a numpy array try: list(your_array)')
        
    def log_scaling(self):
        '''
        changes bool for logscaling
        Default = False
        changes to the other option each time it is called
        '''
        self.logarithmic_scale_bool = not self.logarithmic_scale_bool
    
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
    def set_density(self, points_number):
        '''
        set_density(points_number)
        
        Changes the desnity of points in the same range to the input value
        requires the string equation to be supplied
        Only creates grids with same number of points of each axis.
        
        Parameters:
        ---------------
        points_number -number of points to evaluate on
        
        Returns: none
        '''
        if self.form_0_str == None:
            # Error
            raise TypeError('Error: You need to supply the 0-form equation to do this, look at \'give_eqn\' method')
        else:
            # redefine the grids
            x = np.linspace(self.xg[0,0], self.xg[0,-1], points_number)
            y = np.linspace(self.yg[0,0], self.yg[-1,0], points_number)
            self.xg, self.yg = np.meshgrid(x,y)
            # based on these change other, dependant variables
            self.pt_den_x = len(self.xg[0, :])
            self.pt_den_y = len(self.yg[:, 0])
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
    # Write more useful methods plot, exterior derivative, Hodge etc.
    # #####################################################################
    
    # define a fucntion to plot a zero form pressed.
    def plot(self, axis):
        '''
        Finilises the plotting
        Uses the attribues of the object as set originally and as customised
        with methods to create a plot of the 2-form.
        
        parametes:
        -------------
        axis - matplotlib axis that 0-form will be plotted on
        '''
        
        # Extract L from the x and y grids
        Lx = 0.5*(self.xg[0, -1] - self.xg[0, 0])
        Ly = 0.5*(self.yg[-1, 0] - self.yg[0, 0])
        x0 = self.xg[0, 0] + Lx
        y0 = self.yg[0, 0] + Ly
        
        # reset axis limits
        ax_Lx = Lx + Lx/self.delta_factor
        ax_Ly = Ly + Ly/self.delta_factor
        axis.set_xlim(-ax_Lx + x0, ax_Lx + x0)
        axis.set_ylim(-ax_Ly + y0, ax_Ly + y0)
        
        # cehck requests as to density of lines
        if self.denser != 1:
            if self.form_0_str == None:
                # This cannot be done if a string has not been supplied
                # ERROR
                raise TypeError('Error: You need to supply the 0-form equation to do this, look at \'give_eqn\' method')
            else:
                # get the supplied form as a string
                zero_form_str = str(simplify(self.form_0_str))
                # set up grids for contours
                contour_x, contour_y = np.linspace(self.xg[0,0] , self.xg[0,-1] , self.pt_den_x*self.denser), np.linspace(self.yg[0,0] , self.yg[-1,0], self.pt_den_y*self.denser)
                contour_x_grid, contour_y_grid = np.meshgrid(contour_x, contour_y)
                # format the given ftring
                zero_form_str = zero_form_str.replace('x', 'contour_x_grid')
                zero_form_str = zero_form_str.replace('y', 'contour_y_grid')
                # evaluate bearing in mind zeros
                if zero_form_str.find('contour_x_grid') & zero_form_str.find('contour_y_grid') == -1:
                    form_0_contour = eval(zero_form_str)*np.ones(np.shape(contour_x_grid))
                else:
                    form_0_contour = eval(zero_form_str)
                    
                form_0 = form_0_contour
                xg = contour_x_grid
                yg = contour_y_grid
        else:
            form_0 = self.form_0
            xg = self.xg
            yg = self.yg
        
        
        # set all insignificant values to zero:
        form_0[np.abs(form_0) < 1e-15] = 0
                
        # deal with sinularities that appear on evaluated points
        isnan_arr = np.isnan(form_0)
        for i in range(len(xg[0, :])):
            for j in range(len(yg[:, 0])):
                # set to zero points that are not defined or inf
                if isnan_arr[j, i] or abs(form_0[j, i]) == np.inf or abs(form_0[j, i]) > 1e15:
                    # colour this region as a red dot, not square to
                    # not confuse with high mag 2-forms in stacks. or worse, in
                    # blocks
                    circ = patch.Circle((xg[j, i], yg[j, i]), Lx*0.05/3, color='red')
                    axis.add_patch(circ)
                    form_0[j, i] = 0
        
        if self.logarithmic_scale_bool:
            mag1 = np.abs(form_0) + 1
            form_0_norm = form_0/(mag1)
            logmag = np.log10(mag1)
            form_0 = form_0_norm*logmag

        else:
            pass
        
        CS = axis.contour(xg, yg, form_0, levels=self.lines, cmap=self.cmap)
        axis.clabel(CS, inline=self.inline_bool, fontsize=self.fontsize)
    
    # define a method to compute the exterior derivative
    def ext_d(self):
        '''
        Takes in no argument
        computes the exterior derivative and returns it as the 1-form object
        Returns 1 form object
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
            
            # supply these to the 1-form object function and return object
            result_1_form = form_1(self.xg, self.yg, form_1_x, form_1_y, form_1_x_unformated, form_1_y_unformated)
            return result_1_form
    
    # deifne a method to complete the exterior derivative numerically
    def num_ext_d(self, edge_order=1):
        '''
        Takes in 1 argument:
        -- edge_order: determines order same as in numpy gradient {1 or 2}
        
        Return 1 object - 1-form
        computes the exterior derivative numerically only
        The equations do not need to be given
        If given, they do not get passed onto the 1-form object anyway
        NUMERICAL ONLY
        '''
        
        # from numpy gradient, get the gradient array
        fy, fx = np.gradient(self.form_0, edge_order=edge_order)
        
        # supply these to the 1-form object function
        result_1_form = form_1(self.xg, self.yg, fx, fy)
        
        # return the new object to user
        return result_1_form
    
    # deinfe a method for Hodge of a 0-form
    def num_hodge(self):
        '''
        Takes in no arguments
        
        It calulates the Hodge on R^2 by the standard definition:
        1* = (dx^dy)
        Does so numerically via instance provided arrays
        IF equations were given, this method will lose them
        
        returns a 2-form
        '''
        # check if equations have been given:
        # if they have, doing it only numerically would create
        # a mismatch, avoid that
        if self.form_0_str != None:
            print('Warning: You supplied equations, doing it numerically only will lose these')
        
        # now complete the process numerically
        # pass these in to the object to create a new one and return
        new_object = form_2(self.xg, self.yg, self.form_0)  # N.B no equations to supply
        return new_object
    
    def hodge(self):
        '''
        Takes in no arguments
        
        It calulates the Hodge on R^2 by the standard definition:
        1* = (dx^dy)
        Does so analytically via instance provided equtions
        changes the equations AND the numerical answers
        
        returns a 2-form
        
        '''
        # can only be done if equations have been given, check:
        if self.form_0_str != None:
            # some equations are there, compute the Hodge on these:
            
            # get numerical solutions, evaulated on local strings
            # to relate parameter to the self grids and keep strings, because
            # need to supply these unformatted:
            form_2_str_unformated = self.form_0_str + '' 
            string_2_form = self.form_0_str  # to be formated
            # from these strings, get the numerical 2-form:
            string_2_form = string_2_form.replace('x', '(self.xg)')
            string_2_form = string_2_form.replace('y', '(self.yg)')
            
            if string_2_form.find('x') & string_2_form.find('y') == -1:
                string_2_form = '(' + str(string_2_form) + ')* np.ones(np.shape(self.xg))'
            
            # evaulated numerically
            form_2_result = eval(string_2_form)
            
            # create and return object
            new_object = form_2(self.xg, self.yg, form_2_result, form_2_eq=form_2_str_unformated)
            return new_object
        else:
            # ERROR
            raise TypeError('You need to supply the 2-form equation to do this, look at \'give_eqn\' method')

    # define a fucntion to compute a wedge product
    def wedge(self, form_second, degree=0, keep_object=False):
        '''
        Parameters:
        ----------------
        form_second - the form to wedge the 0-form with.
                    Can be supplied as a DFormPy instance, a tuple of equations,
                    or a single string equation depending on what form is to be
                    wedged.
                    To wedge with 1-form, supply 1-form instance, or tuple of
                    component equations as strings in terms of x and y.
                    To wedge with 0-form or 2-form, supply corresponding
                    instances or a single equation. When using equations,
                    to distinguish between them, provide parmater 'degree'.
        degree - default is 0. Only used when a single string is supplied
                    as form_second, to distinguish betwen 0-form and 2-form
                    for 0-form, degree=0, for 2-form, degree=2.
                    Determines what form is to be wegded with the
                    given 0-form.
        keep_object - bool -default=False - Only needed when 0-form /\ 0-form 
                    If False, a new object is created
                    as a result of the wedge. If True, the 0-form acted on
                    is modified to be the result of the wedge. 
        
        To do so here, strings for the form must be supplied.
        Computes the Wedge product using strings, ANALYTICALLY
        
        Returns:
        --------------
        Wedged with 0-form returns a 0-form object if keep_object is False
                    (default), and returns nothing when it is True
        Wedged with a 1-form, returns a 1-form instance
        Wedged with a 2-form, returns a 2-form instance
        
        '''
        
        # test if equations were given first:
        if self.form_0_str is None:
            raise ValueError('Error: You need to supply the 0-form equation to do this, look at \'give_eqn\' method')
        
        # set up variable to store order of supplied form, initially assume 1-form
        order = 0
        
        # get needed second obejct strings dep. on input
        if isinstance(form_second, tuple):
            # if equations were given here take these, if numerical grids were given - error!
            # check size , should be a 1-form
            if len(form_second) == 2:
                # 0-form/\1-form, check if strings supplied
                if isinstance(form_second[0], str) and isinstance(form_second[1], str):
                    to_wedge_x_2_str = form_second[0]
                    to_wedge_y_2_str = form_second[1]
                    order = 1
                else:
                    raise ValueError('for analytical calulation, supply 1-form equations as strings')
            else:
                raise ValueError('too many or too little equations given in tuple')
        elif isinstance(form_second, str):
            # single string, could be 0-form or 2-form, check given degree:
            if degree == 0:
                to_wedge_0_form_str = form_second
                order = 0
            elif degree == 2:
                to_wedge_2_form_str = form_second
                order = 2
            else:
                raise ValueError('not possible digree given or supplied one string for a 1-form')
        else:
            # object supplied, get numericals checking which object is given:
            if isinstance(form_second, form_1):
                if form_second.form_1_str_x is None or form_second.form_1_str_y is None:
                     raise ValueError('supplied 1-form instance must contain equations for analytical calculation')
                else:
                    to_wedge_x_2_str = form_second.form_1_str_x
                    to_wedge_y_2_str = form_second.form_1_str_y
                    order = 1
            elif isinstance(form_second, form_0):
                if form_second.form_0_str is None:
                    raise ValueError('supplied 0-form instance must contain equations for analytical calculation')
                else:
                    to_wedge_0_form_str = form_second.form_0_str
                    order = 0       
            elif isinstance(form_second, form_2):
                if form_second.form_2_str is None:
                    raise ValueError('supplied 2-form instance must contain equations for analytical calculation')
                else:
                    to_wedge_2_form_str = form_second.form_2_str
                    order = 2
            else:
                raise TypeError('Supplied form to wedge with is not recognised')
        
        # Deal with 0-form/\1-form:
        if order == 1:
            # first, find the result of the 1-form:
            new_str_x = str(simplify('(' + self.form_0_str + ')*(' +  to_wedge_x_2_str + ')'))
            new_str_y = str(simplify('(' + self.form_0_str + ')*(' +  to_wedge_y_2_str + ')'))
            # keep it as it is locally to supply it to object maker later
            form_1_str_x_loc = new_str_x + ''
            form_1_str_y_loc = new_str_y + ''
            # format it to be in terms of grids and:
            # check against constant and zero 1-forms being supplied
            # get the numerical evaluation of it
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

            # return the new one to the user:
            new_object = form_1(self.xg, self.yg, form_1_x, form_1_y, F_x_eqn=form_1_str_x_loc, F_y_eqn=form_1_str_y_loc)
            return new_object
        
        elif order == 0:
            form_0_str = str(simplify( '(' + self.form_0_str + ')*(' +  to_wedge_0_form_str + ')'))
            # keep it as it is locally to supply it to object maker later
            form_0_str_loc = form_0_str + ''
            # format it to be in terms of grids and:
            # check against constant and zero 2-forms being supplied
            # get the numerical evaluation of it
            form_0_str = form_0_str.replace('x', 'self.xg')
            form_0_str = form_0_str.replace('y', 'self.yg')
            if form_0_str.find('x') & form_0_str.find('y') == -1:
                form_0_str = '(' + str(form_0_str) + ')* np.ones(np.shape(self.xg))'
            
            # evaluate it numerically on the grid supplied
            form_0_result = eval(form_0_str)
            
            # depending on keep_object, return:
            if keep_object:
                self.form_0 = form_0_result
                self.form_0_str = form_0_str_loc
            elif not keep_object:
                new_object = form_0(self.xg, self.yg, form_0_result, form_0_str_loc)
                # return the new one to the user:
                return new_object
            else:
                raise ValueError('Error, Invalid input for \'keep_object\'')

        elif order is 2:
            form_2_str = str(simplify( '(' + self.form_0_str + ')*(' +  to_wedge_2_form_str + ')'))
            # keep it as it is locally to supply it to object maker later
            form_2_str_loc = form_2_str + ''
            # format it to be in terms of grids and:
            # check against constant and zero 2-forms being supplied
            # get the numerical evaluation of it
            form_2_str = form_2_str.replace('x', 'self.xg')
            form_2_str = form_2_str.replace('y', 'self.yg')
            if form_2_str.find('x') & form_2_str.find('y') == -1:
                form_2_str = '(' + str(form_2_str) + ')* np.ones(np.shape(self.xg))'
            
            # evaluate it numerically on the grid supplied
            form_2_result = eval(form_2_str)
            
            # create new instance and return to user
            new_object = form_2(self.xg, self.yg, form_2_result, form_2_str_loc)
            return new_object
        else:
            # should never happen, but in case
            raise ValueError('Variable change during code running, look at \'order\' parameter')
    
    # define a method for numerical wedge product
    def num_wedge(self, form_second, degree=0, keep_object=False):
        '''
        Parameters:
        ----------------
        form_second - the form to wedge the 0-form with.
                    Can be supplied as a DFormPy instance, a tuple of grids of
                    same size and dimensions as this 0-form,
                    or a single grid of scaling function values depending on
                    what form is to be wedged.
                    To wedge with 1-form, supply 1-form instance, or tuple of
                    component grids of same size as 1-form acted on.
                    To wedge with 0-form or 2-form, supply corresponding
                    instances or a single grid. When using grids,
                    to distinguish between them, provide parmater 'degree'.
        degree - default is 0. Only used when a single grid is supplied
                    as form_second, to distinguish betwen 0-form and 2-form
                    for 0-form, degree=0, for 2-form, degree=2.
                    Determines what form is to be wegded with the
                    given 0-form.
        keep_object - bool -default=False - only used when 0-form is wedged
                    with a 0-form. If False, a new object is created as 
                    a result of the wedge. If True, the 1-form acted on
                    is modified to be the result of the wedge. 
        
        Computes the Wedge product numerically
        
        Returns:
        --------------
        Wedged with 0-form returns a 0-form object if keep_object is False
                    (default), and returns nothing when it is True
        Wedged with a 1-form, returns a 1-form instance
        Wedged with a 2-form, returns a 2-form instance
        
        '''
        
        # test if equations were given first:
        if self.form_0_str is None:
            pass
        else:
            print('The first 0-form you are completing the wedge with has equations supplied, these will be lost')
        
        # set up variable to store order of supplied form, initially assume 0-form
        order = 0
        
        # get needed second obejct grids dep. on input
        if isinstance(form_second, tuple):
            # check size to see what it is to be wedged with.
            # tuple should only be length 2 --> 1-form/\1-form
            if len(form_second) == 2:
                # 0-form/\1-form, extract components
                # if numerical grids were given, take these, if equations, change to values on grids:
                if isinstance(form_second[0], str) and isinstance(form_second[1], str):
                    new_str_x = form_second[0].replace('x', '(self.xg)')
                    new_str_x = new_str_x.replace('y', '(self.yg)')
                    new_str_y = form_second[1].replace('x', '(self.xg)')
                    new_str_y = new_str_y.replace('y', '(self.yg)')
                    if new_str_x.find('x') & new_str_x.find('y') == -1:
                        new_str_x = '(' + str(new_str_x) + ')* np.ones(np.shape(self.xg))'
                    if new_str_y.find('x') & new_str_y.find('y') == -1:
                        new_str_y = '(' + str(new_str_y) + ')* np.ones(np.shape(self.yg))'
                    f12_x = eval(new_str_x)
                    f12_y = eval(new_str_y)
                    order = 1
                elif isinstance(form_second[0], np.ndarray) and isinstance(form_second[1], np.ndarray):
                    f12_x = form_second[0]
                    f12_y = form_second[1]
                    order = 1
                else:
                    raise ValueError('Not recognised input tuple')
            else:
                raise ValueError('too many or too little equations given in tuple')
        
        elif isinstance(form_second, np.ndarray):
            # check degree:
            if degree == 0:
                to_wedge_0_form = form_second
                order = 0
            elif degree == 1:
                raise ValueError('for degree 1, supply a 1-form, not a single grid')
            elif degree == 2:
                to_wedge_2_form = form_second
                order = 2
        
        elif isinstance(form_second, str):
            # single string, could be 0-form or 2-form, check given degree:
            if degree == 0:
                    str_0_form = form_second.replace('x', '(self.xg)')
                    str_0_form = str_0_form.replace('y', '(self.yg)')
                    if str_0_form.find('x') & str_0_form.find('y') == -1:
                        str_0_form = '(' + str(str_0_form) + ')* np.ones(np.shape(self.xg))'
                    
                    to_wedge_0_form = eval(str_0_form)
                    order = 0
            elif degree == 2:
                str_2_form = form_second.replace('x', '(self.xg)')
                str_2_form = str_2_form.replace('y', '(self.yg)')
                if str_2_form.find('x') & str_2_form.find('y') == -1:
                    str_2_form = '(' + str(str_2_form) + ')* np.ones(np.shape(self.xg))'
                
                to_wedge_2_form = eval(str_2_form)
                order = 2
            else:
                raise ValueError('not possible digree given or supplied one string for a 1-form')
        
        # object supplied, get grids checking which object is given:
        
        elif isinstance(form_second, form_1):
            f12_x = form_second.F_x
            f12_y = form_second.F_y
            order = 1
        elif isinstance(form_second, form_0):
            to_wedge_0_form = form_second.form_0
            order = 0
        elif isinstance(form_second, form_2):
            order = 2
            to_wedge_2_form = form_second.form_2
        else:
            raise TypeError('Supplied form to wedge with is not recognised')
        
        # Use given inputs to evaluate the result:
        
        # Deal with 0-form/\1-form:
        if order == 1:
            # first, find the result of the 1-form
            new_form_1_x = self.form_0 * f12_x
            new_form_1_y = self.form_0 * f12_y
            
            # create instance and return
            new_object = form_1(self.xg, self.yg, new_form_1_x, new_form_1_y)
            return new_object
        
        elif order == 0:
            # from these get the numerical 0-form
            form_0_result = self.form_0 * to_wedge_0_form
            
            # depending on keep_object, return:
            if keep_object:
                self.form_0 = form_0_result
            elif not keep_object:
                new_object = form_0(self.xg, self.yg, form_0_result)
                # return the new one to the user:
                return new_object
            else:
                raise ValueError('Error, Invalid input for \'keep_object\'')

        elif order == 2:
            # from these get the numerical 0-form
            form_2_result = self.form_0 * to_wedge_2_form
            
            # create instance and return
            new_object = form_2(self.xg, self.yg, form_2_result)
            return new_object
        else:
            # should never happen, but in case
            raise ValueError('Variable change during code running, look at \'order\' parameter')
        
# %%

'''

function to create a vector field object and define methods for it

'''

# define a function that will set up a vector field object that can be customised and
# plotted
class vector_field():
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
    # set up all variables
    def __init__(self, xg, yg, F_x, F_y, F_x_eqn=None, F_y_eqn=None):
        self.xg = xg
        self.yg = yg
        self.F_x = F_x
        self.F_y = F_y
        self.pt_den = len(xg[:, 0])  # + 1 , assume square grids
        self.orientation = 'mid'
        self.scale = 1
        self.color = 'black'
        self.logarithmic_scale_bool = 0
        # self.base = 10
        self.scale_bool = True
        self.delta_factor = 10
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
        # self.base = base
    
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
    def set_density(self, points_number):
        '''
        Changes the size of stack in direction perp. to VF
        It is done in in terms of the fraction of plot size
        Note, not strictly needed, can change it by instance.fract(fraction)
        
        Parmaeters:
        ---------------
        fraction - float/int - size of stack in terms of the fraction of plot size
        
        Returns: None
        '''
        if self.str_x == None or self.str_y == None:
            # Error
            raise ValueError('Error: You need to supply the 1-form equation to do this, look at \'give_eqn\' method')
        else:
            # redefine the grids
            x = np.linspace(self.xg[0,0], self.xg[0,-1], points_number)
            y = np.linspace(self.yg[0,0], self.yg[-1,0], points_number)
            self.xg, self.yg = np.meshgrid(x,y)
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
    def plot(self, axis):
        '''
        Finilises the plotting
        Uses the attribues of the object as set originally and as customised
        with methods to create a plot of the VF
        Takes in 1 argument:
        --- axis - matplotlib axes instance, plots on these
        
        No Returns    
        
        '''
        
        # get the lengths of x and y from their grids
        x_len = len(self.xg[:, 0])
        y_len = len(self.yg[0, :])
        
        # Extract L from the x and y grids
        Lx = 0.5*(self.xg[0, -1] - self.xg[0, 0])
        Ly = 0.5*(self.yg[-1, 0] - self.yg[0, 0])
        L = 0.5*(Lx + Ly)
        x0 = self.xg[0, 0] + Lx
        y0 = self.yg[0, 0] + Ly
        
        # reset axis limits
        ax_Lx = Lx + Lx/self.delta_factor
        ax_Ly = Ly + Ly/self.delta_factor
        axis.set_xlim(-ax_Lx + x0, ax_Lx + x0)
        axis.set_ylim(-ax_Ly + y0, ax_Ly + y0)
        
        # for arrows to work, with nan and infs
        # make a local variable of F_x and F_y
        # so that thye don't alter globally
        F_x_local = self.F_x * 1
        F_y_local = self.F_y * 1
        
        # prevent any magnitudes from being inf or nan
        # only here, need to do it to u and v not just mag
        
        # find the distance between neightbouring points on the grid
        dist_points = self.xg[0, 1] - self.xg[0, 0]
        
        # deal with infs and nans in mag
        isnan_arrx = np.isnan(F_x_local)
        isnan_arry = np.isnan(F_y_local)
        for i in range(x_len):
            for j in range(y_len):
                # set to zero points that are not defined or inf
                if isnan_arrx[i, j] or isnan_arry[i, j]:
                    #colour this region as a shaded square
                    rect = patch.Rectangle((self.xg[i, j] - dist_points/2, self.yg[i, j]  - dist_points/2), dist_points, dist_points, color='#B5B5B5')
                    axis.add_patch(rect)
                    F_x_local[i,j] = F_y_local[i,j] = 0
                if abs(F_x_local[i, j]) == np.inf or abs(F_y_local[i, j]) == np.inf or abs(F_y_local[i, j]) > 1e15 or abs(F_x_local[i, j]) > 1e15:
                    # colour this point as a big red dot
                    circ = patch.Circle((self.xg[i, j], self.yg[i, j]), Lx*0.05/3, color='red')
                    axis.add_patch(circ)
                    F_x_local[i,j] = F_y_local[i,j] = 0
#            isnan_arrx = np.isnan(F_x_local)
#            isnan_arry = np.isnan(F_y_local)
#            for i in range(x_len):
#                for j in range(y_len):
#                    if isnan_arrx[i,j] or isnan_arry[i,j] or abs(F_x_local[i, j]) == np.inf or abs(F_y_local[i, j]) == np.inf or abs(F_y_local[i, j]) > 1e15 or abs(F_x_local[i, j]) > 1e15:
#                        
#                        F_x_local[i,j] = F_y_local[i,j] = 0

        # set all insignificant values to zero:
        F_x_local[np.abs(F_x_local) < 1e-15] = 0
        F_y_local[np.abs(F_y_local) < 1e-15] = 0
        
        # find the magnitude corresponding to each point and store in mag array
        mag = np.sqrt(F_x_local**2 + F_y_local**2)
        
        # find the maximum magnitude for scaling
        max_size = np.max(mag)   # careful with singularities, else ---> nan
        
        # Rescale components if log scaling is selected
        if self.logarithmic_scale_bool:
            mag1 = mag + 1
            # min_size = np.min(mag1)
            
            unorm = F_x_local/mag1
            vnorm = F_y_local/mag1
            
            # logsf = np.log10(mag1/min_size)
            logmag = np.log10(mag1)
            F_x_local = unorm*logmag
            F_y_local = vnorm*logmag
            
            mag = np.sqrt(F_x_local**2 + F_y_local**2)
            max_size = np.max(mag)
            
        # deal with requested autoscaling
        if self.scale_bool is False:
            ScaleFactor = self.scale
        elif self.scale_bool is True:
            ScaleFactor = max_size/(0.9*(2*Lx/self.pt_den))
        
        # plot using matplotlib quiver
        axis.quiver(self.xg, self.yg, F_x_local, F_y_local, pivot=self.orientation, scale=ScaleFactor, scale_units='xy', color=self.color) 
    
    
    def zoom(self, target=[0, 0], mag=2, dpd=9, inset=True, axis=None, insize=0.3):
        '''
        Create a new window which displays the field zoomed at a certain point
        User gives arguments
        Target: Determines the zoom location, coordinates
        mag: +ve float, determines zooming amount
        dpd: +int, determines how many points on each axis
        
        inset - bool - if true, zoomed field plotted on given axis
        axis - matplotlib axes instance - axis to plot on if instance it True
        insize - float - size of inset as fraction of total figure
        
        
        Returns:
        --------
        if inset is False:
            zoomed in VF object
        if inset is True, inset axis and zoomed in VF object in this order.
        
        '''
        
        # Requires user to provide eqn of the 1-form they are zooming on.
        
        if self.str_x == None or self.str_y == None:
            # ERROR
            raise TypeError('No equation provided')
        else:
        
            # Zoom must be one or greater
            if mag < 1:
                raise ValueError('Zoom must be greater than one')
            else:
                
                if insize > 1 or insize < 0:
                    raise ValueError('Insize must be +ve and less than one')
                else:
                    
                    # If no inset, set the size of the zoom axis to allow normal plotting
                    if inset == False:
                        insize = 1
                        
                    # Target coordinates
                    x_m = target[0]
                    y_m = target[1]
                    
                    # Get the size of the original VF
                    Lx = 0.5*(self.xg[0,-1] - self.xg[0,0])
                    Ly = 0.5*(self.yg[-1,0] - self.yg[0,0])
                    
                    # Zoom axis range
                    d_range_x = insize*Lx/mag
                    d_range_y = insize*Ly/mag
                    
                    # Set up zoom window grids
                    dx = np.linspace(-d_range_x + x_m, d_range_x + x_m, dpd)
                    dy = np.linspace(-d_range_y + y_m, d_range_y + y_m, dpd)
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
                    
                    # from that create VF instance
                    zoom_vf = vector_field(dxg, dyg, u_zoom, v_zoom, self.str_x, self.str_y)
                    
                    q = 1
                    
                    xi = (q*x_m - self.xg[0,0])/(2*Lx)
                    yi = (q*y_m - self.yg[0,0])/(2*Ly)
                    
                    # depending on preferances, return to user and plot
                    if inset == True:
                        if axis != None:
                            # Create inset axis in the current axis.
                            
                            zoom_inset_ax = axis.inset_axes([(xi - 0.5*insize), (yi - 0.5*insize), insize, insize])
                            zoom_vf.plot(zoom_inset_ax)
                            
                            # return the zoomed on axis
                            # also return zoomed in form in case user wants that.
                            return zoom_inset_ax, zoom_vf
                        else:
                            raise ValueError('Cannot inset without supplied axis')
                    else:
                        # inset is false, just return the new zoomed in instance
                        return zoom_vf
        
    def deriv(self, target=[0, 0], mag=2, dpd=9, inset=False, axis=None, insize=0.3):
        '''
        Creates new vector field object at a target location, showing the derivative field at this point.
        User gives arguments:
        Target - derivative plot location
        mag - Magnification level
        dpd - New plot point density
        
        inset - bool - if true, field deriv is plotted on given axis
        axis - matplotlib axes instance - axis to plot on if instance it True
        
        Returns:
        --------
        if inset is False:
            deriv VF object
        if inset is True, inset axis and deriv VF object in this order.
        
        '''
        if self.str_x == None or self.str_y == None:
            # ERROR
            raise TypeError('No equation provided')
        else:
             # Zoom must be one or greater
            if mag < 1:
                raise ValueError('Zoom must be greater than one')
            else:
                
                if insize > 1 or insize < 0:
                    raise ValueError('Insize must be +ve and less than one')
                else:
                    
                    # If no inset, set the size of the zoom axis to allow normal plotting
                    if inset == False:
                        insize = 1
            
                    # Target coordinates
                    x_m = target[0]
                    y_m = target[1]
                    
                    # Get the size of the original VF
                    Lx = 0.5*(self.xg[0,-1] - self.xg[0,0])
                    Ly = 0.5*(self.yg[-1,0] - self.yg[0,0])
                    
                    # Zoom axis range
                    d_range_x = insize*Lx/mag
                    d_range_y = insize*Ly/mag
                    
                    # Set up zoom window grids
                    dx = np.linspace(-d_range_x + x_m, d_range_x + x_m, dpd)
                    dy = np.linspace(-d_range_y + y_m, d_range_y + y_m, dpd)
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
                    
                    # from that create VF instance
                    deriv_vf = vector_field(dxg, dyg, U, V, self.str_x, self.str_y)
                    
                    q = 1
                    
                    # Coordinates for plotting the inset axis
                    xi = (q*x_m - self.xg[0,0])/(2*Lx)
                    yi = (q*y_m - self.yg[0,0])/(2*Ly)
                    
                    # depending on preferances, return to user and plot
                    if inset == True:
                        if axis != None:
                            # Create inset axis in the current axis.
                            
                            deriv_inset_ax = axis.inset_axes([(xi - 0.5*insize), (yi - 0.5*insize), insize, insize])
                            deriv_vf.plot(deriv_inset_ax)
                            
                            # return the zoomed on axis
                            # also return zoomed in form in case user wants that.
                            return deriv_inset_ax, deriv_vf
                        else:
                            raise ValueError('Cannot inset without supplied axis')
                    else:
                        # inset is false, just return the new zoomed in instance
                        return deriv_vf
    
        
    def div(self, target=[0,0], mag=2, dpd=9, inset=False, axis=None, insize=0.3):
        '''
        Creates new vector field object at a target location, showing the Divergence of the field at this point.
        User gives arguments:
        Target - derivative plot location
        Zoom - Magnification level
        dpd - New plot point density
        
        inset - bool - if true, field div is plotted on given axis
        axis - matplotlib axes instance - axis to plot on if instance it True
        
        Returns:
        --------
        if inset is False:
            div VF object
        if inset is True, inset axis and div VF object in this order.
        
        '''
        if self.str_x == None or self.str_y == None:
            # ERROR
            raise TypeError('No equation provided')
        else:
             # Zoom must be one or greater
            if mag < 1:
                raise ValueError('Zoom must be greater than one')
            else:
                
                if insize > 1 or insize < 0:
                    raise ValueError('Insize must be +ve and less than one')
                else:
                    
                    # If no inset, set the size of the zoom axis to allow normal plotting
                    if inset == False:
                        insize = 1
            
                    # Target coordinates
                    x_m = target[0]
                    y_m = target[1]
                    
                    # Get the size of the original VF# Get the size of the original VF
                    Lx = 0.5*(self.xg[0,-1] - self.xg[0,0])
                    Ly = 0.5*(self.yg[-1,0] - self.yg[0,0])
                    
                    # Zoom axis range
                    d_range_x = insize*Lx/mag
                    d_range_y = insize*Ly/mag
                    
                    # Set up zoom window grids
                    dx = np.linspace(-d_range_x + x_m, d_range_x + x_m, dpd)
                    dy = np.linspace(-d_range_y + y_m, d_range_y + y_m, dpd)
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
                    
                    
                   # from that create VF instance
                    div_vf = vector_field(dxg, dyg, U_div, V_div, self.str_x, self.str_y)
                    
                    q = 1
                    
                    # Coordinates for plotting the inset axis
                    xi = (q*x_m - self.xg[0,0])/(2*Lx)
                    yi = (q*y_m - self.yg[0,0])/(2*Ly)
                    
                    # depending on preferances, return to user and plot
                    if inset == True:
                        if axis != None:
                            # Create inset axis in the current axis.
                            
                            div_inset_ax = axis.inset_axes([(xi - 0.5*insize), (yi - 0.5*insize), insize, insize])
                            div_vf.plot(div_inset_ax)
                            
                            # return the zoomed on axis
                            # also return zoomed in form in case user wants that.
                            return div_inset_ax, div_vf
                        else:
                            raise ValueError('Cannot inset without supplied axis')
                    else:
                        # inset is false, just return the new zoomed in instance
                        return div_vf
        
    def curl(self, target=[0,0], mag=2, dpd=9, inset=False, axis=None, insize=0.3):
        '''
        Creates new vector field object at a target location, showing local rotation (Curl)
        User gives arguments:
        Target - derivative plot location
        Zoom - Magnification level
        dpd - New plot point density
        
        inset - bool - if true, field curl is plotted on given axis
        axis - matplotlib axes instance - axis to plot on if instance it True
        
        Returns:
        --------
        if inset is False:
            div VF object
        if inset is True, inset axis and curl VF object in this order.
        
        '''
        if self.str_x == None or self.str_y == None:
            # ERROR
            raise TypeError('No equation provided')
        else:
             # Zoom must be one or greater
            if mag < 1:
                raise ValueError('Zoom must be greater than one')
            else:
                
                if insize > 1 or insize < 0:
                    raise ValueError('Insize must be +ve and less than one')
                else:
                    
                    # If no inset, set the size of the zoom axis to allow normal plotting
                    if not isinstance(inset, float) and not isinstance(inset, int):
                        insize = 0.4
            
                    # Target coordinates
                    x_m = target[0]
                    y_m = target[1]
                    
                    # Get the size of the original VF# Get the size of the original VF
                    Lx = 0.5*(self.xg[0,-1] - self.xg[0,0])
                    Ly = 0.5*(self.yg[-1,0] - self.yg[0,0])
                    
                    # Zoom axis range
                    d_range_x = insize*Lx/mag
                    d_range_y = insize*Ly/mag
                    
                    # Set up zoom window grids
                    dx = np.linspace(-d_range_x + x_m, d_range_x + x_m, dpd)
                    dy = np.linspace(-d_range_y + y_m, d_range_y + y_m, dpd)
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
                        
                    # from that create VF instance
                    curl_vf = vector_field(dxg, dyg, U_curl, V_curl, self.str_x, self.str_y)
                    
                    q = 1
                    
                    # Coordinates for plotting the inset axis
                    xi = (q*x_m - self.xg[0,0])/(2*Lx)
                    yi = (q*y_m - self.yg[0,0])/(2*Ly)
                    
                    # depending on preferances, return to user and plot
                    if inset == True:
                        if axis != None:
                            # Create inset axis in the current axis.
                            
                            curl_inset_ax = axis.inset_axes([(xi - 0.5*insize), (yi - 0.5*insize), insize, insize])
                            curl_vf.plot(curl_inset_ax)
                            
                            # return the zoomed on axis
                            # also return zoomed in form in case user wants that.
                            return curl_inset_ax, curl_vf
                        else:
                            raise ValueError('Cannot inset without supplied axis')
                    else:
                        # inset is false, just return the new zoomed in instance
                        return curl_vf
    
    # define a method to change a supplied Vector filed to the 1-form
    def covariant(self, g=[['1', '0'], ['0', '1']]):
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
            raise TypeError('Metric components are inconsistent')
        
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

        # based on what was given into the Vector field
        # return a 1-form object with these parameters
        if analytics:
            result_form = form_1(self.xg, self.yg, form_x, form_y, x_str_form, y_str_form)
        elif not analytics:
            result_form = form_1(self.xg, self.yg, form_x, form_y)
    
        # return the found object
        return result_form
