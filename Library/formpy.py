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
    if c == 0:
        return ((2*s + 1)/(2*(n-1)))
    else:
        return (s/(n-1))


# define a function taht will deal with plotting 1-forms:
def form_1_plot(xg, yg, F_x, F_y):
    # define the object that will call stackplot
    class form_1():
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
            self.pt_den = len(xg[:, 0]) + 1  # assume square grids
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
        
        # write some functions that will allow the user to chenge some of the above
        # variables
        
        # define a method to change figure size
        def fig_size(self, n, m):
            self.figure.set_size_inches(n, m)
            # then change figure
        
        # change colour
        def colour(self, color):
            self.color = str(color)
        
        # change arrowsheads
        def arrow_heads(self):
            self.arrowheads = not self.arrowheads
        
        # change w_head
        def head_width(self, wide):
            self.w_head = float(wide)
        
        # change h_head
        def head_height(self, high):
            self.h_head = float(high)
        
        # change orientation:
        def orient(self, string):
            self.orientation = str(string)
        
        # plot stacks boolean
        def stacks(self):
            self.stack_bool = not self.stack_bool
        
        # plot arrows boolean
        def arrows(self):
            self.arrow_bool = not self.arrow_bool
        
        # change boolean that det. if to sclae logarithmically
        def log_scaling(self):
            self.logarithmic_scale_bool = not self.logarithmic_scale_bool
        
        # define a method to be able to change bool that det. if arrows autoscale
        def autoscale(self):
            self.scale_bool = not self.scale_bool
        
        # define methods to change s_max and s_min
        def max_sheets(self, maximum):
            self.s_max = maximum
        
        # define method to change fraction of sheetsize w.r.t graoh size:
        def sheet_size(self, fraction):
            self.fract = fraction
        
        #define a method to change spare spacing around figure
        def surround_space(self, delta_denominator):
            self.delta_factor = delta_denominator
        
        # stcakplot: but it takes the above defined variables:
        def plot(self):
            
            # from self, get axis
            axis = self.axis
            
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
            
            # define length of sheet as a fraction of total graph scale
            sheet_L = L * self.fract
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
            ScaleFactor = max_size/(0.9*(2*L/self.pt_den))
            
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
                A_x = self.xg + (sheet_L/2)*I_sin
                A_y = self.yg - (sheet_L/2)*I_cos
                B_x = self.xg - (sheet_L/2)*I_sin
                B_y = self.yg + (sheet_L/2)*I_cos
                
                # define points of stack arrowheads as arrays for all stacks
                p_sh1x = self.xg + (s_L/2)*I_cos + (sheet_L*self.w_head)*I_sin
                p_sh1y = self.yg + (s_L/2)*I_sin - (sheet_L*self.w_head)*I_cos
                p_sh2x = self.xg + (s_L/2)*I_cos - (sheet_L*self.w_head)*I_sin
                p_sh2y = self.yg + (s_L/2)*I_sin + (sheet_L*self.w_head)*I_cos
                p_sh3x = self.xg + (s_L*0.5 + s_L*self.h_head)*I_cos
                p_sh3y = self.yg + (s_L*0.5 + s_L*self.h_head)*I_sin
                
                # define these for when there is only 1 line in the stack plot:
                P_sh1x = self.xg + (sheet_L*self.w_head)*I_sin
                P_sh1y = self.yg - (sheet_L*self.w_head)*I_cos
                P_sh2x = self.xg - (sheet_L*self.w_head)*I_sin
                P_sh2y = self.yg + (sheet_L*self.w_head)*I_cos
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
    form_1_object = form_1(xg, yg, F_x, F_y)
    #return it to user
    return form_1_object


# %%










