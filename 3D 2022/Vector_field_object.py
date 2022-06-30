import numpy as np
from mayavi import mlab
import matplotlib.pyplot as plt
import sympy
from sympy import sympify
from sympy import diff, simplify
from sympy.parsing.sympy_parser import parse_expr

from numpy import sin, cos, tan, sqrt, log, arctan, arcsin, arccos, tanh
from numpy import sinh, cosh, arcsinh, arccosh, arctanh, exp, pi, e


class vector_field3():
    '''
   
    '''

    # set up all variables
    def __init__(self, xg, yg, zg, F_x, F_y, F_z, F_x_eqn=None, F_y_eqn=None, F_z_eqn=None):

        self.xg = xg
        self.yg = yg
        self.zg = zg
        self.F_x = F_x
        self.F_y = F_y
        self.F_z = F_z
        self.pt_den = len(xg[:, 0])  # + 1 , assume square grids
        self.orientation = 'mid'
        self.scale = 1
        self.color = 'black'
        self.logarithmic_scale_bool = 0
        # self.base = 10
        self.scale_bool = True
        self.delta_factor = 10
        if F_x_eqn is not None:
            self.str_x = str(simplify(F_x_eqn))  # to start with, user must change to access some methods
            # Note, the string must be given with x, y, z as variables
        else:
            self.str_x = None
        
        if F_y_eqn is not None:
            self.str_y = str(simplify(F_y_eqn))
        else:
            self.str_y = None

        if F_z_eqn is not None:
            self.str_z = str(simplify(F_z_eqn))
        else:
            self.str_z = None


    def give_eqn(self, equation_str_x, equation_str_y, equation_str_z):
        '''
        Takes in 1-argument, string
        This must be the equation of the supplied numerical 0-form
        It must be in terms of x and y.
        Has to be given, for some methods to be calculatable.
        '''
        # set equation parameters to simplified inputs
        self.str_x = str(simplify(equation_str_x))
        self.str_y = str(simplify(equation_str_y))
        self.str_z = str(simplify(equation_str_z))
        # make the values match automatically to limit how often mismatch occurs
        # substitute these into the equation:
        # but keep it local
        str_x = self.str_x + ''
        str_y = self.str_y + ''
        str_z = self.str_z + ''

        str_x = str_x.replace('x', '(self.xg)')
        str_x = str_x.replace('y', '(self.yg)')
        str_x = str_x.replace('z', '(self.zg)')

        str_y = str_y.replace('x', '(self.xg)')
        str_y = str_y.replace('y', '(self.yg)')
        str_y = str_y.replace('z', '(self.zg)')

        str_z = str_z.replace('x', '(self.xg)')
        str_z = str_z.replace('y', '(self.yg)')
        str_z = str_z.replace('z', '(self.zg)')

        # check kagainst constant form components:
        if str_x.find('x') & str_x.find('y') & str_x.find('z') == -1:
            str_x = '(' + str(str_x) + ')* np.ones(np.shape(self.xg))'
        if str_y.find('x') & str_y.find('y') & str_y.find('z') == -1:
            str_y = '(' + str(str_y) + ')* np.ones(np.shape(self.yg))'
        if str_z.find('x') & str_z.find('y') & str_z.find('z') == -1:
            str_z = '(' + str(str_z) + ')* np.ones(np.shape(self.zg))'
        
        
        # re-evaluate the 2-form numerically
        self.F_x = eval(str_x)
        self.F_y = eval(str_y)
        self.F_z = eval(str_z)
        
    def return_string(self):
        '''
        Takes in no arguments, returns the unformatted strings back to user
        This is done in case user wants to access strings
        that got here not by input but by ext. alg.
        '''
        return self.str_x, self.str_y, self.str_z


    def plot(self):
        '''
        Finilises the plotting
        Uses the attribues of the object as set originally and as customised
        with methods to create a plot of the VF
        Takes in 1 argument:
        --- axis - matplotlib axes instance, plots on these
        
        No Returns    
        
        '''
        
        # get the lengths of x and y from their grids
        x_len = len(self.xg[:, 0, 0])
        y_len = len(self.yg[0, :, 0])
        z_len = len(self.zg[0, 0, :])
        
        # Extract L from the x and y grids
        Lx = 0.5*(self.xg[0, -1, -1] - self.xg[0, 0, 0])
        Ly = 0.5*(self.yg[-1, 0, -1] - self.yg[0, 0, 0])
        Lz = 0.5*(self.yg[-1, -1, 0] - self.yg[0, 0, 0])
        
        L = 0.5*(Lx + Ly + Lz)
        x0 = self.xg[0, 0, 0] + Lx
        y0 = self.yg[0, 0, 0] + Ly
        z0 = self.zg[0, 0, 0] + Lz
        
        # reset axis limits
        ax_Lx = Lx + Lx/self.delta_factor
        ax_Ly = Ly + Ly/self.delta_factor
        ax_Lz = Lz + Lz/self.delta_factor
        '''
        axis.set_xlim(-ax_Lx + x0, ax_Lx + x0)
        axis.set_ylim(-ax_Ly + y0, ax_Ly + y0)
        axis.set_zlim(-ax_Lz + z0, ax_Lz + z0)
        '''
        # for arrows to work, with nan and infs
        # make a local variable of F_x and F_y
        # so that thye don't alter globally
        F_x_local = self.F_x * 1
        F_y_local = self.F_y * 1
        F_z_local = self.F_z * 1
        
        # prevent any magnitudes from being inf or nan
        # only here, need to do it to u and v not just mag
        
        # find the distance between neightbouring points on the grid
        dist_points = self.xg[0, 1] - self.xg[0, 0]
        
        # deal with infs and nans in mag
        isnan_arrx = np.isnan(F_x_local)
        isnan_arry = np.isnan(F_y_local)
        isnan_arrz = np.isnan(F_z_local)
        """
        for i in range(x_len):
            for j in range(y_len):
                for k in range(z_len):
                    # set to zero points that are not defined or inf
                    if isnan_arrx[i, j, k] or isnan_arry[i, j, k] or isnan_arrz[i, j, k]:
                        #colour this region as a shaded square
                        rect = patch.Rectangle((self.xg[i, j] - dist_points/2, self.yg[i, j]  - dist_points/2), dist_points, dist_points, color='#B5B5B5')
                        axis.add_patch(rect)
                        F_x_local[i,j] = F_y_local[i,j] = 0
                    if abs(F_x_local[i, j]) == np.inf or abs(F_y_local[i, j]) == np.inf or abs(F_y_local[i, j]) > 1e15 or abs(F_x_local[i, j]) > 1e15:
                        # colour this point as a big red dot
                        circ = patch.Circle((self.xg[i, j], self.yg[i, j]), Lx*0.05/3, color='red')
                        axis.add_patch(circ)
                        F_x_local[i,j] = F_y_local[i,j] = 0
        """
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
        F_z_local[np.abs(F_z_local) < 1e-15] = 0
        
        # find the magnitude corresponding to each point and store in mag array
        mag = np.sqrt(F_x_local**2 + F_y_local**2 + F_z_local**2)
        
        # find the maximum magnitude for scaling
        max_size = np.max(mag)   # careful with singularities, else ---> nan
        
        # Rescale components if log scaling is selected
        if self.logarithmic_scale_bool:
            mag1 = mag + 1
            # min_size = np.min(mag1)
            
            unorm = F_x_local/mag1
            vnorm = F_y_local/mag1
            wnorm = F_z_local/mag1
            
            # logsf = np.log10(mag1/min_size)
            logmag = np.log10(mag1)
            F_x_local = unorm*logmag
            F_y_local = vnorm*logmag
            F_z_local = wnorm*logmag
            
            mag = np.sqrt(F_x_local**2 + F_y_local**2 + F_z_local**2)
            max_size = np.max(mag)
            
        # deal with requested autoscaling
        if self.scale_bool is False:
            ScaleFactor = self.scale
        elif self.scale_bool is True:
            ScaleFactor = max_size/(0.9*(2*Lx/self.pt_den))
        
        # plot using matplotlib quiver
        
        mlab.figure(bgcolor = (1,1,1), fgcolor = (0,0,0))
        mlab.quiver3d(self.xg, self.yg, self.zg, F_x_local, F_y_local, F_z_local, colormap='jet')
        cbar = mlab.colorbar(orientation='vertical')
        axes = mlab.axes(color = (0,0,0), nb_labels = 5)
        axes.label_text_property.font_family = 'courier'
        axes.label_text_property.font_size = 1
        axes.title_text_property.font_family = 'times'
        axes.title_text_property.font_size = 3

        cbar.label_text_property.font_family = 'courier'
        cbar.label_text_property.font_size = 1

        mlab.show()
        




















