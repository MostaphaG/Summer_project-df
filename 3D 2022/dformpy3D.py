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
    The class which creates a vector field object. 
    Requires grid (3D), field components.
    in order to calculate  curl and div, it is necessary to
    provide string expressions for the field component to
    the object as well
    .plot() can be used to plot out the vector field and the curl
    (either on the same axes or different)

    Methods: (applied to the object)

    .give_eqn(equation_str_x, equation_str_y, equation_str_z) 
    .curl()
    .div()
    .plot() - applied to either curl of object
    [if applied to object can use add_curl ='yes' argument to plot
    curl and field on the same axes]
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

    def curl(self):

        if self.str_x == None or self.str_y == None or self.str_z == None:
            # ERROR
            raise TypeError('No equation provided')
        else:

            Fx = parse_expr(self.str_x)
            Fy = parse_expr(self.str_y)
            Fz = parse_expr(self.str_z)        

            # differentiate expressions w.r.t x, y, z.
            ddx_Fx = Fx.diff(sympy.symbols('x'))
            ddy_Fx = Fx.diff(sympy.symbols('y'))
            ddz_Fx = Fx.diff(sympy.symbols('z'))

            ddx_Fy = Fy.diff(sympy.symbols('x'))
            ddy_Fy = Fy.diff(sympy.symbols('y'))
            ddz_Fy = Fy.diff(sympy.symbols('z'))

            ddx_Fz = Fz.diff(sympy.symbols('x'))
            ddy_Fz = Fz.diff(sympy.symbols('y'))
            ddz_Fz = Fz.diff(sympy.symbols('z'))
            # Define a string expression for Curl of the input vector field
            CurlX = (str(simplify(ddy_Fz)) + '-' + str(simplify(ddz_Fy)))
            CurlY = (str(simplify(ddz_Fx)) + '-' + str(simplify(ddx_Fz)))
            CurlZ = (str(simplify(ddx_Fy)) + '-' + str(simplify(ddy_Fx)))
        
            # replace input coordinates with coordinate point lists
            CurlX = CurlX.replace('x', '(self.xg)')
            CurlX = CurlX.replace('y', '(self.yg)')
            CurlX = CurlX.replace('z', '(self.zg)')

            CurlY = CurlY.replace('x', '(self.xg)')
            CurlY = CurlY.replace('y', '(self.yg)')
            CurlY = CurlY.replace('z', '(self.zg)')

            CurlZ = CurlZ.replace('x', '(self.xg)')
            CurlZ = CurlZ.replace('y', '(self.yg)')
            CurlZ = CurlZ.replace('z', '(self.zg)')

            # if input scalar, define a scalar list the size of a coordinate obj.
            if CurlX.find('x') & CurlX.find('y') & CurlX.find('z') == -1:
                CurlX = '(' + str(CurlX) + ')* np.ones(np.shape(self.xg))'
            if CurlY.find('x') & CurlY.find('y') & CurlY.find('z') == -1:
                CurlY = '(' + str(CurlY) + ')* np.ones(np.shape(self.yg))'
            if CurlZ.find('x') & CurlZ.find('y') & CurlZ.find('z') == -1:
                CurlZ = '(' + str(CurlZ) + ')* np.ones(np.shape(self.zg))'
            
            # evaluate the curl expression, creating list of curl values
            Curl_X = eval(CurlX)
            Curl_Y = eval(CurlY)
            Curl_Z = eval(CurlZ)

            curl_vf = vector_field3(self.xg, self.yg, self.zg, Curl_X, Curl_Y, Curl_Z, self.str_x, self.str_y, self.str_z)
            return(curl_vf)

    def div(self, at_x=None, at_y=None, at_z=None):
        if self.str_x == None or self.str_y == None or self.str_z == None:
            # ERROR
            raise TypeError('No equation provided')
        else:
            sep = (abs(np.min(self.xg)) + abs(np.max(self.xg)))/100

            New_grid = np.linspace(np.min(self.xg), np.max(self.xg), 100)

            # convert separation to string and find the number of decimals in the separation
            
            #new, denser grid using new coordinates
            xng, yng, zng = np.meshgrid(New_grid, New_grid, New_grid)


            # account for user input mistakes,
            # making sure that input coordinate is within specified range
           
            #Field components into string (simplified expressions)
            Fx = str(simplify(self.str_x))
            Fy = str(simplify(self.str_y))
            Fz = str(simplify(self.str_z))
            
            # convert string field components into sympy expression
            F_x = parse_expr(Fx)
            F_y = parse_expr(Fy)
            F_z = parse_expr(Fz)        

            # take derivatives of the components to define a dot product 
            # between del and F
            ddx_Fx = F_x.diff(sympy.symbols('x'))

            ddy_Fy = F_y.diff(sympy.symbols('y'))
            
            ddz_Fz = F_z.diff(sympy.symbols('z'))

            # string mathematical expression for divergence of the input field
            ddx_Fx = str(simplify(ddx_Fx))
            ddy_Fy = str(simplify(ddy_Fy))
            ddz_Fz = str(simplify(ddz_Fz))


            DivExpr = ('('+ddx_Fx+') + ('
                        +ddy_Fy+') + ('
                        +ddz_Fz+')')

            # print out the expression              
            print(DivExpr)

            if at_x!=None and at_y!=None and at_z!=None:
                if at_x > np.max(self.xg) or at_x < np.min(self.xg):
                    print("Input coordinate has to be within range specified by your grid")
                    exit()
                elif at_y > np.max(self.xg) or at_y < np.min(self.xg):
                    print("Input coordinate has to be within range specified by your grid")
                    exit()
                elif at_y > np.max(self.xg) or at_y < np.min(self.xg):
                    print("Input coordinate has to be within range specified by your grid")
                    exit()

                sep_str = str(sep)
                sep_decimals = sep_str[::-1].find('.')

                # round the input coordinates so that it is easier to relate them to grid points
                at_x = np.round(at_x, sep_decimals)
                at_y = np.round(at_y, sep_decimals)
                at_z = np.round(at_z, sep_decimals)

                # replace x, y, z with corresponding grid components,
                # so that the expression can be evaluated
                DivExpr = DivExpr.replace('x','xng')
                DivExpr = DivExpr.replace('y','yng')
                DivExpr = DivExpr.replace('z','zng')

                
                # if no coordinate specified, create a scalar matrix to define field
                if DivExpr.find('x') & DivExpr.find('y') & DivExpr.find('z') == -1:
                    DivExpr = '(' + str(DivExpr) + ')* np.ones(np.shape(xng))'

                # evaluate the divergence
                Div = eval(DivExpr)

                # define a list of coordinates for the input point
                coord = [at_x,at_y,at_z]

                # find index for point in grid which is closest to the input crd
                coord_idx = np.argwhere((abs((xng)-(coord[0]))<=(sep/2)) &
                                        (abs((yng)-(coord[1]))<=(sep/2)) &
                                        (abs((zng)-(coord[2]))<=(sep/2)))[0]

                # allocate indices to variables
                a = coord_idx[0]
                b = coord_idx[1]
                c = coord_idx[2]

                # define grid point whish is closest to the input
                x = str(np.round(xng[a][b][c],2))
                y = str(np.round(yng[a][b][c],2))
                z = str(np.round(zng[a][b][c],2))

                # print out the point and divergence at that point
                print('Divergence at the grid point closest to the input, ['+x+','+y+','+z+'], = '+str(np.round(Div[a][b][c], 1)))


    def plot(self, add_curl = None):



        # for arrows to work, with nan and infs
        # make a local variable of F_x and F_y
        # so that thye don't alter globally
        F_x_local = self.F_x * 1
        F_y_local = self.F_y * 1
        F_z_local = self.F_z * 1

        isnan_arrx = np.isnan(F_x_local)
        isnan_arry = np.isnan(F_y_local)
        isnan_arrz = np.isnan(F_z_local)

                # set all insignificant values to zero:
        F_x_local[np.abs(F_x_local) < 1e-15] = 0
        F_y_local[np.abs(F_y_local) < 1e-15] = 0
        F_z_local[np.abs(F_z_local) < 1e-15] = 0
        
        # find the magnitude corresponding to each point and store in mag array
        mag = np.sqrt(F_x_local**2 + F_y_local**2 + F_z_local**2)
        
        # find the maximum magnitude for scaling
        max_size = np.max(mag)   # careful with singularities, else ---> nan
        
        # Rescale components if log scaling is selected

        xmin = int(np.min(self.xg))
        ymin = int(np.min(self.yg))
        zmin = int(np.min(self.zg))
        xmax = int(np.max(self.xg))
        ymax = int(np.max(self.yg))
        zmax = int(np.max(self.zg))

        mlab.figure(bgcolor = (1,1,1), fgcolor = (0,0,0))
        mlab.quiver3d(self.xg, self.yg, self.zg, F_x_local, F_y_local, F_z_local, colormap='jet', line_width=3.0)

        mlab.outline(extent = [xmin,xmax,ymin,ymax,zmin,zmax], line_width=1.0)
        axes = mlab.axes(color = (0,0,0), nb_labels = 5, extent = [xmin,xmax,ymin,ymax,zmin,zmax], line_width=3.0)


        if add_curl is not None:
        
            if self.str_x == None or self.str_y == None or self.str_z == None:
            # ERROR
                raise TypeError('No equation provided')
            else:

                Fx = parse_expr(self.str_x)
                Fy = parse_expr(self.str_y)
                Fz = parse_expr(self.str_z)        

                # differentiate expressions w.r.t x, y, z.

                ddy_Fx = Fx.diff(sympy.symbols('y'))
                ddz_Fx = Fx.diff(sympy.symbols('z'))

                ddx_Fy = Fy.diff(sympy.symbols('x'))

                ddz_Fy = Fy.diff(sympy.symbols('z'))

                ddx_Fz = Fz.diff(sympy.symbols('x'))
                ddy_Fz = Fz.diff(sympy.symbols('y'))

                # Define a string expression for Curl of the input vector field
                CurlX = (str(simplify(ddy_Fz)) + '-' + str(simplify(ddz_Fy)))
                CurlY = (str(simplify(ddz_Fx)) + '-' + str(simplify(ddx_Fz)))
                CurlZ = (str(simplify(ddx_Fy)) + '-' + str(simplify(ddy_Fx)))
            
                # replace input coordinates with coordinate point lists
                CurlX = CurlX.replace('x', '(self.xg)')
                CurlX = CurlX.replace('y', '(self.yg)')
                CurlX = CurlX.replace('z', '(self.zg)')

                CurlY = CurlY.replace('x', '(self.xg)')
                CurlY = CurlY.replace('y', '(self.yg)')
                CurlY = CurlY.replace('z', '(self.zg)')

                CurlZ = CurlZ.replace('x', '(self.xg)')
                CurlZ = CurlZ.replace('y', '(self.yg)')
                CurlZ = CurlZ.replace('z', '(self.zg)')

                # if input scalar, define a scalar list the size of a coordinate obj.
                if CurlX.find('x') & CurlX.find('y') & CurlX.find('z') == -1:
                    CurlX = '(' + str(CurlX) + ')* np.ones(np.shape(self.xg))'
                if CurlY.find('x') & CurlY.find('y') & CurlY.find('z') == -1:
                    CurlY = '(' + str(CurlY) + ')* np.ones(np.shape(self.yg))'
                if CurlZ.find('x') & CurlZ.find('y') & CurlZ.find('z') == -1:
                    CurlZ = '(' + str(CurlZ) + ')* np.ones(np.shape(self.zg))'
                
                # evaluate the curl expression, creating list of curl values
                Curl_X = eval(CurlX)
                Curl_Y = eval(CurlY)
                Curl_Z = eval(CurlZ)


                mlab.quiver3d(self.xg, self.yg, self.zg, Curl_X, Curl_Y, Curl_Z, color=(1.0, 0.0, 1.0), opacity=0.4)



        
        for i in range(len(self.xg[:,0,0])):
            for j in range(len(self.yg[0,:,0])):
                for k in range(len(self.zg[0,0,:])):

                    if isnan_arrx[i][j][k] or isnan_arry[i][j][k] or isnan_arrz[i][j][k]:
                        F_x_local[i,j,k] = F_y_local[i,j,k] = F_z_local[i,j,k] = 0
                        sing = [self.xg[i,j,k],self.yg[i,j,k], self.zg[i,j,k]]
                            
                        mlab.points3d(sing[0], sing[1], sing[2], color = (1,0,0))
                            

                    if abs(F_x_local[i,j,k]) == np.inf or abs(F_y_local[i,j,k]) == np.inf or abs(F_z_local[i,j,k]) == np.inf or abs(F_y_local[i,j,k]) > 1e15 or abs(F_x_local[i,j,k]) > 1e15 or abs(F_z_local[i,j,k]) > 1e15:
                        F_x_local[i,j,k] = F_y_local[i,j,k] = F_z_local[i,j,k] = 0
                        sing = [self.xg[i,j,k],self.yg[i,j,k], self.zg[i,j,k]]

                        mlab.points3d(sing[0], sing[1], sing[2], color = (0,0,1))
        mlab.show()



class form_0_3d():

    """
    
    """

    def __init__(self, xg, yg, zg, form_0, form_0_eqn=None):
        self.xg = xg
        self.yg = yg
        self.zg = zg
        self.form_0 = form_0
        self.pt_den_x = len(xg[0, :, :])
        self.pt_den_y = len(xg[:, 0, :])
        self.pt_den_z = len(xg[:, :, 0])
        self.delta_factor = 10
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
        if string.find('x') & string.find('y') & string.find('z') == -1:
            string = '(' + str(string) + ')* np.ones(np.shape(xg))'
        else:
            string = string.replace('x', '(self.xg)')
            string = string.replace('y', '(self.yg)')
            string = string.replace('z', '(self.zg)')
        

        self.form_0 = eval(string)

    def return_string(self):
        '''
        Takes in no arguments, returns the unformatted string back to user
        This is done in case user wants to access strings
        that got here not by input but by ext. alg.
        '''
        return self.form_0_str

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
            x = np.linspace(self.xg[0,0,0], self.xg[-1,-1,0], points_number)
            y = np.linspace(self.yg[0,0,0], self.yg[-1,0,-1], points_number)
            z = np.linspace(self.zg[0,0,0], self.zg[0,-1,-1], points_number)
            self.xg, self.yg, self.zg= np.meshgrid(x,y,z)
            # based on these change other, dependant variables
            self.pt_den_x = len(self.xg[0, :, :])
            self.pt_den_y = len(self.yg[:, 0, :])
            self.pt_den_z = len(self.zg[:, :, 0])
            # substitute these into the equation:
            # but keep it local
            str_0 = self.form_0_str + ''
            str_0 = str_0.replace('x', '(self.xg)')
            str_0 = str_0.replace('y', '(self.yg)')
            str_0 = str_0.replace('z', '(self.zg)')
            # correct for constant forms
            if str_0.find('x') & str_0.find('y') & str_0.find('z') == -1:
                str_0 = '(' + str(str_0) + ')* np.ones(np.shape(self.xg))'
            # re-evaluate the 2-form numerically
            self.form_0 = eval(str_0)

    def plot(self, cross_sec_plane=None):


            
            form_0 = self.form_0
            xg = self.xg
            yg = self.yg
            zg = self.zg
            
            
            # set all insignificant values to zero:
            form_0[np.abs(form_0) < 1e-15] = 0

            xmin = int(np.min(self.xg))
            ymin = int(np.min(self.yg))
            zmin = int(np.min(self.zg))
            xmax = int(np.max(self.xg))
            ymax = int(np.max(self.yg))
            zmax = int(np.max(self.zg))

            fig = mlab.figure(bgcolor = (1,1,1), fgcolor = (0,0,0))

            # deal with sinularities that appear on evaluated points
            isnan_arr = np.isnan(form_0)
            for i in range(len(xg[0, :, :])):
                for j in range(len(yg[:, 0, :])):
                    for k in range(len(zg[:, :, 0])):
                        # set to zero points that are not defined or inf
                        if isnan_arr[k, j, i] or abs(form_0[k, j, i]) == np.inf or abs(form_0[k, j, i]) > 1e15:
                            # colour this region as a red dot, not square to
                            # not confuse with high mag 2-forms in stacks. or worse, in
                            # blocks
                            mlab.points3d((xg[k, j, i], yg[k, j, i]), zg[k, j, i], color=(1,0,0))
                            
                            form_0[k, j, i] = 0
            
            if self.logarithmic_scale_bool:
                mag1 = np.abs(form_0) + 1
                form_0_norm = form_0/(mag1)
                logmag = np.log10(mag1)
                form_0 = form_0_norm*logmag

            else:
                pass
            
            
            opac = 0.5
            colourmap='jet'

            cnt = mlab.contour3d(form_0, colormap=colourmap, opacity = opac, contours=self.lines, figure = fig)
            mlab.colorbar(object = cnt, orientation='vertical')
            mlab.outline(line_width=1.0)
            mlab.axes(color = (0,0,0), ranges = (xmin, xmax, ymin, ymax, zmin, zmax), nb_labels = 5, line_width=3.0)

            if str(cross_sec_plane)=='y':
                opac=0.05
                colour=(0.5,0.5,0.5)
                mlab.clf(fig)
                cnt = mlab.contour3d(form_0, color=colour, opacity = opac, contours=self.lines, figure = fig)
                mlab.colorbar(object = cnt, orientation='vertical')
                mlab.outline(line_width=1.0)
                mlab.axes(color = (0,0,0), ranges = (xmin, xmax, ymin, ymax, zmin, zmax), nb_labels = 5, line_width=3.0)
                mlab.pipeline.scalar_cut_plane(cnt)
            elif str(cross_sec_plane)=='n' or cross_sec_plane==None:
                pass
            else:
                print("please specify cross_sec_plane to be either 'y' or 'n'")
                exit()


            
           
            
            mlab.show()



class form_1_3d():

    def __init__(self, xg, yg, zg, F_x, F_y, F_z, F_x_eqn=None, F_y_eqn=None, F_z_eqn=None):
        self.xg = xg
        self.yg = yg
        self.zg = zg
        self.F_x = F_x
        self.F_y = F_y
        self.F_z = F_z
        self.s_max = 6
        self.s_min = 1
        self.fract = 0.05
        self.scale = 1
        self.w_head = 1/8
        self.h_head = 1/4
        self.arrowheads = True
        self.color = (1,0,1)
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

        if F_z_eqn is not None:
            self.form_1_str_z = str(simplify(F_z_eqn))
        else:
            self.form_1_str_z = None

    def give_eqn(self, equation_str_x, equation_str_y, equation_str_z):
        '''

        '''
        # set equation parameters to simplified inputs
        self.form_1_str_x = str(simplify(equation_str_x))
        self.form_1_str_y = str(simplify(equation_str_y))
        self.form_1_str_z = str(simplify(equation_str_z))
        # make the values match automatically to limit how often mismatch occurs
        # substitute these into the equation, but keep it local: 
        str_x = self.form_1_str_x + ''
        str_y = self.form_1_str_y + ''
        str_z = self.form_1_str_z + ''

        str_x = str_x.replace('x', '(self.xg)')
        str_x = str_x.replace('y', '(self.yg)')
        str_x = str_x.replace('z', '(self.zg)')

        str_y = str_y.replace('x', '(self.xg)')
        str_y = str_y.replace('y', '(self.yg)')
        str_y = str_y.replace('z', '(self.zg)')

        str_z = str_z.replace('x', '(self.xg)')
        str_z = str_z.replace('y', '(self.yg)')
        str_z = str_z.replace('z', '(self.zg)')

        # check against constant forms, to have correct shape
        if str_x.find('x') & str_x.find('y') & str_x.find('z') == -1:
            str_x = '(' + str(str_x) + ')* np.ones(np.shape(self.xg))'
        if str_y.find('x') & str_y.find('y') & str_y.find('z') == -1:
            str_y = '(' + str(str_y) + ')* np.ones(np.shape(self.yg))'
        if str_z.find('x') & str_z.find('y') & str_z.find('z') == -1:
            str_z = '(' + str(str_z) + ')* np.ones(np.shape(self.zg))'
        

        
        # evaluate formatted equations and save
        self.F_x = eval(str_x)
        self.F_y = eval(str_y)
        self.F_z = eval(str_z)

    def return_string(self):

        '''
        Returns unformatted strings for component equations back to user
        Done in case user wants to access strings that got here by ext. alg.
        
        Parmateres: None
        Returns: None
        
        '''
        return self.form_1_str_x, self.form_1_str_y, self.form_1_str_z

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
    





























