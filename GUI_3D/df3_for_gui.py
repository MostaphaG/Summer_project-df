import numpy as np
from mayavi import mlab
import matplotlib.pyplot as plt
import sympy
from sympy import diff, simplify
from sympy.parsing.sympy_parser import parse_expr

from tvtk.api import tvtk
from tvtk.common import configure_input_data, configure_source_data

from numpy import sin, cos, tan, sqrt, log, arctan, arcsin, arccos, tanh
from numpy import sinh, cosh, arcsinh, arccosh, arctanh, exp, pi, e






class vector_field3():

    '''
    The class which creates a vector field object. 
    Requires grid (3D), field components.
    in order to calculate  curl, div, zoomed field, covariant vector object or derivative field,
    it is necessary to
    provide string expressions for the field component to
    the object as well
    .plot() can be used to plot out the vector field and the curl (either on the same axes or different)


    Methods: (applied to the object)

    .give_eqn(equation_str_x, equation_str_y, equation_str_z) - provides equations ; 
    .curl() - creates a separate curl object ; 
    .div(at_x=None, at_y=None, at_z=None) [if coordinates provided, gives divergence at the point] ;
    .log_scaling() - scales the field logarithmycally; 
    .zoom(magnitude, target, point_density) - creates a zoomed in field object ; 
    .autoscale() - automatic scale for better visuals ;
    .covariant() - creates a covariant vector field object (1-form) ;
    .deriv(target, magnitude, point_density) - creates a derivative field object at the target location ;
    .plot() - applied to either curl of object, plots an object of teh .vector_field3() class
    [if applied to object can use add_curl ='y' argument to plot
    curl and field on the same axes] ;
    '''


    def __init__(self, xg, yg, zg, F_x, F_y, F_z, scaling, sing_scl, F_x_eqn=None, F_y_eqn=None, F_z_eqn=None):
        # Define initial variables supplied to the object

        self.xg = xg # coordinate grid
        self.yg = yg
        self.zg = zg
        self.F_x = F_x # field magnitude at the grid points
        self.F_y = F_y
        self.F_z = F_z
        self.pt_den = len(xg[0, :, 0]) # amount of points in a grid
        self.orientation = 'mid'
        self.scale = scaling
        self.sing_scl = sing_scl
        self.logarithmic_scale_bool = 0

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
        Takes in 3 arguments (x, y & z components), string
        It must be in terms of x, y and z.
        Has to be given, for curl, div and some other methods to be calculatable.
        '''
        # set equation parameters to simplified inputs
        self.str_x = str(simplify(equation_str_x))
        self.str_y = str(simplify(equation_str_y))
        self.str_z = str(simplify(equation_str_z))

        # keep the equations local
        str_x = self.str_x + ''
        str_y = self.str_y + ''
        str_z = self.str_z + ''

        # replace coordinates in the strings with coordinate grid arrays
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
        
        
        # re-evaluate the field components numerically
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
        '''
        When applied to a vector field object, creates another vector field
        object which represents curl of the initial vector field.

        --------------------------------------------------------
        Returns:

        separate vector field object

        --------------------------------------------------------
        Arguments:

        None
        '''
        # Check whether the equations are provided
        if self.str_x == None or self.str_y == None or self.str_z == None:
            # ERROR
            raise TypeError('No equation provided')
        else:

            # Create SymPy objects from the input strings
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

            # create a vf object
            curl_vf = vector_field3(self.xg, self.yg, self.zg, Curl_X, Curl_Y, Curl_Z, self.str_x, self.str_y, self.str_z)
            return(curl_vf)



    def div(self, at_x, at_y, at_z):


        # Check whether the string equations are supplied
        if self.str_x == None or self.str_y == None or self.str_z == None:
            # ERROR
            raise TypeError('No equation provided')
        else:


            if at_x==0:
                at_x = 1e-15
            if at_y==0:
                at_y = 1e-15
            if at_z==0:
                at_z = 1e-15
           
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

            DivExpr = DivExpr.replace('x', str(at_x))
            DivExpr = DivExpr.replace('y', str(at_y))
            DivExpr = DivExpr.replace('z', str(at_z))


            
            if DivExpr.find('x') & DivExpr.find('y') & DivExpr.find('z') == -1:
                DivExpr = '(' + str(DivExpr) + ')* 1'

                # evaluate the divergence
            Div = eval(DivExpr)

            if Div >1e+15:
                Div = "+inf"
            elif Div <-1e+15:
                Div = "-inf"
                
            return Div



    def log_scaling(self):
        '''
        Takes no arguments
        Changes the boolean that determines if scaling is logarithmic
        Whenever it is called, it changes that boolean to opposite
        The form object is initialised with this as False
        '''
        self.logarithmic_scale_bool = not self.logarithmic_scale_bool



    def zoom(self, mag, target, dpd):
        
        """
        Redefines a new, zoomed in grid on which a new VF is calculated.
        Takes in three arguments:

        mag - magnitude of zoom ; 
        target - target point [x, y, z] ;
        dpd - points amount in the new grid ; 

        Returns a new, zoomed in vector field 
        """

        if self.str_x == None or self.str_y == None or self.str_z == None:
            # ERROR
            raise TypeError('No equation provided, see \'give_eqn\' method')
        else:
             # Zoom must be one or greater
            if mag < 1:
                raise ValueError('Mag must be greater than one')
            else:

                # Target coordinates
                x_m = target[0]
                y_m = target[1]
                z_m = target[2]
                
                # Get the size of the original VF
                Lx = 0.5*(self.xg[-1, -1, -1] - self.xg[0, 0, 0])
                Ly = 0.5*(self.yg[-1, -1, -1] - self.yg[0, 0, 0])
                Lz = 0.5*(self.zg[-1, -1, -1] - self.zg[0, 0, 0])

                
                # Zoom axis range
                d_range_x = Lx/mag
                d_range_y = Ly/mag
                d_range_z = Lz/mag
                
                # Set up zoom window grids
                dx = np.linspace(-d_range_x + x_m, d_range_x + x_m, dpd)
                dy = np.linspace(-d_range_y + y_m, d_range_y + y_m, dpd)
                dz = np.linspace(-d_range_z + z_m, d_range_z + z_m, dpd)
                dxg, dyg, dzg = np.meshgrid(dx, dy, dz)
                
                # Create variables for the user provided equation strings
                u_str = self.str_x
                v_str = self.str_y
                k_str = self.str_z
                
                # Check if the equations provided contain x and y terms
                if u_str.find('x') & u_str.find('y') & u_str.find('z')== -1:
                    u_str = '(' + str(u_str) + ')* np.ones(np.shape(dxg))'
                else:
                    u_str = u_str.replace('x', 'dxg')
                    u_str = u_str.replace('y', 'dyg')
                    u_str = u_str.replace('z', 'dzg')
            
                if v_str.find('x') & v_str.find('y') & v_str.find('z') == -1:
                    v_str = '(' + str(v_str) + ')* np.ones(np.shape(dyg))'
                else:
                    v_str = v_str.replace('x', 'dxg')
                    v_str = v_str.replace('y', 'dyg')
                    v_str = v_str.replace('z', 'dzg')

                if k_str.find('x') & k_str.find('y') & k_str.find('z') == -1:
                    k_str = '(' + str(k_str) + ')* np.ones(np.shape(dzg))'
                else:
                    k_str = k_str.replace('x', 'dxg')
                    k_str = k_str.replace('y', 'dyg')
                    k_str = k_str.replace('z', 'dzg')
                    
                # Generate arrays for the components of the zoom field
                u_zoom = eval(u_str)
                v_zoom = eval(v_str)
                k_zoom = eval(k_str)
                
                # crate the zoomed in form
                zoom_form = vector_field3(dxg, dyg, dzg, u_zoom, v_zoom, k_zoom, self.str_x, self.str_y, self.str_z)

                
                return zoom_form



    def autoscale(self):
        '''
        Takes no arguments
        Changes the boolean that determines if arrows are autoscaled
        Whenever it is called, it changes that boolean to opposite
        The form object is initialised with this as False
        '''
        self.scale_bool = not self.scale_bool



    def covariant(self, g=[['1', '0', '0'], ['0', '1', '0'], ['0', '0', '1']]):
            '''
            Passes in everything it can (all it has been supplied)
            to the 1-form object.
            Works via the metric on R3
            Can supply the metric in as equations or as evaluated arrays
            Format of the metric is a list of numpy arrays
            0th array is the top row, its 0th component is 00, 1st is 01, 2nd is 02
            1st array is the middle row, its 0th comp is 10, 1st is 11, 2nd is 12.
            2nd array is the bottom row, its 0th comp is 20, 1st is 21, 2nd is 22.
            Note, if it is supplied as arrays, they must come from numpy grids
            via meshgrid, if it is supplied as strings, needs to be in terms of
            x, y, z, and contain no special funtions
            
            Returns a single object (1-form object)
            '''
            
            # extract what is needed form the metric depending on what the user
            # supplied
            # check if its has string components
            if type(g[0][0]) == str and type(g[0][1]) == str and type(g[0][2]) == str \
               and type(g[1][0]) == str and type(g[1][1]) == str and type(g[1][2]) == str \
               and type(g[2][0]) == str and type(g[2][1]) == str and type(g[2][2]) == str:
                # deal with supplied string metric
                # need to format it, correct it for constants and evaluate it's numerical equivalent
                str_comp_00 = g[0][0] + ''
                str_comp_01 = g[0][1] + ''
                str_comp_02 = g[0][2] + ''
                str_comp_10 = g[1][0] + ''
                str_comp_11 = g[1][1] + ''
                str_comp_12 = g[1][2] + ''
                str_comp_20 = g[2][0] + ''
                str_comp_21 = g[2][1] + ''
                str_comp_22 = g[2][2] + ''

                str_comp_00 = str_comp_00.replace('x', '(self.xg)')
                str_comp_00 = str_comp_00.replace('y', '(self.yg)')
                str_comp_00 = str_comp_00.replace('z', '(self.zg)')
                str_comp_01 = str_comp_01.replace('x', '(self.xg)')
                str_comp_01 = str_comp_01.replace('y', '(self.yg)')
                str_comp_01 = str_comp_01.replace('z', '(self.zg)')
                str_comp_02 = str_comp_02.replace('x', '(self.xg)')
                str_comp_02 = str_comp_02.replace('y', '(self.yg)')
                str_comp_02 = str_comp_02.replace('z', '(self.zg)')

                str_comp_10 = str_comp_10.replace('x', '(self.xg)')
                str_comp_10 = str_comp_10.replace('y', '(self.yg)')
                str_comp_10 = str_comp_10.replace('z', '(self.zg)')
                str_comp_11 = str_comp_11.replace('x', '(self.xg)')
                str_comp_11 = str_comp_11.replace('y', '(self.yg)')
                str_comp_11 = str_comp_11.replace('z', '(self.zg)')
                str_comp_12 = str_comp_12.replace('x', '(self.xg)')
                str_comp_12 = str_comp_12.replace('y', '(self.yg)')
                str_comp_12 = str_comp_12.replace('z', '(self.zg)')

                str_comp_20 = str_comp_20.replace('x', '(self.xg)')
                str_comp_20 = str_comp_20.replace('y', '(self.yg)')
                str_comp_20 = str_comp_20.replace('z', '(self.zg)')
                str_comp_21 = str_comp_21.replace('x', '(self.xg)')
                str_comp_21 = str_comp_21.replace('y', '(self.yg)')
                str_comp_21 = str_comp_21.replace('z', '(self.zg)')
                str_comp_22 = str_comp_22.replace('x', '(self.xg)')
                str_comp_22 = str_comp_22.replace('y', '(self.yg)')
                str_comp_22 = str_comp_22.replace('z', '(self.zg)')

                # check against constant form components:
                if str_comp_00.find('x') & str_comp_00.find('y') & str_comp_00.find('z') == -1:
                    str_comp_00 = '(' + str(str_comp_00) + ')* np.ones(np.shape(self.xg))'
                if str_comp_01.find('x') & str_comp_01.find('y') & str_comp_01.find('z') == -1:
                    str_comp_01 = '(' + str(str_comp_01) + ')* np.ones(np.shape(self.xg))'
                if str_comp_02.find('x') & str_comp_02.find('y') & str_comp_02.find('z') == -1:
                    str_comp_02 = '(' + str(str_comp_02) + ')* np.ones(np.shape(self.xg))'
                
                if str_comp_10.find('x') & str_comp_10.find('y') & str_comp_10.find('z') == -1:
                    str_comp_10 = '(' + str(str_comp_10) + ')* np.ones(np.shape(self.xg))'
                if str_comp_11.find('x') & str_comp_11.find('y') & str_comp_11.find('z') == -1:
                    str_comp_11 = '(' + str(str_comp_11) + ')* np.ones(np.shape(self.xg))'
                if str_comp_12.find('x') & str_comp_12.find('y') & str_comp_12.find('z') == -1:
                    str_comp_12 = '(' + str(str_comp_12) + ')* np.ones(np.shape(self.xg))'

                if str_comp_20.find('x') & str_comp_20.find('y') & str_comp_20.find('z') == -1:
                    str_comp_20 = '(' + str(str_comp_20) + ')* np.ones(np.shape(self.xg))'
                if str_comp_21.find('x') & str_comp_21.find('y') & str_comp_21.find('z') == -1:
                    str_comp_21 = '(' + str(str_comp_21) + ')* np.ones(np.shape(self.xg))'
                if str_comp_22.find('x') & str_comp_22.find('y') & str_comp_22.find('z') == -1:
                    str_comp_22 = '(' + str(str_comp_22) + ')* np.ones(np.shape(self.xg))'


                # evaluate the components numerically
                # store numerical metric
                comp_00 = eval(str_comp_00)
                comp_01 = eval(str_comp_01)
                comp_02 = eval(str_comp_02)
                comp_10 = eval(str_comp_10)
                comp_11 = eval(str_comp_11)
                comp_12 = eval(str_comp_12)
                comp_20 = eval(str_comp_20)
                comp_21 = eval(str_comp_21)
                comp_22 = eval(str_comp_22)

                g_num = [[comp_00, comp_01, comp_02], [comp_10, comp_11, comp_12], [comp_20, comp_21, comp_22]]
                
                # set up a dummy variable to store the fact that numericals were given
                # not to check again later
                analytics = True
                
            elif type(g[0][0]) == np.ndarray and type(g[0][1]) == np.ndarray and type(g[0][2]) == np.ndarray\
                 and type(g[1][0]) == np.ndarray and type(g[1][1]) == np.ndarray and type(g[1][2]) == np.ndarray\
                 and type(g[2][0]) == np.ndarray and type(g[2][1]) == np.ndarray and type(g[2][2]) == np.ndarray:
                # deal with the metric being supplied as components
                # if the user has vector field equations, warn that these can't
                # be passed anymore, because we don't have equations for this
                # metric
                if self.str_x == None and self.str_y == None and self.str_z == None:
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
            form_x = self.F_x * g_num[0][0] + self.F_y * g_num[0][1] + self.F_z * g_num[0][2]
            form_y = self.F_x * g_num[1][0] + self.F_y * g_num[1][1] + self.F_z * g_num[1][2]
            form_z = self.F_x * g_num[2][0] + self.F_y * g_num[2][1] + self.F_z * g_num[2][2]
            
            # if the equations were given, evaluate these analytically too:
            # only if vector file doriginally has equations
            if analytics:
                if self.str_x == None and self.str_y == None and self.str_z == None:
                    print('You supplied the metric as equations (or it was default), but did not give VF equations, therefore only numericals will be completed')
                    analytics = False
                else:
                    x_str_form = '(' + self.str_x + ')*(' + g[0][0] + ') + (' + self.str_y + ')*(' + g[0][1] + ') + (' + self.str_z + ')*(' + g[0][2] + ')'
                    y_str_form = '(' + self.str_x + ')*(' + g[1][0] + ') + (' + self.str_y + ')*(' + g[1][1] + ') + (' + self.str_z + ')*(' + g[1][2] + ')'
                    z_str_form = '(' + self.str_x + ')*(' + g[2][0] + ') + (' + self.str_y + ')*(' + g[2][1] + ') + (' + self.str_z + ')*(' + g[2][2] + ')'
                    # simplify them
                    x_str_form = str(simplify(x_str_form))
                    y_str_form = str(simplify(y_str_form))
                    z_str_form = str(simplify(z_str_form))
            else:
                pass

            # based on what was given into the Vector field
            # return a 1-form object with these parameters
            if analytics:
                result_form = form_1_3d(self.xg, self.yg, self.zg, form_x, form_y, form_z, x_str_form, y_str_form, z_str_form)
            elif not analytics:
                result_form = form_1_3d(self.xg, self.yg, self.zg, form_x, form_y, form_z)
        
            # return the found object
            return result_form



    def deriv(self, target, mag, dpd):

        '''
        Creates new vector field object at a target location, showing the derivative field at this point.
        User gives arguments:
        Target - derivative plot location
        mag - Magnification level
        dpd - New plot point density
        
        
        
        Returns:
        --------
    
        deriv VF object

        
        '''
        
        if self.str_x == None or self.str_y == None or self.str_z == None:
            # ERROR
            raise TypeError('No equation provided')
        else:
             # Zoom must be one or greater
            if mag < 1:
                raise ValueError('Zoom must be greater than one')
            else:
                
                
        
                # Target coordinates
                x_m = target[0]
                y_m = target[1]
                z_m = target[2]
                
                # Get the size of the original VF
                Lx = 0.5*(self.xg[-1, -1, -1] - self.xg[0, 0, 0])
                Ly = 0.5*(self.yg[-1, -1, -1] - self.yg[0, 0, 0])
                Lz = 0.5*(self.zg[-1, -1, -1] - self.zg[0, 0, 0])
                
                # Zoom axis range
                d_range_x = Lx/mag
                d_range_y = Ly/mag
                d_range_z = Lz/mag
                
                # Set up zoom window grids
                dx = np.linspace(-d_range_x + x_m, d_range_x + x_m, dpd)
                dy = np.linspace(-d_range_y + y_m, d_range_y + y_m, dpd)
                dz = np.linspace(-d_range_z + z_m, d_range_z + z_m, dpd)
                dxg, dyg, dzg = np.meshgrid(dx, dy, dz)

                dxg[dxg==0]=1e-15
                dyg[dyg==0]=1e-15
                dzg[dzg==0]=1e-15

                if x_m==0:
                    x_m=1e-15
                if y_m==0:
                    y_m=1e-15
                if z_m==0:
                    z_m=1e-15
                
                # Create variables for the user provided equation strings
                u_str = self.str_x
                v_str = self.str_y
                k_str = self.str_z

                # Create string to evaluate the field at the target location
                u_str_point = u_str.replace('x', 'x_m')
                u_str_point = u_str_point.replace('y', 'y_m')
                u_str_point = u_str_point.replace('z', 'z_m')
                
                v_str_point = v_str.replace('x', 'x_m')
                v_str_point = v_str_point.replace('y', 'y_m')
                v_str_point = v_str_point.replace('z', 'z_m')

                k_str_point = k_str.replace('x', 'x_m')
                k_str_point = k_str_point.replace('y', 'y_m')
                k_str_point = k_str_point.replace('z', 'z_m')
                
                # Check if the equations provided contain x, y, z terms
                


                if u_str.find('x') & u_str.find('y') & u_str.find('z')== -1:
                    u_str_grid = '(' + str(u_str) + ')* np.ones(np.shape(dxg))'
                else:
                    u_str_grid = u_str.replace('x', 'dxg')
                    u_str_grid = u_str_grid.replace('y', 'dyg')
                    u_str_grid = u_str_grid.replace('z', 'dzg')
            
                if v_str.find('x') & v_str.find('y') & v_str.find('z') == -1:
                    v_str_grid = '(' + str(v_str) + ')* np.ones(np.shape(dyg))'
                else:
                    v_str_grid = v_str.replace('x', 'dxg')
                    v_str_grid = v_str_grid.replace('y', 'dyg')
                    v_str_grid = v_str_grid.replace('z', 'dzg')

                if k_str.find('x') & k_str.find('y') & k_str.find('z') == -1:
                    k_str_grid = '(' + str(k_str) + ')* np.ones(np.shape(dzg))'
                else:
                    k_str_grid = k_str.replace('x', 'dxg')
                    k_str_grid = k_str_grid.replace('y', 'dyg')
                    k_str_grid = k_str_grid.replace('z', 'dzg') 
                # Generate arrays for the components of the derivative field     
                
                


                U = eval(u_str_grid) - eval(u_str_point)
                V = eval(v_str_grid) - eval(v_str_point)
                K = eval(k_str_grid) - eval(k_str_point)
                




                # from that create VF instance
                deriv_vf = vector_field3(dxg, dyg, dzg, U, V, K, self.str_x, self.str_y, self.str_z)
                
                
                
                # depending on preferances, return to user and plot
                
                return deriv_vf



    def plot(self, cmap, cmap_idx, clr, sing_clr, scaling_fctrs, opacity):

        


        # Check whether the scaling input is float and > 0.0
        scaling = self.scale
        redball_scaling = self.sing_scl

        # Check whether the opacity input is float and 0.0 > opac > 1.0



        # for arrows to work, with nan and infs
        # make a local variable of F_x and F_y
        # so that thye don't alter globally
        F_x_local = self.F_x
        F_y_local = self.F_y
        F_z_local = self.F_z

        

        # set all insignificant values to zero:
        F_x_local[np.abs(F_x_local) < 1e-15] = 0
        F_y_local[np.abs(F_y_local) < 1e-15] = 0
        F_z_local[np.abs(F_z_local) < 1e-15] = 0

        F_x_local[np.abs(F_x_local) > 1e+15] = np.inf
        F_y_local[np.abs(F_y_local) > 1e+15] = np.inf
        F_z_local[np.abs(F_z_local) > 1e+15] = np.inf
        
        # define grid dimension magnitude

        xmin = int((np.min(self.xg))-1)
        ymin = int((np.min(self.yg))-1)
        zmin = int((np.min(self.zg))-1)
        xmax = int((np.max(self.xg))+1)
        ymax = int((np.max(self.yg))+1)
        zmax = int((np.max(self.zg))+1)

        


        # all np.inf -> np.nan for convenience purposes
        F_x_local[np.isinf(F_x_local)] = np.nan
        F_y_local[np.isinf(F_y_local)] = np.nan
        F_z_local[np.isinf(F_z_local)] = np.nan

        
        # Field magnitude at each point
        F = np.sqrt(F_x_local**2+F_y_local**2+F_z_local**2)

        #if self.logarithmic_scale_bool==True:
        #F = np.log(F)

        # Indices where field is not defined (Field. mag. = np.nan)
        Idx = np.argwhere(np.isnan(F))

        # Set all NaN values to zero so that it does not disturb the plotting
        F[np.isnan(F)] = 0
        F_x_local[np.isnan(F_x_local)] = 0
        F_y_local[np.isnan(F_y_local)] = 0
        F_z_local[np.isnan(F_z_local)] = 0


        if self.logarithmic_scale_bool==True:
            mag1 = F + 1
            # min_size = np.min(mag1)
            
            unorm = F_x_local/mag1
            vnorm = F_y_local/mag1
            knorm = F_z_local/mag1
            
            # logsf = np.log10(mag1/min_size)
            logmag = np.log10(mag1)
            F_x_local = unorm*logmag
            F_y_local = vnorm*logmag
            F_z_local = knorm*logmag

        """
        if sng_size==None:
            sng_size=0.5
        else:
            pass

        stck_coords = [self.xg, F_x_local*scaling, self.yg, F_y_local*scaling, self.zg, F_z_local*scaling]
        red_balls_data = [self.xg[Idx[:,0],Idx[:,1],Idx[:,2]] , self.yg[Idx[:,0],Idx[:,1],Idx[:,2]], self.zg[Idx[:,0],Idx[:,1],Idx[:,2]], sng_size*redball_scaling]
        axes_limits = [xmin,xmax,ymin,ymax,zmin,zmax]

        return stck_coords, red_balls_data, axes_limits
        """

        # Plot red spheres at singular points
        mlab.points3d(self.xg[Idx[:,0],Idx[:,1],Idx[:,2]],
                      self.yg[Idx[:,0],Idx[:,1],Idx[:,2]],
                      self.zg[Idx[:,0],Idx[:,1],Idx[:,2]],
                      color = (sing_clr[0]/255, sing_clr[1]/255, sing_clr[2]/255),scale_factor=scaling_fctrs[1], resolution=100)
    
        if cmap_idx==0:
            qivfld=mlab.quiver3d(self.xg, self.yg, self.zg, F_x_local, F_y_local, F_z_local, scale_factor=scaling_fctrs[0], 
                                 colormap=cmap, line_width=1.0, mode='arrow', scale_mode = 'vector', opacity=opacity)
        if cmap_idx==1:
            qivfld=mlab.quiver3d(self.xg, self.yg, self.zg, F_x_local, F_y_local, F_z_local, scale_factor=scaling_fctrs[0],
                                 line_width=1.0, mode='arrow', scale_mode = 'vector', color=(clr[0]/255, clr[1]/255, clr[2]/255), opacity=opacity)
            
        cbar = mlab.vectorbar(object=qivfld)
        cbar.scalar_bar.unconstrained_font_size = True
        cbar.label_text_property.font_size=10
        cbar.scalar_bar_representation.orientation=1
        cbar.scalar_bar_representation.position = [0.05, 0.05]
        cbar.scalar_bar_representation.position2 = [0.05, 0.85]


        # add outline box and axes
        mlab.outline(extent = [xmin,xmax,ymin,ymax,zmin,zmax], line_width=1.0)
        axes = mlab.axes(color = (0,0,0), nb_labels = 5, extent = [xmin,xmax,ymin,ymax,zmin,zmax], line_width=3.0)
        axes.axes.font_factor = 0.8

        # show yhe plot on the figure
        mlab.show()







class form_0_3d():

    """
    The class which creates a 0 differential form object. 
    Requires grid (3D), scalar potential.
    in order to calculate set density of surfaces or change amount of levels, it is necessary to
    provide string expressions for the field component to
    the object as well
    .plot() can be used to plot out the 0 form. Can also provide the cross section plane to see how
    the potential changes in any direction
    

    Methods: (applied to the object)

    .give_eqn(equation_str_x, equation_str_y, equation_str_z) - provides equation for the 0-form; 
    .return_string() - returns string to the user ;
    .levels() - amunt of levels in surface plot ; 
    .log_scaling() - scales 0-form logarithmically ; 
    .set_density() - sets density of grid points ; 
    .zoom() - creates a zoomed in 0-form ; 
    .ext_d() - takes exterior derivative of 0-form, creates 1-form object ; 
    .num_ext_d() - takes numerical exterior derivative of 0-form, creates 1-form object ; 
    .hodge() - applies hodge star operator to the 0-form ; 
    .num_hodge() - applies hodge star operator to the 0-form numerically ; 
    .wedge(...) - wedges 0-form with another p-form ; 
    .num_wedge(...) - wedges 0-form with another p-form numerically ;  
    .plot()
    [if applied to object can use cross_sec_plane ='y' argument to plot
    the cut plane and see the colormesh of the potential]
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
                        It must be in terms of x, y, z.
                        Has to be given, for some methods to be calculatable.
        
        Returns: None
        
        '''
        self.form_0_str = equation_str
        
        # update the numerical values to always match
        string = self.form_0_str + ''
        
        # Check if the equations provided contain x, y, z terms
        # and format them to be evaluated
        if string.find('x') & string.find('y') & string.find('z') == -1:
            string = '(' + str(string) + ')* np.ones(np.shape(xg))'
        else:
            string = string.replace('x', '(self.xg)')
            string = string.replace('y', '(self.yg)')
            string = string.replace('z', '(self.zg)')
        
        # evaluate the string equation
        self.form_0 = eval(string)



    def return_string(self):
        '''
        Takes in no arguments, returns the unformatted string back to user
        This is done in case user wants to access strings
        that got here not by input but by ext. alg.
        '''
        return self.form_0_str



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
            # re-evaluate the 0-form numerically
            self.form_0 = eval(str_0)



    def zoom(self, mag, target, dpd):
        
        """
        Redefines a new, zoomed in grid on which a new 0-form is calculated.
        Takes in three arguments:

        mag - magnitude of zoom ; 
        target - target point [x, y, z] ;
        dpd - points amount in the new grid ; 

        Returns a new, zoomed in 0-form 
        """


        if self.form_0_str == None:
            # ERROR
            raise TypeError('No equation provided, see \'give_eqn\' method')
        else:
             # Zoom must be one or greater
            if mag < 1:
                raise ValueError('Mag must be greater than one')
            else:

                # Target coordinates
                x_m = target[0]
                y_m = target[1]
                z_m = target[2]
                
                # Get the size of the original VF
                Lx = 0.5*(self.xg[-1, -1, -1] - self.xg[0, 0, 0])
                Ly = 0.5*(self.yg[-1, -1, -1] - self.yg[0, 0, 0])
                Lz = 0.5*(self.zg[-1, -1, -1] - self.zg[0, 0, 0])

                
                # Zoom axis range
                d_range_x = Lx/mag
                d_range_y = Ly/mag
                d_range_z = Lz/mag
                
                # Set up zoom window grids
                dx = np.linspace(-d_range_x + x_m, d_range_x + x_m, dpd)
                dy = np.linspace(-d_range_y + y_m, d_range_y + y_m, dpd)
                dz = np.linspace(-d_range_z + z_m, d_range_z + z_m, dpd)
                dxg, dyg, dzg = np.meshgrid(dx, dy, dz)
                
                # Create variables for the user provided equation strings
                f0_str = self.form_0_str

                
                # Check if the equations provided contain x, y, z terms
                if f0_str.find('x') & f0_str.find('y') & f0_str.find('z')== -1:
                    f0_str = '(' + str(f0_str) + ')* np.ones(np.shape(dxg))'
                else:
                    f0_str = f0_str.replace('x', 'dxg')
                    f0_str = f0_str.replace('y', 'dyg')
                    f0_str = f0_str.replace('z', 'dzg')
            

                    
                # Generate arrays for the components of the zoom field
                f0_zoom = eval(f0_str)

                
                # crate the zoomed in form
                zoom_form = form_0_3d(dxg, dyg, dzg, f0_zoom, self.form_0_str)

                
                return zoom_form
 


    def plot(self, cross_sec_plane=None, colormap=None):

            '''
            Plots the 0 diferential form object.

            --------------------------------------------------
            Returns: 

            None

            --------------------------------------------------
            Arguments:
            *optional:
            cross_sec_plane = 'y' - adds scalar cross section slider to clearly view the potential from inside
            '''
            
            # create a local variable of the object within method
            form_0 = self.form_0


            # set all insignificant values to zero:
            form_0[np.abs(form_0) < 1e-15] = 0

            # deal with sinularities that appear on evaluated points
            form_0[np.isinf(form_0)] = np.nan
            form_0[np.isnan(form_0)] = 2*np.nanmax(form_0)


            # set boundaries for plotting axes and outline box
            xmin = float(np.min(self.xg))
            ymin = float(np.min(self.yg))
            zmin = float(np.min(self.zg))
            xmax = float(np.max(self.xg))
            ymax = float(np.max(self.yg))
            zmax = float(np.max(self.zg))



            # take account for logarithmic scaling
            if self.logarithmic_scale_bool:
                mag1 = np.abs(form_0) + 1
                form_0_norm = form_0/(mag1)
                logmag = np.log10(mag1)
                form_0 = form_0_norm*logmag
            else:
                pass
            
            # set initial values
            opac = 0.5
            colourmap='jet'

            if colormap != None:
                colourmap = colormap
            else:
                pass

            
            

                       
            # if the argument is 'y',  make the contours transparent and add the slider to view scalar cross section
            if str(cross_sec_plane)=='y':

                opac=0.05
                colour=(0.5,0.5,0.5)
                # clear figure from previous plot
                cnt = mlab.contour3d(form_0, color=colour, opacity = opac, contours=self.lines)
                mlab.colorbar(object = cnt, orientation='vertical')

                mlab.axes(color = (0,0,0), ranges = (xmin, xmax, ymin, ymax, zmin, zmax), nb_labels = 5, line_width=1.0)
                mlab.pipeline.scalar_cut_plane(cnt)
            elif str(cross_sec_plane)=='n' or cross_sec_plane==None:
                pass
            else:
                print("please specify cross_sec_plane to be either 'y' or 'n'")
                exit()
            

     
            axes_limits = [xmin,xmax,ymin,ymax,zmin,zmax]

            # plot the 0 form contours, add colorbar, axes and an outline
            cnt = mlab.contour3d(form_0, colormap=colourmap, opacity = opac, contours=self.lines)
            mlab.colorbar(object = cnt, orientation='vertical')
            mlab.outline(line_width=1.0)
            mlab.axes(color = (0,0,0), ranges = [xmin, xmax, ymin, ymax, zmin, zmax], nb_labels = 5, line_width=3.0)
            







class form_1_3d():

    """
    The class which creates a differential 1 form object. 
    Requires grid (3D), x y and z components.
    .plot() can be used to plot out the 1 form.
    

    Methods: (applied to the object)

    .give_eqn(equation_str_x, equation_str_y, equation_str_z) - provides equation for the 0-form; 
    .return_string() - returns strings to the user ;
    .log_scaling() - scales 1-form logarithmically ; 
    .set_density() - sets density of grid points ; 
    .zoom() - creates a zoomed in 1-form ; 
    .contravariant() - creates a contravariant vector field object
    .ext_d() - takes exterior derivative of 1-form, creates 2-form object ; 
    .num_ext_d() - takes numerical exterior derivative of 1-form, creates 2-form object ; 
    .hodge() - applies hodge star operator to the 1-form ; 
    .num_hodge() - applies hodge star operator to the 1-form numerically ; 
    .wedge(...) - wedges 1-form with another p-form ; 
    .num_wedge(...) - wedges 1-form with another p-form numerically ; 
    .interior_d() - takes interior derivative of 1-form, creates 0-form object ; 
    .num_interior_d() - takes interior derivative of 1-form numerically, creates 0-form object ; 
    .plot()
    [if applied to object can use cross_sec_plane ='y' argument to plot
    the cut plane and see the colormesh of the potential]
    """

    def __init__(self, xg, yg, zg, F_x, F_y, F_z, F_x_eqn=None, F_y_eqn=None, F_z_eqn=None, scaling=None, sng_scl=None):
        
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

        if scaling is not None:
            self.scaling = float(scaling)
        else:
            self.scaling = 1.0

        if sng_scl is not None:
            self.sng_scl = float(sng_scl)
        else:
            self.sng_scl = 1.0



    def give_eqn(self, equation_str_x, equation_str_y, equation_str_z):
        '''
        Allows user to supply equations to instance, if not initially done so
        
        Parameters:
        ------------
        equation_str_x/y/z - str - equation of the supplied numerical 1-form component x/y/z
                        It must be in terms of x, y, z.
                        Has to be given, for some methods to be calculatable.
        
        Returns: None
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

    

    def zoom(self, mag, target, dpd):
    
        """
        Redefines a new, zoomed in grid on which a new 1-form is calculated.
        Takes in three arguments:

        mag - magnitude of zoom ; 
        target - target point [x, y, z] ;
        dpd - points amount in the new grid ; 

        Returns a new, zoomed in 1-form 
        """

        if self.form_1_str_x == None or self.form_1_str_y == None or self.form_1_str_z == None:
            # ERROR
            raise TypeError('No equation provided, see \'give_eqn\' method')
        else:
             # Zoom must be one or greater
            if mag < 1:
                raise ValueError('Mag must be greater than one')
            else:

                # Target coordinates
                x_m = target[0]
                y_m = target[1]
                z_m = target[2]
                
                # Get the size of the original VF
                Lx = 0.5*(self.xg[-1, -1, -1] - self.xg[0, 0, 0])
                Ly = 0.5*(self.yg[-1, -1, -1] - self.yg[0, 0, 0])
                Lz = 0.5*(self.zg[-1, -1, -1] - self.zg[0, 0, 0])

                
                # Zoom axis range
                d_range_x = Lx/mag
                d_range_y = Ly/mag
                d_range_z = Lz/mag
                
                # Set up zoom window grids
                dx = np.linspace(-d_range_x + x_m, d_range_x + x_m, dpd)
                dy = np.linspace(-d_range_y + y_m, d_range_y + y_m, dpd)
                dz = np.linspace(-d_range_z + z_m, d_range_z + z_m, dpd)
                dxg, dyg, dzg = np.meshgrid(dx, dy, dz)
                
                # Create variables for the user provided equation strings
                u_str = self.form_1_str_x
                v_str = self.form_1_str_y
                k_str = self.form_1_str_z
                
                # Check if the equations provided contain x and y terms
                if u_str.find('x') & u_str.find('y') & u_str.find('z')== -1:
                    u_str = '(' + str(u_str) + ')* np.ones(np.shape(dxg))'
                else:
                    u_str = u_str.replace('x', 'dxg')
                    u_str = u_str.replace('y', 'dyg')
                    u_str = u_str.replace('z', 'dzg')
            
                if v_str.find('x') & v_str.find('y') & v_str.find('z') == -1:
                    v_str = '(' + str(v_str) + ')* np.ones(np.shape(dyg))'
                else:
                    v_str = v_str.replace('x', 'dxg')
                    v_str = v_str.replace('y', 'dyg')
                    v_str = v_str.replace('z', 'dzg')

                if k_str.find('x') & k_str.find('y') & k_str.find('z') == -1:
                    k_str = '(' + str(k_str) + ')* np.ones(np.shape(dzg))'
                else:
                    k_str = k_str.replace('x', 'dxg')
                    k_str = k_str.replace('y', 'dyg')
                    k_str = k_str.replace('z', 'dzg')
                    
                # Generate arrays for the components of the zoom field
                u_zoom = eval(u_str)
                v_zoom = eval(v_str)
                k_zoom = eval(k_str)
                
                # crate the zoomed in form
                zoom_form = form_1_3d(dxg, dyg, dzg, u_zoom, v_zoom, k_zoom, self.form_1_str_x, self.form_1_str_y, self.form_1_str_z)

                
                return zoom_form
 


    def set_density(self, points_number):
        '''
        Changes the desnity of points in the same range to the input value
        Requires the string equation to be supplied to not 'extrapolate'
        
        Only creates grid with same number of points in each
        direction. Cannot be used for any custom grids
        
        Parameters:
        --------------
        points_number - new number of points to use per axis
        
        Returns: None
        '''
        if self.form_1_str_x == None or self.form_1_str_y == None or self.form_1_str_z == None:
            # Error
            raise ValueError('Error: You need to supply the 1-form equation to do this, see \'give_eqn\' method')
        else:
            # redefine the grids
            x = np.linspace(self.xg[0,0,0], self.xg[-1,-1,0], points_number)
            y = np.linspace(self.yg[0,0,0], self.yg[-1,0,-1], points_number)
            z = np.linspace(self.zg[0,0,0], self.zg[0,-1,-1], points_number)
            self.xg, self.yg, self.zg= np.meshgrid(x,y,z)
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
            
            # check against constant forms, to have correct array shape
            if str_x.find('x') & str_x.find('y') & str_x.find('z') == -1:
                str_x = '(' + str(str_x) + ')* np.ones(np.shape(self.xg))'
            if str_y.find('x') & str_y.find('y') & str_y.find('z') == -1:
                str_y = '(' + str(str_y) + ')* np.ones(np.shape(self.yg))'
            if str_z.find('x') & str_z.find('y') & str_z.find('z') == -1:
                str_z = '(' + str(str_z) + ')* np.ones(np.shape(self.zg))'
        
            # re-evaluate the 1-form numerically
            self.F_x = eval(str_x)
            self.F_y = eval(str_y)
            self.F_z = eval(str_z)



    def plot(self, tip_width = None, tip_height = None, stack_side = None, opacity = None, singularity_size = None):

        '''
        plots the 1 form

        ---------------------------------
        Returns:

        None

        ----------------------------------
        Arguments:
        *Optional:

        tip_width - radius of the vector tip
        tip_height - length of the vector tip
        stack_side - dimensions of the stack plane
        opacity - opacity of the vector
        singularity_size - dimension of red spheres placed upon singular points
        '''

        # initial conditions
        tp_wdth = 0.05
        tp_hgth = 0.15
        side = 0.5
        opc = 1.0
        sng_size = 0.5*self.sng_scl

        # check for user input
        if tip_height != None:
            if isinstance(tip_height, float)==True and tip_height>=0.0:
                tp_hgth = tip_height       
            else:             
                print('The vector tip height has to be a float (0.0 < tip_height)')
                exit()

        if tip_width != None:
            if isinstance(tip_width, float)==True and tip_width>=0.0 and tip_width<=side/2:
                tp_wdth = tip_width       
            else:
                print('The vector tip radius has to be a float (0.0 < tip_width < stack_side/2)')
                exit()

        if stack_side != None:
            if isinstance(stack_side, float)==True and stack_side>=0.0:
                side = stack_side       
            else:
                print('The stack side size has to be a float (0.0 < stack_side)')
                exit()

        if opacity != None:
            if isinstance(opacity, float)==True and opacity>=0.0 and opacity<=1.0:
                opc = opacity      
            else:
                print('The opacity has to be a float (0.0 < opacity < 1.0)')
                exit()

        if singularity_size != None:
            if isinstance(singularity_size, float)==True and singularity_size>=0.0:
                sng_size = singularity_size      
            else:
                print('The size of the singularity sphere has to be a float (0.0 < singularity_size)')
                exit()


        # create mlab figure
        #v = self.scene

        # create a 2D array of points (each list contains a point from the meshgrid)
        pts = np.vstack(list(zip(self.xg.ravel(), self.yg.ravel(), self.zg.ravel())))

        # create a 2D array of points (each list contains vector field components corresponding to a point)
        F_list = np.vstack(list(zip(self.F_x.ravel(), self.F_y.ravel(), self.F_z.ravel())))



        ### Getting rid of zero points

        # indices of zero field
        Idx = np.argwhere(np.all(F_list==0,axis=1))

        # delete these points with 0 field since we do not want stacks to be plotted there
        pts_new = np.delete(pts, [Idx[:]], axis=0)
        F_new = np.delete(F_list, [Idx[:]], axis=0)



        ### Getting rid of singular points

        # np.inf -> np.nan
        F_new[np.isinf(F_new)] = np.nan

        # indices & points where field(at least one component) is not defined
        Idx_nan = np.argwhere(np.isnan(F_new))

        Idx_nan = Idx_nan[:,0]

        pts_nan = pts_new[Idx_nan]



        # delete singularities from list of points for stack plot
        pts_new = np.delete(pts_new, [Idx_nan], axis=0)
        F_new = np.delete(F_new, [Idx_nan], axis=0)



        # field_magnitude
        mag = np.sqrt(F_new[:,0]**2 + F_new[:,1]**2 + F_new[:,2]**2)



        if self.logarithmic_scale_bool:
            mag1 = mag + 1
            # min_size = np.min(mag1)
            
            unorm = F_new[:,0]/mag1
            vnorm = F_new[:,1]/mag1
            knorm = F_new[:,2]/mag1
            
            # logsf = np.log10(mag1/min_size)
            logmag = np.log10(mag1)
            F_new[:,0] = unorm*logmag
            F_new[:,1] = vnorm*logmag
            F_new[:,2] = knorm*logmag


        mag = np.sqrt(F_new[:,0]**2 + F_new[:,1]**2 + F_new[:,2]**2)
        mag_lst = np.vstack(list(zip(mag.ravel())))


        mag_max = np.nanmax(mag)
        mag_min = np.nanmin(mag)

        # divide the magnitude range into five parts (to scale stacks between 1-stack and 5-stack)
        sep = (mag_max-mag_min)/5

        # indices where magnitude of the field is
        # larger that (mag_min + (N)*sep)
        Idx1 = np.argwhere(np.all(mag_lst>=(mag_min+sep),axis=1))
        Idx2 = np.argwhere(np.all(mag_lst>=(mag_min+(2*sep)),axis=1))
        Idx3 = np.argwhere(np.all(mag_lst>=(mag_min+(3*sep)),axis=1))
        Idx4 = np.argwhere(np.all(mag_lst>=(mag_min+(4*sep)),axis=1))



        # Create list of points for each stack component. 1-stack will be plotted on each defined point and
        # a plane will be added at each point where the magnitude reaches the separation constant treshold
        pts1 = pts_new[Idx1]
        pts1 = np.vstack(list(zip(pts1[:,:,0].ravel(),pts1[:,:,1].ravel(),pts1[:,:,2].ravel())))
        pts2 = pts_new[Idx2]
        pts2 = np.vstack(list(zip(pts2[:,:,0].ravel(),pts2[:,:,1].ravel(),pts2[:,:,2].ravel())))
        pts3 = pts_new[Idx3]
        pts3 = np.vstack(list(zip(pts3[:,:,0].ravel(),pts3[:,:,1].ravel(),pts3[:,:,2].ravel())))
        pts4 = pts_new[Idx4]
        pts4 = np.vstack(list(zip(pts4[:,:,0].ravel(),pts4[:,:,1].ravel(),pts4[:,:,2].ravel())))



        # direction vector lists
        F_new_1 = F_new[Idx1]
        F_new_1 = np.vstack(list(zip(F_new_1[:,:,0].ravel(),F_new_1[:,:,1].ravel(),F_new_1[:,:,2].ravel())))
        F_new_2 = F_new[Idx2]
        F_new_2 = np.vstack(list(zip(F_new_2[:,:,0].ravel(),F_new_2[:,:,1].ravel(),F_new_2[:,:,2].ravel())))
        F_new_3 = F_new[Idx3]
        F_new_3 = np.vstack(list(zip(F_new_3[:,:,0].ravel(),F_new_3[:,:,1].ravel(),F_new_3[:,:,2].ravel())))
        F_new_4 = F_new[Idx4]
        F_new_4 = np.vstack(list(zip(F_new_4[:,:,0].ravel(),F_new_4[:,:,1].ravel(),F_new_4[:,:,2].ravel())))


        #-------------------------------------------------------------------------
        F = []  # create 1d list to input the directions of glyphs when plotting
                # (2d arrays not allowed for input). That is the only purpose
                # of these lists
        for i in range(len(F_new)):
            F.append(F_new[i])


        F1 = []
        for i in range(len(F_new_1)):
            F1.append(F_new_1[i])

        F2 = []
        for i in range(len(F_new_2)):
            F2.append(F_new_2[i])

        F3 = []
        for i in range(len(F_new_3)):
            F3.append(F_new_3[i])

        F4 = []
        for i in range(len(F_new_4)):
            F4.append(F_new_4[i])
        #-------------------------------------------------------------------------



        ### --------------- PLOTTING --------------------------------------


        # Define glyphs (Cone for tip, shrinked in x direction box for stack plane)
        # Each consequent box is shifted along the unit vector (initially x direction)
        # to not overlap with previous plane => create a stack of increased density
        cyl1 = tvtk.ConeSource(radius = tp_wdth*self.scaling,
                                height = tp_hgth*self.scaling,
                                capping = False,
                                center = (0.075*self.scaling, 0, 0),
                                    )

        box = tvtk.CubeSource(x_length=0.01,
                            y_length = side*self.scaling,
                            z_length = side*self.scaling)

        box1 = tvtk.CubeSource(x_length=0.01,
                            y_length = side*self.scaling,
                            z_length = side*self.scaling,
                            center = (-0.04, 0, 0))

        box2 = tvtk.CubeSource(x_length=0.01,
                            y_length = side*self.scaling,
                            z_length = side*self.scaling,
                            center = (-0.08, 0, 0))

        box3 = tvtk.CubeSource(x_length=0.01,
                            y_length = side*self.scaling,
                            z_length = side*self.scaling,
                            center = (-0.12, 0, 0))

        box4 = tvtk.CubeSource(x_length=0.01,
                            y_length = side*self.scaling,
                            z_length = side*self.scaling,
                            center = (-0.16, 0, 0))


        # Create input for the glyph -- the sources are placed at these input
        # points.
        pd = tvtk.PolyData(points=pts_new)
        pd.point_data.vectors = F
        pd1 = tvtk.PolyData(points=pts1)
        pd1.point_data.vectors = F1
        pd2 = tvtk.PolyData(points=pts2)
        pd2.point_data.vectors = F2
        pd3 = tvtk.PolyData(points=pts3)
        pd3.point_data.vectors = F3
        pd4 = tvtk.PolyData(points=pts4)
        pd4.point_data.vectors = F4

        # define Glyph objects which serve as an indicator that glyph sources will be placed
        # upon the points which are used as input data
        g = tvtk.Glyph3D(scale_mode='data_scaling_off', vector_mode = 'use_vector')
        gb = tvtk.Glyph3D(scale_mode='data_scaling_off', vector_mode = 'use_vector')
        g1 = tvtk.Glyph3D(scale_mode='data_scaling_off', vector_mode = 'use_vector')
        g2 = tvtk.Glyph3D(scale_mode='data_scaling_off', vector_mode = 'use_vector')
        g3 = tvtk.Glyph3D(scale_mode='data_scaling_off', vector_mode = 'use_vector')
        g4 = tvtk.Glyph3D(scale_mode='data_scaling_off', vector_mode = 'use_vector')

        # corresponding points (PolyData) to corresponding Glyph object
        configure_input_data(g, pd)
        configure_input_data(gb, pd)
        configure_input_data(g1, pd1)
        configure_input_data(g2, pd2)
        configure_input_data(g3, pd3)
        configure_input_data(g4, pd4)

        # geometrical sources for glyph object at each PolyData point set
        configure_source_data(g, cyl1.output)
        cyl1.update()
        g.update()
        configure_source_data(gb, box.output)
        box.update()
        gb.update()
        configure_source_data(g1, box1.output)
        box1.update()
        g1.update()
        configure_source_data(g2, box2.output)
        box2.update()
        g2.update()
        configure_source_data(g3, box3.output)
        box3.update()
        g3.update()
        configure_source_data(g4, box4.output)
        box4.update()
        g4.update()


        # mapper which maps the glyphs
        m = tvtk.PolyDataMapper()
        mb = tvtk.PolyDataMapper()
        m1 = tvtk.PolyDataMapper()
        m2 = tvtk.PolyDataMapper()
        m3 = tvtk.PolyDataMapper()
        m4 = tvtk.PolyDataMapper()

        # properties of geometrical objects

        #0.565,0.641,0.46
        pc = tvtk.Property(opacity=opc, color=(0.64,0.008,0.87), edge_visibility='y', edge_color=(0,0,0))
        p1 = tvtk.Property(opacity=opc, color=(0.64,0.008,0.87), edge_visibility='y', edge_color=(0,0,0))

        # map glyphs(which now represent the geometrical objects)
        # to PolyData sets
        configure_input_data(m, g.output)
        configure_input_data(mb, gb.output)
        configure_input_data(m1, g1.output)
        configure_input_data(m2, g2.output)
        configure_input_data(m3, g3.output)
        configure_input_data(m4, g4.output)

        # Create actors (objects created from the Glyphs and PolyData sets)
        # to be represented on a figure
        a = tvtk.Actor(mapper=m, property=pc)
        ab = tvtk.Actor(mapper=mb, property=p1)
        a1 = tvtk.Actor(mapper=m1, property=p1)
        a2 = tvtk.Actor(mapper=m2, property=p1)
        a3 = tvtk.Actor(mapper=m3, property=p1)
        a4 = tvtk.Actor(mapper=m4, property=p1)

        # add actors to the figure
        # only add 1-stack for constant fields


        # axes boundaries
        xmin = float((np.min(self.xg)) - 1)
        ymin = float((np.min(self.yg)) - 1)
        zmin = float((np.min(self.zg)) - 1)
        xmax = float((np.max(self.xg)) + 1)
        ymax = float((np.max(self.yg)) + 1)
        zmax = float((np.max(self.zg)) + 1)

        #----------------------------------------
        """
        v = figure

        if (mag_max-mag_min)<=0.05:
            v.scene.add_actor(a)
            v.scene.add_actor(ab)
        else:
            v.scene.add_actor(a)
            v.scene.add_actor(ab)
            v.scene.add_actor(a1)
            v.scene.add_actor(a2)
            v.scene.add_actor(a3)
            v.scene.add_actor(a4)    


        # plot red spheres at the singularities
        # add axes and outline
        mlab.points3d(pts_nan[:,0],
                    pts_nan[:,1],
                    pts_nan[:,2], color = (1,0,0),scale_factor=sng_size, resolution=36)
        mlab.outline(extent = [xmin,xmax,ymin,ymax,zmin,zmax], line_width=1.0)
        mlab.axes(color = (0,0,0), nb_labels = 5, extent = [xmin,xmax,ymin,ymax,zmin,zmax], line_width=3.0)

        

        mlab.show()
        """
        #--------------------------------------------

        stck_coords = [a, ab, a1, a2, a3, a4]
        red_balls_data = [pts_nan[:,0], pts_nan[:,1], pts_nan[:,2], sng_size]
        axes_limits = [xmin,xmax,ymin,ymax,zmin,zmax]

        return stck_coords, red_balls_data, axes_limits







class form_2_3d():


   


    def __init__(self, xg, yg, zg, Fx=None, Fy=None, Fz=None, Fx_eqn=None, Fy_eqn=None, Fz_eqn=None):
        self.xg = xg
        self.yg = yg
        self.zg = zg
        self.Fx = Fx
        self.Fy = Fy
        self.Fz = Fz
        self.s_max = 6
        self.s_min = 2
        self.pt_den_x = len(xg[0, :, :])
        self.pt_den_y = len(yg[:, 0, :])
        self.pt_den_z = len(zg[:, :, 0])
        self.fract_x = 2/((self.pt_den_x - 1))
        self.fract_y = 2/((self.pt_den_y - 1))
        self.fract_z = 2/((self.pt_den_z - 1))
        self.colour_list = [(1,0,0), (0,0,1)]
        self.logarithmic_scale_bool = 0
        # self.base = 10
        self.delta_factor = 10
        if Fx_eqn is not None:
            self.Fx_eqn = str(simplify(Fx_eqn))  # to start with, user must change to access some methods
            # Note, the string must be given with x and y as variables
        else:
            self.Fx_eqn = None
        if Fy_eqn is not None:
            self.Fy_eqn = str(simplify(Fy_eqn))  # to start with, user must change to access some methods
            # Note, the string must be given with x and y as variables
        else:
            self.Fy_eqn = None
        if Fz_eqn is not None:
            self.Fz_eqn = str(simplify(Fz_eqn))  # to start with, user must change to access some methods
            # Note, the string must be given with x and y as variables
        else:
            self.Fz_eqn = None
    
    

    def give_eqn(self, equation_str_x, equation_str_y, equation_str_z ):
        '''
        Takes in 1-argument, string
        This must be the equation of the supplied numerical 0-form
        It must be in terms of x and y.
        Has to be given, for some methods to be calculatable.
        '''
        self.Fx_eqn = str(simplify(equation_str_x))
        self.Fy_eqn = str(simplify(equation_str_y))
        self.Fz_eqn = str(simplify(equation_str_z))
        
        # update the numerical values to always match
        string_x = self.Fx_eqn + ''
        string_x = string_x.replace('x', '(self.xg)')
        string_x = string_x.replace('y', '(self.yg)')
        string_x = string_x.replace('z', '(self.zg)')

        string_y = self.Fy_eqn + ''
        string_y = string_y.replace('x', '(self.xg)')
        string_y = string_y.replace('y', '(self.yg)')
        string_y = string_y.replace('z', '(self.zg)')

        string_z = self.Fz_eqn + ''
        string_z = string_z.replace('x', '(self.xg)')
        string_z = string_z.replace('y', '(self.yg)')
        string_z = string_z.replace('z', '(self.zg)')
        
        # correct for consatnt form before evaluating
        if string_x.find('x') & string_x.find('y') & string_x.find('z') == -1:
            string_x = '(' + str(string_x) + ')* np.ones(np.shape(self.xg))'
        else:
            pass

        if string_y.find('x') & string_y.find('y') & string_y.find('z') == -1:
            string_y = '(' + str(string_y) + ')* np.ones(np.shape(self.yg))'
        else:
            pass
        
        if string_z.find('x') & string_z.find('y') & string_z.find('z') == -1:
            string_z = '(' + str(string_z) + ')* np.ones(np.shape(self.zg))'
        else:
            pass

        # re-evaluate the 2-form numerically
        self.Fx = eval(string_x)
        self.Fy = eval(string_y)
        self.Fz = eval(string_z)
    
    
    
    def return_string(self):
        '''
        Takes in no arguments, returns the unformatted string back to user
        This is done in case user wants to access strings
        that got here not by input but by ext. alg.
        '''
        return self.Fx_eqn, self.Fy_eqn, self.Fz_eqn
    
    

    def log_scaling(self):
        '''
        Takes no arguments
        Changes the boolean that determines if scaling is logarithmic
        Whenever it is called, it changes that boolean to opposite
        The form object is initialised with this as False (as 0)
        '''
        self.logarithmic_scale_bool = not self.logarithmic_scale_bool

    
    
    def zoom(self, mag, target, dpd):
    
        """
        Redefines a new, zoomed in grid on which a new 2-form is calculated.
        Takes in three arguments:

        mag - magnitude of zoom ; 
        target - target point [x, y, z] ;
        dpd - points amount in the new grid ; 

        Returns a new, zoomed in 2-form 
        """

        if self.Fx_eqn == None or self.Fy_eqn == None or self.Fz_eqn == None:
            # ERROR
            raise TypeError('No equation provided, see \'give_eqn\' method')
        else:
             # Zoom must be one or greater
            if mag < 1:
                raise ValueError('Mag must be greater than one')
            else:

                # Target coordinates
                x_m = target[0]
                y_m = target[1]
                z_m = target[2]
                
                # Get the size of the original VF
                Lx = 0.5*(self.xg[-1, -1, -1] - self.xg[0, 0, 0])
                Ly = 0.5*(self.yg[-1, -1, -1] - self.yg[0, 0, 0])
                Lz = 0.5*(self.zg[-1, -1, -1] - self.zg[0, 0, 0])

                
                # Zoom axis range
                d_range_x = Lx/mag
                d_range_y = Ly/mag
                d_range_z = Lz/mag
                
                # Set up zoom window grids
                dx = np.linspace(-d_range_x + x_m, d_range_x + x_m, dpd)
                dy = np.linspace(-d_range_y + y_m, d_range_y + y_m, dpd)
                dz = np.linspace(-d_range_z + z_m, d_range_z + z_m, dpd)
                dxg, dyg, dzg = np.meshgrid(dx, dy, dz)
                
                # Create variables for the user provided equation strings
                u_str = self.Fx_eqn
                v_str = self.Fy_eqn
                k_str = self.Fz_eqn
                
                # Check if the equations provided contain x and y terms
                if u_str.find('x') & u_str.find('y') & u_str.find('z')== -1:
                    u_str = '(' + str(u_str) + ')* np.ones(np.shape(dxg))'
                else:
                    u_str = u_str.replace('x', 'dxg')
                    u_str = u_str.replace('y', 'dyg')
                    u_str = u_str.replace('z', 'dzg')
            
                if v_str.find('x') & v_str.find('y') & v_str.find('z') == -1:
                    v_str = '(' + str(v_str) + ')* np.ones(np.shape(dyg))'
                else:
                    v_str = v_str.replace('x', 'dxg')
                    v_str = v_str.replace('y', 'dyg')
                    v_str = v_str.replace('z', 'dzg')

                if k_str.find('x') & k_str.find('y') & k_str.find('z') == -1:
                    k_str = '(' + str(k_str) + ')* np.ones(np.shape(dzg))'
                else:
                    k_str = k_str.replace('x', 'dxg')
                    k_str = k_str.replace('y', 'dyg')
                    k_str = k_str.replace('z', 'dzg')
                    
                # Generate arrays for the components of the zoom field
                u_zoom = eval(u_str)
                v_zoom = eval(v_str)
                k_zoom = eval(k_str)
                
                # crate the zoomed in form
                zoom_form = form_2_3d(dxg, dyg, dzg, u_zoom, v_zoom, k_zoom, self.Fx_eqn, self.Fy_eqn, self.Fz_eqn)

                
                return zoom_form
 
    
    
    def set_density2(self, points_number_x, points_number_y, points_number_z):
        '''
        
        Changes number of points on grids to given, if equations have been given
        
        Parameters:
        -------------
        points_number_x - int - number of points to put along the x axis
        points_number_y - int - number of points to put along the y axis
        points_number_z - int - number of points to put along the z axis
        
        Returns: None
        '''
        if self.form_2_str == None:
            # Error
            raise TypeError('Error: You need to supply the 2-form equation to do this, look at \'give_eqn\' method')
        else:
            # redefine the grids
            x = np.linspace(self.xg[0,0,0], self.xg[0,-1,0], points_number_x)
            y = np.linspace(self.yg[0,0,0], self.yg[-1,0,0], points_number_y)
            z = np.linspace(self.zg[0,0,0], self.zg[0,0,-1], points_number_z)
            self.xg, self.yg, self.zg = np.meshgrid(x, y, z)
            # based on these change other, dependant variables
            self.pt_den_x = len(self.xg[0, :, :])
            self.pt_den_y = len(self.yg[:, 0, :])
            self.pt_den_z = len(self.yg[:, :, 0])
            self.fract_x = 2/(self.pt_den_x - 1)
            self.fract_y = 2/(self.pt_den_y - 1)
            self.fract_z = 2/(self.pt_den_z - 1)
            # substitute these into the equation:
            # but keep it local
            str_2 = self.form_2_str + ''
            str_2 = str_2.replace('x', '(self.xg)')
            str_2 = str_2.replace('y', '(self.yg)')
            str_2 = str_2.replace('z', '(self.zg)')
            
            # correct for consatnt form before evaluating
            if str_2.find('x') & str_2.find('y') & str_2.find('z') == -1:
                str_2 = '(' + str(str_2) + ')* np.ones(np.shape(self.xg))'
            else:
                pass
            # re-evaluate the 2-form numerically
            self.form_2 = eval(str_2)
    
    

    def plot(self, scn):
        
        """
        Plots the 2-form. Retunrs nothing, takes in no arguments
        """

        if self.Fx is None and self.Fy is None and self.Fz is None:
            print('Please, provide at least one component of the field')
        
        Fz = self.Fz
        Fx = self.Fx
        Fy = self.Fy

        Fmag = np.sqrt(Fx**2 + Fy**2 + Fz**2)



        def form_2(F, direction):
            


            gr_sep = abs(self.xg[1,1,1]-self.xg[0,0,0])

            pts = np.vstack(list(zip(self.xg.ravel(), self.yg.ravel(), self.zg.ravel())))

            mag_lst = np.vstack(list(zip(F.ravel())))
            mag_lst[np.isinf(mag_lst)] = np.nan
            Idx_nan = np.argwhere(np.isnan(mag_lst))

            Idx_nan = Idx_nan[:,0]

            pts_nan = pts[Idx_nan]

            pts = np.delete(pts, [Idx_nan], axis=0)
            mag_lst = np.delete(mag_lst, [Idx_nan], axis=0)

            f_max = np.nanmax(abs(F))

            sep = (f_max)/5

            

            ### Find points where field has different magnitudes

            Idx1 = np.argwhere(np.all(mag_lst>(0),axis=1) & np.all(mag_lst<=(0+sep),axis=1))
            Idx2 = np.argwhere(np.all(mag_lst>(0+sep),axis=1) & np.all(mag_lst<=(0+2*sep),axis=1))
            Idx3 = np.argwhere(np.all(mag_lst>(0+2*sep),axis=1) & np.all(mag_lst<=(0+3*sep),axis=1))
            Idx4 = np.argwhere(np.all(mag_lst>(0+3*sep),axis=1) & np.all(mag_lst<=(0+4*sep),axis=1))
            Idx5 = np.argwhere(np.all(mag_lst>(0+4*sep),axis=1))

            Idx_1 = np.argwhere(np.all(mag_lst<(0),axis=1) & np.all(mag_lst>=(0-sep),axis=1))
            Idx_2 = np.argwhere(np.all(mag_lst<(0-sep),axis=1) & np.all(mag_lst>=(0-2*sep),axis=1))
            Idx_3 = np.argwhere(np.all(mag_lst<(0-2*sep),axis=1) & np.all(mag_lst>=(0-3*sep),axis=1))
            Idx_4 = np.argwhere(np.all(mag_lst<(0-3*sep),axis=1) & np.all(mag_lst>=(0-4*sep),axis=1))
            Idx_5 = np.argwhere(np.all(mag_lst<(0-4*sep),axis=1))



            pts1 = pts[Idx1[:,0]]
            pts_1 = pts[Idx_1[:,0]]
            pts2 = pts[Idx2[:,0]]
            pts_2 = pts[Idx_2[:,0]]
            pts3 = pts[Idx3[:,0]]
            pts_3 = pts[Idx_3[:,0]]
            pts4 = pts[Idx4[:,0]]
            pts_4 = pts[Idx_4[:,0]]
            pts5 = pts[Idx5[:,0]]
            pts_5 = pts[Idx_5[:,0]]


            #----points for lines------------------


            if direction=='z':
            #------1 mag-----------

                lnx1 = pts1+[0,gr_sep/2,gr_sep/2], pts1+[0,gr_sep/2,-gr_sep/2], pts1+[0,-gr_sep/2,gr_sep/2], pts1+[0,-gr_sep/2,-gr_sep/2]
                lnx1 = np.concatenate(lnx1, axis=0)

                lny1 = pts1+[gr_sep/2,0,gr_sep/2], pts1+[gr_sep/2,0,-gr_sep/2], pts1+[-gr_sep/2,0,gr_sep/2], pts1+[-gr_sep/2,0,-gr_sep/2]
                lny1 = np.concatenate(lny1, axis=0)

                lnz1 = pts1+[gr_sep/2,gr_sep/2,0], pts1+[gr_sep/2,-gr_sep/2,0], pts1+[-gr_sep/2,gr_sep/2,0], pts1+[-gr_sep/2,-gr_sep/2,0]
                lnz1 = np.concatenate(lnz1, axis=0)


                lnx_1 = pts_1+[0,gr_sep/2,gr_sep/2], pts_1+[0,gr_sep/2,-gr_sep/2], pts_1+[0,-gr_sep/2,gr_sep/2], pts_1+[0,-gr_sep/2,-gr_sep/2]
                lnx_1 = np.concatenate(lnx_1, axis=0)

                lny_1 = pts_1+[gr_sep/2,0,gr_sep/2], pts_1+[gr_sep/2,0,-gr_sep/2], pts_1+[-gr_sep/2,0,gr_sep/2], pts_1+[-gr_sep/2,0,-gr_sep/2]
                lny_1 = np.concatenate(lny_1, axis=0)

                lnz_1 = pts_1+[gr_sep/2,gr_sep/2,0], pts_1+[gr_sep/2,-gr_sep/2,0], pts_1+[-gr_sep/2,gr_sep/2,0], pts_1+[-gr_sep/2,-gr_sep/2,0]
                lnz_1 = np.concatenate(lnz_1, axis=0)


                #--------2 mag--------------



                lnx2 = pts2+[0,gr_sep/2,gr_sep/2], pts2+[0,gr_sep/2,-gr_sep/2], pts2+[0,-gr_sep/2,gr_sep/2], pts2+[0,-gr_sep/2,-gr_sep/2],\
                    pts2+[0,0,gr_sep/2], pts2+[0,0,-gr_sep/2]
                lnx2 = np.concatenate(lnx2, axis=0)

                lny2 = pts2+[gr_sep/2,0,gr_sep/2], pts2+[gr_sep/2,0,-gr_sep/2], pts2+[-gr_sep/2,0,gr_sep/2], pts2+[-gr_sep/2,0,-gr_sep/2],\
                    pts2+[0,0,gr_sep/2], pts2+[0,0,-gr_sep/2]
                lny2 = np.concatenate(lny2, axis=0)

                lnz2 = pts2+[gr_sep/2,gr_sep/2,0], pts2+[gr_sep/2,-gr_sep/2,0], pts2+[-gr_sep/2,gr_sep/2,0], pts2+[-gr_sep/2,-gr_sep/2,0],\
                    pts2+[gr_sep/2,0,0], pts2+[-gr_sep/2,0,0], pts2+[0,gr_sep/2,0], pts2+[0,-gr_sep/2,0], pts2+[0,0,0]
                lnz2 = np.concatenate(lnz2, axis=0)


                lnx_2 = pts_2+[0,gr_sep/2,gr_sep/2], pts_2+[0,gr_sep/2,-gr_sep/2], pts_2+[0,-gr_sep/2,gr_sep/2], pts_2+[0,-gr_sep/2,-gr_sep/2],\
                        pts_2+[0,0,gr_sep/2], pts_2+[0,0,-gr_sep/2]
                lnx_2 = np.concatenate(lnx_2, axis=0)

                lny_2 = pts_2+[gr_sep/2,0,gr_sep/2], pts_2+[gr_sep/2,0,-gr_sep/2], pts_2+[-gr_sep/2,0,gr_sep/2], pts_2+[-gr_sep/2,0,-gr_sep/2],\
                        pts_2+[0,0,gr_sep/2], pts_2+[0,0,-gr_sep/2]
                lny_2 = np.concatenate(lny_2, axis=0)

                lnz_2 = pts_2+[gr_sep/2,gr_sep/2,0], pts_2+[gr_sep/2,-gr_sep/2,0], pts_2+[-gr_sep/2,gr_sep/2,0], pts_2+[-gr_sep/2,-gr_sep/2,0],\
                        pts_2+[gr_sep/2,0,0], pts_2+[-gr_sep/2,0,0], pts_2+[0,gr_sep/2,0], pts_2+[0,-gr_sep/2,0],\
                        pts_2+[0,0,0]
                lnz_2 = np.concatenate(lnz_2, axis=0)


                #----------3 mag----------------


                lnx3 = pts3+[0,gr_sep/2,gr_sep/2], pts3+[0,gr_sep/2,-gr_sep/2], pts3+[0,-gr_sep/2,gr_sep/2], pts3+[0,-gr_sep/2,-gr_sep/2],\
                    pts3+[0,gr_sep/6,gr_sep/2], pts3+[0,gr_sep/6,-gr_sep/2], pts3+[0,-gr_sep/6,gr_sep/2], pts3+[0,-gr_sep/6,-gr_sep/2]
                lnx3 = np.concatenate(lnx3, axis=0)

                lny3 = pts3+[gr_sep/2,0,gr_sep/2], pts3+[gr_sep/2,0,-gr_sep/2], pts3+[-gr_sep/2,0,gr_sep/2], pts3+[-gr_sep/2,0,-gr_sep/2],\
                    pts3+[gr_sep/6,0,gr_sep/2], pts3+[gr_sep/6,0,-gr_sep/2], pts3+[-gr_sep/6,0,gr_sep/2], pts3+[-gr_sep/6,0,-gr_sep/2]
                lny3 = np.concatenate(lny3, axis=0)

                lnz3 = pts3+[gr_sep/2,gr_sep/2,0], pts3+[gr_sep/2,-gr_sep/2,0], pts3+[-gr_sep/2,gr_sep/2,0], pts3+[-gr_sep/2,-gr_sep/2,0],\
                    pts3+[-gr_sep/2,gr_sep/6,0],pts3+[-gr_sep/2,-gr_sep/6,0],pts3+[gr_sep/2,gr_sep/6,0],pts3+[gr_sep/2,-gr_sep/6,0],\
                    pts3+[-gr_sep/6,-gr_sep/2,0],pts3+[-gr_sep/6,gr_sep/2,0],pts3+[gr_sep/6,gr_sep/2,0],pts3+[gr_sep/6,-gr_sep/2,0],\
                    pts3+[-gr_sep/6,-gr_sep/6,0],pts3+[-gr_sep/6,gr_sep/6,0],pts3+[gr_sep/6,gr_sep/6,0],pts3+[gr_sep/6,-gr_sep/6,0]
                lnz3 = np.concatenate(lnz3, axis=0)


                lnx_3 = pts_3+[0,gr_sep/2,gr_sep/2], pts_3+[0,gr_sep/2,-gr_sep/2], pts_3+[0,-gr_sep/2,gr_sep/2], pts_3+[0,-gr_sep/2,-gr_sep/2],\
                        pts_3+[0,gr_sep/6,gr_sep/2], pts_3+[0,gr_sep/6,-gr_sep/2], pts_3+[0,-gr_sep/6,gr_sep/2], pts_3+[0,-gr_sep/6,-gr_sep/2]
                lnx_3 = np.concatenate(lnx_3, axis=0)

                lny_3 = pts_3+[gr_sep/2,0,gr_sep/2], pts_3+[gr_sep/2,0,-gr_sep/2], pts_3+[-gr_sep/2,0,gr_sep/2], pts_3+[-gr_sep/2,0,-gr_sep/2],\
                        pts_3+[gr_sep/6,0,gr_sep/2], pts_3+[gr_sep/6,0,-gr_sep/2], pts_3+[-gr_sep/6,0,gr_sep/2], pts_3+[-gr_sep/6,0,-gr_sep/2]
                lny_3 = np.concatenate(lny_3, axis=0)

                lnz_3 = pts_3+[gr_sep/2,gr_sep/2,0], pts_3+[gr_sep/2,-gr_sep/2,0], pts_3+[-gr_sep/2,gr_sep/2,0], pts_3+[-gr_sep/2,-gr_sep/2,0],\
                        pts_3+[-gr_sep/2,gr_sep/6,0],pts_3+[-gr_sep/2,-gr_sep/6,0],pts_3+[gr_sep/2,gr_sep/6,0],pts_3+[gr_sep/2,-gr_sep/6,0],\
                        pts_3+[-gr_sep/6,-gr_sep/2,0],pts_3+[-gr_sep/6,gr_sep/2,0],pts_3+[gr_sep/6,gr_sep/2,0],pts_3+[gr_sep/6,-gr_sep/2,0],\
                        pts_3+[-gr_sep/6,-gr_sep/6,0],pts_3+[-gr_sep/6,gr_sep/6,0],pts_3+[gr_sep/6,gr_sep/6,0],pts_3+[gr_sep/6,-gr_sep/6,0]
                lnz_3 = np.concatenate(lnz_3, axis=0)



                lnx4 = pts4+[0,gr_sep/2,gr_sep/2], pts4+[0,gr_sep/2,-gr_sep/2], pts4+[0,-gr_sep/2,gr_sep/2], pts4+[0,-gr_sep/2,-gr_sep/2],\
                       pts4+[0,0,-gr_sep/2], pts4+[0,-gr_sep/4,-gr_sep/2], pts4+[0,gr_sep/4,-gr_sep/2],\
                       pts4+[0,0,gr_sep/2], pts4+[0,-gr_sep/4,gr_sep/2], pts4+[0,gr_sep/4,gr_sep/2]
                lnx4 = np.concatenate(lnx4, axis=0)

                lny4 = pts4+[gr_sep/2,0,gr_sep/2], pts4+[gr_sep/2,0,-gr_sep/2], pts4+[-gr_sep/2,0,gr_sep/2], pts4+[-gr_sep/2,0,-gr_sep/2],\
                       pts4+[0,0,-gr_sep/2], pts4+[-gr_sep/4,0,-gr_sep/2], pts4+[gr_sep/4,0,-gr_sep/2],\
                       pts4+[0,0,gr_sep/2], pts4+[-gr_sep/4,0,gr_sep/2], pts4+[gr_sep/4,0,gr_sep/2]
                lny4 = np.concatenate(lny4, axis=0)

                lnz4 = pts4+[gr_sep/2,gr_sep/2,0], pts4+[gr_sep/2,-gr_sep/2,0], pts4+[-gr_sep/2,gr_sep/2,0], pts4+[-gr_sep/2,-gr_sep/2,0],\
                       pts4+[0,0,0],\
                       pts4+[-gr_sep/4,0,0], pts4+[gr_sep/4,0,0], pts4+[0,-gr_sep/4,0], pts4+[0,gr_sep/4,0],\
                       pts4+[-gr_sep/4,-gr_sep/4,0], pts4+[-gr_sep/4,gr_sep/4,0], pts4+[gr_sep/4,-gr_sep/4,0], pts4+[gr_sep/4,gr_sep/4,0],\
                       pts4+[-gr_sep/2,0,0], pts4+[gr_sep/2,0,0], pts4+[0,-gr_sep/2,0], pts4+[0,gr_sep/2,0],\
                       pts4+[-gr_sep/2,gr_sep/4,0], pts4+[-gr_sep/2,-gr_sep/4,0], pts4+[gr_sep/2,-gr_sep/4,0], pts4+[gr_sep/2,gr_sep/4,0],\
                       pts4+[-gr_sep/4,gr_sep/2,0], pts4+[-gr_sep/4,-gr_sep/2,0], pts4+[gr_sep/4,-gr_sep/2,0], pts4+[gr_sep/4,gr_sep/2,0]
                lnz4 = np.concatenate(lnz4, axis=0)


                lnx_4 = pts_4+[0,gr_sep/2,gr_sep/2], pts_4+[0,gr_sep/2,-gr_sep/2], pts_4+[0,-gr_sep/2,gr_sep/2], pts_4+[0,-gr_sep/2,-gr_sep/2],\
                        pts_4+[0,0,-gr_sep/2], pts_4+[0,-gr_sep/4,-gr_sep/2], pts_4+[0,gr_sep/4,-gr_sep/2],\
                        pts_4+[0,0,gr_sep/2], pts_4+[0,-gr_sep/4,gr_sep/2], pts_4+[0,gr_sep/4,gr_sep/2]
                lnx_4 = np.concatenate(lnx_4, axis=0)

                lny_4 = pts_4+[gr_sep/2,0,gr_sep/2], pts_4+[gr_sep/2,0,-gr_sep/2], pts_4+[-gr_sep/2,0,gr_sep/2], pts_4+[-gr_sep/2,0,-gr_sep/2],\
                        pts_4+[0,0,-gr_sep/2], pts_4+[-gr_sep/4,0,-gr_sep/2], pts_4+[gr_sep/4,0,-gr_sep/2],\
                        pts_4+[0,0,gr_sep/2], pts_4+[-gr_sep/4,0,gr_sep/2], pts_4+[gr_sep/4,0,gr_sep/2]
                lny_4 = np.concatenate(lny_4, axis=0)

                lnz_4 = pts_4+[gr_sep/2,gr_sep/2,0], pts_4+[gr_sep/2,-gr_sep/2,0], pts_4+[-gr_sep/2,gr_sep/2,0], pts_4+[-gr_sep/2,-gr_sep/2,0],\
                        pts_4+[0,0,0],\
                        pts_4+[-gr_sep/4,0,0], pts_4+[gr_sep/4,0,0], pts_4+[0,-gr_sep/4,0], pts_4+[0,gr_sep/4,0],\
                        pts_4+[-gr_sep/4,-gr_sep/4,0], pts_4+[-gr_sep/4,gr_sep/4,0], pts_4+[gr_sep/4,-gr_sep/4,0], pts_4+[gr_sep/4,gr_sep/4,0],\
                        pts_4+[-gr_sep/2,0,0], pts_4+[gr_sep/2,0,0], pts_4+[0,-gr_sep/2,0], pts_4+[0,gr_sep/2,0],\
                        pts_4+[-gr_sep/2,gr_sep/4,0], pts_4+[-gr_sep/2,-gr_sep/4,0], pts_4+[gr_sep/2,-gr_sep/4,0], pts_4+[gr_sep/2,gr_sep/4,0],\
                        pts_4+[-gr_sep/4,gr_sep/2,0], pts_4+[-gr_sep/4,-gr_sep/2,0], pts_4+[gr_sep/4,-gr_sep/2,0], pts_4+[gr_sep/4,gr_sep/2,0]
                lnz_4 = np.concatenate(lnz_4, axis=0)



                #--------5 mag-------------------


                lnx5 = pts5+[0,gr_sep/2,gr_sep/2], pts5+[0,gr_sep/2,-gr_sep/2], pts5+[0,-gr_sep/2,gr_sep/2], pts5+[0,-gr_sep/2,-gr_sep/2],\
                       pts5+[0,-gr_sep/10,-gr_sep/2], pts5+[0,gr_sep/10,-gr_sep/2], pts5+[0,-(3*gr_sep/10),-gr_sep/2], pts5+[0,(3*gr_sep/10),-gr_sep/2],\
                       pts5+[0,-gr_sep/10,gr_sep/2], pts5+[0,gr_sep/10,gr_sep/2], pts5+[0,-(3*gr_sep/10),gr_sep/2], pts5+[0,(3*gr_sep/10),gr_sep/2]
                lnx5 = np.concatenate(lnx5, axis=0)

                lny5 = pts5+[gr_sep/2,0,gr_sep/2], pts5+[gr_sep/2,0,-gr_sep/2], pts5+[-gr_sep/2,0,gr_sep/2], pts5+[-gr_sep/2,0,-gr_sep/2],\
                       pts5+[-gr_sep/10,0,-gr_sep/2], pts5+[gr_sep/10,0,-gr_sep/2], pts5+[-(3*gr_sep/10),0,-gr_sep/2], pts5+[(3*gr_sep/10),0,-gr_sep/2],\
                       pts5+[-gr_sep/10,0,gr_sep/2], pts5+[gr_sep/10,0,gr_sep/2], pts5+[-(3*gr_sep/10),0,gr_sep/2], pts5+[(3*gr_sep/10),0,gr_sep/2]  
                lny5 = np.concatenate(lny5, axis=0)

                lnz5 = pts5+[gr_sep/2,gr_sep/2,0], pts5+[gr_sep/2,-gr_sep/2,0], pts5+[-gr_sep/2,gr_sep/2,0], pts5+[-gr_sep/2,-gr_sep/2,0],\
                       pts5+[-gr_sep/10,gr_sep/10,0], pts5+[-gr_sep/10,-gr_sep/10,0], pts5+[gr_sep/10,gr_sep/10,0], pts5+[gr_sep/10,-gr_sep/10,0],\
                       pts5+[-(3*gr_sep/10),gr_sep/10,0],pts5+[-(3*gr_sep/10),-gr_sep/10,0], pts5+[-(3*gr_sep/10),(3*gr_sep/10),0], pts5+[-(3*gr_sep/10),-(3*gr_sep/10),0],\
                       pts5+[(3*gr_sep/10),gr_sep/10,0],pts5+[(3*gr_sep/10),-gr_sep/10,0], pts5+[(3*gr_sep/10),(3*gr_sep/10),0], pts5+[(3*gr_sep/10),-(3*gr_sep/10),0],\
                       pts5+[-(gr_sep/10),(3*gr_sep/10),0], pts5+[-(gr_sep/10),-(3*gr_sep/10),0], pts5+[(gr_sep/10),(3*gr_sep/10),0], pts5+[(gr_sep/10),-(3*gr_sep/10),0],\
                       pts5+[-(gr_sep/2),gr_sep/10,0],pts5+[-(gr_sep/2),-gr_sep/10,0], pts5+[-(gr_sep/2),(3*gr_sep/10),0], pts5+[-(gr_sep/2),-(3*gr_sep/10),0],\
                       pts5+[(gr_sep/2),gr_sep/10,0],pts5+[(gr_sep/2),-gr_sep/10,0], pts5+[(gr_sep/2),(3*gr_sep/10),0], pts5+[(gr_sep/2),-(3*gr_sep/10),0],\
                       pts5+[gr_sep/10,-(gr_sep/2),0],pts5+[-gr_sep/10,-(gr_sep/2),0], pts5+[(3*gr_sep/10),-(gr_sep/2),0], pts5+[-(3*gr_sep/10),-(gr_sep/2),0],\
                       pts5+[gr_sep/10,(gr_sep/2),0],pts5+[-gr_sep/10,(gr_sep/2),0], pts5+[(3*gr_sep/10),(gr_sep/2),0], pts5+[-(3*gr_sep/10),(gr_sep/2),0]
                lnz5 = np.concatenate(lnz5, axis=0)




                lnx_5 = pts_5+[0,gr_sep/2,gr_sep/2], pts_5+[0,gr_sep/2,-gr_sep/2], pts_5+[0,-gr_sep/2,gr_sep/2], pts_5+[0,-gr_sep/2,-gr_sep/2],\
                        pts_5+[0,-gr_sep/10,-gr_sep/2], pts_5+[0,gr_sep/10,-gr_sep/2], pts_5+[0,-(3*gr_sep/10),-gr_sep/2], pts_5+[0,(3*gr_sep/10),-gr_sep/2],\
                        pts_5+[0,-gr_sep/10,gr_sep/2], pts_5+[0,gr_sep/10,gr_sep/2], pts_5+[0,-(3*gr_sep/10),gr_sep/2], pts_5+[0,(3*gr_sep/10),gr_sep/2]
                lnx_5 = np.concatenate(lnx_5, axis=0)

                lny_5 = pts_5+[gr_sep/2,0,gr_sep/2], pts_5+[gr_sep/2,0,-gr_sep/2], pts_5+[-gr_sep/2,0,gr_sep/2], pts_5+[-gr_sep/2,0,-gr_sep/2],\
                        pts_5+[-gr_sep/10,0,-gr_sep/2], pts_5+[gr_sep/10,0,-gr_sep/2], pts_5+[-(3*gr_sep/10),0,-gr_sep/2], pts_5+[(3*gr_sep/10),0,-gr_sep/2],\
                        pts_5+[-gr_sep/10,0,gr_sep/2], pts_5+[gr_sep/10,0,gr_sep/2], pts_5+[-(3*gr_sep/10),0,gr_sep/2], pts_5+[(3*gr_sep/10),0,gr_sep/2]  
                lny_5 = np.concatenate(lny_5, axis=0)

                lnz_5 = pts_5+[gr_sep/2,gr_sep/2,0], pts_5+[gr_sep/2,-gr_sep/2,0], pts_5+[-gr_sep/2,gr_sep/2,0], pts_5+[-gr_sep/2,-gr_sep/2,0],\
                        pts_5+[-gr_sep/10,gr_sep/10,0], pts_5+[-gr_sep/10,-gr_sep/10,0], pts_5+[gr_sep/10,gr_sep/10,0], pts_5+[gr_sep/10,-gr_sep/10,0],\
                        pts_5+[-(3*gr_sep/10),gr_sep/10,0],pts_5+[-(3*gr_sep/10),-gr_sep/10,0], pts_5+[-(3*gr_sep/10),(3*gr_sep/10),0], pts_5+[-(3*gr_sep/10),-(3*gr_sep/10),0],\
                        pts_5+[(3*gr_sep/10),gr_sep/10,0],pts_5+[(3*gr_sep/10),-gr_sep/10,0], pts_5+[(3*gr_sep/10),(3*gr_sep/10),0], pts_5+[(3*gr_sep/10),-(3*gr_sep/10),0],\
                        pts_5+[-(gr_sep/10),(3*gr_sep/10),0], pts_5+[-(gr_sep/10),-(3*gr_sep/10),0], pts_5+[(gr_sep/10),(3*gr_sep/10),0], pts_5+[(gr_sep/10),-(3*gr_sep/10),0],\
                        pts_5+[-(gr_sep/2),gr_sep/10,0],pts_5+[-(gr_sep/2),-gr_sep/10,0], pts_5+[-(gr_sep/2),(3*gr_sep/10),0], pts_5+[-(gr_sep/2),-(3*gr_sep/10),0],\
                        pts_5+[(gr_sep/2),gr_sep/10,0],pts_5+[(gr_sep/2),-gr_sep/10,0], pts_5+[(gr_sep/2),(3*gr_sep/10),0], pts_5+[(gr_sep/2),-(3*gr_sep/10),0],\
                        pts_5+[gr_sep/10,-(gr_sep/2),0],pts_5+[-gr_sep/10,-(gr_sep/2),0], pts_5+[(3*gr_sep/10),-(gr_sep/2),0], pts_5+[-(3*gr_sep/10),-(gr_sep/2),0],\
                        pts_5+[gr_sep/10,(gr_sep/2),0],pts_5+[-gr_sep/10,(gr_sep/2),0], pts_5+[(3*gr_sep/10),(gr_sep/2),0], pts_5+[-(3*gr_sep/10),(gr_sep/2),0]
                lnz_5 = np.concatenate(lnz_5, axis=0)




            if direction=='x':
            #------1 mag-----------

                lnx1 = pts1+[0,gr_sep/2,gr_sep/2], pts1+[0,gr_sep/2,-gr_sep/2], pts1+[0,-gr_sep/2,gr_sep/2], pts1+[0,-gr_sep/2,-gr_sep/2]
                lnx1 = np.concatenate(lnx1, axis=0)

                lny1 = pts1+[gr_sep/2,0,gr_sep/2], pts1+[gr_sep/2,0,-gr_sep/2], pts1+[-gr_sep/2,0,gr_sep/2], pts1+[-gr_sep/2,0,-gr_sep/2]
                lny1 = np.concatenate(lny1, axis=0)

                lnz1 = pts1+[gr_sep/2,gr_sep/2,0], pts1+[gr_sep/2,-gr_sep/2,0], pts1+[-gr_sep/2,gr_sep/2,0], pts1+[-gr_sep/2,-gr_sep/2,0]
                lnz1 = np.concatenate(lnz1, axis=0)


                lnx_1 = pts_1+[0,gr_sep/2,gr_sep/2], pts_1+[0,gr_sep/2,-gr_sep/2], pts_1+[0,-gr_sep/2,gr_sep/2], pts_1+[0,-gr_sep/2,-gr_sep/2]
                lnx_1 = np.concatenate(lnx_1, axis=0)

                lny_1 = pts_1+[gr_sep/2,0,gr_sep/2], pts_1+[gr_sep/2,0,-gr_sep/2], pts_1+[-gr_sep/2,0,gr_sep/2], pts_1+[-gr_sep/2,0,-gr_sep/2]
                lny_1 = np.concatenate(lny_1, axis=0)

                lnz_1 = pts_1+[gr_sep/2,gr_sep/2,0], pts_1+[gr_sep/2,-gr_sep/2,0], pts_1+[-gr_sep/2,gr_sep/2,0], pts_1+[-gr_sep/2,-gr_sep/2,0]
                lnz_1 = np.concatenate(lnz_1, axis=0)


                #--------2 mag--------------



                lnz2 = pts2+[gr_sep/2,gr_sep/2,0], pts2+[gr_sep/2,-gr_sep/2,0], pts2+[-gr_sep/2,gr_sep/2,0], pts2+[-gr_sep/2,-gr_sep/2,0],\
                       pts2+[gr_sep/2,0,0], pts2+[-gr_sep/2,0,0]
                lnz2 = np.concatenate(lnz2, axis=0)

                lny2 = pts2+[gr_sep/2,0,gr_sep/2], pts2+[gr_sep/2,0,-gr_sep/2], pts2+[-gr_sep/2,0,gr_sep/2], pts2+[-gr_sep/2,0,-gr_sep/2],\
                       pts2+[gr_sep/2,0,0], pts2+[-gr_sep/2,0,0]
                lny2 = np.concatenate(lny2, axis=0)

                lnx2 = pts2+[0,gr_sep/2,gr_sep/2], pts2+[0,-gr_sep/2,gr_sep/2], pts2+[0,gr_sep/2,-gr_sep/2], pts2+[0,-gr_sep/2,-gr_sep/2],\
                       pts2+[0,0,gr_sep/2], pts2+[0,0,-gr_sep/2], pts2+[0,-gr_sep/2,0], pts2+[0,gr_sep/2,0],\
                       pts2+[0,0,0]
                lnx2 = np.concatenate(lnx2, axis=0)


                lnz_2 = pts_2+[gr_sep/2,gr_sep/2,0], pts_2+[gr_sep/2,-gr_sep/2,0], pts_2+[-gr_sep/2,gr_sep/2,0], pts_2+[-gr_sep/2,-gr_sep/2,0],\
                        pts_2+[gr_sep/2,0,0], pts_2+[-gr_sep/2,0,0]
                lnz_2 = np.concatenate(lnz_2, axis=0)

                lny_2 = pts_2+[gr_sep/2,0,gr_sep/2], pts_2+[gr_sep/2,0,-gr_sep/2], pts_2+[-gr_sep/2,0,gr_sep/2], pts_2+[-gr_sep/2,0,-gr_sep/2],\
                        pts_2+[gr_sep/2,0,0], pts_2+[-gr_sep/2,0,0]
                lny_2 = np.concatenate(lny_2, axis=0)

                lnx_2 = pts_2+[0,gr_sep/2,gr_sep/2], pts_2+[0,-gr_sep/2,gr_sep/2], pts_2+[0,gr_sep/2,-gr_sep/2], pts_2+[0,-gr_sep/2,-gr_sep/2],\
                        pts_2+[0,0,gr_sep/2], pts_2+[0,0,-gr_sep/2], pts_2+[0,-gr_sep/2,0], pts_2+[0,gr_sep/2,0],\
                        pts_2+[0,0,0]
                lnx_2 = np.concatenate(lnx_2, axis=0)


                #----------3 mag----------------


                lnz3 = pts3+[gr_sep/2,gr_sep/2,0], pts3+[gr_sep/2,-gr_sep/2,0], pts3+[-gr_sep/2,gr_sep/2,0], pts3+[-gr_sep/2,-gr_sep/2,0],\
                        pts3+[gr_sep/2,gr_sep/6,0], pts3+[gr_sep/2,-gr_sep/6,0], pts3+[-gr_sep/2,gr_sep/6,0], pts3+[-gr_sep/2,-gr_sep/6,0]
                lnz3 = np.concatenate(lnz3, axis=0)

                lny3 = pts3+[gr_sep/2,0,gr_sep/2], pts3+[gr_sep/2,0,-gr_sep/2], pts3+[-gr_sep/2,0,gr_sep/2], pts3+[-gr_sep/2,0,-gr_sep/2],\
                        pts3+[gr_sep/2,0,gr_sep/6], pts3+[gr_sep/2,0,-gr_sep/6], pts3+[-gr_sep/2,0,gr_sep/6], pts3+[-gr_sep/2,0,-gr_sep/6]
                lny3 = np.concatenate(lny3, axis=0)

                lnx3 = pts3+[0,gr_sep/2,gr_sep/2], pts3+[0,gr_sep/2,-gr_sep/2], pts3+[0,-gr_sep/2,gr_sep/2], pts3+[0,-gr_sep/2,-gr_sep/2],\
                        pts3+[0,-gr_sep/2,gr_sep/6],pts3+[0,-gr_sep/2,-gr_sep/6],pts3+[0,gr_sep/2,gr_sep/6],pts3+[0,gr_sep/2,-gr_sep/6],\
                        pts3+[0,-gr_sep/6,-gr_sep/2],pts3+[0,-gr_sep/6,gr_sep/2],pts3+[0,gr_sep/6,gr_sep/2],pts3+[0,gr_sep/6,-gr_sep/2],\
                        pts3+[0,-gr_sep/6,-gr_sep/6],pts3+[0,-gr_sep/6,gr_sep/6],pts3+[0,gr_sep/6,gr_sep/6],pts3+[0,gr_sep/6,-gr_sep/6]
                lnx3 = np.concatenate(lnx3, axis=0)


                lnz_3 = pts_3+[gr_sep/2,gr_sep/2,0], pts_3+[gr_sep/2,-gr_sep/2,0], pts_3+[-gr_sep/2,gr_sep/2,0], pts_3+[-gr_sep/2,-gr_sep/2,0],\
                        pts_3+[gr_sep/2,gr_sep/6,0], pts_3+[gr_sep/2,-gr_sep/6,0], pts_3+[-gr_sep/2,gr_sep/6,0], pts_3+[-gr_sep/2,-gr_sep/6,0]
                lnz_3 = np.concatenate(lnz_3, axis=0)

                lny_3 = pts_3+[gr_sep/2,0,gr_sep/2], pts_3+[gr_sep/2,0,-gr_sep/2], pts_3+[-gr_sep/2,0,gr_sep/2], pts_3+[-gr_sep/2,0,-gr_sep/2],\
                        pts_3+[gr_sep/2,0,gr_sep/6], pts_3+[gr_sep/2,0,-gr_sep/6], pts_3+[-gr_sep/2,0,gr_sep/6], pts_3+[-gr_sep/2,0,-gr_sep/6]
                lny_3 = np.concatenate(lny_3, axis=0)

                lnx_3 = pts_3+[0,gr_sep/2,gr_sep/2], pts_3+[0,gr_sep/2,-gr_sep/2], pts_3+[0,-gr_sep/2,gr_sep/2], pts_3+[0,-gr_sep/2,-gr_sep/2],\
                        pts_3+[0,-gr_sep/2,gr_sep/6],pts_3+[0,-gr_sep/2,-gr_sep/6],pts_3+[0,gr_sep/2,gr_sep/6],pts_3+[0,gr_sep/2,-gr_sep/6],\
                        pts_3+[0,-gr_sep/6,-gr_sep/2],pts_3+[0,-gr_sep/6,gr_sep/2],pts_3+[0,gr_sep/6,gr_sep/2],pts_3+[0,gr_sep/6,-gr_sep/2],\
                        pts_3+[0,-gr_sep/6,-gr_sep/6],pts_3+[0,-gr_sep/6,gr_sep/6],pts_3+[0,gr_sep/6,gr_sep/6],pts_3+[0,gr_sep/6,-gr_sep/6]
                lnx_3 = np.concatenate(lnx_3, axis=0)


            #----------4 mag----------------


                lnz4 = pts4+[gr_sep/2,gr_sep/2,0], pts4+[gr_sep/2,-gr_sep/2,0], pts4+[-gr_sep/2,gr_sep/2,0], pts4+[-gr_sep/2,-gr_sep/2,0],\
                    pts4+[-gr_sep/2,0,0], pts4+[-gr_sep/2,-gr_sep/4,0], pts4+[-gr_sep/2,gr_sep/4,0],\
                    pts4+[gr_sep/2,0,0], pts4+[gr_sep/2,-gr_sep/4,0], pts4+[gr_sep/2,gr_sep/4,0]
                lnz4 = np.concatenate(lnz4, axis=0)

                lny4 = pts4+[gr_sep/2,0,gr_sep/2], pts4+[gr_sep/2,0,-gr_sep/2], pts4+[-gr_sep/2,0,gr_sep/2], pts4+[-gr_sep/2,0,-gr_sep/2],\
                    pts4+[-gr_sep/2,0,0], pts4+[-gr_sep/2,0,-gr_sep/4], pts4+[gr_sep/2,0,-gr_sep/4],\
                    pts4+[gr_sep/2,0,0], pts4+[-gr_sep/2,0,gr_sep/4], pts4+[gr_sep/2,0,gr_sep/4]
                lny4 = np.concatenate(lny4, axis=0)

                lnx4 = pts4+[0,gr_sep/2,gr_sep/2], pts4+[0,-gr_sep/2,gr_sep/2], pts4+[0,gr_sep/2,-gr_sep/2], pts4+[0,-gr_sep/2,-gr_sep/2],\
                    pts4+[0,0,0],pts4+[0,0,gr_sep/4],pts4+[0,0,-gr_sep/4],\
                    pts4+[0,-gr_sep/4,0], pts4+[0,gr_sep/4,0],\
                    pts4+[0,-gr_sep/4,-gr_sep/4], pts4+[0,gr_sep/4,-gr_sep/4], pts4+[0,-gr_sep/4,gr_sep/4], pts4+[0,gr_sep/4,gr_sep/4],\
                    pts4+[0,0,-gr_sep/2], pts4+[0,0,gr_sep/2], pts4+[0,-gr_sep/2,0], pts4+[0,gr_sep/2,0],\
                    pts4+[0,gr_sep/4,-gr_sep/2], pts4+[0,-gr_sep/4,-gr_sep/2], pts4+[0,-gr_sep/4,gr_sep/2], pts4+[0,gr_sep/4,gr_sep/2],\
                    pts4+[0,gr_sep/2,-gr_sep/4], pts4+[0,-gr_sep/2,-gr_sep/4], pts4+[0,-gr_sep/2,gr_sep/4], pts4+[0,gr_sep/2,gr_sep/4]
                lnx4 = np.concatenate(lnx4, axis=0)


                lnz_4 = pts_4+[gr_sep/2,gr_sep/2,0], pts_4+[gr_sep/2,-gr_sep/2,0], pts_4+[-gr_sep/2,gr_sep/2,0], pts_4+[-gr_sep/2,-gr_sep/2,0],\
                    pts_4+[-gr_sep/2,0,0], pts_4+[-gr_sep/2,-gr_sep/4,0], pts_4+[-gr_sep/2,gr_sep/4,0],\
                    pts_4+[gr_sep/2,0,0], pts_4+[gr_sep/2,-gr_sep/4,0], pts_4+[gr_sep/2,gr_sep/4,0]
                lnz_4 = np.concatenate(lnz_4, axis=0)

                lny_4 = pts_4+[gr_sep/2,0,gr_sep/2], pts_4+[gr_sep/2,0,-gr_sep/2], pts_4+[-gr_sep/2,0,gr_sep/2], pts_4+[-gr_sep/2,0,-gr_sep/2],\
                    pts_4+[-gr_sep/2,0,0], pts_4+[-gr_sep/2,0,-gr_sep/4], pts_4+[gr_sep/2,0,-gr_sep/4],\
                    pts_4+[gr_sep/2,0,0], pts_4+[-gr_sep/2,0,gr_sep/4], pts_4+[gr_sep/2,0,gr_sep/4]
                lny_4 = np.concatenate(lny_4, axis=0)

                lnx_4 = pts_4+[0,gr_sep/2,gr_sep/2], pts_4+[0,-gr_sep/2,gr_sep/2], pts_4+[0,gr_sep/2,-gr_sep/2], pts_4+[0,-gr_sep/2,-gr_sep/2],\
                    pts_4+[0,0,0],pts_4+[0,0,gr_sep/4],pts_4+[0,0,-gr_sep/4],\
                    pts_4+[0,-gr_sep/4,0], pts_4+[0,gr_sep/4,0],\
                    pts_4+[0,-gr_sep/4,-gr_sep/4], pts_4+[0,gr_sep/4,-gr_sep/4], pts_4+[0,-gr_sep/4,gr_sep/4], pts_4+[0,gr_sep/4,gr_sep/4],\
                    pts_4+[0,0,-gr_sep/2], pts_4+[0,0,gr_sep/2], pts_4+[0,-gr_sep/2,0], pts_4+[0,gr_sep/2,0],\
                    pts_4+[0,gr_sep/4,-gr_sep/2], pts_4+[0,-gr_sep/4,-gr_sep/2], pts_4+[0,-gr_sep/4,gr_sep/2], pts_4+[0,gr_sep/4,gr_sep/2],\
                    pts_4+[0,gr_sep/2,-gr_sep/4], pts_4+[0,-gr_sep/2,-gr_sep/4], pts_4+[0,-gr_sep/2,gr_sep/4], pts_4+[0,gr_sep/2,gr_sep/4]
                lnx_4 = np.concatenate(lnx_4, axis=0)



                #--------5 mag-------------------


                lnz5 = pts5+[gr_sep/2,gr_sep/2,0], pts5+[gr_sep/2,-gr_sep/2,0], pts5+[-gr_sep/2,gr_sep/2,0], pts5+[-gr_sep/2,-gr_sep/2,0],\
                        pts5+[-gr_sep/2,-gr_sep/10,0], pts5+[-gr_sep/2,gr_sep/10,0], pts5+[-gr_sep/2,-(3*gr_sep/10),0], pts5+[-gr_sep/2,(3*gr_sep/10),0],\
                        pts5+[gr_sep/2,-gr_sep/10,0], pts5+[gr_sep/2,gr_sep/10,0], pts5+[gr_sep/2,-(3*gr_sep/10),0], pts5+[gr_sep/2,(3*gr_sep/10),0]
                lnz5 = np.concatenate(lnz5, axis=0)

                lny5 = pts5+[gr_sep/2,0,gr_sep/2], pts5+[gr_sep/2,0,-gr_sep/2], pts5+[-gr_sep/2,0,gr_sep/2], pts5+[-gr_sep/2,0,-gr_sep/2],\
                        pts5+[-gr_sep/2,0,-gr_sep/10], pts5+[gr_sep/2,0,-gr_sep/10], pts5+[-gr_sep/2,0,-(3*gr_sep/10)], pts5+[-gr_sep/2,0,(3*gr_sep/10)],\
                        pts5+[-gr_sep/2,0,gr_sep/10], pts5+[gr_sep/2,0,gr_sep/10], pts5+[gr_sep/2,0,-(3*gr_sep/10)], pts5+[gr_sep/2,0,(3*gr_sep/10)]  
                lny5 = np.concatenate(lny5, axis=0)

                lnx5 = pts5+[0,gr_sep/2,gr_sep/2], pts_5+[0,gr_sep/2,-gr_sep/2], pts5+[0,-gr_sep/2,gr_sep/2], pts5+[0,-gr_sep/2,-gr_sep/2],\
                        pts5+[0,-gr_sep/10,gr_sep/10], pts5+[0,-gr_sep/10,-gr_sep/10], pts5+[0,gr_sep/10,gr_sep/10], pts5+[0,gr_sep/10,-gr_sep/10],\
                        pts5+[0,-(3*gr_sep/10),gr_sep/10],pts5+[0,-(3*gr_sep/10),-gr_sep/10], pts5+[0,-(3*gr_sep/10),(3*gr_sep/10)], pts5+[0,-(3*gr_sep/10),-(3*gr_sep/10)],\
                        pts5+[0,(3*gr_sep/10),gr_sep/10],pts5+[0,(3*gr_sep/10),-gr_sep/10], pts5+[0,(3*gr_sep/10),(3*gr_sep/10)], pts5+[0,(3*gr_sep/10),-(3*gr_sep/10)],\
                        pts5+[0,-(gr_sep/10),(3*gr_sep/10)], pts5+[0,-(gr_sep/10),-(3*gr_sep/10)], pts5+[0,(gr_sep/10),(3*gr_sep/10)], pts5+[0,(gr_sep/10),-(3*gr_sep/10)],\
                        pts5+[0,-(gr_sep/2),gr_sep/10],pts5+[0,-(gr_sep/2),-gr_sep/10], pts5+[0,-(gr_sep/2),(3*gr_sep/10)], pts5+[0,-(gr_sep/2),-(3*gr_sep/10)],\
                        pts5+[0,(gr_sep/2),gr_sep/10],pts5+[0,(gr_sep/2),-gr_sep/10], pts5+[0,(gr_sep/2),(3*gr_sep/10)], pts5+[0,(gr_sep/2),-(3*gr_sep/10)],\
                        pts5+[0,gr_sep/10,-(gr_sep/2)],pts5+[0,-gr_sep/10,-(gr_sep/2)], pts5+[0,(3*gr_sep/10),-(gr_sep/2)], pts5+[0,-(3*gr_sep/10),-(gr_sep/2)],\
                        pts5+[0,gr_sep/10,(gr_sep/2)],pts5+[0,-gr_sep/10,(gr_sep/2)], pts5+[0,(3*gr_sep/10),(gr_sep/2)], pts5+[0,-(3*gr_sep/10),(gr_sep/2)]
                lnx5 = np.concatenate(lnx5, axis=0)




                lnz_5 = pts_5+[gr_sep/2,gr_sep/2,0], pts_5+[gr_sep/2,-gr_sep/2,0], pts_5+[-gr_sep/2,gr_sep/2,0], pts_5+[-gr_sep/2,-gr_sep/2,0],\
                        pts_5+[-gr_sep/2,-gr_sep/10,0], pts_5+[-gr_sep/2,gr_sep/10,0], pts_5+[-gr_sep/2,-(3*gr_sep/10),0], pts_5+[-gr_sep/2,(3*gr_sep/10),0],\
                        pts_5+[gr_sep/2,-gr_sep/10,0], pts_5+[gr_sep/2,gr_sep/10,0], pts_5+[gr_sep/2,-(3*gr_sep/10),0], pts_5+[gr_sep/2,(3*gr_sep/10),0]
                lnz_5 = np.concatenate(lnz_5, axis=0)

                lny_5 = pts_5+[gr_sep/2,0,gr_sep/2], pts_5+[gr_sep/2,0,-gr_sep/2], pts_5+[-gr_sep/2,0,gr_sep/2], pts_5+[-gr_sep/2,0,-gr_sep/2],\
                        pts_5+[-gr_sep/2,0,-gr_sep/10], pts_5+[gr_sep/2,0,-gr_sep/10], pts_5+[-gr_sep/2,0,-(3*gr_sep/10)], pts_5+[-gr_sep/2,0,(3*gr_sep/10)],\
                        pts_5+[-gr_sep/2,0,gr_sep/10], pts_5+[gr_sep/2,0,gr_sep/10], pts_5+[gr_sep/2,0,-(3*gr_sep/10)], pts_5+[gr_sep/2,0,(3*gr_sep/10)]  
                lny_5 = np.concatenate(lny_5, axis=0)

                lnx_5 = pts_5+[0,gr_sep/2,gr_sep/2], pts_5+[0,gr_sep/2,-gr_sep/2], pts_5+[0,-gr_sep/2,gr_sep/2], pts_5+[0,-gr_sep/2,-gr_sep/2],\
                        pts_5+[0,-gr_sep/10,gr_sep/10], pts_5+[0,-gr_sep/10,-gr_sep/10], pts_5+[0,gr_sep/10,gr_sep/10], pts_5+[0,gr_sep/10,-gr_sep/10],\
                        pts_5+[0,-(3*gr_sep/10),gr_sep/10],pts_5+[0,-(3*gr_sep/10),-gr_sep/10], pts_5+[0,-(3*gr_sep/10),(3*gr_sep/10)], pts_5+[0,-(3*gr_sep/10),-(3*gr_sep/10)],\
                        pts_5+[0,(3*gr_sep/10),gr_sep/10],pts_5+[0,(3*gr_sep/10),-gr_sep/10], pts_5+[0,(3*gr_sep/10),(3*gr_sep/10)], pts_5+[0,(3*gr_sep/10),-(3*gr_sep/10)],\
                        pts_5+[0,-(gr_sep/10),(3*gr_sep/10)], pts_5+[0,-(gr_sep/10),-(3*gr_sep/10)], pts_5+[0,(gr_sep/10),(3*gr_sep/10)], pts_5+[0,(gr_sep/10),-(3*gr_sep/10)],\
                        pts_5+[0,-(gr_sep/2),gr_sep/10],pts_5+[0,-(gr_sep/2),-gr_sep/10], pts_5+[0,-(gr_sep/2),(3*gr_sep/10)], pts_5+[0,-(gr_sep/2),-(3*gr_sep/10)],\
                        pts_5+[0,(gr_sep/2),gr_sep/10],pts_5+[0,(gr_sep/2),-gr_sep/10], pts_5+[0,(gr_sep/2),(3*gr_sep/10)], pts_5+[0,(gr_sep/2),-(3*gr_sep/10)],\
                        pts_5+[0,gr_sep/10,-(gr_sep/2)],pts_5+[0,-gr_sep/10,-(gr_sep/2)], pts_5+[0,(3*gr_sep/10),-(gr_sep/2)], pts_5+[0,-(3*gr_sep/10),-(gr_sep/2)],\
                        pts_5+[0,gr_sep/10,(gr_sep/2)],pts_5+[0,-gr_sep/10,(gr_sep/2)], pts_5+[0,(3*gr_sep/10),(gr_sep/2)], pts_5+[0,-(3*gr_sep/10),(gr_sep/2)]
                lnx_5 = np.concatenate(lnx_5, axis=0)





            if direction=='y':
            #------1 mag-----------

                lnx1 = pts1+[0,gr_sep/2,gr_sep/2], pts1+[0,gr_sep/2,-gr_sep/2], pts1+[0,-gr_sep/2,gr_sep/2], pts1+[0,-gr_sep/2,-gr_sep/2]
                lnx1 = np.concatenate(lnx1, axis=0)

                lny1 = pts1+[gr_sep/2,0,gr_sep/2], pts1+[gr_sep/2,0,-gr_sep/2], pts1+[-gr_sep/2,0,gr_sep/2], pts1+[-gr_sep/2,0,-gr_sep/2]
                lny1 = np.concatenate(lny1, axis=0)

                lnz1 = pts1+[gr_sep/2,gr_sep/2,0], pts1+[gr_sep/2,-gr_sep/2,0], pts1+[-gr_sep/2,gr_sep/2,0], pts1+[-gr_sep/2,-gr_sep/2,0]
                lnz1 = np.concatenate(lnz1, axis=0)


                lnx_1 = pts_1+[0,gr_sep/2,gr_sep/2], pts_1+[0,gr_sep/2,-gr_sep/2], pts_1+[0,-gr_sep/2,gr_sep/2], pts_1+[0,-gr_sep/2,-gr_sep/2]
                lnx_1 = np.concatenate(lnx_1, axis=0)

                lny_1 = pts_1+[gr_sep/2,0,gr_sep/2], pts_1+[gr_sep/2,0,-gr_sep/2], pts_1+[-gr_sep/2,0,gr_sep/2], pts_1+[-gr_sep/2,0,-gr_sep/2]
                lny_1 = np.concatenate(lny_1, axis=0)

                lnz_1 = pts_1+[gr_sep/2,gr_sep/2,0], pts_1+[gr_sep/2,-gr_sep/2,0], pts_1+[-gr_sep/2,gr_sep/2,0], pts_1+[-gr_sep/2,-gr_sep/2,0]
                lnz_1 = np.concatenate(lnz_1, axis=0)


                #--------2 mag--------------



                lnx2 = pts2+[0,gr_sep/2,gr_sep/2], pts2+[0,gr_sep/2,-gr_sep/2], pts2+[0,-gr_sep/2,gr_sep/2], pts2+[0,-gr_sep/2,-gr_sep/2],\
                       pts2+[0,gr_sep/2,0], pts2+[0,-gr_sep/2,0]
                lnx2 = np.concatenate(lnx2, axis=0)

                lnz2 = pts2+[gr_sep/2,gr_sep/2,0], pts2+[gr_sep/2,-gr_sep/2,0], pts2+[-gr_sep/2,gr_sep/2,0], pts2+[-gr_sep/2,-gr_sep/2,0],\
                       pts2+[0,gr_sep/2,0], pts2+[0,-gr_sep/2,0]
                lnz2 = np.concatenate(lnz2, axis=0)

                lny2 = pts2+[gr_sep/2,0,gr_sep/2], pts2+[gr_sep/2,0,-gr_sep/2], pts2+[-gr_sep/2,0,gr_sep/2], pts2+[-gr_sep/2,0,-gr_sep/2],\
                        pts2+[gr_sep/2,0,0], pts2+[-gr_sep/2,0,0], pts2+[0,0,gr_sep/2], pts2+[0,0,-gr_sep/2],\
                        pts2+[0,0,0]
                lny2 = np.concatenate(lny2, axis=0)


                lnx_2 = pts_2+[0,gr_sep/2,gr_sep/2], pts_2+[0,gr_sep/2,-gr_sep/2], pts_2+[0,-gr_sep/2,gr_sep/2], pts_2+[0,-gr_sep/2,-gr_sep/2],\
                        pts_2+[0,gr_sep/2,0], pts_2+[0,-gr_sep/2,0]
                lnx_2 = np.concatenate(lnx_2, axis=0)

                lnz_2 = pts_2+[gr_sep/2,gr_sep/2,0], pts_2+[gr_sep/2,-gr_sep/2,0], pts_2+[-gr_sep/2,gr_sep/2,0], pts_2+[-gr_sep/2,-gr_sep/2,0],\
                        pts_2+[0,gr_sep/2,0], pts_2+[0,-gr_sep/2,0]
                lnz_2 = np.concatenate(lnz_2, axis=0)

                lny_2 = pts_2+[gr_sep/2,0,gr_sep/2], pts_2+[gr_sep/2,0,-gr_sep/2], pts_2+[-gr_sep/2,0,gr_sep/2], pts_2+[-gr_sep/2,0,-gr_sep/2],\
                        pts_2+[gr_sep/2,0,0], pts_2+[-gr_sep/2,0,0], pts_2+[0,0,gr_sep/2], pts_2+[0,0,-gr_sep/2],\
                        pts_2+[0,0,0]
                lny_2 = np.concatenate(lny_2, axis=0)


                #----------3 mag----------------


                lnx3 = pts3+[0,gr_sep/2,gr_sep/2], pts3+[0,gr_sep/2,-gr_sep/2], pts3+[0,-gr_sep/2,gr_sep/2], pts3+[0,-gr_sep/2,-gr_sep/2],\
                       pts3+[0,gr_sep/2,gr_sep/6], pts3+[0,gr_sep/2,-gr_sep/6], pts3+[0,-gr_sep/2,gr_sep/6], pts3+[0,-gr_sep/2,-gr_sep/6]
                lnx3 = np.concatenate(lnx3, axis=0)

                lnz3 = pts3+[gr_sep/2,gr_sep/2,0], pts3+[gr_sep/2,-gr_sep/2,0], pts3+[-gr_sep/2,gr_sep/2,0], pts3+[-gr_sep/2,-gr_sep/2,0],\
                       pts3+[gr_sep/6,gr_sep/2,0], pts3+[gr_sep/6,-gr_sep/2,0], pts3+[-gr_sep/6,gr_sep/2,0], pts3+[-gr_sep/6,-gr_sep/2,0]
                lnz3 = np.concatenate(lnz3, axis=0)

                lny3 = pts3+[gr_sep/2,0,gr_sep/2], pts3+[gr_sep/2,0,-gr_sep/2], pts3+[-gr_sep/2,0,gr_sep/2], pts3+[-gr_sep/2,0,-gr_sep/2],\
                       pts3+[-gr_sep/2,0,gr_sep/6],pts3+[-gr_sep/2,0,-gr_sep/6],pts3+[gr_sep/2,0,gr_sep/6],pts3+[gr_sep/2,0,-gr_sep/6],\
                       pts3+[-gr_sep/6,0,-gr_sep/2],pts3+[-gr_sep/6,0,gr_sep/2],pts3+[gr_sep/6,0,gr_sep/2],pts3+[gr_sep/6,0,-gr_sep/2],\
                       pts3+[-gr_sep/6,0,-gr_sep/6],pts3+[-gr_sep/6,0,gr_sep/6],pts3+[gr_sep/6,0,gr_sep/6],pts3+[gr_sep/6,0,-gr_sep/6]
                lny3 = np.concatenate(lny3, axis=0)



                lnx_3 = pts_3+[0,gr_sep/2,gr_sep/2], pts_3+[0,gr_sep/2,-gr_sep/2], pts_3+[0,-gr_sep/2,gr_sep/2], pts_3+[0,-gr_sep/2,-gr_sep/2],\
                    pts_3+[0,gr_sep/2,gr_sep/6], pts_3+[0,gr_sep/2,-gr_sep/6], pts_3+[0,-gr_sep/2,gr_sep/6], pts_3+[0,-gr_sep/2,-gr_sep/6]
                lnx_3 = np.concatenate(lnx_3, axis=0)

                lnz_3 = pts_3+[gr_sep/2,gr_sep/2,0], pts_3+[gr_sep/2,-gr_sep/2,0], pts_3+[-gr_sep/2,gr_sep/2,0], pts_3+[-gr_sep/2,-gr_sep/2,0],\
                    pts_3+[gr_sep/6,gr_sep/2,0], pts_3+[gr_sep/6,-gr_sep/2,0], pts_3+[-gr_sep/6,gr_sep/2,0], pts_3+[-gr_sep/6,-gr_sep/2,0]
                lnz_3 = np.concatenate(lnz_3, axis=0)

                lny_3 = pts_3+[gr_sep/2,0,gr_sep/2], pts_3+[gr_sep/2,0,-gr_sep/2], pts_3+[-gr_sep/2,0,gr_sep/2], pts_3+[-gr_sep/2,0,-gr_sep/2],\
                    pts_3+[-gr_sep/2,0,gr_sep/6],pts_3+[-gr_sep/2,0,-gr_sep/6],pts_3+[gr_sep/2,0,gr_sep/6],pts_3+[gr_sep/2,0,-gr_sep/6],\
                    pts_3+[-gr_sep/6,0,-gr_sep/2],pts_3+[-gr_sep/6,0,gr_sep/2],pts_3+[gr_sep/6,0,gr_sep/2],pts_3+[gr_sep/6,0,-gr_sep/2],\
                    pts_3+[-gr_sep/6,0,-gr_sep/6],pts_3+[-gr_sep/6,0,gr_sep/6],pts_3+[gr_sep/6,0,gr_sep/6],pts_3+[gr_sep/6,0,-gr_sep/6]
                lny_3 = np.concatenate(lny_3, axis=0)




                lnx4 = pts4+[0,gr_sep/2,gr_sep/2], pts4+[0,gr_sep/2,-gr_sep/2], pts4+[0,-gr_sep/2,gr_sep/2], pts4+[0,-gr_sep/2,-gr_sep/2],\
                    pts4+[0,-gr_sep/2,0], pts4+[0,-gr_sep/2,-gr_sep/4], pts4+[0,gr_sep/2,-gr_sep/4],\
                    pts4+[0,gr_sep/2,0], pts4+[0,-gr_sep/2,gr_sep/4], pts4+[0,gr_sep/2,gr_sep/4]
                lnx4 = np.concatenate(lnx4, axis=0)

                lnz4 = pts4+[gr_sep/2,gr_sep/2,0], pts4+[gr_sep/2,-gr_sep/2,0], pts4+[-gr_sep/2,gr_sep/2,0], pts4+[-gr_sep/2,-gr_sep/2,0],\
                    pts4+[0,-gr_sep/2,0], pts4+[-gr_sep/4,-gr_sep/2,0], pts4+[gr_sep/4,-gr_sep/2,0],\
                    pts4+[0,gr_sep/2,0], pts4+[-gr_sep/4,gr_sep/2,0], pts4+[gr_sep/4,gr_sep/2,0]
                lnz4 = np.concatenate(lnz4, axis=0)

                lny4 = pts4+[gr_sep/2,0,gr_sep/2], pts4+[gr_sep/2,0,-gr_sep/2], pts4+[-gr_sep/2,0,gr_sep/2], pts4+[-gr_sep/2,0,-gr_sep/2],\
                    pts4+[0,0,0],\
                    pts4+[-gr_sep/4,0,0], pts4+[gr_sep/4,0,0], pts4+[0,0,-gr_sep/4], pts4+[0,0,gr_sep/4],\
                    pts4+[-gr_sep/4,0,-gr_sep/4], pts4+[-gr_sep/4,0,gr_sep/4], pts4+[gr_sep/4,0,-gr_sep/4], pts4+[gr_sep/4,0,gr_sep/4],\
                    pts4+[-gr_sep/2,0,0], pts4+[gr_sep/2,0,0], pts4+[0,0,-gr_sep/2], pts4+[0,0,gr_sep/2],\
                    pts4+[-gr_sep/2,0,gr_sep/4], pts4+[-gr_sep/2,0,-gr_sep/4], pts4+[gr_sep/2,0,-gr_sep/4], pts4+[gr_sep/2,0,gr_sep/4],\
                    pts4+[-gr_sep/4,0,gr_sep/2], pts4+[-gr_sep/4,0,-gr_sep/2], pts4+[gr_sep/4,0,-gr_sep/2], pts4+[gr_sep/4,0,gr_sep/2]
                lny4 = np.concatenate(lny4, axis=0)


                lnx_4 = pts_4+[0,gr_sep/2,gr_sep/2], pts_4+[0,gr_sep/2,-gr_sep/2], pts_4+[0,-gr_sep/2,gr_sep/2], pts_4+[0,-gr_sep/2,-gr_sep/2],\
                    pts_4+[0,-gr_sep/2,0], pts_4+[0,-gr_sep/2,-gr_sep/4], pts_4+[0,gr_sep/2,-gr_sep/4],\
                    pts_4+[0,gr_sep/2,0], pts_4+[0,-gr_sep/2,gr_sep/4], pts_4+[0,gr_sep/2,gr_sep/4]
                lnx_4 = np.concatenate(lnx_4, axis=0)

                lnz_4 = pts_4+[gr_sep/2,gr_sep/2,0], pts_4+[gr_sep/2,-gr_sep/2,0], pts_4+[-gr_sep/2,gr_sep/2,0], pts_4+[-gr_sep/2,-gr_sep/2,0],\
                    pts_4+[0,-gr_sep/2,0], pts_4+[-gr_sep/4,-gr_sep/2,0], pts_4+[gr_sep/4,-gr_sep/2,0],\
                    pts_4+[0,gr_sep/2,0], pts_4+[-gr_sep/4,gr_sep/2,0], pts_4+[gr_sep/4,gr_sep/2,0]
                lnz_4 = np.concatenate(lnz_4, axis=0)

                lny_4 = pts_4+[gr_sep/2,0,gr_sep/2], pts_4+[gr_sep/2,0,-gr_sep/2], pts_4+[-gr_sep/2,0,gr_sep/2], pts_4+[-gr_sep/2,0,-gr_sep/2],\
                    pts_4+[0,0,0],\
                    pts_4+[-gr_sep/4,0,0], pts_4+[gr_sep/4,0,0], pts_4+[0,0,-gr_sep/4], pts_4+[0,0,gr_sep/4],\
                    pts_4+[-gr_sep/4,0,-gr_sep/4], pts_4+[-gr_sep/4,0,gr_sep/4], pts_4+[gr_sep/4,0,-gr_sep/4], pts_4+[gr_sep/4,0,gr_sep/4],\
                    pts_4+[-gr_sep/2,0,0], pts_4+[gr_sep/2,0,0], pts_4+[0,0,-gr_sep/2], pts_4+[0,0,gr_sep/2],\
                    pts_4+[-gr_sep/2,0,gr_sep/4], pts_4+[-gr_sep/2,0,-gr_sep/4], pts_4+[gr_sep/2,0,-gr_sep/4], pts_4+[gr_sep/2,0,gr_sep/4],\
                    pts_4+[-gr_sep/4,0,gr_sep/2], pts_4+[-gr_sep/4,0,-gr_sep/2], pts_4+[gr_sep/4,0,-gr_sep/2], pts_4+[gr_sep/4,0,gr_sep/2]
                lny_4 = np.concatenate(lny_4, axis=0)



                #--------5 mag-------------------


                lnx5 = pts5+[0,gr_sep/2,gr_sep/2], pts5+[0,gr_sep/2,-gr_sep/2], pts5+[0,-gr_sep/2,gr_sep/2], pts5+[0,-gr_sep/2,-gr_sep/2],\
                        pts5+[0,-gr_sep/2,-gr_sep/10], pts5+[0,gr_sep/2,-gr_sep/10], pts5+[0,-gr_sep/2,-(3*gr_sep/10)], pts5+[0,-gr_sep/2,(3*gr_sep/10)],\
                        pts5+[0,-gr_sep/2,gr_sep/10], pts5+[0,gr_sep/2,gr_sep/10], pts5+[0,gr_sep/2,-(3*gr_sep/10)], pts5+[0,gr_sep/2,(3*gr_sep/10)]
                lnx5 = np.concatenate(lnx5, axis=0)

                lnz5 = pts5+[gr_sep/2,gr_sep/2,0], pts5+[gr_sep/2,-gr_sep/2,0], pts5+[-gr_sep/2,gr_sep/2,0], pts5+[-gr_sep/2,-gr_sep/2,0],\
                        pts5+[-gr_sep/10,-gr_sep/2,0], pts5+[gr_sep/10,-gr_sep/2,0], pts5+[-(3*gr_sep/10),-gr_sep/2,0], pts5+[(3*gr_sep/10),-gr_sep/2,0],\
                        pts5+[-gr_sep/10,gr_sep/2,0], pts5+[gr_sep/10,gr_sep/2,0], pts5+[-(3*gr_sep/10),gr_sep/2,0], pts5+[(3*gr_sep/10),gr_sep/2,0]  
                lnz5 = np.concatenate(lnz5, axis=0)

                lny5 = pts5+[gr_sep/2,0,gr_sep/2], pts5+[gr_sep/2,0,-gr_sep/2], pts5+[-gr_sep/2,0,gr_sep/2], pts5+[-gr_sep/2,0,-gr_sep/2],\
                        pts5+[-gr_sep/10,0,gr_sep/10], pts5+[-gr_sep/10,0,-gr_sep/10], pts5+[gr_sep/10,0,gr_sep/10], pts5+[gr_sep/10,0,-gr_sep/10],\
                        pts5+[-(3*gr_sep/10),0,gr_sep/10],pts5+[-(3*gr_sep/10),0,-gr_sep/10], pts5+[-(3*gr_sep/10),0,(3*gr_sep/10)], pts5+[-(3*gr_sep/10),0,-(3*gr_sep/10)],\
                        pts5+[(3*gr_sep/10),0,gr_sep/10],pts5+[(3*gr_sep/10),0,-gr_sep/10], pts5+[(3*gr_sep/10),0,(3*gr_sep/10)], pts5+[(3*gr_sep/10),0,-(3*gr_sep/10)],\
                        pts5+[-(gr_sep/10),0,(3*gr_sep/10)], pts5+[-(gr_sep/10),0,-(3*gr_sep/10)], pts5+[(gr_sep/10),0,(3*gr_sep/10)], pts5+[(gr_sep/10),0,-(3*gr_sep/10)],\
                        pts5+[-(gr_sep/2),0,gr_sep/10],pts5+[-(gr_sep/2),0,-gr_sep/10], pts5+[-(gr_sep/2),0,(3*gr_sep/10)], pts5+[-(gr_sep/2),0,-(3*gr_sep/10)],\
                        pts5+[(gr_sep/2),0,gr_sep/10],pts5+[(gr_sep/2),0,-gr_sep/10], pts5+[(gr_sep/2),0,(3*gr_sep/10)], pts5+[(gr_sep/2),0,-(3*gr_sep/10)],\
                        pts5+[gr_sep/10,0,-(gr_sep/2)],pts5+[-gr_sep/10,0,-(gr_sep/2)], pts5+[(3*gr_sep/10),0,-(gr_sep/2)], pts5+[-(3*gr_sep/10),0,-(gr_sep/2)],\
                        pts5+[gr_sep/10,0,(gr_sep/2)],pts5+[-gr_sep/10,0,(gr_sep/2)], pts5+[(3*gr_sep/10),0,(gr_sep/2)], pts5+[-(3*gr_sep/10),0,(gr_sep/2)]
                lny5 = np.concatenate(lny5, axis=0)




                lnx_5 = pts_5+[0,gr_sep/2,gr_sep/2], pts_5+[0,gr_sep/2,-gr_sep/2], pts_5+[0,-gr_sep/2,gr_sep/2], pts_5+[0,-gr_sep/2,-gr_sep/2],\
                        pts_5+[0,-gr_sep/2,-gr_sep/10], pts_5+[0,gr_sep/2,-gr_sep/10], pts_5+[0,-gr_sep/2,-(3*gr_sep/10)], pts_5+[0,-gr_sep/2,(3*gr_sep/10)],\
                        pts_5+[0,-gr_sep/2,gr_sep/10], pts_5+[0,gr_sep/2,gr_sep/10], pts_5+[0,gr_sep/2,-(3*gr_sep/10)], pts_5+[0,gr_sep/2,(3*gr_sep/10)]
                lnx_5 = np.concatenate(lnx_5, axis=0)

                lnz_5 = pts_5+[gr_sep/2,gr_sep/2,0], pts_5+[gr_sep/2,-gr_sep/2,0], pts_5+[-gr_sep/2,gr_sep/2,0], pts_5+[-gr_sep/2,-gr_sep/2,0],\
                        pts_5+[-gr_sep/10,-gr_sep/2,0], pts_5+[gr_sep/10,-gr_sep/2,0], pts_5+[-(3*gr_sep/10),-gr_sep/2,0], pts_5+[(3*gr_sep/10),-gr_sep/2,0],\
                        pts_5+[-gr_sep/10,gr_sep/2,0], pts_5+[gr_sep/10,gr_sep/2,0], pts_5+[-(3*gr_sep/10),gr_sep/2,0], pts_5+[(3*gr_sep/10),gr_sep/2,0]  
                lnz_5 = np.concatenate(lnz_5, axis=0)

                lny_5 = pts_5+[gr_sep/2,0,gr_sep/2], pts_5+[gr_sep/2,0,-gr_sep/2], pts_5+[-gr_sep/2,0,gr_sep/2], pts_5+[-gr_sep/2,0,-gr_sep/2],\
                        pts_5+[-gr_sep/10,0,gr_sep/10], pts_5+[-gr_sep/10,0,-gr_sep/10], pts_5+[gr_sep/10,0,gr_sep/10], pts_5+[gr_sep/10,0,-gr_sep/10],\
                        pts_5+[-(3*gr_sep/10),0,gr_sep/10],pts_5+[-(3*gr_sep/10),0,-gr_sep/10], pts_5+[-(3*gr_sep/10),0,(3*gr_sep/10)], pts_5+[-(3*gr_sep/10),0,-(3*gr_sep/10)],\
                        pts_5+[(3*gr_sep/10),0,gr_sep/10],pts_5+[(3*gr_sep/10),0,-gr_sep/10], pts_5+[(3*gr_sep/10),0,(3*gr_sep/10)], pts_5+[(3*gr_sep/10),0,-(3*gr_sep/10)],\
                        pts_5+[-(gr_sep/10),0,(3*gr_sep/10)], pts_5+[-(gr_sep/10),0,-(3*gr_sep/10)], pts_5+[(gr_sep/10),0,(3*gr_sep/10)], pts_5+[(gr_sep/10),0,-(3*gr_sep/10)],\
                        pts_5+[-(gr_sep/2),0,gr_sep/10],pts_5+[-(gr_sep/2),0,-gr_sep/10], pts_5+[-(gr_sep/2),0,(3*gr_sep/10)], pts_5+[-(gr_sep/2),0,-(3*gr_sep/10)],\
                        pts_5+[(gr_sep/2),0,gr_sep/10],pts_5+[(gr_sep/2),0,-gr_sep/10], pts_5+[(gr_sep/2),0,(3*gr_sep/10)], pts_5+[(gr_sep/2),0,-(3*gr_sep/10)],\
                        pts_5+[gr_sep/10,0,-(gr_sep/2)],pts_5+[-gr_sep/10,0,-(gr_sep/2)], pts_5+[(3*gr_sep/10),0,-(gr_sep/2)], pts_5+[-(3*gr_sep/10),0,-(gr_sep/2)],\
                        pts_5+[gr_sep/10,0,(gr_sep/2)],pts_5+[-gr_sep/10,0,(gr_sep/2)], pts_5+[(3*gr_sep/10),0,(gr_sep/2)], pts_5+[-(3*gr_sep/10),0,(gr_sep/2)]
                lny_5 = np.concatenate(lny_5, axis=0)



            #===============================PLOTTING========================================================================





            ###----coord for line centres around each point-----------

            ###-------------------------------------------------------


            def plotter(lnx, lny, lnz, clr):



                line1 = tvtk.LineSource(point1=(-gr_sep/2,0,0), point2=(gr_sep/2,0,0))
                line2 = tvtk.LineSource(point1=(0,-gr_sep/2,0), point2=(0,gr_sep/2,0))
                line3 = tvtk.LineSource(point1=(0,0,-gr_sep/2), point2=(0,0,gr_sep/2))

                pd1 = tvtk.PolyData(points=lnx)
                pd2 = tvtk.PolyData(points=lny)
                pd3 = tvtk.PolyData(points=lnz)

                g1 = tvtk.Glyph3D(scale_mode='data_scaling_off', vector_mode = 'use_vector')
                g2 = tvtk.Glyph3D(scale_mode='data_scaling_off', vector_mode = 'use_vector')
                g3 = tvtk.Glyph3D(scale_mode='data_scaling_off', vector_mode = 'use_vector')



                configure_input_data(g1, pd1)
                configure_input_data(g2, pd2)
                configure_input_data(g3, pd3)

                if direction=='z':
                    configure_source_data(g1, line1.output)
                    line1.update()
                    g1.update()
                    configure_source_data(g2, line2.output)
                    line2.update()
                    g2.update()
                    configure_source_data(g3, line3.output)
                    line3.update()
                    g3.update()


                if direction=='x':
                    configure_source_data(g1, line1.output)
                    line1.update()
                    g1.update()
                    configure_source_data(g2, line2.output)
                    line2.update()
                    g2.update()
                    configure_source_data(g3, line3.output)
                    line3.update()
                    g3.update()

                if direction=='y':
                    configure_source_data(g1, line1.output)
                    line1.update()
                    g1.update()
                    configure_source_data(g2, line2.output)
                    line2.update()
                    g2.update()
                    configure_source_data(g3, line3.output)
                    line3.update()
                    g3.update()

                m1 = tvtk.PolyDataMapper()
                m2 = tvtk.PolyDataMapper()
                m3 = tvtk.PolyDataMapper()

                pc1 = tvtk.Property(opacity=0.9, color=clr)


                configure_input_data(m1, g1.output)
                configure_input_data(m2, g2.output)
                configure_input_data(m3, g3.output)

                a1 = tvtk.Actor(mapper=m1, property=pc1)
                a2 = tvtk.Actor(mapper=m2, property=pc1)
                a3 = tvtk.Actor(mapper=m3, property=pc1)

                scn.add_actor(a1)
                scn.add_actor(a2)
                scn.add_actor(a3)


            plotter(lnx1, lny1, lnz1, (0.9,0,0))
            plotter(lnx_1, lny_1, lnz_1, (0,0,0.9))

            plotter(lnx2, lny2, lnz2, (0.9,0,0))
            plotter(lnx_2, lny_2, lnz_2, (0,0,0.9))

            plotter(lnx3, lny3, lnz3, (0.9,0,0))
            plotter(lnx_3, lny_3, lnz_3, (0,0,0.9))

            plotter(lnx4, lny4, lnz4, (0.9,0,0))
            plotter(lnx_4, lny_4, lnz_4, (0,0,0.9))

            plotter(lnx5, lny5, lnz5, (0.9,0,0))
            plotter(lnx_5, lny_5, lnz_5, (0,0,0.9))
            

            xmin = float((np.min(self.xg)) - gr_sep)
            ymin = float((np.min(self.yg)) - gr_sep)
            zmin = float((np.min(self.zg)) - gr_sep)
            xmax = float((np.max(self.xg)) + gr_sep)
            ymax = float((np.max(self.yg)) + gr_sep)
            zmax = float((np.max(self.zg)) + gr_sep)

            if direction=='z':
                xlab='dx'
                ylab='dy'
                zlab='dx/\dy'
            else:
                pass

            if direction=='x':
                xlab='dy/\dz '
                ylab='dy'
                zlab='dz'
            else:
                pass


            if direction=='y':
                xlab='dx'
                ylab='dx/\dz'
                zlab='dz'
            else:
                pass

            scn.mlab.points3d(([xmin,ymin,zmin],[xmin,ymin,zmax],[xmin,ymax,zmin],[xmin,ymax,zmax],[xmax,ymin,zmin],[xmax,ymin,zmax],[xmax,ymax,zmin],[xmax,ymax,zmax]),opacity=0.0)
            scn.mlab.points3d(pts_nan[:,0],
                                pts_nan[:,1],
                                pts_nan[:,2], color = (1,0,0),scale_factor=gr_sep, resolution=36)
            scn.mlab.axes(extent = [xmin,xmax,ymin,ymax,zmin,zmax], nb_labels = 5, line_width=3.0, color = (0,0,0), xlabel=xlab, ylabel=ylab, zlabel=zlab)

            xfoc = (xmin+xmax)/2
            yfoc = (ymin+ymax)/2
            zfoc = (zmin+zmax)/2


            scn.mlab.view(focalpoint=[xfoc,yfoc,zfoc])
            scn.mlab.show()




        if Fz is not None and np.count_nonzero(Fz) != 0:
            form_2(Fz,'z')
        
        if Fx is not None and np.count_nonzero(Fx) != 0:
            form_2(Fx,'x')

        if Fy is not None and np.count_nonzero(Fy) != 0:
            form_2(Fy,'y')
       



   


class form_3_3d():

    '''
    defines a 3-form object and returns it to user
    Takes 4 arguments basic, these are the 3 grids in 3D, which muse be square
    and of equal sizes. Then 1 argument for the dx^dy^dz component
    based on the same grids. Also takes in equation which is needed for some
    operaions.

    Methods: (applied to the object)

    .give_eqn(equation_str_x, equation_str_y, equation_str_z) - provides equations ; 
    .return_string() - returns equations as strings of a 2-form. ;
    .log_scaling() - scales the field logarithmycally ; 
    .zoom(magnitude, target, point_density) - creates a zoomed in field object ; 
    .plot() - applied to the 3-form object, gives out a plot of non-zero components ;
    .hodge() - applies hodge star operator to the 3-form ; 
    .wedge() - wedges the 3-form with another p-form object ;
    .interior_d() - takes the interior derivative of the form

    '''

    def __init__(self, xg, yg, zg, form_3, form_3_eqn=None):
        self.xg = xg
        self.yg = yg
        self.zg = zg
        self.form_3 = form_3
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
        
        if form_3_eqn is not None:
            self.form_3_str = str(simplify(form_3_eqn))  # user must change to access some methods
        else:
            self.form_3_str = None



    def give_eqn(self, equation_str):
        '''
        
        Allows user to supply equation to instance, if not initially done so
        
        Parameters:
        ------------
        equation_str - str - equation of the supplied numerical 0-form
                        It must be in terms of x, y, z.
                        Has to be given, for some methods to be calculatable.
        
        Returns: None
        
        '''
        self.form_3_str = str(simplify(equation_str))
        
        # update the numerical values to always match
        string = self.form_3_str + ''
        
        # Check if the equations provided contain x and y terms
        # and format them to be evaluated
        if string.find('x') & string.find('y') & string.find('z') == -1:
            string = '(' + str(string) + ')* np.ones(np.shape(xg))'
        else:
            string = string.replace('x', '(self.xg)')
            string = string.replace('y', '(self.yg)')
            string = string.replace('z', '(self.zg)')
        
        # evaluate the string equation
        self.form_3 = eval(string)



    def return_string(self):
        '''
        Takes in no arguments, returns the unformatted string back to user
        This is done in case user wants to access strings
        that got here not by input but by ext. alg.
        '''
        return self.form_3_str



    def log_scaling(self):
        '''
        changes bool for logscaling
        Default = False
        changes to the other option each time it is called
        '''
        self.logarithmic_scale_bool = not self.logarithmic_scale_bool



    def zoom(self, mag, target, dpd):
        
        """
        Redefines a new, zoomed in grid on which a new 2-form is calculated.
        Takes in three arguments:

        mag - magnitude of zoom ; 
        target - target point [x, y, z] ;
        dpd - points amount in the new grid ; 

        Returns a new, zoomed in 3-form 
        """


        if self.form_3_str == None:
            # ERROR
            raise TypeError('No equation provided, see \'give_eqn\' method')
        else:
             # Zoom must be one or greater
            if mag < 1:
                raise ValueError('Mag must be greater than one')
            else:

                # Target coordinates
                x_m = target[0]
                y_m = target[1]
                z_m = target[2]
                
                # Get the size of the original VF
                Lx = 0.5*(self.xg[-1, -1, -1] - self.xg[0, 0, 0])
                Ly = 0.5*(self.yg[-1, -1, -1] - self.yg[0, 0, 0])
                Lz = 0.5*(self.zg[-1, -1, -1] - self.zg[0, 0, 0])

                
                # Zoom axis range
                d_range_x = Lx/mag
                d_range_y = Ly/mag
                d_range_z = Lz/mag
                
                # Set up zoom window grids
                dx = np.linspace(-d_range_x + x_m, d_range_x + x_m, dpd)
                dy = np.linspace(-d_range_y + y_m, d_range_y + y_m, dpd)
                dz = np.linspace(-d_range_z + z_m, d_range_z + z_m, dpd)
                dxg, dyg, dzg = np.meshgrid(dx, dy, dz)
                
                # Create variables for the user provided equation strings
                f3_str = self.form_3_str

                
                # Check if the equations provided contain x and y terms
                if f3_str.find('x') & f3_str.find('y') & f3_str.find('z')== -1:
                    f3_str = '(' + str(f3_str) + ')* np.ones(np.shape(dxg))'
                else:
                    f3_str = f3_str.replace('x', 'dxg')
                    f3_str = f3_str.replace('y', 'dyg')
                    f3_str = f3_str.replace('z', 'dzg')
            

                    
                # Generate arrays for the components of the zoom field
                f3_zoom = eval(f3_str)

                
                # crate the zoomed in form
                zoom_form = form_3_3d(dxg, dyg, dzg, f3_zoom, self.form_3_str)

                
                return zoom_form
 


    def plot(self, scn):

        """
        Plots the 3-form. Takes in no arguments, returns nothing.
        """

        xg = self.xg
        yg = self.yg
        zg = self.zg
        form_3 = self.form_3

        gr_sep = abs(self.xg[1,1,1]-self.xg[0,0,0])

        pts = np.vstack(list(zip(xg.ravel(), yg.ravel(), zg.ravel())))

        if self.logarithmic_scale_bool:
            mag1 = np.abs(form_3) + 1
            form_3_norm = form_3/(mag1)
            logmag = np.log10(mag1)
            form_3 = form_3_norm*logmag
        else:
            pass

        mag_lst = np.vstack(list(zip(form_3.ravel())))
        mag_lst[np.isinf(mag_lst)] = np.nan
        Idx_nan = np.argwhere(np.isnan(mag_lst))

        Idx_nan = Idx_nan[:,0]

        pts_nan = pts[Idx_nan]

        pts = np.delete(pts, [Idx_nan], axis=0)
        mag_lst = np.delete(mag_lst, [Idx_nan], axis=0)

        f_max = np.nanmax(abs(self.form_3))

        sep = (f_max)/5

        

        ### Find points where field has different magnitudes

        Idx1 = np.argwhere(np.all(mag_lst>(0),axis=1) & np.all(mag_lst<=(0+sep),axis=1))
        Idx2 = np.argwhere(np.all(mag_lst>(0+sep),axis=1) & np.all(mag_lst<=(0+2*sep),axis=1))
        Idx3 = np.argwhere(np.all(mag_lst>(0+2*sep),axis=1) & np.all(mag_lst<=(0+3*sep),axis=1))
        Idx4 = np.argwhere(np.all(mag_lst>(0+3*sep),axis=1) & np.all(mag_lst<=(0+4*sep),axis=1))
        Idx5 = np.argwhere(np.all(mag_lst>(0+4*sep),axis=1))

        Idx_1 = np.argwhere(np.all(mag_lst<(0),axis=1) & np.all(mag_lst>=(0-sep),axis=1))
        Idx_2 = np.argwhere(np.all(mag_lst<(0-sep),axis=1) & np.all(mag_lst>=(0-2*sep),axis=1))
        Idx_3 = np.argwhere(np.all(mag_lst<(0-2*sep),axis=1) & np.all(mag_lst>=(0-3*sep),axis=1))
        Idx_4 = np.argwhere(np.all(mag_lst<(0-3*sep),axis=1) & np.all(mag_lst>=(0-4*sep),axis=1))
        Idx_5 = np.argwhere(np.all(mag_lst<(0-4*sep),axis=1))



        pts1 = pts[Idx1[:,0]]
        pts_1 = pts[Idx_1[:,0]]
        pts2 = pts[Idx2[:,0]]
        pts_2 = pts[Idx_2[:,0]]
        pts3 = pts[Idx3[:,0]]
        pts_3 = pts[Idx_3[:,0]]
        pts4 = pts[Idx4[:,0]]
        pts_4 = pts[Idx_4[:,0]]
        pts5 = pts[Idx5[:,0]]
        pts_5 = pts[Idx_5[:,0]]


        #----points for lines------------------


        
        #------1 mag-----------

        lnx1 = pts1+[0,gr_sep/2,gr_sep/2], pts1+[0,gr_sep/2,-gr_sep/2], pts1+[0,-gr_sep/2,gr_sep/2], pts1+[0,-gr_sep/2,-gr_sep/2]
        lnx1 = np.concatenate(lnx1, axis=0)

        lny1 = pts1+[gr_sep/2,0,gr_sep/2], pts1+[gr_sep/2,0,-gr_sep/2], pts1+[-gr_sep/2,0,gr_sep/2], pts1+[-gr_sep/2,0,-gr_sep/2]
        lny1 = np.concatenate(lny1, axis=0)

        lnz1 = pts1+[gr_sep/2,gr_sep/2,0], pts1+[gr_sep/2,-gr_sep/2,0], pts1+[-gr_sep/2,gr_sep/2,0], pts1+[-gr_sep/2,-gr_sep/2,0]
        lnz1 = np.concatenate(lnz1, axis=0)


        lnx_1 = pts_1+[0,gr_sep/2,gr_sep/2], pts_1+[0,gr_sep/2,-gr_sep/2], pts_1+[0,-gr_sep/2,gr_sep/2], pts_1+[0,-gr_sep/2,-gr_sep/2]
        lnx_1 = np.concatenate(lnx_1, axis=0)

        lny_1 = pts_1+[gr_sep/2,0,gr_sep/2], pts_1+[gr_sep/2,0,-gr_sep/2], pts_1+[-gr_sep/2,0,gr_sep/2], pts_1+[-gr_sep/2,0,-gr_sep/2]
        lny_1 = np.concatenate(lny_1, axis=0)

        lnz_1 = pts_1+[gr_sep/2,gr_sep/2,0], pts_1+[gr_sep/2,-gr_sep/2,0], pts_1+[-gr_sep/2,gr_sep/2,0], pts_1+[-gr_sep/2,-gr_sep/2,0]
        lnz_1 = np.concatenate(lnz_1, axis=0)


        #--------2 mag--------------



        lnx2 = pts2+[0,gr_sep/2,gr_sep/2], pts2+[0,-gr_sep/2,gr_sep/2], pts2+[0,gr_sep/2,-gr_sep/2], pts2+[0,-gr_sep/2,-gr_sep/2],\
                pts2+[0,0,gr_sep/2], pts2+[0,0,-gr_sep/2], pts2+[0,gr_sep/2,0], pts2+[0,-gr_sep/2,0],\
                pts2+[0,0,0]
        lnx2 = np.concatenate(lnx2, axis=0)

        lny2 = pts2+[gr_sep/2,0,gr_sep/2], pts2+[gr_sep/2,0,-gr_sep/2], pts2+[-gr_sep/2,0,gr_sep/2], pts2+[-gr_sep/2,0,-gr_sep/2],\
                pts2+[gr_sep/2,0,0], pts2+[-gr_sep/2,0,0], pts2+[0,0,gr_sep/2], pts2+[0,0,-gr_sep/2],\
                pts2+[0,0,0]
        lny2 = np.concatenate(lny2, axis=0)

        lnz2 = pts2+[gr_sep/2,gr_sep/2,0], pts2+[gr_sep/2,-gr_sep/2,0], pts2+[-gr_sep/2,gr_sep/2,0], pts2+[-gr_sep/2,-gr_sep/2,0],\
                pts2+[gr_sep/2,0,0], pts2+[-gr_sep/2,0,0], pts2+[0,gr_sep/2,0], pts2+[0,-gr_sep/2,0],\
                pts2+[0,0,0]
        lnz2 = np.concatenate(lnz2, axis=0)


        lnx_2 = pts_2+[0,gr_sep/2,gr_sep/2], pts_2+[0,-gr_sep/2,gr_sep/2], pts_2+[0,gr_sep/2,-gr_sep/2], pts_2+[0,-gr_sep/2,-gr_sep/2],\
                pts_2+[0,0,gr_sep/2], pts_2+[0,0,-gr_sep/2], pts_2+[0,gr_sep/2,0], pts_2+[0,-gr_sep/2,0],\
                pts_2+[0,0,0]
        lnx_2 = np.concatenate(lnx_2, axis=0)

        lny_2 = pts_2+[gr_sep/2,0,gr_sep/2], pts_2+[gr_sep/2,0,-gr_sep/2], pts_2+[-gr_sep/2,0,gr_sep/2], pts_2+[-gr_sep/2,0,-gr_sep/2],\
                pts_2+[gr_sep/2,0,0], pts_2+[-gr_sep/2,0,0], pts_2+[0,0,gr_sep/2], pts_2+[0,0,-gr_sep/2],\
                pts_2+[0,0,0]
        lny_2 = np.concatenate(lny_2, axis=0)

        lnz_2 = pts_2+[gr_sep/2,gr_sep/2,0], pts_2+[gr_sep/2,-gr_sep/2,0], pts_2+[-gr_sep/2,gr_sep/2,0], pts_2+[-gr_sep/2,-gr_sep/2,0],\
                pts_2+[gr_sep/2,0,0], pts_2+[-gr_sep/2,0,0], pts_2+[0,gr_sep/2,0], pts_2+[0,-gr_sep/2,0],\
                pts_2+[0,0,0]
        lnz_2 = np.concatenate(lnz_2, axis=0)


        #----------3 mag----------------


        lnx3 = pts3+[0,gr_sep/2,gr_sep/2], pts3+[0,-gr_sep/2,gr_sep/2], pts3+[0,gr_sep/2,-gr_sep/2], pts3+[0,-gr_sep/2,-gr_sep/2],\
                pts3+[0,gr_sep/6,-gr_sep/2],pts3+[0,-gr_sep/6,-gr_sep/2],pts3+[0,gr_sep/6,gr_sep/2],pts3+[0,-gr_sep/6,gr_sep/2],\
                pts3+[0,-gr_sep/2,-gr_sep/6],pts3+[0,gr_sep/2,-gr_sep/6],pts3+[0,gr_sep/2,gr_sep/6],pts3+[0,-gr_sep/2,gr_sep/6],\
                pts3+[0,-gr_sep/6,-gr_sep/6],pts3+[0,gr_sep/6,-gr_sep/6],pts3+[0,gr_sep/6,gr_sep/6],pts3+[0,-gr_sep/6,gr_sep/6]
        lnx3 = np.concatenate(lnx3, axis=0)

        lny3 = pts3+[gr_sep/2,0,gr_sep/2], pts3+[gr_sep/2,0,-gr_sep/2], pts3+[-gr_sep/2,0,gr_sep/2], pts3+[-gr_sep/2,0,-gr_sep/2],\
                pts3+[-gr_sep/2,0,gr_sep/6],pts3+[-gr_sep/2,0,-gr_sep/6],pts3+[gr_sep/2,0,gr_sep/6],pts3+[gr_sep/2,0,-gr_sep/6],\
                pts3+[-gr_sep/6,0,-gr_sep/2],pts3+[-gr_sep/6,0,gr_sep/2],pts3+[gr_sep/6,0,gr_sep/2],pts3+[gr_sep/6,0,-gr_sep/2],\
                pts3+[-gr_sep/6,0,-gr_sep/6],pts3+[-gr_sep/6,0,gr_sep/6],pts3+[gr_sep/6,0,gr_sep/6],pts3+[gr_sep/6,0,-gr_sep/6]
        lny3 = np.concatenate(lny3, axis=0)

        lnz3 = pts3+[gr_sep/2,gr_sep/2,0], pts3+[gr_sep/2,-gr_sep/2,0], pts3+[-gr_sep/2,gr_sep/2,0], pts3+[-gr_sep/2,-gr_sep/2,0],\
                pts3+[-gr_sep/2,gr_sep/6,0],pts3+[-gr_sep/2,-gr_sep/6,0],pts3+[gr_sep/2,gr_sep/6,0],pts3+[gr_sep/2,-gr_sep/6,0],\
                pts3+[-gr_sep/6,-gr_sep/2,0],pts3+[-gr_sep/6,gr_sep/2,0],pts3+[gr_sep/6,gr_sep/2,0],pts3+[gr_sep/6,-gr_sep/2,0],\
                pts3+[-gr_sep/6,-gr_sep/6,0],pts3+[-gr_sep/6,gr_sep/6,0],pts3+[gr_sep/6,gr_sep/6,0],pts3+[gr_sep/6,-gr_sep/6,0]
        lnz3 = np.concatenate(lnz3, axis=0)


        lnx_3 = pts_3+[0,gr_sep/2,gr_sep/2], pts_3+[0,-gr_sep/2,gr_sep/2], pts_3+[0,gr_sep/2,-gr_sep/2], pts_3+[0,-gr_sep/2,-gr_sep/2],\
                pts_3+[0,gr_sep/6,-gr_sep/2],pts_3+[0,-gr_sep/6,-gr_sep/2],pts_3+[0,gr_sep/6,gr_sep/2],pts_3+[0,-gr_sep/6,gr_sep/2],\
                pts_3+[0,-gr_sep/2,-gr_sep/6],pts_3+[0,gr_sep/2,-gr_sep/6],pts_3+[0,gr_sep/2,gr_sep/6],pts_3+[0,-gr_sep/2,gr_sep/6],\
                pts_3+[0,-gr_sep/6,-gr_sep/6],pts_3+[0,gr_sep/6,-gr_sep/6],pts_3+[0,gr_sep/6,gr_sep/6],pts_3+[0,-gr_sep/6,gr_sep/6]
        lnx_3 = np.concatenate(lnx_3, axis=0)

        lny_3 = pts_3+[gr_sep/2,0,gr_sep/2], pts_3+[gr_sep/2,0,-gr_sep/2], pts_3+[-gr_sep/2,0,gr_sep/2], pts_3+[-gr_sep/2,0,-gr_sep/2],\
                pts_3+[-gr_sep/2,0,gr_sep/6],pts_3+[-gr_sep/2,0,-gr_sep/6],pts_3+[gr_sep/2,0,gr_sep/6],pts_3+[gr_sep/2,0,-gr_sep/6],\
                pts_3+[-gr_sep/6,0,-gr_sep/2],pts_3+[-gr_sep/6,0,gr_sep/2],pts_3+[gr_sep/6,0,gr_sep/2],pts_3+[gr_sep/6,0,-gr_sep/2],\
                pts_3+[-gr_sep/6,0,-gr_sep/6],pts_3+[-gr_sep/6,0,gr_sep/6],pts_3+[gr_sep/6,0,gr_sep/6],pts_3+[gr_sep/6,0,-gr_sep/6]
        lny_3 = np.concatenate(lny_3, axis=0)

        lnz_3 = pts_3+[gr_sep/2,gr_sep/2,0], pts_3+[gr_sep/2,-gr_sep/2,0], pts_3+[-gr_sep/2,gr_sep/2,0], pts_3+[-gr_sep/2,-gr_sep/2,0],\
                pts_3+[-gr_sep/2,gr_sep/6,0],pts_3+[-gr_sep/2,-gr_sep/6,0],pts_3+[gr_sep/2,gr_sep/6,0],pts_3+[gr_sep/2,-gr_sep/6,0],\
                pts_3+[-gr_sep/6,-gr_sep/2,0],pts_3+[-gr_sep/6,gr_sep/2,0],pts_3+[gr_sep/6,gr_sep/2,0],pts_3+[gr_sep/6,-gr_sep/2,0],\
                pts_3+[-gr_sep/6,-gr_sep/6,0],pts_3+[-gr_sep/6,gr_sep/6,0],pts_3+[gr_sep/6,gr_sep/6,0],pts_3+[gr_sep/6,-gr_sep/6,0]
        lnz_3 = np.concatenate(lnz_3, axis=0)



        lnx4 = pts4+[0,gr_sep/2,gr_sep/2], pts4+[0,-gr_sep/2,gr_sep/2], pts4+[0,gr_sep/2,-gr_sep/2], pts4+[0,-gr_sep/2,-gr_sep/2],\
                pts4+[0,0,0],\
                pts4+[0,0,-gr_sep/4], pts4+[0,0,gr_sep/4], pts4+[0,-gr_sep/4,0], pts4+[0,gr_sep/4,0],\
                pts4+[0,-gr_sep/4,-gr_sep/4], pts4+[0,gr_sep/4,-gr_sep/4], pts4+[0,-gr_sep/4,gr_sep/4], pts4+[0,gr_sep/4,gr_sep/4],\
                pts4+[0,0,-gr_sep/2], pts4+[0,0,gr_sep/2], pts4+[0,-gr_sep/2,0], pts4+[0,gr_sep/2,0],\
                pts4+[0,gr_sep/4,-gr_sep/2], pts4+[0,-gr_sep/4,-gr_sep/2], pts4+[0,-gr_sep/4,gr_sep/2], pts4+[0,gr_sep/4,gr_sep/2],\
                pts4+[0,gr_sep/2,-gr_sep/4], pts4+[0,-gr_sep/2,-gr_sep/4], pts4+[0,-gr_sep/2,gr_sep/4], pts4+[0,gr_sep/2,gr_sep/4]
        lnx4 = np.concatenate(lnx4, axis=0)

        lny4 = pts4+[gr_sep/2,0,gr_sep/2], pts4+[gr_sep/2,0,-gr_sep/2], pts4+[-gr_sep/2,0,gr_sep/2], pts4+[-gr_sep/2,0,-gr_sep/2],\
                pts4+[0,0,0],\
                pts4+[-gr_sep/4,0,0], pts4+[gr_sep/4,0,0], pts4+[0,0,-gr_sep/4], pts4+[0,0,gr_sep/4],\
                pts4+[-gr_sep/4,0,-gr_sep/4], pts4+[-gr_sep/4,0,gr_sep/4], pts4+[gr_sep/4,0,-gr_sep/4], pts4+[gr_sep/4,0,gr_sep/4],\
                pts4+[-gr_sep/2,0,0], pts4+[gr_sep/2,0,0], pts4+[0,0,-gr_sep/2], pts4+[0,0,gr_sep/2],\
                pts4+[-gr_sep/2,0,gr_sep/4], pts4+[-gr_sep/2,0,-gr_sep/4], pts4+[gr_sep/2,0,-gr_sep/4], pts4+[gr_sep/2,0,gr_sep/4],\
                pts4+[-gr_sep/4,0,gr_sep/2], pts4+[-gr_sep/4,0,-gr_sep/2], pts4+[gr_sep/4,0,-gr_sep/2], pts4+[gr_sep/4,0,gr_sep/2]
        lny4 = np.concatenate(lny4, axis=0)

        lnz4 = pts4+[gr_sep/2,gr_sep/2,0], pts4+[gr_sep/2,-gr_sep/2,0], pts4+[-gr_sep/2,gr_sep/2,0], pts4+[-gr_sep/2,-gr_sep/2,0],\
                pts4+[0,0,0],\
                pts4+[-gr_sep/4,0,0], pts4+[gr_sep/4,0,0], pts4+[0,-gr_sep/4,0], pts4+[0,gr_sep/4,0],\
                pts4+[-gr_sep/4,-gr_sep/4,0], pts4+[-gr_sep/4,gr_sep/4,0], pts4+[gr_sep/4,-gr_sep/4,0], pts4+[gr_sep/4,gr_sep/4,0],\
                pts4+[-gr_sep/2,0,0], pts4+[gr_sep/2,0,0], pts4+[0,-gr_sep/2,0], pts4+[0,gr_sep/2,0],\
                pts4+[-gr_sep/2,gr_sep/4,0], pts4+[-gr_sep/2,-gr_sep/4,0], pts4+[gr_sep/2,-gr_sep/4,0], pts4+[gr_sep/2,gr_sep/4,0],\
                pts4+[-gr_sep/4,gr_sep/2,0], pts4+[-gr_sep/4,-gr_sep/2,0], pts4+[gr_sep/4,-gr_sep/2,0], pts4+[gr_sep/4,gr_sep/2,0]
        lnz4 = np.concatenate(lnz4, axis=0)


        lnx_4 = pts_4+[0,gr_sep/2,gr_sep/2], pts_4+[0,-gr_sep/2,gr_sep/2], pts_4+[0,gr_sep/2,-gr_sep/2], pts_4+[0,-gr_sep/2,-gr_sep/2],\
                pts_4+[0,0,0],\
                pts_4+[0,0,-gr_sep/4], pts_4+[0,0,gr_sep/4], pts_4+[0,-gr_sep/4,0], pts_4+[0,gr_sep/4,0],\
                pts_4+[0,-gr_sep/4,-gr_sep/4], pts_4+[0,gr_sep/4,-gr_sep/4], pts_4+[0,-gr_sep/4,gr_sep/4], pts_4+[0,gr_sep/4,gr_sep/4],\
                pts_4+[0,0,-gr_sep/2], pts_4+[0,0,gr_sep/2], pts_4+[0,-gr_sep/2,0], pts_4+[0,gr_sep/2,0],\
                pts_4+[0,gr_sep/4,-gr_sep/2], pts_4+[0,-gr_sep/4,-gr_sep/2], pts_4+[0,-gr_sep/4,gr_sep/2], pts_4+[0,gr_sep/4,gr_sep/2],\
                pts_4+[0,gr_sep/2,-gr_sep/4], pts_4+[0,-gr_sep/2,-gr_sep/4], pts_4+[0,-gr_sep/2,gr_sep/4], pts_4+[0,gr_sep/2,gr_sep/4]
        lnx_4 = np.concatenate(lnx_4, axis=0)

        lny_4 = pts_4+[gr_sep/2,0,gr_sep/2], pts_4+[gr_sep/2,0,-gr_sep/2], pts_4+[-gr_sep/2,0,gr_sep/2], pts_4+[-gr_sep/2,0,-gr_sep/2],\
                pts_4+[0,0,0],\
                pts_4+[-gr_sep/4,0,0], pts_4+[gr_sep/4,0,0], pts_4+[0,0,-gr_sep/4], pts_4+[0,0,gr_sep/4],\
                pts_4+[-gr_sep/4,0,-gr_sep/4], pts_4+[-gr_sep/4,0,gr_sep/4], pts_4+[gr_sep/4,0,-gr_sep/4], pts_4+[gr_sep/4,0,gr_sep/4],\
                pts_4+[-gr_sep/2,0,0], pts_4+[gr_sep/2,0,0], pts_4+[0,0,-gr_sep/2], pts_4+[0,0,gr_sep/2],\
                pts_4+[-gr_sep/2,0,gr_sep/4], pts_4+[-gr_sep/2,0,-gr_sep/4], pts_4+[gr_sep/2,0,-gr_sep/4], pts_4+[gr_sep/2,0,gr_sep/4],\
                pts_4+[-gr_sep/4,0,gr_sep/2], pts_4+[-gr_sep/4,0,-gr_sep/2], pts_4+[gr_sep/4,0,-gr_sep/2], pts_4+[gr_sep/4,0,gr_sep/2]
        lny_4 = np.concatenate(lny_4, axis=0)

        lnz_4 = pts_4+[gr_sep/2,gr_sep/2,0], pts_4+[gr_sep/2,-gr_sep/2,0], pts_4+[-gr_sep/2,gr_sep/2,0], pts_4+[-gr_sep/2,-gr_sep/2,0],\
                pts_4+[0,0,0],\
                pts_4+[-gr_sep/4,0,0], pts_4+[gr_sep/4,0,0], pts_4+[0,-gr_sep/4,0], pts_4+[0,gr_sep/4,0],\
                pts_4+[-gr_sep/4,-gr_sep/4,0], pts_4+[-gr_sep/4,gr_sep/4,0], pts_4+[gr_sep/4,-gr_sep/4,0], pts_4+[gr_sep/4,gr_sep/4,0],\
                pts_4+[-gr_sep/2,0,0], pts_4+[gr_sep/2,0,0], pts_4+[0,-gr_sep/2,0], pts_4+[0,gr_sep/2,0],\
                pts_4+[-gr_sep/2,gr_sep/4,0], pts_4+[-gr_sep/2,-gr_sep/4,0], pts_4+[gr_sep/2,-gr_sep/4,0], pts_4+[gr_sep/2,gr_sep/4,0],\
                pts_4+[-gr_sep/4,gr_sep/2,0], pts_4+[-gr_sep/4,-gr_sep/2,0], pts_4+[gr_sep/4,-gr_sep/2,0], pts_4+[gr_sep/4,gr_sep/2,0]
        lnz_4 = np.concatenate(lnz_4, axis=0)



        #--------5 mag-------------------


        lnx5 = pts5+[0,gr_sep/2,gr_sep/2], pts5+[0,gr_sep/2,-gr_sep/2], pts5+[0,-gr_sep/2,gr_sep/2], pts5+[0,-gr_sep/2,-gr_sep/2],\
                pts5+[0,-gr_sep/10,gr_sep/10], pts5+[0,-gr_sep/10,-gr_sep/10], pts5+[0,gr_sep/10,gr_sep/10], pts5+[0,gr_sep/10,-gr_sep/10],\
                pts5+[0,-(3*gr_sep/10),gr_sep/10],pts5+[0,-(3*gr_sep/10),-gr_sep/10], pts5+[0,-(3*gr_sep/10),(3*gr_sep/10)], pts5+[0,-(3*gr_sep/10),-(3*gr_sep/10)],\
                pts5+[0,(3*gr_sep/10),gr_sep/10],pts5+[0,(3*gr_sep/10),-gr_sep/10], pts5+[0,(3*gr_sep/10),(3*gr_sep/10)], pts5+[0,(3*gr_sep/10),-(3*gr_sep/10)],\
                pts5+[0,-(gr_sep/10),(3*gr_sep/10)], pts5+[0,-(gr_sep/10),-(3*gr_sep/10)], pts5+[0,(gr_sep/10),(3*gr_sep/10)], pts5+[0,(gr_sep/10),-(3*gr_sep/10)],\
                pts5+[0,-(gr_sep/2),gr_sep/10],pts5+[0,-(gr_sep/2),-gr_sep/10], pts5+[0,-(gr_sep/2),(3*gr_sep/10)], pts5+[0,-(gr_sep/2),-(3*gr_sep/10)],\
                pts5+[0,(gr_sep/2),gr_sep/10],pts5+[0,(gr_sep/2),-gr_sep/10], pts5+[0,(gr_sep/2),(3*gr_sep/10)], pts5+[0,(gr_sep/2),-(3*gr_sep/10)],\
                pts5+[0,gr_sep/10,-(gr_sep/2)],pts5+[0,-gr_sep/10,-(gr_sep/2)], pts5+[0,(3*gr_sep/10),-(gr_sep/2)], pts5+[0,-(3*gr_sep/10),-(gr_sep/2)],\
                pts5+[0,gr_sep/10,(gr_sep/2)],pts5+[0,-gr_sep/10,(gr_sep/2)], pts5+[0,(3*gr_sep/10),(gr_sep/2)], pts5+[0,-(3*gr_sep/10),(gr_sep/2)]
        lnx5 = np.concatenate(lnx5, axis=0)

        lny5 = pts5+[gr_sep/2,0,gr_sep/2], pts5+[gr_sep/2,0,-gr_sep/2], pts5+[-gr_sep/2,0,gr_sep/2], pts5+[-gr_sep/2,0,-gr_sep/2],\
                pts5+[-gr_sep/10,0,gr_sep/10], pts5+[-gr_sep/10,0,-gr_sep/10], pts5+[gr_sep/10,0,gr_sep/10], pts5+[gr_sep/10,0,-gr_sep/10],\
                pts5+[-(3*gr_sep/10),0,gr_sep/10],pts5+[-(3*gr_sep/10),0,-gr_sep/10], pts5+[-(3*gr_sep/10),0,(3*gr_sep/10)], pts5+[-(3*gr_sep/10),0,-(3*gr_sep/10)],\
                pts5+[(3*gr_sep/10),0,gr_sep/10],pts5+[(3*gr_sep/10),0,-gr_sep/10], pts5+[(3*gr_sep/10),0,(3*gr_sep/10)], pts5+[(3*gr_sep/10),0,-(3*gr_sep/10)],\
                pts5+[-(gr_sep/10),0,(3*gr_sep/10)], pts5+[-(gr_sep/10),0,-(3*gr_sep/10)], pts5+[(gr_sep/10),0,(3*gr_sep/10)], pts5+[(gr_sep/10),0,-(3*gr_sep/10)],\
                pts5+[-(gr_sep/2),0,gr_sep/10],pts5+[-(gr_sep/2),0,-gr_sep/10], pts5+[-(gr_sep/2),0,(3*gr_sep/10)], pts5+[-(gr_sep/2),0,-(3*gr_sep/10)],\
                pts5+[(gr_sep/2),0,gr_sep/10],pts5+[(gr_sep/2),0,-gr_sep/10], pts5+[(gr_sep/2),0,(3*gr_sep/10)], pts5+[(gr_sep/2),0,-(3*gr_sep/10)],\
                pts5+[gr_sep/10,0,-(gr_sep/2)],pts5+[-gr_sep/10,0,-(gr_sep/2)], pts5+[(3*gr_sep/10),0,-(gr_sep/2)], pts5+[-(3*gr_sep/10),0,-(gr_sep/2)],\
                pts5+[gr_sep/10,0,(gr_sep/2)],pts5+[-gr_sep/10,0,(gr_sep/2)], pts5+[(3*gr_sep/10),0,(gr_sep/2)], pts5+[-(3*gr_sep/10),0,(gr_sep/2)]
        lny5 = np.concatenate(lny5, axis=0)

        lnz5 = pts5+[gr_sep/2,gr_sep/2,0], pts5+[gr_sep/2,-gr_sep/2,0], pts5+[-gr_sep/2,gr_sep/2,0], pts5+[-gr_sep/2,-gr_sep/2,0],\
                pts5+[-gr_sep/10,gr_sep/10,0], pts5+[-gr_sep/10,-gr_sep/10,0], pts5+[gr_sep/10,gr_sep/10,0], pts5+[gr_sep/10,-gr_sep/10,0],\
                pts5+[-(3*gr_sep/10),gr_sep/10,0],pts5+[-(3*gr_sep/10),-gr_sep/10,0], pts5+[-(3*gr_sep/10),(3*gr_sep/10),0], pts5+[-(3*gr_sep/10),-(3*gr_sep/10),0],\
                pts5+[(3*gr_sep/10),gr_sep/10,0],pts5+[(3*gr_sep/10),-gr_sep/10,0], pts5+[(3*gr_sep/10),(3*gr_sep/10),0], pts5+[(3*gr_sep/10),-(3*gr_sep/10),0],\
                pts5+[-(gr_sep/10),(3*gr_sep/10),0], pts5+[-(gr_sep/10),-(3*gr_sep/10),0], pts5+[(gr_sep/10),(3*gr_sep/10),0], pts5+[(gr_sep/10),-(3*gr_sep/10),0],\
                pts5+[-(gr_sep/2),gr_sep/10,0],pts5+[-(gr_sep/2),-gr_sep/10,0], pts5+[-(gr_sep/2),(3*gr_sep/10),0], pts5+[-(gr_sep/2),-(3*gr_sep/10),0],\
                pts5+[(gr_sep/2),gr_sep/10,0],pts5+[(gr_sep/2),-gr_sep/10,0], pts5+[(gr_sep/2),(3*gr_sep/10),0], pts5+[(gr_sep/2),-(3*gr_sep/10),0],\
                pts5+[gr_sep/10,-(gr_sep/2),0],pts5+[-gr_sep/10,-(gr_sep/2),0], pts5+[(3*gr_sep/10),-(gr_sep/2),0], pts5+[-(3*gr_sep/10),-(gr_sep/2),0],\
                pts5+[gr_sep/10,(gr_sep/2),0],pts5+[-gr_sep/10,(gr_sep/2),0], pts5+[(3*gr_sep/10),(gr_sep/2),0], pts5+[-(3*gr_sep/10),(gr_sep/2),0]
        lnz5 = np.concatenate(lnz5, axis=0)




        lnx_5 = pts_5+[0,gr_sep/2,gr_sep/2], pts_5+[0,gr_sep/2,-gr_sep/2], pts_5+[0,-gr_sep/2,gr_sep/2], pts_5+[0,-gr_sep/2,-gr_sep/2],\
                pts_5+[0,-gr_sep/10,gr_sep/10], pts_5+[0,-gr_sep/10,-gr_sep/10], pts_5+[0,gr_sep/10,gr_sep/10], pts_5+[0,gr_sep/10,-gr_sep/10],\
                pts_5+[0,-(3*gr_sep/10),gr_sep/10],pts_5+[0,-(3*gr_sep/10),-gr_sep/10], pts_5+[0,-(3*gr_sep/10),(3*gr_sep/10)], pts_5+[0,-(3*gr_sep/10),-(3*gr_sep/10)],\
                pts_5+[0,(3*gr_sep/10),gr_sep/10],pts_5+[0,(3*gr_sep/10),-gr_sep/10], pts_5+[0,(3*gr_sep/10),(3*gr_sep/10)], pts_5+[0,(3*gr_sep/10),-(3*gr_sep/10)],\
                pts_5+[0,-(gr_sep/10),(3*gr_sep/10)], pts_5+[0,-(gr_sep/10),-(3*gr_sep/10)], pts_5+[0,(gr_sep/10),(3*gr_sep/10)], pts_5+[0,(gr_sep/10),-(3*gr_sep/10)],\
                pts_5+[0,-(gr_sep/2),gr_sep/10],pts_5+[0,-(gr_sep/2),-gr_sep/10], pts_5+[0,-(gr_sep/2),(3*gr_sep/10)], pts_5+[0,-(gr_sep/2),-(3*gr_sep/10)],\
                pts_5+[0,(gr_sep/2),gr_sep/10],pts_5+[0,(gr_sep/2),-gr_sep/10], pts_5+[0,(gr_sep/2),(3*gr_sep/10)], pts_5+[0,(gr_sep/2),-(3*gr_sep/10)],\
                pts_5+[0,gr_sep/10,-(gr_sep/2)],pts_5+[0,-gr_sep/10,-(gr_sep/2)], pts_5+[0,(3*gr_sep/10),-(gr_sep/2)], pts_5+[0,-(3*gr_sep/10),-(gr_sep/2)],\
                pts_5+[0,gr_sep/10,(gr_sep/2)],pts_5+[0,-gr_sep/10,(gr_sep/2)], pts_5+[0,(3*gr_sep/10),(gr_sep/2)], pts_5+[0,-(3*gr_sep/10),(gr_sep/2)]
        lnx_5 = np.concatenate(lnx_5, axis=0)

        lny_5 = pts_5+[gr_sep/2,0,gr_sep/2], pts_5+[gr_sep/2,0,-gr_sep/2], pts_5+[-gr_sep/2,0,gr_sep/2], pts_5+[-gr_sep/2,0,-gr_sep/2],\
                pts_5+[-gr_sep/10,0,gr_sep/10], pts_5+[-gr_sep/10,0,-gr_sep/10], pts_5+[gr_sep/10,0,gr_sep/10], pts_5+[gr_sep/10,0,-gr_sep/10],\
                pts_5+[-(3*gr_sep/10),0,gr_sep/10],pts_5+[-(3*gr_sep/10),0,-gr_sep/10], pts_5+[-(3*gr_sep/10),0,(3*gr_sep/10)], pts_5+[-(3*gr_sep/10),0,-(3*gr_sep/10)],\
                pts_5+[(3*gr_sep/10),0,gr_sep/10],pts_5+[(3*gr_sep/10),0,-gr_sep/10], pts_5+[(3*gr_sep/10),0,(3*gr_sep/10)], pts_5+[(3*gr_sep/10),0,-(3*gr_sep/10)],\
                pts_5+[-(gr_sep/10),0,(3*gr_sep/10)], pts_5+[-(gr_sep/10),0,-(3*gr_sep/10)], pts_5+[(gr_sep/10),0,(3*gr_sep/10)], pts_5+[(gr_sep/10),0,-(3*gr_sep/10)],\
                pts_5+[-(gr_sep/2),0,gr_sep/10],pts_5+[-(gr_sep/2),0,-gr_sep/10], pts_5+[-(gr_sep/2),0,(3*gr_sep/10)], pts_5+[-(gr_sep/2),0,-(3*gr_sep/10)],\
                pts_5+[(gr_sep/2),0,gr_sep/10],pts_5+[(gr_sep/2),0,-gr_sep/10], pts_5+[(gr_sep/2),0,(3*gr_sep/10)], pts_5+[(gr_sep/2),0,-(3*gr_sep/10)],\
                pts_5+[gr_sep/10,0,-(gr_sep/2)],pts_5+[-gr_sep/10,0,-(gr_sep/2)], pts_5+[(3*gr_sep/10),0,-(gr_sep/2)], pts_5+[-(3*gr_sep/10),0,-(gr_sep/2)],\
                pts_5+[gr_sep/10,0,(gr_sep/2)],pts_5+[-gr_sep/10,0,(gr_sep/2)], pts_5+[(3*gr_sep/10),0,(gr_sep/2)], pts_5+[-(3*gr_sep/10),0,(gr_sep/2)]
        lny_5 = np.concatenate(lny_5, axis=0)

        lnz_5 = pts_5+[gr_sep/2,gr_sep/2,0], pts_5+[gr_sep/2,-gr_sep/2,0], pts_5+[-gr_sep/2,gr_sep/2,0], pts_5+[-gr_sep/2,-gr_sep/2,0],\
                pts_5+[-gr_sep/10,gr_sep/10,0], pts_5+[-gr_sep/10,-gr_sep/10,0], pts_5+[gr_sep/10,gr_sep/10,0], pts_5+[gr_sep/10,-gr_sep/10,0],\
                pts_5+[-(3*gr_sep/10),gr_sep/10,0],pts_5+[-(3*gr_sep/10),-gr_sep/10,0], pts_5+[-(3*gr_sep/10),(3*gr_sep/10),0], pts_5+[-(3*gr_sep/10),-(3*gr_sep/10),0],\
                pts_5+[(3*gr_sep/10),gr_sep/10,0],pts_5+[(3*gr_sep/10),-gr_sep/10,0], pts_5+[(3*gr_sep/10),(3*gr_sep/10),0], pts_5+[(3*gr_sep/10),-(3*gr_sep/10),0],\
                pts_5+[-(gr_sep/10),(3*gr_sep/10),0], pts_5+[-(gr_sep/10),-(3*gr_sep/10),0], pts_5+[(gr_sep/10),(3*gr_sep/10),0], pts_5+[(gr_sep/10),-(3*gr_sep/10),0],\
                pts_5+[-(gr_sep/2),gr_sep/10,0],pts_5+[-(gr_sep/2),-gr_sep/10,0], pts_5+[-(gr_sep/2),(3*gr_sep/10),0], pts_5+[-(gr_sep/2),-(3*gr_sep/10),0],\
                pts_5+[(gr_sep/2),gr_sep/10,0],pts_5+[(gr_sep/2),-gr_sep/10,0], pts_5+[(gr_sep/2),(3*gr_sep/10),0], pts_5+[(gr_sep/2),-(3*gr_sep/10),0],\
                pts_5+[gr_sep/10,-(gr_sep/2),0],pts_5+[-gr_sep/10,-(gr_sep/2),0], pts_5+[(3*gr_sep/10),-(gr_sep/2),0], pts_5+[-(3*gr_sep/10),-(gr_sep/2),0],\
                pts_5+[gr_sep/10,(gr_sep/2),0],pts_5+[-gr_sep/10,(gr_sep/2),0], pts_5+[(3*gr_sep/10),(gr_sep/2),0], pts_5+[-(3*gr_sep/10),(gr_sep/2),0]
        lnz_5 = np.concatenate(lnz_5, axis=0)



        #===============================PLOTTING========================================================================




        ###----coord for line centres around each point-----------

        ###-------------------------------------------------------


        def plotter(lnx, lny, lnz, clr):



            line1 = tvtk.LineSource(point1=(-gr_sep/2,0,0), point2=(gr_sep/2,0,0))
            line2 = tvtk.LineSource(point1=(0,-gr_sep/2,0), point2=(0,gr_sep/2,0))
            line3 = tvtk.LineSource(point1=(0,0,-gr_sep/2), point2=(0,0,gr_sep/2))

            pd1 = tvtk.PolyData(points=lnx)
            pd2 = tvtk.PolyData(points=lny)
            pd3 = tvtk.PolyData(points=lnz)

            g1 = tvtk.Glyph3D(scale_mode='data_scaling_off', vector_mode = 'use_vector')
            g2 = tvtk.Glyph3D(scale_mode='data_scaling_off', vector_mode = 'use_vector')
            g3 = tvtk.Glyph3D(scale_mode='data_scaling_off', vector_mode = 'use_vector')



            configure_input_data(g1, pd1)
            configure_input_data(g2, pd2)
            configure_input_data(g3, pd3)

            configure_source_data(g1, line1.output)
            line1.update()
            g1.update()
            configure_source_data(g2, line2.output)
            line2.update()
            g2.update()
            configure_source_data(g3, line3.output)
            line3.update()
            g3.update()



            m1 = tvtk.PolyDataMapper()
            m2 = tvtk.PolyDataMapper()
            m3 = tvtk.PolyDataMapper()

            pc1 = tvtk.Property(opacity=0.9, color=clr)


            configure_input_data(m1, g1.output)
            configure_input_data(m2, g2.output)
            configure_input_data(m3, g3.output)

            a1 = tvtk.Actor(mapper=m1, property=pc1)
            a2 = tvtk.Actor(mapper=m2, property=pc1)
            a3 = tvtk.Actor(mapper=m3, property=pc1)

            scn.add_actor(a1)
            scn.add_actor(a2)
            scn.add_actor(a3)


        plotter(lnx1, lny1, lnz1, (0.9,0,0))
        plotter(lnx_1, lny_1, lnz_1, (0,0,0.9))

        plotter(lnx2, lny2, lnz2, (0.9,0,0))
        plotter(lnx_2, lny_2, lnz_2, (0,0,0.9))

        plotter(lnx3, lny3, lnz3, (0.9,0,0))
        plotter(lnx_3, lny_3, lnz_3, (0,0,0.9))

        plotter(lnx4, lny4, lnz4, (0.9,0,0))
        plotter(lnx_4, lny_4, lnz_4, (0,0,0.9))

        plotter(lnx5, lny5, lnz5, (0.9,0,0))
        plotter(lnx_5, lny_5, lnz_5, (0,0,0.9))
        

        xmin = float((np.min(self.xg)) - gr_sep)
        ymin = float((np.min(self.yg)) - gr_sep)
        zmin = float((np.min(self.zg)) - gr_sep)
        xmax = float((np.max(self.xg)) + gr_sep)
        ymax = float((np.max(self.yg)) + gr_sep)
        zmax = float((np.max(self.zg)) + gr_sep)


        xlab='dx'
        ylab='dy'
        zlab='dz'
        

        scn.mlab.points3d(([xmin,ymin,zmin],[xmin,ymin,zmax],[xmin,ymax,zmin],[xmin,ymax,zmax],[xmax,ymin,zmin],[xmax,ymin,zmax],[xmax,ymax,zmin],[xmax,ymax,zmax]),opacity=0.0)
        scn.mlab.points3d(pts_nan[:,0],
                            pts_nan[:,1],
                            pts_nan[:,2], color = (1,0,0),scale_factor=gr_sep, resolution=36)
        scn.mlab.axes(extent = [xmin,xmax,ymin,ymax,zmin,zmax], nb_labels = 5, line_width=3.0, color = (0,0,0), xlabel=xlab, ylabel=ylab, zlabel=zlab)

        xfoc = (xmin+xmax)/2
        yfoc = (ymin+ymax)/2
        zfoc = (zmin+zmax)/2


        scn.mlab.view(focalpoint=[xfoc,yfoc,zfoc])
        scn.mlab.show()


























