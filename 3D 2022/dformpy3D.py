import numpy as np
from mayavi import mlab
import matplotlib.pyplot as plt
import sympy
from sympy import sympify
from sympy import diff, simplify
from sympy.parsing.sympy_parser import parse_expr

from tvtk.api import tvtk
from tvtk.common import configure_input_data, configure_source_data

from numpy import matrix, sin, cos, tan, sqrt, log, arctan, arcsin, arccos, tanh
from numpy import sinh, cosh, arcsinh, arccosh, arctanh, exp, pi, e


def find_2_form(expressions, coords, xg, yg, zg=None, m=3):
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
    [if applied to object can use add_curl ='y' argument to plot
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
        It must be in terms of x, y and z.
        Has to be given, for curl and div to be calculatable.
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

            curl_vf = vector_field3(self.xg, self.yg, self.zg, Curl_X, Curl_Y, Curl_Z, self.str_x, self.str_y, self.str_z)
            return(curl_vf)

    def div(self, at_x=None, at_y=None, at_z=None):
        '''
        When applied to a vector field object, prints out
        the total divergence equation for the field.

        If three coordinates for a point (x, y, z) are supplied,
        also prints out approximate divergence of the field
        at that points

        --------------------------------------------------------
        Returns:

        Nothing

        --------------------------------------------------------
        Arguments:

        *Optional:
            at_x - x coordinate of the point
            at_y -  ...
            at_z -  ...
        '''

        # Check whether the string equations are supplied
        if self.str_x == None or self.str_y == None or self.str_z == None:
            # ERROR
            raise TypeError('No equation provided')
        else:

            # define separation between points to define a new grid with boundaries of the input
            # grid but with pre-defined point density
            sep = (abs(np.min(self.xg)) + abs(np.max(self.xg)))/100

            # define the new grid spacing
            New_grid = np.linspace(np.min(self.xg), np.max(self.xg), 100)
            
            #new, denser grid using new coordinates
            xng, yng, zng = np.meshgrid(New_grid, New_grid, New_grid)

            # convert the x, y, g coordinate list into a 2D array containing
            # x, y and z coordinates within one list (i.e. list defining each point of the grid)
            pts = np.vstack(list(zip(xng.ravel(), yng.ravel(), zng.ravel())))


           
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


            # Check whether supplied coordinates are within the grid boundaries
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

                # turn the separation into a string to find decimal points
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

                # convert the divergence value list into a 2D array containing
                # div at each point within one list
                Div_lst = np.vstack(list(zip(Div.ravel())))

                # convert all np.inf to np.nan
                Div_lst[np.isinf(Div_lst)] = np.nan

                

                # calculate magnitude of difference between each point component and the input coordinates
                # to find points on the grid which are the closest to the input one
                
                pts_diff = (pts[:]-[at_x,at_y,at_z])
                diffmag = np.sqrt(pts_diff[:,0]**2 + pts_diff[:,1]**2 + pts_diff[:,2]**2)
                diffmag = np.vstack(list(zip(diffmag.ravel())))



                # index of the point closest to the input
                a = np.argwhere(diffmag == (np.nanmin(diffmag)))
                a = a[:,0]

                print('Divergence at the grid point/points closest to the input, '+ str(pts[a]) +', = '+ str(Div_lst[a]))
                

    def log_scaling(self):
        '''
        Takes no arguments
        Changes the boolean that determines if scaling is logarithmic
        Whenever it is called, it changes that boolean to opposite
        The form object is initialised with this as False
        '''
        self.logarithmic_scale_bool = not self.logarithmic_scale_bool


    def zoom(self, mag, target, dpd):
        

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


                # evaluate the components numerically, inputting them into a
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
        
        inset - bool - if true, field deriv is plotted on given axis
        axis - matplotlib axes instance - axis to plot on if instance it True
        
        Returns:
        --------
        if inset is False:
            deriv VF object
        if inset is True, inset axis and deriv VF object in this order.
        
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
                
                # Check if the equations provided contain x and y terms
                


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

    def plot(self, add_curl = None, arrow_colour = None, arrow_cmap=None, scaling=None, opacity=None, curl_opacity=None):

        '''
        Plots a vector_field3 object.

        ------------------------------------------
        Returns:

        None

        ------------------------------------------
        Arguments:
        *Optional:
        add_curl = 'y' - adds another vector field to the plot, which represents curl of the initial vector field
        arrow_colour - takes input of either red, green, blue, magenta, black, yellow or cyan colours. (for const. magnitude/overlay curl plots)
        arrow-cmap - takes input of possible Mayavi colormaps for plots of varying magnitude. (default = viridis)
        scaling - scale arrows. Has to be > 0.
        opacity - changes opacity of the initial quiver plot 0>=opac>=1.0
        curl_opacity - changes opacity of the overlay curl plot. Same boundaries as for regular opacity argument
        '''

        # set initial conditions for when the arguments are not supplied
        cmap = 'viridis'
        clr = (0.2,1,0)

        scl = 1.0
        
        opc_crl = 1.0
        opc = 1.0

        # check for color input and set the corresponding RGB triplet
        if arrow_colour=='red' or arrow_colour=='r' or arrow_colour=='Red':
            clr = (1,0,0)
        elif arrow_colour=='black' or arrow_colour=='k' or arrow_colour=='Black':
            clr = (0,0,0)
        elif arrow_colour=='green' or arrow_colour=='g' or arrow_colour=='Green':
            clr = (0,1,0)
        elif arrow_colour=='blue' or arrow_colour=='b' or arrow_colour=='Blue':
            clr = (0,0,1)
        elif arrow_colour=='magenta' or arrow_colour=='m' or arrow_colour=='Magenta':
            clr = (1,0,1)
        elif arrow_colour=='yellow' or arrow_colour=='y' or arrow_colour=='Yellow':
            clr = (1,1,0)
        elif arrow_colour=='cyan' or arrow_colour=='c' or arrow_colour=='Cyan':
            clr = (0,1,1)
        else:
            pass
        
        # Set colormap to the user input (no need to check for errors since)
        # the Mayavi package will give out an error if invalid cmap name will be supplied
        if arrow_cmap != None:
            cmap = arrow_cmap
        else:
            pass

        # Check whether the scaling input is float and > 0.0
        if scaling != None:
            if isinstance(scaling, float)==True and scaling>=0.0:
                scl = scaling       
            else:
                scl=1.0
                print('The scaling factor has to be a float (0.0 < scaling)')
                exit()

        # Check whether the opacity input is float and 0.0 > opac > 1.0
        if opacity != None:
            if isinstance(opacity, float)==True and opacity<=1.0 and opacity>=0.0:   
                opc = opacity
            else:

                print('The opacity has to be a float (0.0 < opacity < 1.0)')
                exit()

        # Check whether the curl opacity input is float and 0.0 > opac > 1.0
        if curl_opacity != None:
            if isinstance(curl_opacity, float)==True and curl_opacity<=1.0 and curl_opacity>=0.0:   
                opc_crl = curl_opacity
            else:

                print('The opacity has to be a float (0.0 < opacity < 1.0)')
                exit()

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
        
        # define grid dimension magnitude

        Lx = (self.xg[-1, -1, -1] - self.xg[0, 0, 0])
        Ly = (self.yg[-1, -1, -1] - self.yg[0, 0, 0])
        Lz = (self.zg[-1, -1, -1] - self.zg[0, 0, 0])

        # use the magnitude to define min and max values of x, y and z directions
        # these to be used for mlab.axes() plot
        xmin = int((np.min(self.xg)))
        ymin = int((np.min(self.yg)))
        zmin = int((np.min(self.zg)))
        xmax = int((np.max(self.xg)))
        ymax = int((np.max(self.yg)))
        zmax = int((np.max(self.zg)))

        # define mlab figure upon which the field will be plotted
        fig = mlab.figure(bgcolor = (1,1,1), fgcolor = (0,0,0))

        # all np.inf -> np.nan for convenience purposes
        F_x_local[np.isinf(F_x_local)] = np.nan
        F_y_local[np.isinf(F_y_local)] = np.nan
        F_z_local[np.isinf(F_z_local)] = np.nan

        # Field magnitude at each point
        F = np.sqrt(F_x_local**2+F_y_local**2+F_z_local**2)



        # Indices where field is not defined (Field. mag. = np.nan)
        Idx = np.argwhere(np.isnan(F))

        # Set all NaN values to zero so that it does not disturb the plotting
        F_x_local[np.isnan(F_x_local)] = 0
        F_y_local[np.isnan(F_y_local)] = 0
        F_z_local[np.isnan(F_z_local)] = 0


        # convert to log scale if method applied to the vf object 
        if self.logarithmic_scale_bool:
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


        # Field magnitude at each point
        F = np.sqrt(F_x_local**2+F_y_local**2+F_z_local**2)
        max_size = np.max(F)

        if self.scale_bool is False:
            ScaleFactor = scl
        elif self.scale_bool is True:
            ScaleFactor = max_size/(0.9*(2*Lx/self.pt_den))

        # if magnitude is constant, use color instead of colormap for the quiver
        if abs((np.nanmax(F)-np.nanmin(F)))>=0.001:
            mlab.quiver3d(self.xg, self.yg, self.zg, F_x_local, F_y_local, F_z_local, colormap=cmap, line_width=3.0, mode='arrow', scale_factor=ScaleFactor, scale_mode = 'vector', opacity = opc)
        else:
            mlab.quiver3d(self.xg, self.yg, self.zg, F_x_local, F_y_local, F_z_local, color=clr, line_width=3.0, mode='arrow', scale_factor=ScaleFactor, scale_mode = 'vector', opacity = opc)

        
        # adds curl plot if the argument is supplied
        # Same code as in .curl() method
        if add_curl != None:
            if str(add_curl)=='y':
                # check whether equations are supplied
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

                    Curl_X[np.isinf(Curl_X)] = 0
                    Curl_Y[np.isinf(Curl_Y)] = 0
                    Curl_Z[np.isinf(Curl_Z)] = 0

                    Curl_X[np.isnan(Curl_X)] = 0
                    Curl_Y[np.isnan(Curl_Y)] = 0
                    Curl_Z[np.isnan(Curl_Z)] = 0

                    Curl_X[np.abs(F_x_local) < 1e-15] = 0
                    Curl_Y[np.abs(F_y_local) < 1e-15] = 0
                    Curl_Z[np.abs(F_z_local) < 1e-15] = 0

                    mlab.quiver3d(self.xg, self.yg, self.zg, Curl_X, Curl_Y, Curl_Z, color=(1.0, 0.0, 1.0), opacity=opc_crl, mode='arrow',scale_factor=ScaleFactor/10, scale_mode = 'vector')
            elif str(add_curl)=='n':
                pass
            else:
                print('add curl argument has to be either [y] or [n]')
                exit()
        

        
        # Plot red spheres at singular points
        mlab.points3d(self.xg[Idx[:,0],Idx[:,1],Idx[:,2]],
                      self.yg[Idx[:,0],Idx[:,1],Idx[:,2]],
                      self.zg[Idx[:,0],Idx[:,1],Idx[:,2]], color = (1,0,0),scale_factor=0.4)


        # add outline box and axes
        mlab.outline(extent = [xmin,xmax,ymin,ymax,zmin,zmax], line_width=1.0)
        axes = mlab.axes(color = (0,0,0), nb_labels = 5, extent = [xmin,xmax,ymin,ymax,zmin,zmax], line_width=3.0)

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

    .give_eqn(equation_str_x, equation_str_y, equation_str_z) 
    .levels()
    .set_density()
    .log_scaling()
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
        
        # Check if the equations provided contain x and y terms
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
            # re-evaluate the 2-form numerically
            self.form_0 = eval(str_0)



    def zoom(self, mag, target, dpd):
        



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

                
                # Check if the equations provided contain x and y terms
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
            xmin = int(np.min(self.xg))
            ymin = int(np.min(self.yg))
            zmin = int(np.min(self.zg))
            xmax = int(np.max(self.xg))
            ymax = int(np.max(self.yg))
            zmax = int(np.max(self.zg))

            # create figure to host the plot
            fig = mlab.figure(bgcolor = (1,1,1), fgcolor = (0,0,0))

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

            # plot the 0 form contours, add colorbar, axes and an outline
            cnt = mlab.contour3d(form_0, colormap=colourmap, opacity = opac, contours=self.lines, figure = fig)
            mlab.colorbar(object = cnt, orientation='vertical')
            mlab.outline(line_width=1.0)
            mlab.axes(color = (0,0,0), ranges = (xmin, xmax, ymin, ymax, zmin, zmax), nb_labels = 5, line_width=3.0)

            
            # if the argument is 'y',  make the contours transparent and add the slider to view scalar cross section
            if str(cross_sec_plane)=='y':

                opac=0.05
                colour=(0.5,0.5,0.5)
                # clear figure from previous plot
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
            coords = ['x', 'y', 'z']
            # from these, find the derivatives
            form_1_x_str = str(diff(sympy_expr_form_0, coords[0]))
            form_1_y_str = str(diff(sympy_expr_form_0, coords[1]))
            form_1_z_str = str(diff(sympy_expr_form_0, coords[2]))
            # need to uspply these unformatted, so save those:
            form_1_x_unformated, form_1_y_unformated, form_1_z_unformated = form_1_x_str*1, form_1_y_str*1, form_1_z_str*1
            # from these strings, get the numerical 1-form:
            form_1_x_str = form_1_x_str.replace('x', '(self.xg)')
            form_1_x_str = form_1_x_str.replace('y', '(self.yg)')
            form_1_x_str = form_1_x_str.replace('z', '(self.zg)')

            form_1_y_str = form_1_y_str.replace('x', '(self.xg)')
            form_1_y_str = form_1_y_str.replace('y', '(self.yg)')
            form_1_y_str = form_1_y_str.replace('z', '(self.zg)')

            form_1_z_str = form_1_z_str.replace('x', '(self.xg)')
            form_1_z_str = form_1_z_str.replace('y', '(self.yg)')
            form_1_z_str = form_1_z_str.replace('z', '(self.zg)')


            if form_1_x_str.find('x') & form_1_x_str.find('y') & form_1_x_str.find('z') == -1:
                form_1_x_str = '(' + str(form_1_x_str) + ')* np.ones(np.shape(self.xg))'
            if form_1_y_str.find('x') & form_1_y_str.find('y') & form_1_y_str.find('z') == -1:
                form_1_y_str = '(' + str(form_1_y_str) + ')* np.ones(np.shape(self.yg))'
            if form_1_z_str.find('x') & form_1_z_str.find('y') & form_1_z_str.find('z') == -1:
                form_1_z_str = '(' + str(form_1_z_str) + ')* np.ones(np.shape(self.yg))'


            form_1_x = eval(form_1_x_str)
            form_1_y = eval(form_1_y_str)
            form_1_z = eval(form_1_z_str)
            
            # supply these to the 1-form object function and return object
            result_1_form = form_1_3d(self.xg, self.yg, self.zg, form_1_x, form_1_y, form_1_z, form_1_x_unformated, form_1_y_unformated, form_1_z_unformated)

            return result_1_form



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
        fy, fx, fz= np.gradient(self.form_0, edge_order=edge_order)
        
        # supply these to the 1-form object function
        result_1_form = form_1_3d(self.xg, self.yg, self.zg, fx, fy, fz)
        
        # return the new object to user
        return result_1_form



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
            string_2_form = string_2_form.replace('z', '(self.zg)')
            
            if string_2_form.find('x') & string_2_form.find('y') & string_2_form.find('z') == -1:
                string_2_form = '(' + str(string_2_form) + ')* np.ones(np.shape(self.xg))'
            
            # evaulated numerically
            form_2_result = eval(string_2_form)
            
            # create and return object
            new_object = form_2_3d(self.xg, self.yg, self.zg, Fx = form_2_result, Fy = form_2_result, Fz = form_2_result, Fx_eqn=form_2_str_unformated, Fy_eqn=form_2_str_unformated, Fz_eqn=form_2_str_unformated)
            return new_object
        else:
            # ERROR
            raise TypeError('You need to supply the 2-form equations to do this, look at \'give_eqn\' method')



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
        new_object = form_2_3d(self.xg, self.yg, self.zg, Fx = self.form_0, Fy = self.form_0, Fz = self.form_0)  # N.B no equations to supply
        return new_object



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
            if len(form_second) == 3:
                # 0-form/\1-form, check if strings supplied
                if isinstance(form_second[0], str) and isinstance(form_second[1], str) and isinstance(form_second[2], str):
                    to_wedge_x_2_str = form_second[0]
                    to_wedge_y_2_str = form_second[1]
                    to_wedge_z_2_str = form_second[2]
                    if degree == 1:
                        order = 1
                    elif degree == 2:
                        order = 2
                else:
                    raise ValueError('for analytical calulation, supply 1-form equations as strings')
            else:
                raise ValueError('too many or too little equations given in tuple')
        elif isinstance(form_second, str):
            # single string, could be 0-form or 2-form, check given degree:
            if degree == 0:
                to_wedge_0_form_str = form_second
                order = 0

            elif degree == 3:
                to_wedge_3_form_str = form_second
                order = 3
            else:
                raise ValueError('not possible digree given or supplied one string for a 1-form')
        else:
            # object supplied, get numericals checking which object is given:
            if isinstance(form_second, form_1_3d):
                if form_second.form_1_str_x is None or form_second.form_1_str_y is None or form_second.form_1_str_z is None:
                     raise ValueError('supplied 1-form instance must contain equations for analytical calculation')
                else:
                    to_wedge_x_2_str = form_second.form_1_str_x
                    to_wedge_y_2_str = form_second.form_1_str_y
                    to_wedge_z_2_str = form_second.form_1_str_z
                    order = 1
            elif isinstance(form_second, form_0_3d):
                if form_second.form_0_str is None:
                    raise ValueError('supplied 0-form instance must contain equations for analytical calculation')
                else:
                    to_wedge_0_form_str = form_second.form_0_str
                    order = 0       
            elif isinstance(form_second, form_2_3d):
                if form_second.Fx_eqn is None and form_second.Fy_eqn is None and form_second.Fz_eqn is None:
                    raise ValueError('supplied 2-form instance must contain equations for analytical calculation')
                else:
                    to_wedge_2_form_str_x = form_second.Fx_eqn
                    to_wedge_2_form_str_y = form_second.Fy_eqn
                    to_wedge_2_form_str_z = form_second.Fz_eqn
                    order = 2
            elif isinstance(form_second, form_3_3d):
                if form_second.form_3_str is None:
                    raise ValueError('supplied 3-form instance must contain equations for analytical calculation')
                else:
                    to_wedge_3_form_str = form_second.form_3_str
                    order = 3
            else:
                raise TypeError('Supplied form to wedge with is not recognised')
        
        # Deal with 0-form/\1-form:
        if order == 1:
            # first, find the result of the 1-form:
            new_str_x = str(simplify('(' + self.form_0_str + ')*(' +  to_wedge_x_2_str + ')'))
            new_str_y = str(simplify('(' + self.form_0_str + ')*(' +  to_wedge_y_2_str + ')'))
            new_str_z = str(simplify('(' + self.form_0_str + ')*(' +  to_wedge_z_2_str + ')'))
            # keep it as it is locally to supply it to object maker later
            form_1_str_x_loc = new_str_x + ''
            form_1_str_y_loc = new_str_y + ''
            form_1_str_z_loc = new_str_z + ''
            # format it to be in terms of grids and:
            # check against constant and zero 1-forms being supplied
            # get the numerical evaluation of it
            new_str_x = new_str_x.replace('x', '(self.xg)')
            new_str_x = new_str_x.replace('y', '(self.yg)')
            new_str_x = new_str_x.replace('z', '(self.zg)')

            new_str_y = new_str_y.replace('x', '(self.xg)')
            new_str_y = new_str_y.replace('y', '(self.yg)')
            new_str_y = new_str_y.replace('z', '(self.zg)')

            new_str_z = new_str_z.replace('x', '(self.xg)')
            new_str_z = new_str_z.replace('y', '(self.yg)')
            new_str_z = new_str_z.replace('z', '(self.zg)')
            
            if new_str_x.find('x') & new_str_x.find('y') & new_str_x.find('z')== -1:
                new_str_x = '(' + str(new_str_x) + ')* np.ones(np.shape(self.xg))'
            if new_str_y.find('x') & new_str_y.find('y') & new_str_y.find('z') == -1:
                new_str_y = '(' + str(new_str_y) + ')* np.ones(np.shape(self.yg))'
            if new_str_z.find('x') & new_str_z.find('y') & new_str_z.find('z') == -1:
                new_str_z = '(' + str(new_str_z) + ')* np.ones(np.shape(self.zg))'
            
            form_1_x = eval(new_str_x)
            form_1_y = eval(new_str_y)
            form_1_z = eval(new_str_z)

            # return the new one to the user:
            new_object = form_1_3d(self.xg, self.yg, self.zg, form_1_x, form_1_y, form_1_z, F_x_eqn=form_1_str_x_loc, F_y_eqn=form_1_str_y_loc, F_z_eqn = form_1_str_z_loc)
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
            form_0_str = form_0_str.replace('z', 'self.zg')
            if form_0_str.find('x') & form_0_str.find('y') & form_0_str.find('z') == -1:
                form_0_str = '(' + str(form_0_str) + ')* np.ones(np.shape(self.xg))'
            
            # evaluate it numerically on the grid supplied
            form_0_result = eval(form_0_str)
            
            # depending on keep_object, return:
            if keep_object:
                self.form_0 = form_0_result
                self.form_0_str = form_0_str_loc
            elif not keep_object:
                new_object = form_0_3d(self.xg, self.yg, self.zg, form_0_result, form_0_str_loc)
                # return the new one to the user:
                return new_object
            else:
                raise ValueError('Error, Invalid input for \'keep_object\'')

        elif order == 2:
            form_2_str_x = str(simplify( '(' + self.form_0_str + ')*(' +  to_wedge_2_form_str_x + ')'))
            form_2_str_y = str(simplify( '(' + self.form_0_str + ')*(' +  to_wedge_2_form_str_y + ')'))
            form_2_str_z = str(simplify( '(' + self.form_0_str + ')*(' +  to_wedge_2_form_str_z + ')'))
            # keep it as it is locally to supply it to object maker later
            form_2_str_x_loc = form_2_str_x + ''
            form_2_str_y_loc = form_2_str_y + ''
            form_2_str_z_loc = form_2_str_z + ''
            # format it to be in terms of grids and:
            # check against constant and zero 2-forms being supplied
            # get the numerical evaluation of it
            form_2_str_x = form_2_str_x.replace('x', 'self.xg')
            form_2_str_x = form_2_str_x.replace('y', 'self.yg')
            form_2_str_x = form_2_str_x.replace('z', 'self.zg')

            form_2_str_y = form_2_str_y.replace('x', 'self.xg')
            form_2_str_y = form_2_str_y.replace('y', 'self.yg')
            form_2_str_y = form_2_str_y.replace('z', 'self.zg')

            form_2_str_z = form_2_str_z.replace('x', 'self.xg')
            form_2_str_z = form_2_str_z.replace('y', 'self.yg')
            form_2_str_z = form_2_str_z.replace('z', 'self.zg')

            if form_2_str_x.find('x') & form_2_str_x.find('y') & form_2_str_x.find('z') == -1:
                form_2_str_x = '(' + str(form_2_str_x) + ')* np.ones(np.shape(self.xg))'
            if form_2_str_y.find('x') & form_2_str_y.find('y') & form_2_str_y.find('z') == -1:
                form_2_str_y = '(' + str(form_2_str_y) + ')* np.ones(np.shape(self.yg))'
            if form_2_str_z.find('x') & form_2_str_z.find('y') & form_2_str_z.find('z') == -1:
                form_2_str_z = '(' + str(form_2_str_z) + ')* np.ones(np.shape(self.zg))'
            
            # evaluate it numerically on the grid supplied
            form_2_x_result = eval(form_2_str_x)
            form_2_y_result = eval(form_2_str_y)
            form_2_z_result = eval(form_2_str_z)
            
            # create new instance and return to user
            new_object = form_2_3d(self.xg, self.yg, self.zg, form_2_x_result, form_2_y_result, form_2_z_result, form_2_str_x_loc, form_2_str_y_loc, form_2_str_z_loc)
            return new_object

        elif order == 3:
            form_3_str = str(simplify( '(' + self.form_0_str + ')*(' +  to_wedge_3_form_str + ')'))
            # keep it as it is locally to supply it to object maker later
            form_3_str_loc = form_3_str + ''
            # format it to be in terms of grids and:
            # check against constant and zero 2-forms being supplied
            # get the numerical evaluation of it
            form_3_str = form_3_str.replace('x', 'self.xg')
            form_3_str = form_3_str.replace('y', 'self.yg')
            form_3_str = form_3_str.replace('z', 'self.zg')
            if form_3_str.find('x') & form_3_str.find('y') & form_3_str.find('z') == -1:
                form_3_str = '(' + str(form_3_str) + ')* np.ones(np.shape(self.xg))'
            
            # evaluate it numerically on the grid supplied
            form_3_result = eval(form_3_str)
            
            # depending on keep_object, return:
            if keep_object:
                self.form_3 = form_3_result
                self.form_3_str = form_3_str_loc
            elif not keep_object:
                new_object = form_3_3d(self.xg, self.yg, self.zg, form_3_result, form_3_str_loc)
                # return the new one to the user:
                return new_object
            else:
                raise ValueError('Error, Invalid input for \'keep_object\'')

        else:
            # should never happen, but in case
            raise ValueError('Variable change during code running, look at \'order\' parameter')
    


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
            if len(form_second) == 3:
                # 0-form/\1-form, extract components
                # if numerical grids were given, take these, if equations, change to values on grids:
                if isinstance(form_second[0], str) and isinstance(form_second[1], str) and isinstance(form_second[2], str):

                    new_str_x = form_second[0].replace('x', '(self.xg)')
                    new_str_x = new_str_x.replace('y', '(self.yg)')
                    new_str_x = new_str_x.replace('z', '(self.zg)')

                    new_str_y = form_second[1].replace('x', '(self.xg)')
                    new_str_y = new_str_y.replace('y', '(self.yg)')
                    new_str_y = new_str_y.replace('z', '(self.zg)')

                    new_str_z = form_second[2].replace('x', '(self.xg)')
                    new_str_z = new_str_z.replace('y', '(self.yg)')
                    new_str_z = new_str_z.replace('z', '(self.zg)')


                    if new_str_x.find('x') & new_str_x.find('y') & new_str_x.find('z')== -1:
                        new_str_x = '(' + str(new_str_x) + ')* np.ones(np.shape(self.xg))'
                    if new_str_y.find('x') & new_str_y.find('y') & new_str_y.find('z') == -1:
                        new_str_y = '(' + str(new_str_y) + ')* np.ones(np.shape(self.yg))'
                    if new_str_z.find('x') & new_str_z.find('y') & new_str_z.find('z') == -1:
                        new_str_z = '(' + str(new_str_z) + ')* np.ones(np.shape(self.zg))'
                    
                    f12_x = eval(new_str_x)
                    f12_y = eval(new_str_y)
                    f12_z = eval(new_str_z)
                    order = 1
                elif isinstance(form_second[0], np.ndarray) and isinstance(form_second[1], np.ndarray) and isinstance(form_second[2], np.ndarray):
                    f12_x = form_second[0]
                    f12_y = form_second[1]
                    f12_z = form_second[2]
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
                raise ValueError('for degree 2, supply a 1-form, not a single grid')
            elif degree == 3:
                to_wedge_3_form = form_second
                order = 3
        
        elif isinstance(form_second, str):
            # single string, could be 0-form or 2-form, check given degree:
            if degree == 0:
                    str_0_form = form_second.replace('x', '(self.xg)')
                    str_0_form = str_0_form.replace('y', '(self.yg)')
                    str_0_form = str_0_form.replace('z', '(self.zg)')
                    if str_0_form.find('x') & str_0_form.find('y') & str_0_form.find('z') == -1:
                        str_0_form = '(' + str(str_0_form) + ')* np.ones(np.shape(self.xg))'
                    
                    to_wedge_0_form = eval(str_0_form)
                    order = 0
            elif degree == 3:
                str_3_form = form_second.replace('x', '(self.xg)')
                str_3_form = str_3_form.replace('y', '(self.yg)')
                str_3_form = str_3_form.replace('z', '(self.zg)')
                if str_3_form.find('x') & str_3_form.find('y') & str_3_form.find('z') == -1:
                    str_3_form = '(' + str(str_3_form) + ')* np.ones(np.shape(self.xg))'
                
                to_wedge_2_form = eval(str_3_form)
                order = 3
            else:
                raise ValueError('not possible digree given or supplied one string for a 1-form and/or 2-dorm')
        
        # object supplied, get grids checking which object is given:
        
        elif isinstance(form_second, form_1_3d):
            f12_x = form_second.F_x
            f12_y = form_second.F_y
            f12_z = form_second.F_z
            order = 1
        elif isinstance(form_second, form_0_3d):
            to_wedge_0_form = form_second.form_0
            order = 0
        elif isinstance(form_second, form_2_3d):
            order = 2
            to_wedge_2_form_x = form_second.Fx
            to_wedge_2_form_y = form_second.Fy
            to_wedge_2_form_z = form_second.Fz
        elif isinstance(form_second, form_3_3d):
            to_wedge_3_form = form_second.form_3
            order = 3
        else:
            raise TypeError('Supplied form to wedge with is not recognised')
        
        # Use given inputs to evaluate the result:
        
        # Deal with 0-form/\1-form:
        if order == 1:
            # first, find the result of the 1-form
            new_form_1_x = self.form_0 * f12_x
            new_form_1_y = self.form_0 * f12_y
            new_form_1_z = self.form_0 * f12_z
            
            # create instance and return
            new_object = form_1_3d(self.xg, self.yg, self.zg, new_form_1_x, new_form_1_y, new_form_1_z)
            return new_object
        
        elif order == 0:
            # from these get the numerical 0-form
            form_0_result = self.form_0 * to_wedge_0_form
            
            # depending on keep_object, return:
            if keep_object:
                self.form_0 = form_0_result
            elif not keep_object:
                new_object = form_0_3d(self.xg, self.yg, self.zg, form_0_result)
                # return the new one to the user:
                return new_object
            else:
                raise ValueError('Error, Invalid input for \'keep_object\'')

        elif order == 2:
            # from these get the numerical 0-form
            form_2_result_x = self.form_0 * to_wedge_2_form_x
            form_2_result_y = self.form_0 * to_wedge_2_form_y
            form_2_result_z = self.form_0 * to_wedge_2_form_z
            
            # create instance and return
            new_object = form_2_3d(self.xg, self.yg, self.zg, form_2_result_x, form_2_result_y, form_2_result_z)
            return new_object
        elif order == 3:
            # from these get the numerical 0-form
            form_3_result = self.form_0 * to_wedge_3_form
            
            # depending on keep_object, return:
            if keep_object:
                self.form_0 = form_3_result
            elif not keep_object:
                new_object = form_3_3d(self.xg, self.yg, self.zg, form_3_result)
                # return the new one to the user:
                return new_object
            else:
                raise ValueError('Error, Invalid input for \'keep_object\'')

        else:
            # should never happen, but in case
            raise ValueError('Variable change during code running, look at \'order\' parameter')



class form_1_3d():

    """
    The class which creates a differential 1 form object. 
    Requires grid (3D), x y and z components.
    .plot() can be used to plot out the 1 form.
    

    Methods: (applied to the object)

    .give_eqn(equation_str_x, equation_str_y, equation_str_z) 
    .log_scaling()
    .plot()
    """

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
        
        Only creates 2 axis with same number of points each
        cannot be used for any custom grids
        
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



    def contravariant(self, g=[['1', '0', '0'], ['0', '1', '0'], ['0', '0', '1']]):
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


                # evaluate the components numerically, inputting them into a
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
                if self.form_1_str_x == None and self.form_1_str_y == None and self.form_1_str_z == None:
                    print('You supplied the metric as equations (or it was default), but did not give VF equations, therefore only numericals will be completed')
                    analytics = False
                else:
                    x_str_form = '(' + self.form_1_str_x + ')*(' + g[0][0] + ') + (' + self.form_1_str_y + ')*(' + g[0][1] + ') + (' + self.form_1_str_z + ')*(' + g[0][2] + ')'
                    y_str_form = '(' + self.form_1_str_x + ')*(' + g[1][0] + ') + (' + self.form_1_str_y + ')*(' + g[1][1] + ') + (' + self.form_1_str_z + ')*(' + g[1][2] + ')'
                    z_str_form = '(' + self.form_1_str_x + ')*(' + g[2][0] + ') + (' + self.form_1_str_y + ')*(' + g[2][1] + ') + (' + self.form_1_str_z + ')*(' + g[2][2] + ')'
                    # simplify them
                    x_str_form = str(simplify(x_str_form))
                    y_str_form = str(simplify(y_str_form))
                    z_str_form = str(simplify(z_str_form))
            else:
                pass

            # based on what was given into the Vector field
            # return a 1-form object with these parameters
            if analytics:
                result_form = vector_field3(self.xg, self.yg, self.zg, form_x, form_y, form_z, x_str_form, y_str_form, z_str_form)
            elif not analytics:
                result_form = vector_field3(self.xg, self.yg, self.zg, form_x, form_y, form_z)
        
            # return the found object
            return result_form



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
        sng_size = 0.5

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
        v = mlab.figure(bgcolor=(1,1,1), fgcolor=(0,0,0))

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
        cyl1 = tvtk.ConeSource(radius = tp_wdth,
                                height = tp_hgth,
                                capping = False,
                                center = (0.075, 0, 0),
                                
                                    )

        box = tvtk.CubeSource(x_length=0.01,
                            y_length = side,
                            z_length = side)

        box1 = tvtk.CubeSource(x_length=0.01,
                            y_length = side,
                            z_length = side,
                            center = (-0.04, 0, 0))

        box2 = tvtk.CubeSource(x_length=0.01,
                            y_length = side,
                            z_length = side,
                            center = (-0.08, 0, 0))

        box3 = tvtk.CubeSource(x_length=0.01,
                            y_length = side,
                            z_length = side,
                            center = (-0.12, 0, 0))

        box4 = tvtk.CubeSource(x_length=0.01,
                            y_length = side,
                            z_length = side,
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



        # axes boundaries
        xmin = int((np.min(self.xg)) - 1)
        ymin = int((np.min(self.yg)) - 1)
        zmin = int((np.min(self.zg)) - 1)
        xmax = int((np.max(self.xg)) + 1)
        ymax = int((np.max(self.yg)) + 1)
        zmax = int((np.max(self.zg)) + 1)

        # plot red spheres at the singularities
        # add axes and outline
        mlab.points3d(pts_nan[:,0],
                    pts_nan[:,1],
                    pts_nan[:,2], color = (1,0,0),scale_factor=sng_size, resolution=36)
        mlab.outline(extent = [xmin,xmax,ymin,ymax,zmin,zmax], line_width=1.0)
        mlab.axes(color = (0,0,0), nb_labels = 5, extent = [xmin,xmax,ymin,ymax,zmin,zmax], line_width=3.0)


        # show the plot
        mlab.show()



    def ext_d(self):
        '''
        
        ext_d()
        
        Computes the exterior derivative and returns it
        as the 2-form object
        '''
        if self.form_1_str_x == None or self.form_1_str_y == None or self.form_1_str_z == None:
                # ERROR
                raise ValueError('Error: You need to supply the 1-form equations to do this, look at \'give_eqn\' method')
        else:
            # the strings have been correctly given, compute the
            # exterior derivative
            # get the inpus from fields of x and u components
            x_comp_str = self.form_1_str_x
            y_comp_str = self.form_1_str_y
            z_comp_str = self.form_1_str_z
            # from found u and v in the interior derivative, set up sympy components
            sympy_expr_x = parse_expr(x_comp_str, evaluate=False)
            sympy_expr_y = parse_expr(y_comp_str, evaluate=False)
            sympy_expr_z = parse_expr(z_comp_str, evaluate=False)
            # combine the 2 into a list:
            
            coords = ['x', 'y', 'z']
            # from these, find the derivatives
            ddx_Fx = str(diff(sympy_expr_x, coords[0]))
            ddy_Fx = str(diff(sympy_expr_x, coords[1]))
            ddz_Fx = str(diff(sympy_expr_x, coords[2]))

            ddx_Fy = str(diff(sympy_expr_y, coords[0]))
            ddy_Fy = str(diff(sympy_expr_y, coords[1]))
            ddz_Fy = str(diff(sympy_expr_y, coords[2]))

            ddx_Fz = str(diff(sympy_expr_z, coords[0]))
            ddy_Fz = str(diff(sympy_expr_z, coords[1]))
            ddz_Fz = str(diff(sympy_expr_z, coords[2]))

            Ex_deriv_x_str = str(simplify('('+ddy_Fz+')-'+'('+ddz_Fy+')'))
            Ex_deriv_y_str = str(simplify('('+ddz_Fx+')-'+'('+ddx_Fz+')'))
            Ex_deriv_z_str = str(simplify('('+ddx_Fy+')-'+'('+ddy_Fx+')'))
            # need to uspply these unformatted, so save those:
            form_2_x_unformated, form_2_y_unformated, form_2_z_unformated = Ex_deriv_x_str*1, Ex_deriv_y_str*1, Ex_deriv_z_str*1
            # from these strings, get the numerical 1-form:
            form_2_x_str = Ex_deriv_x_str.replace('x', '(self.xg)')
            form_2_x_str = form_2_x_str.replace('y', '(self.yg)')
            form_2_x_str = form_2_x_str.replace('z', '(self.zg)')

            form_2_y_str = Ex_deriv_y_str.replace('x', '(self.xg)')
            form_2_y_str = form_2_y_str.replace('y', '(self.yg)')
            form_2_y_str = form_2_y_str.replace('z', '(self.zg)')

            form_2_z_str = Ex_deriv_z_str.replace('x', '(self.xg)')
            form_2_z_str = form_2_z_str.replace('y', '(self.yg)')
            form_2_z_str = form_2_z_str.replace('z', '(self.zg)')


            if form_2_x_str.find('x') & form_2_x_str.find('y') & form_2_x_str.find('z') == -1:
                form_2_x_str = '(' + str(form_2_x_str) + ')* np.ones(np.shape(self.xg))'
            if form_2_y_str.find('x') & form_2_y_str.find('y') & form_2_y_str.find('z') == -1:
                form_2_y_str = '(' + str(form_2_y_str) + ')* np.ones(np.shape(self.yg))'
            if form_2_z_str.find('x') & form_2_z_str.find('y') & form_2_z_str.find('z') == -1:
                form_2_z_str = '(' + str(form_2_z_str) + ')* np.ones(np.shape(self.yg))'


            form_2_x = eval(form_2_x_str)
            form_2_y = eval(form_2_y_str)
            form_2_z = eval(form_2_z_str)
            
            # supply these to the 1-form object function and return object
            result_2_form = form_2_3d(self.xg, self.yg, self.zg, form_2_x, form_2_y, form_2_z, form_2_x_unformated, form_2_y_unformated, form_2_z_unformated)

            return result_2_form



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
        dx = self.xg[0, :, 0]
        dy = self.yg[:, 0, 0]
        dz = self.zg[0, 0, :]

        


        fx = self.F_x 
        fy = self.F_y
        fz = self.F_z
        
        
        
        # Calculate deirvatvies as needed, using numpy gradient.
        _, dy_Fx, _= np.gradient(fx, dx, dy, dz)
        _, _, dz_Fx= np.gradient(fx, dx, dy, dz)

        dx_Fy, _, _ = np.gradient(fy, dx, dy, dz)
        _, _ ,dz_Fy= np.gradient(fy, dx, dy, dz)

        dx_Fz, _, _= np.gradient(fz, dx, dy, dz)
        _, dy_Fz,  _= np.gradient(fz, dx, dy, dz)

        


        # from these, get the 2-form
        form_2_result_x = dy_Fz - dz_Fy
        form_2_result_y = dz_Fx - dx_Fz
        form_2_result_z = dx_Fy - dy_Fx

        """print(form_2_result_x[0,0])
        print(form_2_result_y[0,0])
        print(form_2_result_z[0,0])"""
        


        # return 2-form object to user
        result_form = form_2_3d(self.xg, self.yg, self.zg, Fx = form_2_result_x, Fy = form_2_result_y, Fz = form_2_result_z)
        
        # return it to the user
        return result_form



    def hodge(self, keep_object=False):
        '''
        
        hodge(keep_object=False)
        
        Parameters:
        -------------
        keep_object - determines if the result should be returned as a new
                      1-form or if current one need to be changed.
                      Default is False. When False, a new object is created
                      When true, the acted on is modified.
        
        It calulates the Hodge on R^3 by the standard definition:
        dx -> dy/\dz ; dy -> -dx/\dz ; dz -> dx/\dy
        Does so analytically using the equations provided in the instance
        
        returns: 1-form if keep_object is False, else returns nothing
        '''
        
        # check for equations:
        if self.form_1_str_x == None or self.form_1_str_y == None or self.form_1_str_z == None:
            # ERROR
            raise TypeError('Error: You need to supply the 1-form equations to do this, look at \'give_eqn\' method')
        else:
            # some equations are there, compute the Hodge on these:
            new_str_x = self.form_1_str_x
            new_str_y = '-(' + self.form_1_str_y + ')'
            new_str_z = self.form_1_str_z
            # from these, get numerical solutions, evaulated on local
            # strings changed to relate to the self grids
            # need to supply these unformatted, so save those:
            form_1_x_unformated, form_1_y_unformated, form_1_z_unformated = new_str_x*1, new_str_y*1, new_str_z*1
            # from these strings, get the numerical 1-form:
            new_str_x = new_str_x.replace('x', '(self.xg)')
            new_str_x = new_str_x.replace('y', '(self.yg)')
            new_str_x = new_str_x.replace('z', '(self.zg)')

            new_str_y = new_str_y.replace('x', '(self.xg)')
            new_str_y = new_str_y.replace('y', '(self.yg)')
            new_str_y = new_str_y.replace('z', '(self.zg)')

            new_str_z = new_str_z.replace('x', '(self.xg)')
            new_str_z = new_str_z.replace('y', '(self.yg)')
            new_str_z = new_str_z.replace('z', '(self.zg)')
            # correct for constants
            if new_str_x.find('x') & new_str_x.find('y') & new_str_x.find('z') == -1:
                new_str_x = '(' + str(new_str_x) + ')* np.ones(np.shape(self.xg))'
            if new_str_y.find('x') & new_str_y.find('y') & new_str_y.find('z') == -1:
                new_str_y = '(' + str(new_str_y) + ')* np.ones(np.shape(self.yg))'
            if new_str_z.find('x') & new_str_z.find('y') & new_str_z.find('z') == -1:
                new_str_z = '(' + str(new_str_z) + ')* np.ones(np.shape(self.zg))'
            
            # evaluate
            form_1_x = eval(new_str_x)
            form_1_y = eval(new_str_y)
            form_1_z = eval(new_str_z)
            
            # depending on keep_object, return:
            if keep_object:
                self.F_x = form_1_x
                self.F_y = form_1_y
                self.F_z = form_1_z
                self.form_1_str_x = form_1_x_unformated
                self.form_1_str_y = form_1_y_unformated
                self.form_1_str_z = form_1_z_unformated
            elif not keep_object:
                new_object = form_2_3d(self.xg, self.yg, self.zg, form_1_x, form_1_y, form_1_z, Fx_eqn=form_1_x_unformated, Fy_eqn=form_1_y_unformated, Fz_eqn=form_1_z_unformated)
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
        
        It calulates the Hodge on R^3 by the standard definition:
        dx -> dy/\dz ; dy -> -dx/\dz ; dz -> dx/\dy
        
        Does no numerically using only component arrays.
        If equations have been previously provided, this method will
        loose them
        
        returns: 1-form if keep_object is False, else returns nothing
        '''
        # check if equations have been given:
        # if they have, doing it only numerically would create
        # a mismatch, warn user
        if self.form_1_str_x != None or self.form_1_str_y != None or self.form_1_str_z != None:
            print('Warning: You supplied equations, doing it numerically only will result in a mismacth between numerical values and equations')
        
        # now complete the process numerically save as instructed
        # check keep_object:
        if keep_object:
            # change the object self properties accoringly
            new_x = self.F_x
            new_y = -self.F_y
            new_z = self.F_z
            self.F_x = new_x
            self.F_y = new_y
            self.F_z = new_z
        elif not keep_object:
            # pass these in to the object to create a new one:
            # N.B no equations to supply
            new_object = form_2_3d(self.xg, self.yg, self.zg, self.F_x, -self.F_y, self.F_z)
            # return the new one to the user:
            return new_object
        else:
            raise ValueError('Error, Invalid input for \'keep_object\'')



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
        if self.form_1_str_x == None or self.form_1_str_y == None or self.form_1_str_z == None:
            raise ValueError('Error: You need to supply the 1-form equations to do this, look at \'give_eqn\' method')
        
        # set up variable to store order of supplied form, initially assume 1-form
        order = 1
        
        # get needed second obejct strings dep. on input
        if isinstance(form_second, tuple):
            # if equations were given here take these, if numerical grids were given - error!
            # check size , should be a 1-form
            if len(form_second) == 3:
                # 1-form/\1-form or 1-form/\2-form check if strings supplied
                if isinstance(form_second[0], str) and isinstance(form_second[1], str) and isinstance(form_second[2], str):
                    to_wedge_x_2_str = form_second[0]
                    to_wedge_y_2_str = form_second[1]
                    to_wedge_z_2_str = form_second[2]

                    if degree == 1:
                        order = 1
                    elif degree == 2:
                        order = 2

                else:
                    raise ValueError('for analytical calulation, supply 1-form equations as strings')
            else:
                raise ValueError('too many or too little equations given in tuple')
        elif isinstance(form_second, str):
            # single string, could be 0-form or 3-form, check given degree:
            if degree == 0:
                to_wedge_0_form_str = form_second
                order = 0
            elif degree == 3:
                to_wedge_3_form_str = form_second
                order = 3
            else:
                raise ValueError('not possible digree given or supplied one string for a 1-form')
        else:
            # object supplied, get numericals checking which object is given:
            if isinstance(form_second, form_1_3d):
                if form_second.form_1_str_x is None or form_second.form_1_str_y is None or form_second.form_1_str_z is None:
                     raise ValueError('supplied 1-form instance must contain equations for analytical calculation')
                else:
                    to_wedge_x_2_str = form_second.form_1_str_x
                    to_wedge_y_2_str = form_second.form_1_str_y
                    to_wedge_z_2_str = form_second.form_1_str_z
                    order = 1
            elif isinstance(form_second, form_0_3d):
                if form_second.form_0_str is None:
                    raise ValueError('supplied 0-form instance must contain equations for analytical calculation')
                else:
                    to_wedge_0_form_str = form_second.form_0_str
                    order = 0       
            elif isinstance(form_second, form_2_3d):
                if form_second.Fx_eqn is None and form_second.Fy_eqn is None and form_second.Fz_eqn is None:
                     raise ValueError('supplied 2-form instance must contain at least one equation for analytical calculation')
                else:
                    to_wedge_x_2_str = form_second.Fx_eqn
                    to_wedge_y_2_str = form_second.Fy_eqn
                    to_wedge_z_2_str = form_second.Fz_eqn
                    order = 2
            elif isinstance(form_second, form_3_3d):
                order=None
                print('This operation makes a 4-form, which dformpy cannot handle (yet)')
            else:
                raise TypeError('Supplied form to wedge with is not recognised')
        
        # Deal with 1-form/\1-form:
        if order == 1:
            # first, mathematically:  2-form = f*m - g*h
            form_2_str_x = str(simplify( '(' + self.form_1_str_y + ')*(' +  to_wedge_z_2_str + ')' ))
            form_2_str_y = str(simplify( '(' + self.form_1_str_z + ')*(' +  to_wedge_x_2_str + ')' ))
            form_2_str_z = str(simplify( '(' + self.form_1_str_x + ')*(' +  to_wedge_y_2_str + ')' ))
            # keep it as it is locally to supply it to object maker later
            form_2_str_x_loc = form_2_str_x + ''
            form_2_str_y_loc = form_2_str_y + ''
            form_2_str_z_loc = form_2_str_z + ''
            # format it to be in terms of grids and:
            # check against constant and zero 2-forms being supplied
            # get the numerical evaluation of it
            form_2_str_x = form_2_str_x.replace('x', 'self.xg')
            form_2_str_x = form_2_str_x.replace('y', 'self.yg')
            form_2_str_x = form_2_str_x.replace('z', 'self.zg')

            form_2_str_y = form_2_str_y.replace('x', 'self.xg')
            form_2_str_y = form_2_str_y.replace('y', 'self.yg')
            form_2_str_y = form_2_str_y.replace('z', 'self.zg')

            form_2_str_z = form_2_str_z.replace('x', 'self.xg')
            form_2_str_z = form_2_str_z.replace('y', 'self.yg')
            form_2_str_z = form_2_str_z.replace('z', 'self.zg')

            if form_2_str_x.find('x') & form_2_str_x.find('y') & form_2_str_x.find('z') == -1:
                form_2_str_x = '(' + str(form_2_str_x) + ')* np.ones(np.shape(self.xg))'
            if form_2_str_y.find('x') & form_2_str_y.find('y') & form_2_str_y.find('z') == -1:
                form_2_str_y = '(' + str(form_2_str_y) + ')* np.ones(np.shape(self.yg))'
            if form_2_str_z.find('x') & form_2_str_z.find('y') & form_2_str_z.find('z') == -1:
                form_2_str_z = '(' + str(form_2_str_z) + ')* np.ones(np.shape(self.zg))'
            
            # evaluate it numerically on the grid supplied
            form_2_x_result = eval(form_2_str_x)
            form_2_y_result = eval(form_2_str_y)
            form_2_z_result = eval(form_2_str_z)
            
            # create a 2-form object from this; to return and do so
            ret_object = form_2_3d(self.xg, self.yg, self.zg, form_2_x_result, form_2_y_result, form_2_z_result, form_2_str_x_loc, form_2_str_y_loc, form_2_str_z_loc)
            return ret_object
        
        elif order == 0:
            # form-1 /\ form-0
            # first, find the result of the 1-form:
            new_str_x = str(simplify('(' + self.form_1_str_x + ')*(' +  to_wedge_0_form_str + ')'))
            new_str_y = str(simplify('(' + self.form_1_str_y + ')*(' +  to_wedge_0_form_str + ')'))
            new_str_z = str(simplify('(' + self.form_1_str_z + ')*(' +  to_wedge_0_form_str + ')'))
            # keep it as it is locally to supply it to object maker later
            form_1_str_x_loc = new_str_x + ''
            form_1_str_y_loc = new_str_y + ''
            form_1_str_z_loc = new_str_z + ''
            # format it to be in terms of grids and:
            # check against constant and zero 1-forms being supplied
            # get the numerical evaluation of it
            new_str_x = new_str_x.replace('x', '(self.xg)')
            new_str_x = new_str_x.replace('y', '(self.yg)')
            new_str_x = new_str_x.replace('z', '(self.zg)')

            new_str_y = new_str_y.replace('x', '(self.xg)')
            new_str_y = new_str_y.replace('y', '(self.yg)')
            new_str_y = new_str_y.replace('z', '(self.zg)')

            new_str_z = new_str_z.replace('x', '(self.xg)')
            new_str_z = new_str_z.replace('y', '(self.yg)')
            new_str_z = new_str_z.replace('z', '(self.zg)')
            
            if new_str_x.find('x') & new_str_x.find('y') & new_str_x.find('z') == -1:
                new_str_x = '(' + str(new_str_x) + ')* np.ones(np.shape(self.xg))'
            if new_str_y.find('x') & new_str_y.find('y') & new_str_y.find('z') == -1:
                new_str_y = '(' + str(new_str_y) + ')* np.ones(np.shape(self.yg))'
            if new_str_z.find('x') & new_str_z.find('y') & new_str_z.find('z') == -1:
                new_str_z = '(' + str(new_str_z) + ')* np.ones(np.shape(self.zg))'
            
            form_1_x = eval(new_str_x)
            form_1_y = eval(new_str_y)
            form_1_z = eval(new_str_z)
            
            # depending on keep_object, return:
            if keep_object:
                self.F_x = form_1_x
                self.F_y = form_1_y
                self.F_z = form_1_z
                self.form_1_str_x = form_1_str_x_loc
                self.form_1_str_y = form_1_str_y_loc
                self.form_1_str_z = form_1_str_z_loc
            elif not keep_object:
                new_object = form_1_3d(self.xg, self.yg, self.zg, form_1_x, form_1_y, form_1_z, F_x_eqn=form_1_str_x_loc, F_y_eqn=form_1_str_y_loc, F_z_eqn=form_1_str_z_loc)
                # return the new one to the user:
                return new_object
            else:
                raise ValueError('Error, Invalid input for \'keep_object\'')

        elif order == 2:

            # form-1 /\ form-2
            form_3_str = str(simplify( '(((' + self.form_1_str_x + ')*(' +  to_wedge_x_2_str + ')) + (('+self.form_1_str_y+')*('+ to_wedge_y_2_str+')) + (('+self.form_1_str_z+')*('+ to_wedge_z_2_str+')))'))
            
            # keep it as it is locally to supply it to object maker later
            form_3_str_loc = form_3_str + ''

            # format it to be in terms of grids and:
            # check against constant and zero 2-forms being supplied
            # get the numerical evaluation of it
            form_3_str = form_3_str.replace('x', 'self.xg')
            form_3_str = form_3_str.replace('y', 'self.yg')
            form_3_str = form_3_str.replace('z', 'self.zg')


            if form_3_str.find('x') & form_3_str.find('y') & form_3_str.find('z') == -1:
                form_3_str = '(' + str(form_3_str) + ')* np.ones(np.shape(self.xg))'

            
            # evaluate it numerically on the grid supplied
            form_3_result = eval(form_3_str)

            
            # create a 2-form object from this; to return and do so
            ret_object = form_3_3d(self.xg, self.yg, self.zg, form_3_result, form_3_str_loc)
            return ret_object

        elif order is None:
            # made a form that is always zero on R2, no need to make it
            # Warning already shown, when degree was set
            pass
        else:
            # should never happen, but in case
            raise ValueError('Variable change during code running, look at \'order\' parameter')







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
    
    # #####################################################################
    # Define basic methods to customise this object
    # #####################################################################
    
    # define a mehtod to allow user to supply the string equation
    # of the 2-form
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
    
    # deifne a function to return the string equation to the user
    def return_string(self):
        '''
        Takes in no arguments, returns the unformatted string back to user
        This is done in case user wants to access strings
        that got here not by input but by ext. alg.
        '''
        return self.Fx_eqn_str, self.Fy_eqn_str, self.Fz_eqn_str
    
    # change colour list
   
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
    
   
    
    def zoom(self, mag, target, dpd):
    


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
 
    
    # define a method to change the density of grids in same range
    # requires string input of 1-form:
    def set_density2(self, points_number_x, points_number_y, points_number_z):
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
    

    

    def plot(self):
        
  

        if self.Fx is None and self.Fy is None and self.Fz is None:
            print('Please, provide at least one component of the field')
        
        Fz = self.Fz
        Fx = self.Fx
        Fy = self.Fy

        Fmag = np.sqrt(Fx**2 + Fy**2 + Fz**2)

        if self.logarithmic_scale_bool:
            mag1 = Fmag + 1
            # min_size = np.min(mag1)
            
            unorm = Fx/mag1
            vnorm = Fy/mag1
            knorm = Fz/mag1
            
            # logsf = np.log10(mag1/min_size)
            logmag = np.log10(mag1)
            Fx = unorm*logmag
            Fy = vnorm*logmag
            Fz = knorm*logmag



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


            v = mlab.figure(bgcolor=(1,1,1), fgcolor=(0,0,0))


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

                v.scene.add_actor(a1)
                v.scene.add_actor(a2)
                v.scene.add_actor(a3)


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

            mlab.points3d(([xmin,ymin,zmin],[xmin,ymin,zmax],[xmin,ymax,zmin],[xmin,ymax,zmax],[xmax,ymin,zmin],[xmax,ymin,zmax],[xmax,ymax,zmin],[xmax,ymax,zmax]),opacity=0.0)
            mlab.points3d(pts_nan[:,0],
                                pts_nan[:,1],
                                pts_nan[:,2], color = (1,0,0),scale_factor=gr_sep, resolution=36)
            mlab.axes(extent = [xmin,xmax,ymin,ymax,zmin,zmax], nb_labels = 5, line_width=3.0, color = (0,0,0), xlabel=xlab, ylabel=ylab, zlabel=zlab)

            xfoc = (xmin+xmax)/2
            yfoc = (ymin+ymax)/2
            zfoc = (zmin+zmax)/2


            mlab.view(focalpoint=[xfoc,yfoc,zfoc])
            mlab.show()


        if Fz is not None and np.all(Fz) != 0:
            form_2(Fz,'z')
        
        if Fx is not None and np.all(Fx) != 0:
            form_2(Fx,'x')

        if Fy is not None and np.all(Fy) != 0:
            form_2(Fy,'y')
       
                                

class form_3_3d():

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
        self.form_3_str = equation_str
        
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
 

    def plot(self):



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


        v = mlab.figure(bgcolor=(1,1,1), fgcolor=(0,0,0))


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

            v.scene.add_actor(a1)
            v.scene.add_actor(a2)
            v.scene.add_actor(a3)


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
        

        mlab.points3d(([xmin,ymin,zmin],[xmin,ymin,zmax],[xmin,ymax,zmin],[xmin,ymax,zmax],[xmax,ymin,zmin],[xmax,ymin,zmax],[xmax,ymax,zmin],[xmax,ymax,zmax]),opacity=0.0)
        mlab.points3d(pts_nan[:,0],
                            pts_nan[:,1],
                            pts_nan[:,2], color = (1,0,0),scale_factor=gr_sep, resolution=36)
        mlab.axes(extent = [xmin,xmax,ymin,ymax,zmin,zmax], nb_labels = 5, line_width=3.0, color = (0,0,0), xlabel=xlab, ylabel=ylab, zlabel=zlab)

        xfoc = (xmin+xmax)/2
        yfoc = (ymin+ymax)/2
        zfoc = (zmin+zmax)/2


        mlab.view(focalpoint=[xfoc,yfoc,zfoc])
        mlab.show()




    








