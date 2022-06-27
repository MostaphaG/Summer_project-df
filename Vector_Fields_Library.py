import numpy as np
import mayavi.mlab as mlb
import matplotlib.pyplot as plt
import sympy
from sympy import sympify
from sympy import diff, simplify
from sympy.parsing.sympy_parser import parse_expr

from numpy import sin, cos, tan, sqrt, log, arctan, arcsin, arccos, tanh
from numpy import sinh, cosh, arcsinh, arccosh, arctanh, exp, pi, e



class vector_field():



    """
    Defines a vector field object

    Methods:
    --------------

    .coordinate_list_cartesian(l_lim, h_lim, spacing)

    .coordinate_list_spherical(R, Phi_sep, Theta_sep, R_sep)

    """

    # create class-wide lists into which the functions will return variables

    CartCoords = []
    SphCoords = []
    CylCoords = []
    F = []
    Curl = []
    Div = []
    
#________COORDINATE-SYSTEMS_____________________________________________________________


    def coordinate_list_cartesian(l_lim, h_lim, spacing):
        
        
        
        '''
        defines a list of coordinates in x, y and z directions

        .coordinate_list_cartesian(l_lim, h_lim, spacing)
        
        Parameters:
        --------
        l_lim - lower limit of the coordinate plane __;__ 
        h_lim - upper limit of the coordinate plane (if l_lim=-5 and h_lim is 5,
        ...the plotted figure will have x,y and z dimensions from -5 to 5 included) __;__ 
        spacing - amount of units between each point in the grid
        
        Returns:
        --------
        set of coordinate lists. coordinates[0] = list of x coords __;__ 
        coordinates[1] = list of y coords __;__ 
        coordinates[2] = list of z coords
        
        '''
        

        # define grid to be usedd for the coordinate system
        grid = np.arange(l_lim, h_lim + (spacing*0.5), spacing)

        # define the system itself i.e. meshgrid of points based on
        # grid parameter
        xg, yg, zg = np.meshgrid(grid, grid, grid)  

        # append the corresponding class-wide set of coordinates with three
        # direction based lists (x, y, z)
        vector_field.CartCoords.append(xg)
        vector_field.CartCoords.append(yg)
        vector_field.CartCoords.append(zg)

        # return the nev class-wide set to the class with new values
        return(vector_field.CartCoords)
        


    #------------------------------------------------------------------


    def coordinate_list_spherical(R, Phi_sep, Theta_sep, R_sep):
        
        
        
        '''
        defines a list of coordinates in R, Phi and Theta directions

        .coordinate_list_spherical(R, Phi_sep, Theta_sep, R_sep)
        
        Parameters:
        --------
        R - two-sided limit of the coordinate plane. Equal to radius of the spherical mesh + extra unit for visual purposes __;__
        Phi_sep - Azimuthal separation of the points. Better to choose this value as a fraction of Pi 
                  i.e. fraction of a full circle __;__
        Theta_sep - Longitudal separation of the points. number N from 0 to inf, represent 1/N fraction of np.pi
        R_sep - separation in radial direction (separation between "shells")
        
        Returns:
        --------
        set of coordinate lists. coordinates[0] = list of x coords (sphere)
        coordinates[1] = list of y coords (sphere)
        coordinates[2] = list of z coords (sphere)
        
        '''


        # define the coordinates
        R = np.arange(0, R + (R_sep*0.5), R_sep)
        Th = np.arange(0, np.pi, ((np.pi)/(np.abs(Theta_sep))))
        Phi = np.arange(0, 2*np.pi, ((2*np.pi)/(np.abs(Phi_sep))))

        # define the grid based on the coordinates
        Rg, Thg, Phig = np.meshgrid(R, Th, Phi)

        # convert the grid to cartesian points
        xg = Rg*np.cos(Phig)*np.sin(Thg)
        yg = Rg*np.sin(Phig)*np.sin(Thg)
        zg = Rg*np.cos(Thg)

        # provide the cartesian collection of points to the class-wide vector
        vector_field.SphCoords.append(xg)
        vector_field.SphCoords.append(yg)
        vector_field.SphCoords.append(zg)

        # return the vector
        return(vector_field.SphCoords)


    #------------------------------------------------------------------


    def coordinate_list_cylindrical(Rho, Rho_sep, Phi_sep, l_lim_z, h_lim_z, z_sep):
        
        
        
        '''
        defines a list of coordinates in Rho, z and Phi directions

        .coordinate_list_cylindrical(R, Phi_sep, Theta_sep, R_sep)
        
        Parameters:
        --------
        
        Rho - cross section radius of the cylinder
        Rho_sep - separation of grid points in Rho direction __;__
        Phi_sep - Azimuthal separation of the points. Fraction of Pi 
                  i.e. fraction of a full circle __;__
        l_lim_z - lower limit for z coordinate __;__
        h_lim_z - upper limit for z coordinate __;__
        z_sep - separation of grid points in z direction

        
        Returns:
        --------
        set of coordinate lists. 
        coordinates[0] = list of x coords (cylinder)
        coordinates[1] = list of y coords (cylinder)
        coordinates[2] = list of z coords (cylinder)
        
        '''
        # define the coordinates
        grid = np.arange(l_lim_z, h_lim_z + (z_sep*0.5), z_sep)
        Rho = np.arange(0, Rho + (Rho_sep*0.5), Rho_sep)
        Phi = np.arange(0, 2*np.pi, ((2*np.pi)/(np.abs(Phi_sep))))

        # define the grid based on the coordinates
        Rhog, Phig, zg= np.meshgrid(Rho, Phi, grid)

        # convert the grid to cartesian points
        xg = Rhog*np.cos(Phig)
        yg = Rhog*np.sin(Phig)

        # provide the cartesian collection of points to the class-wide vector
        vector_field.CylCoords.append(xg)
        vector_field.CylCoords.append(yg)
        vector_field.CylCoords.append(zg)

        # return the vector
        return(vector_field.CylCoords)


#_______________________________________________________________________________________


### The functions above define a set of cylindrical, spherical and cartesian
#   coordinates converted to x, y, z ; which will serve as starting points
#   for arrow vectors in the vector field plots


#________VECTOR_FIELD_OBJECT____________________________________________________________


    def vfield(Fx_str, Fy_str, Fz_str, Coord_list):      

        '''
        defines the vector field object. At each point of a pre-defined grid
        this function defines a vector in direction of inputs.

        .vfield(F_x_string, F_y_string, F_z_string, Coordinate_list)
        
        Parameters:
        --------
        
        Fx_str - string of x component of the field __;__
        Fy_str - string of y component of the field __;__
        Fz_str - string of z component of the field __;__
        Coord_list - either Cartesian, Spherical or Cylindrical points grid.
         (defined by functions of this class as well)

        
        Returns:
        --------
        A vector (three dimensions) list F which incldes Fx,
        Fy and Fz components of the field
        '''
        # before every new operation clear the vector from previous results to
        # allow the user to create more than one field objects
        vector_field.F = []

        # simplify the input expressions and convert them to string
        Fx = str(simplify(Fx_str))
        Fy = str(simplify(Fy_str))
        Fz = str(simplify(Fz_str))

        # replace coordinate variables with pre-defined sets of coordinate points
        Fx = Fx.replace('x', '(Coord_list[0])')
        Fx = Fx.replace('y', '(Coord_list[1])')
        Fx = Fx.replace('z', '(Coord_list[2])')

        Fy = Fy.replace('x', '(Coord_list[0])')
        Fy = Fy.replace('y', '(Coord_list[1])')
        Fy = Fy.replace('z', '(Coord_list[2])')

        Fz = Fz.replace('x', '(Coord_list[0])')
        Fz = Fz.replace('y', '(Coord_list[1])')
        Fz = Fz.replace('z', '(Coord_list[2])')

        # if input contains no x, y or z variables convert the input into
        # list of shape of the coordinate list containing only the input
        # scalar value
        if Fx.find('x') & Fx.find('y') & Fx.find('z') == -1:
            Fx = '(' + str(Fx) + ')* np.ones(np.shape(Coord_list[0]))'
        if Fy.find('x') & Fy.find('y') & Fy.find('z') == -1:
            Fy = '(' + str(Fy) + ')* np.ones(np.shape(Coord_list[1]))'
        if Fz.find('x') & Fz.find('y') & Fz.find('z') == -1:
            Fz = '(' + str(Fz) + ')* np.ones(np.shape(Coord_list[2]))'

        # evaluate the expressions, creating a list of vector field direction
        # values for x, y and z directions
        F_x = eval(Fx)
        F_y = eval(Fy)
        F_z = eval(Fz)

        # append the three vector components to the class-wide list
        vector_field.F.append(F_x)
        vector_field.F.append(F_y)
        vector_field.F.append(F_z)

        # return the vector
        return(vector_field.F)

### The function above creates a vector field object, which is a set of three lists
#   containing vector field component for three directions (x, y, z). Combined with
#   the coordinate set it allows to produce a vector field plot


#_______________CURL____________________________________________________________________


    def curl(Fx_str, Fy_str, Fz_str, Coord_list):
        '''
        defines curl of a vector field object. At each point of a pre-defined grid
        this function defines a vector in direction of inputs.

        .curl(Fx_str, Fy_str, Fz_str, Coordinate_list)
        
        Parameters:
        --------
        
        Fx_str - string of x component of the field __;__
        Fy_str - string of y component of the field __;__
        Fz_str - string of z component of the field __;__
        Coord_list - either Cartesian, Spherical or Cylindrical points grid.
         (defined by functions of this class as well)

        
        Returns:
        --------
        A vector (three dimensions) list Curl which incldes Curl_X,
        Curl_Y and Curl_Z components of the curl for supplied field
        '''
        # before every new operation clear the vector from previous results to
        # allow the user to create more than one field objects
        vector_field.Curl =[]

        # simplify the input expressions and convert them to string
        Fx = str(simplify(Fx_str))
        Fy = str(simplify(Fy_str))
        Fz = str(simplify(Fz_str))
        
        # parse expressions into sympy object to allow manipulations
        F_x = parse_expr(Fx)
        F_y = parse_expr(Fy)
        F_z = parse_expr(Fz)        

        # differentiate expressions w.r.t x, y, z.
        ddx_Fx = F_x.diff(sympy.symbols('x'))
        ddy_Fx = F_x.diff(sympy.symbols('y'))
        ddz_Fx = F_x.diff(sympy.symbols('z'))

        ddx_Fy = F_y.diff(sympy.symbols('x'))
        ddy_Fy = F_y.diff(sympy.symbols('y'))
        ddz_Fy = F_y.diff(sympy.symbols('z'))

        ddx_Fz = F_z.diff(sympy.symbols('x'))
        ddy_Fz = F_z.diff(sympy.symbols('y'))
        ddz_Fz = F_z.diff(sympy.symbols('z'))

        # Define a string expression for Curl of the input vector field
        CurlX = (str(simplify(ddy_Fz)) + '-' + str(simplify(ddz_Fy)))
        CurlY = (str(simplify(ddz_Fx)) + '-' + str(simplify(ddx_Fz)))
        CurlZ = (str(simplify(ddx_Fy)) + '-' + str(simplify(ddy_Fx)))
        
        # replace input coordinates with coordinate point lists
        CurlX = CurlX.replace('x', '(Coord_list[0])')
        CurlX = CurlX.replace('y', '(Coord_list[1])')
        CurlX = CurlX.replace('z', '(Coord_list[2])')

        CurlY = CurlY.replace('x', '(Coord_list[0])')
        CurlY = CurlY.replace('y', '(Coord_list[1])')
        CurlY = CurlY.replace('z', '(Coord_list[2])')

        CurlZ = CurlZ.replace('x', '(Coord_list[0])')
        CurlZ = CurlZ.replace('y', '(Coord_list[1])')
        CurlZ = CurlZ.replace('z', '(Coord_list[2])')

        # if input scalar, define a scalar list the size of a coordinate obj.
        if CurlX.find('x') & CurlX.find('y') & CurlX.find('z') == -1:
            CurlX = '(' + str(CurlX) + ')* np.ones(np.shape(Coord_list[0]))'
        if CurlY.find('x') & CurlY.find('y') & CurlY.find('z') == -1:
            CurlY = '(' + str(CurlY) + ')* np.ones(np.shape(Coord_list[1]))'
        if CurlZ.find('x') & CurlZ.find('y') & CurlZ.find('z') == -1:
            CurlZ = '(' + str(CurlZ) + ')* np.ones(np.shape(Coord_list[2]))'
        
        # evaluate the curl expression, creating list of curl values
        Curl_X = eval(CurlX)
        Curl_Y = eval(CurlY)
        Curl_Z = eval(CurlZ)

        # append the class-wide curl list with the curl value 
        # lists for each coordinate
        vector_field.Curl.append(Curl_X)
        vector_field.Curl.append(Curl_Y)
        vector_field.Curl.append(Curl_Z)

        # return the vector
        return(vector_field.Curl)

### The function above defines curl of the vector field of users choice. It
#   creates the three curl vector component lists with values for each coordinate
#   point on the grid. Combined with a set of coordinates can be used to plot the curl
#   vector of a vector field


#_______________DIVERGENCE______________________________________________________________


    def div(Fx_str, Fy_str, Fz_str, Coord_list, at_x, at_y, at_z):

        New_grid = np.arange(np.min(Coord_list[0]), np.max(Coord_list[0]), 0.01)

        xng, yng, zng = np.meshgrid(New_grid, New_grid, New_grid)

        Fx = str(simplify(Fx_str))
        Fy = str(simplify(Fy_str))
        Fz = str(simplify(Fz_str))
        
        F_x = parse_expr(Fx)
        F_y = parse_expr(Fy)
        F_z = parse_expr(Fz)        

        ddx_Fx = F_x.diff(sympy.symbols('x'))
        ddy_Fx = F_x.diff(sympy.symbols('y'))
        ddz_Fx = F_x.diff(sympy.symbols('z'))

        ddx_Fy = F_y.diff(sympy.symbols('x'))
        ddy_Fy = F_y.diff(sympy.symbols('y'))
        ddz_Fy = F_y.diff(sympy.symbols('z'))

        ddx_Fz = F_z.diff(sympy.symbols('x'))
        ddy_Fz = F_z.diff(sympy.symbols('y'))
        ddz_Fz = F_z.diff(sympy.symbols('z'))

        DivExpr = str('{'+str(simplify(ddx_Fx))+'}i + {'
                      +str(simplify(ddy_Fy))+'}j + {'
                      +str(simplify(ddz_Fz))+'}k')
                      
        print(DivExpr)

        Div = str(simplify(ddx_Fx)) + '+' + str(simplify(ddy_Fy)) + '+' + str(simplify(ddz_Fz))
        
        Div = Div.replace('x', 'xng')
        Div = Div.replace('y', 'yng')
        Div = Div.replace('z', 'zng')

        if Div.find('x') & Div.find('y') & Div.find('z') == -1:
            Div = '(' + str(Div) + ')* np.ones(np.shape(xng))'

        Div_val = eval(Div)

        at_x = np.round(at_x, 2)
        at_y = np.round(at_y, 2)
        at_z = np.round(at_z, 2)

        for i in range(np.shape(Div_val)):
            if abs(at_x-xng) & abs(at_y-yng) & abs(at_z-zng) <= 0.01:
                vector_field.Div.append(Div_val[i])

        return(vector_field.Div)
        


        


















