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


        if add_curl is None:
        

            mlab.figure(bgcolor = (1,1,1), fgcolor = (0,0,0))
            mlab.quiver3d(self.xg, self.yg, self.zg, F_x_local, F_y_local, F_z_local, colormap='jet')

            mlab.outline(extent = [xmin,xmax,ymin,ymax,zmin,zmax], line_width=3.0)
            axes = mlab.axes(color = (0,0,0), nb_labels = 5, extent = [xmin,xmax,ymin,ymax,zmin,zmax], line_width=3.0)




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
        else:
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

                mlab.figure(bgcolor = (1,1,1), fgcolor = (0,0,0))
                mlab.quiver3d(self.xg, self.yg, self.zg, F_x_local, F_y_local, F_z_local, colormap='Blues')
                mlab.quiver3d(self.xg, self.yg, self.zg, Curl_X, Curl_Y, Curl_Z, color=(1.0, 0.0, 1.0))

                mlab.outline(extent = [xmin,xmax,ymin,ymax,zmin,zmax], line_width=3.0)
                axes = mlab.axes(color = (0,0,0), nb_labels = 5, extent = [xmin,xmax,ymin,ymax,zmin,zmax], line_width=3.0)




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






















