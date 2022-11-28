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


        # axes boundaries
        xmin = int((np.min(self.xg)) - 1)
        ymin = int((np.min(self.yg)) - 1)
        zmin = int((np.min(self.zg)) - 1)
        xmax = int((np.max(self.xg)) + 1)
        ymax = int((np.max(self.yg)) + 1)
        zmax = int((np.max(self.zg)) + 1)
        """if (mag_max-mag_min)<=0.05:
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
        v.mlab.points3d(pts_nan[:,0],
                    pts_nan[:,1],
                    pts_nan[:,2], color = (1,0,0),scale_factor=sng_size, resolution=36)
        v.mlab.outline(extent = [xmin,xmax,ymin,ymax,zmin,zmax], line_width=1.0)
        v.mlab.axes(color = (0,0,0), nb_labels = 5, extent = [xmin,xmax,ymin,ymax,zmin,zmax], line_width=3.0)

        mlab.show()
        """

        stck_coords = [a, ab, a1, a2, a3, a4]
        red_balls_data = [pts_nan[:,0], pts_nan[:,1], pts_nan[:,2], sng_size]
        axes_limits = [xmin,xmax,ymin,ymax,zmin,zmax]

        return stck_coords, red_balls_data, axes_limits

