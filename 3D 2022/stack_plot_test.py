from tvtk.api import tvtk
from tvtk.common import configure_input_data, configure_source_data
from mayavi import mlab
import numpy as np

v = mlab.figure(bgcolor=(1,1,1), fgcolor=(0,0,0))

grid = np.linspace(-5,5,7)
xg, yg, zg= np.meshgrid(grid, grid, grid)
pts = np.vstack(list(zip(xg.ravel(), yg.ravel(), zg.ravel())))

Fx = -yg/np.sqrt(xg**2+yg**2-zg**2)
Fy = xg/np.sqrt(xg**2+yg**2-zg**2)
Fz = zg/np.sqrt(xg**2+yg**2-zg**2)
F_list = np.vstack(list(zip(Fx.ravel(), Fy.ravel(), Fz.ravel())))




### Getting rid of zero points

Idx = np.argwhere(np.all(F_list==0,axis=1))


pts_new = np.delete(pts, [Idx[:]], axis=0)

F_new = np.delete(F_list, [Idx[:]], axis=0)



# getting rid of singular points

F_new[np.isinf(F_new)] = np.nan

Idx_nan = np.argwhere(np.all(np.isnan(F_new),axis=1))

pts_nan = pts_new[Idx_nan]


pts_new = np.delete(pts_new, [Idx_nan], axis=0)

F_new = np.delete(F_new, [Idx_nan], axis=0)




# field_magnitude
mag = np.sqrt(F_new[:,0]**2 + F_new[:,1]**2 + F_new[:,2]**2)
mag_lst = np.vstack(list(zip(mag.ravel())))


mag_max = np.max(mag)
mag_min = np.min(mag)

sep = (mag_max-mag_min)/5

# indices where magnitude of the field is
# larger that (mag_min + (N)*sep)
Idx1 = np.argwhere(np.all(mag_lst>=(mag_min+sep),axis=1))
Idx2 = np.argwhere(np.all(mag_lst>=(mag_min+(2*sep)),axis=1))
Idx3 = np.argwhere(np.all(mag_lst>=(mag_min+(3*sep)),axis=1))
Idx4 = np.argwhere(np.all(mag_lst>=(mag_min+(4*sep)),axis=1))
Idx5 = np.argwhere(np.all(mag_lst>=(mag_min+(5*sep)),axis=1))




pts1 = pts_new[Idx1]
pts1 = np.vstack(list(zip(pts1[:,:,0].ravel(),pts1[:,:,1].ravel(),pts1[:,:,2].ravel())))
pts2 = pts_new[Idx2]
pts2 = np.vstack(list(zip(pts2[:,:,0].ravel(),pts2[:,:,1].ravel(),pts2[:,:,2].ravel())))
pts3 = pts_new[Idx3]
pts3 = np.vstack(list(zip(pts3[:,:,0].ravel(),pts3[:,:,1].ravel(),pts3[:,:,2].ravel())))
pts4 = pts_new[Idx4]
pts4 = np.vstack(list(zip(pts4[:,:,0].ravel(),pts4[:,:,1].ravel(),pts4[:,:,2].ravel())))
pts5 = pts_new[Idx5]
pts5 = np.vstack(list(zip(pts5[:,:,0].ravel(),pts5[:,:,1].ravel(),pts5[:,:,2].ravel())))


# direction vectors

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



cyl1 = tvtk.ConeSource(radius = 0.05,
                        height = 0.15,
                        capping = False,
                        center = (0.075, 0, 0)
                            )

box = tvtk.CubeSource(x_length=0.01,
                     y_length=0.5,
                     z_length = 0.5)

box1 = tvtk.CubeSource(x_length=0.01,
                     y_length=0.5,
                     z_length = 0.5,
                     center = (-0.04, 0, 0))

box2 = tvtk.CubeSource(x_length=0.01,
                     y_length=0.5,
                     z_length = 0.5,
                     center = (-0.08, 0, 0))

box3 = tvtk.CubeSource(x_length=0.01,
                     y_length=0.5,
                     z_length = 0.5,
                     center = (-0.12, 0, 0))

box4 = tvtk.CubeSource(x_length=0.01,
                     y_length=0.5,
                     z_length = 0.5,
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

g = tvtk.Glyph3D(scale_mode='data_scaling_off', vector_mode = 'use_vector')
gb = tvtk.Glyph3D(scale_mode='data_scaling_off', vector_mode = 'use_vector')
g1 = tvtk.Glyph3D(scale_mode='data_scaling_off', vector_mode = 'use_vector')
g2 = tvtk.Glyph3D(scale_mode='data_scaling_off', vector_mode = 'use_vector')
g3 = tvtk.Glyph3D(scale_mode='data_scaling_off', vector_mode = 'use_vector')
g4 = tvtk.Glyph3D(scale_mode='data_scaling_off', vector_mode = 'use_vector')

configure_input_data(g, pd)
configure_input_data(gb, pd)
configure_input_data(g1, pd1)
configure_input_data(g2, pd2)
configure_input_data(g3, pd3)
configure_input_data(g4, pd4)


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


m = tvtk.PolyDataMapper()
mb = tvtk.PolyDataMapper()
m1 = tvtk.PolyDataMapper()
m2 = tvtk.PolyDataMapper()
m3 = tvtk.PolyDataMapper()
m4 = tvtk.PolyDataMapper()


p1 = tvtk.Property(opacity=1.0, color=(1, 0, 1))

configure_input_data(m, g.output)
configure_input_data(mb, gb.output)
configure_input_data(m1, g1.output)
configure_input_data(m2, g2.output)
configure_input_data(m3, g3.output)
configure_input_data(m4, g4.output)

a = tvtk.Actor(mapper=m, property=p1)
ab = tvtk.Actor(mapper=mb, property=p1)
a1 = tvtk.Actor(mapper=m1, property=p1)
a2 = tvtk.Actor(mapper=m2, property=p1)
a3 = tvtk.Actor(mapper=m3, property=p1)
a4 = tvtk.Actor(mapper=m4, property=p1)

v.scene.add_actor(a)
v.scene.add_actor(ab)
v.scene.add_actor(a1)
v.scene.add_actor(a2)
v.scene.add_actor(a3)
v.scene.add_actor(a4)




xmin = int((np.min(xg)) - 1)
ymin = int((np.min(yg)) - 1)
zmin = int((np.min(zg)) - 1)
xmax = int((np.max(xg)) + 1)
ymax = int((np.max(yg)) + 1)
zmax = int((np.max(zg)) + 1)


mlab.points3d(pts_nan[:,:,0],
              pts_nan[:,:,1],
              pts_nan[:,:,2], color = (1,0,0),scale_factor=1.0, resolution=36)
mlab.outline(extent = [xmin,xmax,ymin,ymax,zmin,zmax], line_width=1.0)
mlab.axes(color = (0,0,0), nb_labels = 5, extent = [xmin,xmax,ymin,ymax,zmin,zmax], line_width=3.0)



mlab.show()
