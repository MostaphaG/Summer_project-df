from mayavi import mlab
import numpy as np
from tvtk.api import tvtk
from tvtk.common import configure_input_data
from tvtk.tools import visual


grid = np.linspace(-3,3,6)
xg, yg, zg = np.meshgrid(grid,grid,grid)

Fx = yg
Fy = -xg
Fz = 0*zg

xi = Fx/np.sqrt(Fx**2+Fy**2+Fz**2)
yj = Fy/np.sqrt(Fx**2+Fy**2+Fz**2)
zk = Fz/np.sqrt(Fx**2+Fy**2+Fz**2)


f = mlab.figure(bgcolor=(1,1,1), fgcolor=(0,0,0))

grid_ax = np.linspace(np.min(grid), np.max(grid), 2)

x_ax, y_ax, z_ax = np.meshgrid(grid_ax,grid_ax,grid_ax)

visual.set_viewer(f)

mag_max = np.max(np.sqrt(Fx**2 + Fy**2 +Fz**2))
mag_min = np.min(np.sqrt(Fx**2 + Fy**2 +Fz**2))

sep = (mag_max-mag_min)/6


plane_thck = 0.01
plane_side = 0.25

pl = 0.1
pr = 0.05


for i in range(len(xg[0, :, :])):
    for j in range(len(yg[:, 0, :])):
        for k in range(len(zg[:, :, 0])):

            mag = np.sqrt(Fx[k, j, i]**2 + Fy[k, j, i]**2 + Fz[k, j, i]**2)
            
            if mag<=(mag_min+sep):

                cone = visual.cone(x=xg[k,j,i]+((pl/2)*xi[k,j,i]),
                                y=yg[k,j,i]+((pl/2)*yj[k,j,i]),
                                z=zg[k,j,i]+((pl/2)*zk[k,j,i]), height=pl, radius=pr, color=(0.9,0,0.9))
                cone.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])

                box = visual.box(x=xg[k,j,i], y=yg[k,j,i], z=zg[k,j,i], size=(plane_thck, plane_side, plane_side),color=(0.9,0,0.9))
                box.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])

            if mag>(mag_min+sep) and mag<=(mag_min+(2*sep)):

                cone = visual.cone(x=xg[k,j,i]+((pl/2)*xi[k,j,i]),
                                   y=yg[k,j,i]+((pl/2)*yj[k,j,i]),
                                   z=zg[k,j,i]+((pl/2)*zk[k,j,i]), height=pl, radius=pr, color=(0.9,0,0.9))
                cone.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])

                box = visual.box(x=xg[k,j,i], y=yg[k,j,i], z=zg[k,j,i], size=(plane_thck, plane_side, plane_side),color=(0.9,0,0.9))
                box.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])

                box1 = visual.box(x=xg[k,j,i]-((pl/2)*xi[k,j,i]),
                                  y=yg[k,j,i]-((pl/2)*yj[k,j,i]),
                                  z=zg[k,j,i]-((pl/2)*zk[k,j,i]), size=(plane_thck, plane_side, plane_side),color=(0.9,0,0.9))
                box1.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])

            if mag>(mag_min+(2*sep)) and mag<=(mag_min+(3*sep)):

                cone = visual.cone(x=xg[k,j,i]+((2*(pl/2))*xi[k,j,i]),
                                   y=yg[k,j,i]+((2*(pl/2))*yj[k,j,i]),
                                   z=zg[k,j,i]+((2*(pl/2))*zk[k,j,i]), height=pl, radius=pr, color=(0.9,0,0.9))
                cone.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])

                box = visual.box(x=xg[k,j,i]+(((pl/2))*xi[k,j,i]),
                                 y=yg[k,j,i]+(((pl/2))*yj[k,j,i]),
                                 z=zg[k,j,i]+(((pl/2))*zk[k,j,i]), size=(plane_thck, plane_side, plane_side),color=(0.9,0,0.9))
                box.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])

                box1 = visual.box(x=xg[k,j,i],
                                  y=yg[k,j,i],
                                  z=zg[k,j,i], size=(plane_thck, plane_side, plane_side),color=(0.9,0,0.9))
                box1.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])

                box2 = visual.box(x=xg[k,j,i]-(((pl/2))*xi[k,j,i]),
                                  y=yg[k,j,i]-(((pl/2))*yj[k,j,i]),
                                  z=zg[k,j,i]-(((pl/2))*zk[k,j,i]), size=(plane_thck, plane_side, plane_side),color=(0.9,0,0.9))
                box2.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])

            if mag>(mag_min+(3*sep)) and mag<=(mag_min+(4*sep)):

                cone = visual.cone(x=xg[k,j,i]+((2*(pl/2))*xi[k,j,i]),
                                   y=yg[k,j,i]+((2*(pl/2))*yj[k,j,i]),
                                   z=zg[k,j,i]+((2*(pl/2))*zk[k,j,i]), height=pl, radius=pr, color=(0.9,0,0.9))
                cone.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])

                box = visual.box(x=xg[k,j,i]+((pl/2)*xi[k,j,i]),
                                 y=yg[k,j,i]+((pl/2)*yj[k,j,i]),
                                 z=zg[k,j,i]+((pl/2)*zk[k,j,i]), size=(plane_thck, plane_side, plane_side),color=(0.9,0,0.9))
                box.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])

                box1 = visual.box(x=xg[k,j,i],
                                  y=yg[k,j,i],
                                  z=zg[k,j,i], size=(plane_thck, plane_side, plane_side),color=(0.9,0,0.9))
                box1.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])

                box2 = visual.box(x=xg[k,j,i]-((pl/2)*xi[k,j,i]),
                                  y=yg[k,j,i]-((pl/2)*yj[k,j,i]),
                                  z=zg[k,j,i]-((pl/2)*zk[k,j,i]), size=(plane_thck, plane_side, plane_side),color=(0.9,0,0.9))
                box2.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])

                box3 = visual.box(x=xg[k,j,i]-((2*(pl/2))*xi[k,j,i]),
                                  y=yg[k,j,i]-((2*(pl/2))*yj[k,j,i]),
                                  z=zg[k,j,i]-((2*(pl/2))*zk[k,j,i]), size=(plane_thck, plane_side, plane_side),color=(0.9,0,0.9))
                box3.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])

            if mag>(mag_min+(4*sep)) and mag<=(mag_min+(5*sep)):

                cone = visual.cone(x=xg[k,j,i]+((3*(pl/2))*xi[k,j,i]),
                                   y=yg[k,j,i]+((3*(pl/2))*yj[k,j,i]),
                                   z=zg[k,j,i]+((3*(pl/2))*zk[k,j,i]), height=pl, radius=pr, color=(0.9,0,0.9))
                cone.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])

                box = visual.box(x=xg[k,j,i]+((2*(pl/2))*xi[k,j,i]),
                                 y=yg[k,j,i]+((2*(pl/2))*yj[k,j,i]),
                                 z=zg[k,j,i]+((2*(pl/2))*zk[k,j,i]), size=(plane_thck, plane_side, plane_side),color=(0.9,0,0.9))
                box.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])

                box1 = visual.box(x=xg[k,j,i]+((pl/2)*xi[k,j,i]),
                                  y=yg[k,j,i]+((pl/2)*yj[k,j,i]),
                                  z=zg[k,j,i]+((pl/2)*zk[k,j,i]), size=(plane_thck, plane_side, plane_side),color=(0.9,0,0.9))
                box1.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])

                box2 = visual.box(x=xg[k,j,i],
                                  y=yg[k,j,i],
                                  z=zg[k,j,i], size=(plane_thck, plane_side, plane_side),color=(0.9,0,0.9))
                box2.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])

                box3 = visual.box(x=xg[k,j,i]-((pl/2)*xi[k,j,i]),
                                  y=yg[k,j,i]-((pl/2)*yj[k,j,i]),
                                  z=zg[k,j,i]-((pl/2)*zk[k,j,i]), size=(plane_thck, plane_side, plane_side),color=(0.9,0,0.9))
                box3.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])

                box4 = visual.box(x=xg[k,j,i]-((2*(pl/2))*xi[k,j,i]),
                                  y=yg[k,j,i]-((2*(pl/2))*yj[k,j,i]),
                                  z=zg[k,j,i]-((2*(pl/2))*zk[k,j,i]), size=(plane_thck, plane_side, plane_side),color=(0.9,0,0.9))
                box4.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])

            if mag>(mag_min+(5*sep)):
                cone = visual.cone(x=xg[k,j,i]+((3*(pl/2))*xi[k,j,i]),
                                   y=yg[k,j,i]+((3*(pl/2))*yj[k,j,i]),
                                   z=zg[k,j,i]+((3*(pl/2))*zk[k,j,i]), height=pl, radius=pr, color=(0.9,0,0.9))
                cone.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])

                box = visual.box(x=xg[k,j,i]+((2*(pl/2))*xi[k,j,i]),
                                 y=yg[k,j,i]+((2*(pl/2))*yj[k,j,i]),
                                 z=zg[k,j,i]+((2*(pl/2))*zk[k,j,i]), size=(plane_thck, plane_side, plane_side),color=(0.9,0,0.9))
                box.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])

                box1 = visual.box(x=xg[k,j,i]+((pl/2)*xi[k,j,i]),
                                  y=yg[k,j,i]+((pl/2)*yj[k,j,i]),
                                  z=zg[k,j,i]+((pl/2)*zk[k,j,i]), size=(plane_thck, plane_side, plane_side),color=(0.9,0,0.9))
                box1.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])

                box2 = visual.box(x=xg[k,j,i],
                                  y=yg[k,j,i],
                                  z=zg[k,j,i], size=(plane_thck, plane_side, plane_side),color=(0.9,0,0.9))
                box2.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])

                box3 = visual.box(x=xg[k,j,i]-((pl/2)*xi[k,j,i]),
                                  y=yg[k,j,i]-((pl/2)*yj[k,j,i]),
                                  z=zg[k,j,i]-((pl/2)*zk[k,j,i]), size=(plane_thck, plane_side, plane_side),color=(0.9,0,0.9))
                box3.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])

                box4 = visual.box(x=xg[k,j,i]-((2*(pl/2))*xi[k,j,i]),
                                  y=yg[k,j,i]-((2*(pl/2))*yj[k,j,i]),
                                  z=zg[k,j,i]-((2*(pl/2))*zk[k,j,i]), size=(plane_thck, plane_side, plane_side),color=(0.9,0,0.9))
                box4.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])

                box5 = visual.box(x=xg[k,j,i]-((3*(pl/2))*xi[k,j,i]),
                                  y=yg[k,j,i]-((3*(pl/2))*yj[k,j,i]),
                                  z=zg[k,j,i]-((3*(pl/2))*zk[k,j,i]), size=(plane_thck, plane_side, plane_side),color=(0.9,0,0.9))
                box5.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])




mlab.points3d(x_ax,y_ax,z_ax, opacity=0.0)
mlab.axes(extent=[np.min(grid_ax)-0.5, np.max(grid_ax)+0.5,np.min(grid_ax)-0.5, np.max(grid_ax)+0.5,np.min(grid_ax)-0.5, np.max(grid_ax)+0.5],
          nb_labels = 5)
mlab.show()














