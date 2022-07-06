from mayavi import mlab
import numpy as np
from tvtk.api import tvtk
from tvtk.common import configure_input_data
from tvtk.tools import visual
import vpython as vp

grid = np.linspace(1,5,6)
xg, yg, zg = np.meshgrid(grid,grid,grid)

Fx = xg/np.sqrt(xg**2+yg**2+zg**2)
Fy = yg/np.sqrt(xg**2+yg**2+zg**2)
Fz = zg/np.sqrt(xg**2+yg**2+zg**2)

xi = Fx/np.sqrt(Fx**2+Fy**2+Fz**2)
yj = Fy/np.sqrt(Fx**2+Fy**2+Fz**2)
zk = Fz/np.sqrt(Fx**2+Fy**2+Fz**2)


f = mlab.figure(bgcolor=(1,1,1), fgcolor=(0,0,0))

grid_ax = np.linspace(np.min(grid), np.max(grid), 2)

x_ax, y_ax, z_ax = np.meshgrid(grid_ax,grid_ax,grid_ax)

visual.set_viewer(f)

for i in range(len(xg[0, :, :])):
    for j in range(len(yg[:, 0, :])):
        for k in range(len(zg[:, :, 0])):


            cone = visual.cone(x=xg[k,j,i]+(0.05*xi[k,j,i]),
                               y=yg[k,j,i]+(0.05*yj[k,j,i]),
                               z=zg[k,j,i]+(0.05*zk[k,j,i]), height=0.1, radius=0.05, color=(0.9,0,0.9))
            box = visual.box(x=xg[k,j,i], y=yg[k,j,i], z=zg[k,j,i], size=(0.01, 0.2,0.2),color=(0.9,0,0.9))
            cone.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])
            box.axis = (Fx[k, j, i], Fy[k, j, i], Fz[k, j, i])


mlab.points3d(x_ax,y_ax,z_ax, opacity=0.0)
mlab.axes(extent=[np.min(grid_ax)-0.5, np.max(grid_ax)+0.5,np.min(grid_ax)-0.5, np.max(grid_ax)+0.5,np.min(grid_ax)-0.5, np.max(grid_ax)+0.5],
          nb_labels = 5)
mlab.show()