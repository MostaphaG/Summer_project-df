from tvtk.api import tvtk
import numpy as np
from mayavi import mlab
import dformpy3D as df3
 


grid = np.linspace(-10, 10, 8)
xg, yg, zg = np.meshgrid(grid, grid, grid)

grid_sep = abs(grid[1]-grid[0])

pts = np.vstack(list(zip(xg.ravel(), yg.ravel(), zg.ravel())))

Fx = 0*xg
Fy = 0*yg
Fz = zg

zmag_lst = np.vstack(list(zip(Fz.ravel())))
z_max = np.nanmax(abs(Fz))
sep = z_max/5



### Find points where field has different magnitudes

Idx1 = np.argwhere(np.all(zmag_lst>(0),axis=1) & np.all(zmag_lst<=(0+sep),axis=1))
Idx2 = np.argwhere(np.all(zmag_lst>(0+sep),axis=1) & np.all(zmag_lst<=(0+2*sep),axis=1))
Idx3 = np.argwhere(np.all(zmag_lst>(0+2*sep),axis=1) & np.all(zmag_lst<=(0+3*sep),axis=1))
Idx4 = np.argwhere(np.all(zmag_lst>(0+3*sep),axis=1) & np.all(zmag_lst<=(0+4*sep),axis=1))
Idx5 = np.argwhere(np.all(zmag_lst>(0+4*sep),axis=1))

Idx_1 = np.argwhere(np.all(zmag_lst<(0),axis=1) & np.all(zmag_lst>=(0-sep),axis=1))
Idx_2 = np.argwhere(np.all(zmag_lst<(0-sep),axis=1) & np.all(zmag_lst>=(0-2*sep),axis=1))
Idx_3 = np.argwhere(np.all(zmag_lst<(0-2*sep),axis=1) & np.all(zmag_lst>=(0-3*sep),axis=1))
Idx_4 = np.argwhere(np.all(zmag_lst<(0-3*sep),axis=1) & np.all(zmag_lst>=(0-4*sep),axis=1))
Idx_5 = np.argwhere(np.all(zmag_lst<(0-4*sep),axis=1))


x1 = pts[Idx1,0]
x_1  = pts[Idx_1,0]
x2 = pts[Idx2,0]
x_2  = pts[Idx_2,0]
x3 = pts[Idx3,0]
x_3  = pts[Idx_3,0]
x4 = pts[Idx4,0]
x_4  = pts[Idx_4,0]
x5 = pts[Idx5,0]
x_5  = pts[Idx_5,0]

y1 = pts[Idx1,1]
y_1  = pts[Idx_1,1]
y2 = pts[Idx2,1]
y_2  = pts[Idx_2,1]
y3 = pts[Idx3,1]
y_3  = pts[Idx_3,1]
y4 = pts[Idx4,1]
y_4  = pts[Idx_4,1]
y5 = pts[Idx5,1]
y_5  = pts[Idx_5,1]

z1 = pts[Idx1,2]
z_1  = pts[Idx_1,2]
z2 = pts[Idx2,2]
z_2  = pts[Idx_2,2]
z3 = pts[Idx3,2]
z_3  = pts[Idx_3,2]
z4 = pts[Idx4,2]
z_4  = pts[Idx_4,2]
z5 = pts[Idx5,2]
z_5  = pts[Idx_5,2]


f = mlab.figure(bgcolor=(1,1,1), fgcolor=(0,0,0))

for i in range(len(x1)):
    x_arr = []
    y_arr = []
    z_arr = []

    x_arr.append(float(x1[i]-(grid_sep/2)))
    x_arr.append(float(x1[i]+(grid_sep/2)))
    y_arr.append(float(y1[i]-(grid_sep/2)))
    y_arr.append(float(y1[i]+(grid_sep/2)))
    z_arr.append(float(z1[i]-(grid_sep/2)))
    z_arr.append(float(z1[i]+(grid_sep/2)))

    x_lst =[]
    x_lst.append(x_arr)

    y_lst =[]
    y_lst.append(y_arr)

    z_lst =[]
    z_lst.append(z_arr)

    for j in range(len(x_lst)):
        r = tvtk.RectilinearGrid()
        r.x_coordinates = x_lst[j]
        r.y_coordinates = y_lst[j]
        r.z_coordinates = z_lst[j]
        r.dimensions = len(x_lst[j]),len(y_lst[j]),len(z_lst[j])


        surf = mlab.pipeline.surface(r, opacity=0)
        mlab.pipeline.surface(mlab.pipeline.extract_edges(surf),
                                color=(1, 0, 0), opacity=0.9)



for i in range(len(x_1)):
    x_arr = []
    y_arr = []
    z_arr = []

    x_arr.append(float(x_1[i]-(grid_sep/2)))
    x_arr.append(float(x_1[i]+(grid_sep/2)))
    y_arr.append(float(y_1[i]-(grid_sep/2)))
    y_arr.append(float(y_1[i]+(grid_sep/2)))
    z_arr.append(float(z_1[i]-(grid_sep/2)))
    z_arr.append(float(z_1[i]+(grid_sep/2)))

    x_lst =[]
    x_lst.append(x_arr)

    y_lst =[]
    y_lst.append(y_arr)

    z_lst =[]
    z_lst.append(z_arr)

    for j in range(len(x_lst)):
        r = tvtk.RectilinearGrid()
        r.x_coordinates = x_lst[j]
        r.y_coordinates = y_lst[j]
        r.z_coordinates = z_lst[j]
        r.dimensions = len(x_lst[j]),len(y_lst[j]),len(z_lst[j])


        surf = mlab.pipeline.surface(r, opacity=0)
        mlab.pipeline.surface(mlab.pipeline.extract_edges(surf),
                                color=(0, 0, 1), opacity=0.9)




for i in range(len(x2)):
    x_arr = []
    y_arr = []
    z_arr = []

    x_arr.append(float(x2[i]-(grid_sep/2)))
    x_arr.append(float(x2[i]))
    x_arr.append(float(x2[i]+(grid_sep/2)))
    y_arr.append(float(y2[i]-(grid_sep/2)))
    y_arr.append(float(y2[i]))
    y_arr.append(float(y2[i]+(grid_sep/2)))
    z_arr.append(float(z2[i]-(grid_sep/2)))

    z_arr.append(float(z2[i]+(grid_sep/2)))

    x_lst =[]
    x_lst.append(x_arr)

    y_lst =[]
    y_lst.append(y_arr)

    z_lst =[]
    z_lst.append(z_arr)

    for j in range(len(x_lst)):
        r = tvtk.RectilinearGrid()
        r.x_coordinates = x_lst[j]
        r.y_coordinates = y_lst[j]
        r.z_coordinates = z_lst[j]
        r.dimensions = len(x_lst[j]),len(y_lst[j]),len(z_lst[j])


        surf = mlab.pipeline.surface(r, opacity=0)
        mlab.pipeline.surface(mlab.pipeline.extract_edges(surf),
                                color=(1, 0, 0), opacity=0.9)



for i in range(len(x_2)):
    x_arr = []
    y_arr = []
    z_arr = []

    x_arr.append(float(x_2[i]-(grid_sep/2)))
    x_arr.append(float(x_2[i]))
    x_arr.append(float(x_2[i]+(grid_sep/2)))
    y_arr.append(float(y_2[i]-(grid_sep/2)))
    y_arr.append(float(y_2[i]))
    y_arr.append(float(y_2[i]+(grid_sep/2)))
    z_arr.append(float(z_2[i]-(grid_sep/2)))

    z_arr.append(float(z_2[i]+(grid_sep/2)))

    x_lst =[]
    x_lst.append(x_arr)

    y_lst =[]
    y_lst.append(y_arr)

    z_lst =[]
    z_lst.append(z_arr)

    for j in range(len(x_lst)):
        r = tvtk.RectilinearGrid()
        r.x_coordinates = x_lst[j]
        r.y_coordinates = y_lst[j]
        r.z_coordinates = z_lst[j]
        r.dimensions = len(x_lst[j]),len(y_lst[j]),len(z_lst[j])


        surf = mlab.pipeline.surface(r, opacity=0)
        mlab.pipeline.surface(mlab.pipeline.extract_edges(surf),
                                color=(0, 0, 1), opacity=0.9)




for i in range(len(x3)):
    x_arr = []
    y_arr = []
    z_arr = []

    x_arr.append(float(x3[i]-(grid_sep/2)))
    x_arr.append(float(x3[i]+(grid_sep/6)))
    x_arr.append(float(x3[i]-(grid_sep/6)))
    x_arr.append(float(x3[i]+(grid_sep/2)))
    y_arr.append(float(y3[i]-(grid_sep/2)))
    y_arr.append(float(y3[i]+(grid_sep/6)))
    y_arr.append(float(y3[i]-(grid_sep/6)))
    y_arr.append(float(y3[i]+(grid_sep/2)))
    z_arr.append(float(z3[i]-(grid_sep/2)))
    z_arr.append(float(z3[i]+(grid_sep/2)))

    x_lst =[]
    x_lst.append(x_arr)

    y_lst =[]
    y_lst.append(y_arr)

    z_lst =[]
    z_lst.append(z_arr)

    for j in range(len(x_lst)):
        r = tvtk.RectilinearGrid()
        r.x_coordinates = x_lst[j]
        r.y_coordinates = y_lst[j]
        r.z_coordinates = z_lst[j]
        r.dimensions = len(x_lst[j]),len(y_lst[j]),len(z_lst[j])


        surf = mlab.pipeline.surface(r, opacity=0)
        mlab.pipeline.surface(mlab.pipeline.extract_edges(surf),
                                color=(1, 0, 0), opacity=0.9)



for i in range(len(x_3)):
    x_arr = []
    y_arr = []
    z_arr = []

    x_arr.append(float(x_3[i]-(grid_sep/2)))
    x_arr.append(float(x_3[i]+(grid_sep/6)))
    x_arr.append(float(x_3[i]-(grid_sep/6)))
    x_arr.append(float(x_3[i]+(grid_sep/2)))
    y_arr.append(float(y_3[i]-(grid_sep/2)))
    y_arr.append(float(y_3[i]+(grid_sep/6)))
    y_arr.append(float(y_3[i]-(grid_sep/6)))
    y_arr.append(float(y_3[i]+(grid_sep/2)))
    z_arr.append(float(z_3[i]-(grid_sep/2)))
    z_arr.append(float(z_3[i]+(grid_sep/2)))


    x_lst =[]
    x_lst.append(x_arr)

    y_lst =[]
    y_lst.append(y_arr)

    z_lst =[]
    z_lst.append(z_arr)

    for j in range(len(x_lst)):
        r = tvtk.RectilinearGrid()
        r.x_coordinates = x_lst[j]
        r.y_coordinates = y_lst[j]
        r.z_coordinates = z_lst[j]
        r.dimensions = len(x_lst[j]),len(y_lst[j]),len(z_lst[j])


        surf = mlab.pipeline.surface(r, opacity=0)
        mlab.pipeline.surface(mlab.pipeline.extract_edges(surf),
                                color=(0, 0, 1), opacity=0.9)




for i in range(len(x4)):
    x_arr = []
    y_arr = []
    z_arr = []

    x_arr.append(float(x4[i]-(grid_sep/2)))
    x_arr.append(float(x4[i]-(grid_sep/4)))
    x_arr.append(float(x4[i]))
    x_arr.append(float(x4[i]+(grid_sep/4)))
    x_arr.append(float(x4[i]+(grid_sep/2)))

    y_arr.append(float(y4[i]-(grid_sep/2)))
    y_arr.append(float(y4[i]+(grid_sep/4)))
    y_arr.append(float(y4[i]))
    y_arr.append(float(y4[i]-(grid_sep/4)))
    y_arr.append(float(y4[i]+(grid_sep/2)))

    z_arr.append(float(z4[i]-(grid_sep/2)))
    z_arr.append(float(z4[i]+(grid_sep/2)))

    x_lst =[]
    x_lst.append(x_arr)

    y_lst =[]
    y_lst.append(y_arr)

    z_lst =[]
    z_lst.append(z_arr)

    for j in range(len(x_lst)):
        r = tvtk.RectilinearGrid()
        r.x_coordinates = x_lst[j]
        r.y_coordinates = y_lst[j]
        r.z_coordinates = z_lst[j]
        r.dimensions = len(x_lst[j]),len(y_lst[j]),len(z_lst[j])


        surf = mlab.pipeline.surface(r, opacity=0)
        mlab.pipeline.surface(mlab.pipeline.extract_edges(surf),
                                color=(1, 0, 0), opacity=0.9)



for i in range(len(x_4)):
    x_arr = []
    y_arr = []
    z_arr = []

    x_arr.append(float(x_4[i]-(grid_sep/2)))
    x_arr.append(float(x_4[i]+(grid_sep/4)))
    x_arr.append(float(x_4[i]))
    x_arr.append(float(x_4[i]-(grid_sep/4)))
    x_arr.append(float(x_4[i]+(grid_sep/2)))
    y_arr.append(float(y_4[i]-(grid_sep/2)))
    y_arr.append(float(y_4[i]+(grid_sep/4)))
    y_arr.append(float(y_4[i]))
    y_arr.append(float(y_4[i]-(grid_sep/4)))
    y_arr.append(float(y_4[i]+(grid_sep/2)))
    z_arr.append(float(z_4[i]-(grid_sep/2)))
    z_arr.append(float(z_4[i]+(grid_sep/2)))


    x_lst =[]
    x_lst.append(x_arr)

    y_lst =[]
    y_lst.append(y_arr)

    z_lst =[]
    z_lst.append(z_arr)

    for j in range(len(x_lst)):
        r = tvtk.RectilinearGrid()
        r.x_coordinates = x_lst[j]
        r.y_coordinates = y_lst[j]
        r.z_coordinates = z_lst[j]
        r.dimensions = len(x_lst[j]),len(y_lst[j]),len(z_lst[j])


        surf = mlab.pipeline.surface(r, opacity=0)
        mlab.pipeline.surface(mlab.pipeline.extract_edges(surf),
                                color=(0, 0, 1), opacity=0.9)





for i in range(len(x5)):
    x_arr = []
    y_arr = []
    z_arr = []

    x_arr.append(float(x5[i]-(grid_sep/2)))
    x_arr.append(float(x5[i]-(grid_sep/4)))
    x_arr.append(float(x5[i]-(grid_sep/12))) 
    x_arr.append(float(x5[i]+(grid_sep/12)))
    x_arr.append(float(x5[i]+(grid_sep/4)))
    x_arr.append(float(x5[i]+(grid_sep/2)))

    y_arr.append(float(y5[i]-(grid_sep/2)))
    y_arr.append(float(y5[i]-(grid_sep/4)))
    y_arr.append(float(y5[i]-(grid_sep/12))) 
    y_arr.append(float(y5[i]+(grid_sep/12)))
    y_arr.append(float(y5[i]+(grid_sep/4)))
    y_arr.append(float(y5[i]+(grid_sep/2)))

    z_arr.append(float(z5[i]-(grid_sep/2)))
    z_arr.append(float(z5[i]+(grid_sep/2)))

    x_lst =[]
    x_lst.append(x_arr)

    y_lst =[]
    y_lst.append(y_arr)

    z_lst =[]
    z_lst.append(z_arr)

    for j in range(len(x_lst)):
        r = tvtk.RectilinearGrid()
        r.x_coordinates = x_lst[j]
        r.y_coordinates = y_lst[j]
        r.z_coordinates = z_lst[j]
        r.dimensions = len(x_lst[j]),len(y_lst[j]),len(z_lst[j])


        surf = mlab.pipeline.surface(r, opacity=0)
        mlab.pipeline.surface(mlab.pipeline.extract_edges(surf),
                                color=(1, 0, 0), opacity=0.9)



for i in range(len(x_5)):
    x_arr = []
    y_arr = []
    z_arr = []

    x_arr.append(float(x_5[i]-(grid_sep/2)))
    x_arr.append(float(x_5[i]-(grid_sep/4)))
    x_arr.append(float(x_5[i]-(grid_sep/12))) 
    x_arr.append(float(x_5[i]+(grid_sep/12)))
    x_arr.append(float(x_5[i]+(grid_sep/4)))
    x_arr.append(float(x_5[i]+(grid_sep/2)))

    y_arr.append(float(y_5[i]-(grid_sep/2)))
    y_arr.append(float(y_5[i]-(grid_sep/4)))
    y_arr.append(float(y_5[i]-(grid_sep/12))) 
    y_arr.append(float(y_5[i]+(grid_sep/12)))
    y_arr.append(float(y_5[i]+(grid_sep/4)))
    y_arr.append(float(y_5[i]+(grid_sep/2)))

    z_arr.append(float(z_5[i]-(grid_sep/2)))
    z_arr.append(float(z_5[i]+(grid_sep/2)))


    x_lst =[]
    x_lst.append(x_arr)

    y_lst =[]
    y_lst.append(y_arr)

    z_lst =[]
    z_lst.append(z_arr)

    for j in range(len(x_lst)):
        r = tvtk.RectilinearGrid()
        r.x_coordinates = x_lst[j]
        r.y_coordinates = y_lst[j]
        r.z_coordinates = z_lst[j]
        r.dimensions = len(x_lst[j]),len(y_lst[j]),len(z_lst[j])


        surf = mlab.pipeline.surface(r, opacity=0)
        mlab.pipeline.surface(mlab.pipeline.extract_edges(surf),
                                color=(0, 0, 1), opacity=0.9)







axmax = np.nanmax(grid+grid_sep/2)
axmin = np.nanmin(grid-grid_sep/2)


mlab.axes(extent=[axmin,axmax,axmin,axmax,axmin,axmax],nb_labels = 9)
mlab.show()


"""f1 = df3.form_1_3d(xg, yg, zg, Fx, Fy, Fz)
f1.plot()

x = np.array([-1.25,-1.25,-1.25,-1.25,-0.75,-0.75,-0.75,-0.75,-1])
y = np.array([-1.25,-1.25,-1.25,-1.25,-0.75,-0.75,-0.75,-0.75,-1])
z = np.array([-1.25,-1.25,-1.25,-1.25,-1,-1,-1,-1,-1])


x1 = np.array([1.25,1.25,1.25,1.25,0.75,0.75,0.75,0.75,1])
y1 = np.array([1.25,1.25,1.25,1.25,0.75,0.75,0.75,0.75,1])
z1 = np.array([1.25,1.25,1.25,1.25,1,1,1,1,1])

x_lst = [x1, x1]
y_lst = [y, y1]
z_lst = [z, z1]


f = mlab.figure(bgcolor=(1,1,1), fgcolor=(0,0,0))


for i in range(len(x_lst)):
    r = tvtk.RectilinearGrid()
    r.x_coordinates = x_lst[i]
    r.y_coordinates = y_lst[i]
    r.z_coordinates = z_lst[i]
    r.dimensions = len(x_lst[i]),len(y_lst[i]),len(z_lst[i])


    surf = mlab.pipeline.surface(r, opacity=0)
    mlab.pipeline.surface(mlab.pipeline.extract_edges(surf),
                            color=(1, 0, 0), )


mlab.axes(extent=[-5,5,-5,5,-5,5])

mlab.show()
"""