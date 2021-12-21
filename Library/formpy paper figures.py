import formpy as fp
import numpy as np
import matplotlib.pyplot as plt

# %%

# Figure presenting funcitonality of a basic Vector field

# set up the axis grids, figure and axis 
r = np.linspace(-5, 5, 15)
xg, yg = np.meshgrid(r, r)
fig = plt.figure(figsize=(7, 7))
ax = fig.gca()
ax.set_aspect('equal')

# set up the vecotr field e^x and e^y components
u = xg*np.cos(yg)
v = yg*np.sin(xg)

# set up the vector filed object
field = fp.vector_field(xg, yg, u, v, fig=fig)
field.give_eqn('x*cos(y)','y*sin(x)')

field.plot()

# The method zoom_inset creates the zoomed field, creates the inset axis and plots the field on the inset axis.
# Target, zoom and dpd carry over from original zoom.
field.zoom_inset((3, -3), 10, 9)



# %%

# Figure for basic 0-form


# %%

# Figure for basic 1-form


# %%

# Figure for basic 2-form


# %%

# Figure to show Hodge operations 


# %%

# Figure to show exterior derivative 

# %%

# Figure to show wedge


# %%

# Figure to show the interior derivative (possibly with metirc inc.)

