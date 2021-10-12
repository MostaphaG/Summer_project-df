# formpy examples
import formpy as fp
import numpy as np
import matplotlib.pyplot as plt

# %%

# plotting a 1-form

# set up needed parameters
v = np.linspace(-6, 6, 31)
xg, yg = np.meshgrid(v, v)
F_x = yg*np.sin(xg)
F_y = xg*np.cos(yg)

# PLOT, note, it will create a figure for user
# we probably don't want that, otherwise we would have to make this
# a method to the matplotlib object, which might mean we need to play with
# their library, which I suppose we can't.
form_obj = fp.form_1(xg, yg, F_x, F_y)

form_obj.axis.set_xlabel('x')
form_obj.axis.set_ylabel('y')
form_obj.arrows()
form_obj.colour('blue')
form_obj.fig_size(8, 9)
form_obj.head_width(0.3)
form_obj.orient('tail')
form_obj.max_sheets(8)
form_obj.sheet_size(0.04)
form_obj.surround_space(6)

form_obj.plot()

# %%

# plotting a 2-form

v = np.linspace(-6, 6, 31)
xg, yg = np.meshgrid(v, v)

form_2 = xg*yg


form_obj = fp.form_2(xg, yg, form_2)
form_obj.plot()


# %%

# plotting a 0-form

v = np.linspace(-4.5, 4.5, 11)
xg, yg = np.meshgrid(v, v)


form_0 = np.cos(xg*yg)

form_obj = fp.form_0(xg, yg, form_0)
form_obj.lines_number(4)
# form_obj.density_increase(20)  # demonstation of an error
form_obj.plot()

plt.pause(5)

# customising grids with an equation:

form_obj.form_0_give_eqn('cos(x*y)')

# change the density:
form_obj.density_increase(20)
form_obj.lines_number(15)
form_obj.plot(keep=False)


