# GUI using dformpy

#%% Import Modules
import timeit
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
import matplotlib as mpl
from sympy import diff, simplify
from sympy.parsing.sympy_parser import parse_expr
from PIL import Image, ImageTk
from math import isnan
from matplotlib import patches as patch
import matplotlib.path as mplPath
from matplotlib import animation
from scipy.integrate import odeint
import matplotlib.path as mPath

import dformpy as fp

# input many numpy functions to deal with user input
from numpy import sin, cos, tan, sqrt, log, arctan, arcsin, arccos, tanh
from numpy import sinh, cosh, arcsinh, arccosh, arctanh, exp, pi, e

# %%
start = timeit.default_timer()

# =============================================================================
# Setup the GUI
# =============================================================================

root = tk.Tk()
root.title('Dformpy GUI')

# Create canvas for plotting
width = root.winfo_screenwidth()
height = root.winfo_screenheight()
root.geometry(str(width) + 'x' + str(height))
# get a toggle on and off switch images and scale the images to wanted size
toggle_image_on = Image.open('toggle_on_image.png')
toggle_image_off = Image.open('toggle_off_image.png')

toggle_image_on = toggle_image_on.resize((65, 25))
toggle_image_off = toggle_image_off.resize((65, 25))

toggle_image_on = ImageTk.PhotoImage(toggle_image_on)
toggle_image_off = ImageTk.PhotoImage(toggle_image_off)

# # Set up notebook
# style_notebook = ttk.Style()
# style_notebook.configure('TNotebook.Tab', font=('URW Gothic L','11','bold') )

# Set up the main frames
right_frame_frame = tk.LabelFrame(root, text='', padx=0, pady=0, bd=0)
right_frame_frame.grid(row=1, column=1, rowspan=100, sticky='N')
# bot frame:
bot_frame_frame = tk.LabelFrame(root, text='', padx=0, pady=0, bd=0)
bot_frame_frame.grid(row=2, column=0)
# plot frame:
plot_frame = tk.LabelFrame(root, text='', padx=0, pady=0, bd=0)
plot_frame.grid(row=1, column=0)

# define notebook for tabs
notebook_right = ttk.Notebook(right_frame_frame)
notebook_right.grid(row=0, column=0)
# notebook for bottom field input, when needed to disappear.
notebook_bot = ttk.Notebook(bot_frame_frame)
notebook_bot.grid(row=0, column=0)
# singularities notebook
notebook_singular = ttk.Notebook(right_frame_frame)
notebook_singular.grid(row=1, column=0)
# plotting options notebook
notebook_small = ttk.Notebook(bot_frame_frame)
notebook_small.grid(row=0, column=1)
# labels for hover over buttons notebook
notebook_instruct = ttk.Notebook(right_frame_frame)
notebook_instruct.grid(row=2, column=0)

# singularities:
singular_frame = tk.LabelFrame(notebook_singular)
singular_frame.grid(row=0, column=1)
# main options:
right_frame = tk.LabelFrame(notebook_right)
right_frame.grid(row=0, column=0)
# field input
bot_frame = tk.LabelFrame(notebook_bot)
bot_frame.grid(row=0, column=0)
# Line integrals
LI_frame = tk.LabelFrame(notebook_right)
LI_frame.grid(row=0, column=1)
# calculus
calculus_frame = tk.LabelFrame(notebook_right)
calculus_frame.grid(row=0, column=3)
# R3
r3_frame = tk.LabelFrame(notebook_right)
r3_frame.grid(row=0, column=4)
# dynamics
dynamics_frame = tk.LabelFrame(notebook_right)
dynamics_frame.grid(row=0, column=2)
# plotting options
small_frame = tk.LabelFrame(notebook_small)
small_frame.grid(row=0, column=0)
# labels for hovering
instruct_frame = tk.LabelFrame(notebook_instruct)
instruct_frame.grid(row=0, column=2)

notebook_right.add(right_frame, text='VF')
notebook_right.add(LI_frame, text='Line Integrals')
notebook_right.add(dynamics_frame, text='Dynamics')
notebook_right.add(calculus_frame, text='Ext. Alegebra')
notebook_right.add(r3_frame, text='R^3')
notebook_bot.add(bot_frame, text='1-Forms')
notebook_singular.add(singular_frame, text='singularities')
notebook_small.add(small_frame, text='Plotting')
# notebook_instruct.add(instruct_frame, text='Instructions')

# make an initial in instructions frame too
# instruct_frame_label = tk.Label(instruct_frame, text='', wraplength=400)
# instruct_frame_label.pack()

# Response to the users tab selection
def tab_selection(event):
    selected_tab = event.widget.select()
    tab_text = event.widget.tab(selected_tab, "text")
    main_axis.clear()
    
    if tab_text == 'VF':
        update_grids()
        tensor.set(0)
        plot_response()

    elif tab_text == 'Ext. Alegebra':
        F2 = fp.form_2(xg, yg, w, str(form_2_entry.get()))
        F2.plot(main_axis)

    print('tab')
    canvas.draw

# Bind the clicks on tabs to a function
notebook_right.bind_all('<<NotebookTabChanged>>', tab_selection)

# =============================================================================
# Set up figure
# =============================================================================

my_dpi = 100
fig = plt.figure(figsize=(730/my_dpi, 573/my_dpi), dpi=my_dpi)
main_axis = fig.gca()
main_axis.set_aspect('equal')
# delta = 10
# ax_L = L + L/10
# main_axis.set_xlim(-ax_L, ax_L)
# main_axis.set_ylim(-ax_L, ax_L)

fig.tight_layout()
canvas = FigureCanvasTkAgg(fig, master=plot_frame)
canvas.draw()
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
toolbar = NavigationToolbar2Tk(canvas, plot_frame)
toolbar.update()  # allow the plot to update based on the toolbar
canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

#  Set up initial plot
L = 5
origin_x = 0
origin_y = 0
pt_den = 15
gpts = np.linspace(origin_x - L, origin_y + L, pt_den)
xg, yg = np.meshgrid(gpts, gpts)
u = yg*np.sin(xg)
v = -xg*np.cos(yg)
w = xg*yg**2
s_max = 5

F1 = fp.form_1(xg, yg, u, v, 'y*sin(x)', '-x*cos(y)')
F1.max_sheets(s_max)
F1.plot(main_axis)
plt.close()

# =============================================================================
# Plot customisation functions and widgets 
# =============================================================================

def format_eq(string):#, LI=0, singular_ty=0, x_to_x_m =0, i=0, j=0):
    # replace all the x and y with xg and yg:
    # if LI == 1 :
    #     string = string.replace('x', 'intervals[0,:]')
    #     string = string.replace('y', 'intervals[1,:]')
    # elif singular_ty == 1:
    #     string = string.replace('x', 'points[' + str(0) + ', ' + str(i) + ']')
    #     string = string.replace('y', 'points[' + str(1) + ', ' + str(j) + ']')
    # elif x_to_x_m == 1:
    #     string = string.replace('x', 'x_m')
    #     string = string.replace('y', 'y_m')
    # else:
        
    string = string.replace('x', 'xg')
    string = string.replace('y', 'yg')
    string = string.replace('z', 'zg')
        
    # where there are special functions, replace them with library directions
    string = string.replace('^', '**')
    string = string.replace('ln', 'log')
    return string

# Equations to component grids
def eq_to_comps(string_x, string_y, xg, yg, check_2_frm=0):
    global equation_x, equation_y
    # use the format_eq fucntion to make given string readable by python
    equation_x = format_eq(string_x)
    equation_y = format_eq(string_y)
    
    # use these to define the field:
    # also: checking if equation equals zero, to then replace it with an array and not just 0:
    u = eval(equation_x)
    v = eval(equation_y)
    
    if equation_x.find('x') & equation_x.find('y') == -1:
        u = eval(equation_x)*np.ones(np.shape(xg))
    if equation_y.find('x') & equation_y.find('y') == -1:
        v = eval(equation_y)*np.ones(np.shape(yg))
        
    # deal properly with zero and constant fields too:
    # check for when the derivative is zero, do not plot it as nothing
    # if the other component is not zero.
    
    if check_2_frm == 1:
        if equation_x == '0' and equation_y != '0':
            u = np.ones(np.shape(xg))
        if equation_y == '0' and equation_x != '0':
            v = np.ones(np.shape(yg))
    else:
        pass
    # return these
    return u, v

def update_grids():
    global xg, yg, u, v
    xpts = np.linspace(eval(origin_x_entry.get()) - eval(L_entry.get()), eval(origin_x_entry.get()) + eval(L_entry.get()), int(pt_den_entry.get()))
    ypts = np.linspace(eval(origin_y_entry.get()) - eval(L_entry.get()), eval(origin_y_entry.get()) + eval(L_entry.get()), int(pt_den_entry.get()))
    xg, yg = np.meshgrid(xpts, ypts)
    u, v = eq_to_comps(str(x_comp_entry.get()), str(y_comp_entry.get()), xg, yg)

# Plot button 
def plot_response():
    global F1
    
    update_grids()
    main_axis.clear()
    
    if tensor.get() == 1:
        F1 = fp.vector_field(xg, yg, u, v, str(x_comp_entry.get()), str(y_comp_entry.get()))
        
    elif tensor.get() == 0:
        F1 = fp.form_1(xg, yg, u, v, str(x_comp_entry.get()), str(y_comp_entry.get()))
        F1.max_sheets(int(s_max_entry.get()))

    F1.logarithmic_scale_bool = logarithmic_scale_int.get()
    
    F1.plot(main_axis)
    canvas.draw()
    
plot_btn = tk.Button(small_frame, text='PLOT', padx=40, pady=20, command=plot_response)
plot_btn.grid(row=0, column=2, rowspan=2)

# Component entry boxes 
x_comp_entry_label = tk.Label(bot_frame, text='dx component')
x_comp_entry_label.grid(row=0, column=0)
x_comp_entry = tk.Entry(bot_frame, width=20, borderwidth=2)
x_comp_entry.grid(row=1, column=0)
x_comp_entry.insert(0, 'y*sin(x)')

y_comp_entry_label = tk.Label(bot_frame, text='dy component')
y_comp_entry_label.grid(row=0, column=1)
y_comp_entry = tk.Entry(bot_frame, width=20, borderwidth=2)
y_comp_entry.grid(row=1, column=1)
y_comp_entry.insert(0, '-x*cos(y)')

# Field drop down selection
def field_selection_response(event):
    # clear the x and y component boxes
    x_comp_entry.delete(0, 'end')
    y_comp_entry.delete(0, 'end')
    # get the index at which this entry is
    selected_index = field_name_list.index(str(field_select_drop.get()))
    # using that index, get the x and y components from their lists
    # and insert these into x and y comp. entry boxes
    x_comp_selected = field_x_list[selected_index]
    y_comp_selected = field_y_list[selected_index]
    x_comp_entry.insert(0, x_comp_selected)
    y_comp_entry.insert(0, y_comp_selected)
    # colour code to be able to distinguish what is being plotted
    x_comp_entry.configure(bg='#C0F6BB')
    y_comp_entry.configure(bg='#C0F6BB')
    plot_response()

field_name_list = ['Default: y*sin(x)dx - x*cos(y)dy',
                   'Simple pendulum: ydx  - sin(x)dy',
                   'Harmonic oscillator: ydx -xdy',
                   'Linear example 1: (x + 2*y)dx + (3*x - 4*y)dy',
                   'Linear example 2: xdx',
                   'Constant: 6dx + 3dy',
                   'Falling cat (Planar 3 link robot)',
                   'Electric Point Charge: -x/((x**2+y**2)**(1.5))dx + -y/((x**2+y**2)**(1.5))dy',
                   'H field of Current Carrying Wire: -y/((x**2+y**2)**(1.5))dx + x/((x**2+y**2)**(1.5))dy',
                   'Flamms paraboloid',
                   'BLACK HOLE!'
                   ]

field_x_list = ['y*sin(x)',
                'y',
                'y',
                'x + 2*y',
                'x',
                '6',
                '(3*cos(y) + 4)/(15 + 6*cos(x) + 6*cos(y))',
                '-x/((x**2+y**2)**(1.5))',
                '-y/((x**2+y**2)**(1.5))',
                'x/(sqrt(x**2 + y**2)*(1-2/(sqrt(x**2 + y**2)))) - y',
                '-2*x*((x^2+y^2)^(-1.5))*(1-(2/sqrt(x^2+y^2)))^(-2)'
                ]

field_y_list = ['- x*cos(y)',
                '-sin(x)',
                '-x',
                '3*x - 4*y',
                '0',
                '3',
                '-(3*cos(x) + 4)/(15 + 6*cos(x) + 6*cos(y))',
                '-y/((x**2+y**2)**(1.5))',
                'x/((x**2+y**2)**(1.5))',
                'y/(sqrt(x**2 + y**2)*(1-2/(sqrt(x**2 + y**2)))) + x',
                '-2*y*((x^2+y^2)^(-1.5))*(1-(2/sqrt(x^2+y^2)))^(-2)'
                ]
 
field_select_drop_label = tk.Label(bot_frame, text='Select Pre-Defined 1-Form:')
field_select_drop_label.grid(row=2, column=0, columnspan=2)
field_select_drop = ttk.Combobox(bot_frame, value=field_name_list, width=40)
field_select_drop.current(0)
field_select_drop.grid(row=3, column=0, columnspan=2)
field_select_drop.bind("<<ComboboxSelected>>", field_selection_response)

tk.Label(small_frame, text='Size').grid(row=2, column=0)
L_entry = tk.Entry(small_frame, width=5, borderwidth=1)
L_entry.grid(row=3, column=0, padx=2)
L_entry.insert(0, L)

tk.Label(small_frame, text='Grid').grid(row=2, column=1)
pt_den_entry = tk.Entry(small_frame, width=5, borderwidth=1)
pt_den_entry.grid(row=3, column=1, padx=2)
pt_den_entry.insert(0, pt_den)

tk.Label(small_frame, text='Max sheets').grid(row=2, column=2)
s_max_entry = tk.Entry(small_frame, width=5, borderwidth=1)
s_max_entry.grid(row=3, column=2, padx=2)
s_max_entry.insert(0, s_max)

tk.Label(small_frame, text='Origin x').grid(row=2, column=3)
origin_x_entry = tk.Entry(small_frame, width=5, borderwidth=1)
origin_x_entry.grid(row=3, column=3, padx=2)
origin_x_entry.insert(0, origin_x)

tk.Label(small_frame, text='Origin y').grid(row=2, column=4)
origin_y_entry = tk.Entry(small_frame, width=5, borderwidth=1)
origin_y_entry.grid(row=3, column=4, padx=2)
origin_y_entry.insert(0, origin_y)

def log_scale_toggle_response():
    # global logartmic_scale_bool
    if logarithmic_scale_int.get() == 0:
        # the burron is off, and has been clicked therefore change the
        # variable to an and the image to on
        logarithmic_scale_int.set(1)
        # logartmic_scale_bool = 1
        logarithmic_scale_toggle.configure(image=toggle_image_on)
    else:
        # the button is on and has been clicked
        # set it to off and change image
        logarithmic_scale_int.set(0)
        # logartmic_scale_bool = 0
        logarithmic_scale_toggle.configure(image=toggle_image_off)
    plot_response()

tk.Label(small_frame, text='Log Scale').grid(row=0, column=1)
logarithmic_scale_int = tk.IntVar()
logarithmic_scale_int.set(0)
logarithmic_scale_toggle = tk.Button(small_frame, image=toggle_image_off, bd=0, command=log_scale_toggle_response)
logarithmic_scale_toggle.grid(row=1, column=1)

# logartmic_scale_toggle.bind('<Enter>', lambda x: hover_instruction_response(0, 1))
# logartmic_scale_toggle.bind('<Leave>', lambda x: hover_instruction_response(0, 0))

# =============================================================================
# VF tab widgets
# =============================================================================

# Track radioubutton selection
tensor = tk.IntVar()
tensor.set(0)

tensor_label = tk.Label(right_frame, text='Arrows/Stacks:')
tensor_label.grid(row=8, column=0)

#  Stack: 0
#  Arrow: 1
def vect_type_response(tensor):
    # main_axis.clear()
    # canvas.draw()
    plot_response()
    canvas.draw()

stack_btn = tk.Radiobutton(right_frame, text='Stack', variable=tensor, value=0, command=lambda: vect_type_response(tensor.get())).grid(row=8, column=1)
arrow_btn = tk.Radiobutton(right_frame, text='Arrow', variable=tensor, value=1, command=lambda: vect_type_response(tensor.get())).grid(row=8, column=2)
# arrow_stack_btn = tk.Radiobutton(right_frame, text='both', variable=tensor, value=2, command=lambda: vect_type_response(tensor.get())).grid(row=8, column=3)


def click_option_handler():
    pass

# Radiobuttons to select what happens when clicking the plot
click_option = tk.IntVar()
click_option.set(0)
click_option_Tools_btn = tk.Radiobutton(right_frame, text='Tools', variable=click_option, value=0, command=lambda: click_option_handler(click_option.get()))
click_option_Zoom_btn = tk.Radiobutton(right_frame, text='Zoom', variable=click_option, value=1, command=lambda: click_option_handler(click_option.get()))
click_option_Deriv_btn = tk.Radiobutton(right_frame, text='Deriv.', variable=click_option, value=2, command=lambda: click_option_handler(click_option.get()))
click_option_Div_btn = tk.Radiobutton(right_frame, text='Div.', variable=click_option, value=3, command=lambda: click_option_handler(click_option.get()))
click_option_Curl_btn = tk.Radiobutton(right_frame, text='Curl', variable=click_option, value=4, command=lambda: click_option_handler(click_option.get()))
click_option_Deriv_btn.configure(state=tk.DISABLED)
click_option_Div_btn.configure(state=tk.DISABLED)
click_option_Curl_btn.configure(state=tk.DISABLED)
click_option_Tools_btn.grid(row=1, column=0)
click_option_Zoom_btn.grid(row=1, column=1)
click_option_Deriv_btn.grid(row=1, column=2)
click_option_Div_btn.grid(row=2, column=0)
click_option_Curl_btn.grid(row=2, column=1)

# click_option_Tools_btn.bind('<Enter>', lambda x: hover_instruction_response(5, 1))
# click_option_Tools_btn.bind('<Leave>', lambda x: hover_instruction_response(5, 0))

# click_option_Zoom_btn.bind('<Enter>', lambda x: hover_instruction_response(6, 1))
# click_option_Zoom_btn.bind('<Leave>', lambda x: hover_instruction_response(6, 0))

# click_option_Deriv_btn.bind('<Enter>', lambda x: hover_instruction_response(7, 1))
# click_option_Deriv_btn.bind('<Leave>', lambda x: hover_instruction_response(7, 0))

# click_option_Div_btn.bind('<Enter>', lambda x: hover_instruction_response(8, 1))
# click_option_Div_btn.bind('<Leave>', lambda x: hover_instruction_response(8, 0))

# click_option_Curl_btn.bind('<Enter>', lambda x: hover_instruction_response(9, 1))
# click_option_Curl_btn.bind('<Leave>', lambda x: hover_instruction_response(9, 0))

def update_deriv():
    pass

# Zoom slider
tk.Label(right_frame, text='Zoom').grid(row=3, column=0)
zoom_slider = tk.Scale(right_frame, from_=1, to=100, orient=tk.HORIZONTAL, resolution=1)
zoom_slider.bind("<ButtonRelease-1>", update_deriv)
zoom_slider.grid(row=3, column=1)

# Zoom pt density
dpd_select = tk.IntVar()
dpd_select.set(5)
dpd_list = [5, 7, 9]
tk.Label(right_frame, text='Inset Plot Point Density:').grid(row=4, column=0)
dpd_drop = tk.OptionMenu(right_frame, dpd_select, *dpd_list, command = update_deriv)
dpd_drop.grid(row=4, column=1)

# Insize
d_length_select = tk.DoubleVar()
d_length_list = [0.1, 0.2, 0.3, 0.4, 0.5]
d_length_select.set(d_length_list[2])
tk.Label(right_frame, text='Inset Fractional Size:').grid(row=5, column=0)
d_length_drop = tk.OptionMenu(right_frame, d_length_select, *d_length_list, command=update_deriv)
d_length_drop.grid(row=5, column=1)

def scale_toggle_response():
    pass

# Autoscale Toggle
ascale_label = tk.Label(right_frame, text='Autoscale arrows:')
ascale_label.grid(row=7, column=0)
ascale_toggle = tk.Button(right_frame, image=toggle_image_off, bd=0, command=scale_toggle_response)
ascale_toggle.grid(row=7, column=1, pady=5)

# ascale_toggle.bind('<Enter>', lambda x: hover_instruction_response(1, 1))
# ascale_toggle.bind('<Leave>', lambda x: hover_instruction_response(1, 0))

def set_inset_target():
    pass

# define entry boxes to allow user to input x_m and y_m
x_m_entry = tk.Entry(right_frame, width=12)
y_m_entry = tk.Entry(right_frame, width=12)
x_m_entry.grid(row=6, column=0)
y_m_entry.grid(row=6, column=1)
# and a button to submit these:
Set_target_btn = tk.Button(right_frame, text='Set Target', command=set_inset_target)
Set_target_btn.grid(row=6, column=2, padx=20)

# Set_target_btn.bind('<Enter>', lambda x: hover_instruction_response(2, 1))
# Set_target_btn.bind('<Leave>', lambda x: hover_instruction_response(2, 0))

def analytic_toggle_response():
    pass

analytic_select = tk.IntVar()
analytic_select.set(0)
tk.Label(right_frame, text= 'Toggle Analytic Label:').grid(row=9, column=0)
analytic_toggle = tk.Button(right_frame, image=toggle_image_off, bd=0, command=analytic_toggle_response)
analytic_toggle.grid(row=9, column=1)

# analytic_toggle.bind('<Enter>', lambda x: hover_instruction_response(4, 1))
# analytic_toggle.bind('<Leave>', lambda x: hover_instruction_response(4, 0))

# =============================================================================
# Exterior algebra tab widgets 
# =============================================================================

# 2-form entry
tk.Label(calculus_frame, text='2-form on R2').grid(row=0, column=1)
form_2_entry = tk.Entry(calculus_frame, width=15, borderwidth=2)
form_2_entry.grid(row=0, column=0)
form_2_entry.insert(0, 'x*y**2')
form_2_entry.configure(bg='#C0F6BB')

def form_2_response():
    pass

# 2-form plot button
form_2_btn = tk.Button(calculus_frame, text='2-form plot', padx=3, pady=5, command=form_2_response)
form_2_btn.grid(row=3, column=1)

# 0-form entry
tk.Label(calculus_frame, text='Zero form:').grid(row=4, column=0)
form_0_entry = tk.Entry(calculus_frame, width=15, borderwidth=2)
form_0_entry.grid(row=4, column=1)
form_0_entry.insert(0, '')

def form_0_response():
    pass

# 0-form plot button
form_0_btn = tk.Button(calculus_frame, text='0-form plot', padx=3, pady=5, command=form_0_response)
form_0_btn.grid(row=4, column=2)

def form_1_stacks_response():
    pass

# 1-form plot button
form_1_btn = tk.Button(calculus_frame, text='1-form plot', padx=3, pady=5, command=form_1_stacks_response)
form_1_btn.grid(row=3, column=0)

def Int_deriv_response():
    pass

# Interior derivative button
INT_btn = tk.Button(calculus_frame, text='Int Deriv', padx=0, pady=2, command=Int_deriv_response)
INT_btn.grid(row=5, column=0)

def Ext_deriv_response():
    pass

# Exterior derivative button
EXT_int_btn = tk.Button(calculus_frame, text='Ext Deriv', padx=0, pady=2, command=Ext_deriv_response)
EXT_int_btn.grid(row=5, column=1)

def wedge_2_response():
    pass

# Wedge product button
wedge_btn = tk.Button(calculus_frame, text='Wedge', padx=0, pady=2, command=wedge_2_response)
wedge_btn.grid(row=6, column=0)


def Hodge_full():
    pass

# hodge button
Hodge_btn = tk.Button(calculus_frame, text='Hodge', padx=5, pady=2, command=Hodge_full)
Hodge_btn.grid(row=7, column=0)

def R2_tools_handler():
    pass

# Radiobuttons
R2_tools_opt = tk.IntVar()
R2_tools_opt.set(0)
R2_tools_Tools_btn = tk.Radiobutton(calculus_frame, text='Tools', variable=R2_tools_opt, value=0, command=lambda: R2_tools_handler(R2_tools_opt.get()))
R2_tools_Zoom_btn = tk.Radiobutton(calculus_frame, text='Zoom', variable=R2_tools_opt, value=1, command=lambda: R2_tools_handler(R2_tools_opt.get()))
R2_tools_int_btn = tk.Radiobutton(calculus_frame, text='Area Int', variable=R2_tools_opt, value=2, command=lambda: R2_tools_handler(R2_tools_opt.get()))
R2_tools_Tools_btn.grid(row=8, column=0)
R2_tools_Zoom_btn.grid(row=8, column=1)
R2_tools_int_btn.grid(row=8, column=2)

def update_2_form_zoom():
    pass

# Zoom slider 
tk.Label(calculus_frame, text='Zoom').grid(row=9, column=0)
zoom_slider_R2 = tk.Scale(calculus_frame, from_=1, to=20, orient=tk.HORIZONTAL)
zoom_slider_R2.bind("<ButtonRelease-1>", update_2_form_zoom)
zoom_slider_R2.grid(row=9, column=1)
zoom_slider_R2.configure(state=tk.DISABLED)

# 2-form zoom pt density
zoomR2pd_select = tk.IntVar()
zoomR2pd_select.set(11)
zoomR2pd_list = [5, 6, 10, 11, 15, 16, 20, 21]
tk.Label(calculus_frame, text='Inset Plot Point Density:').grid(row=10, column=0)
zoomR2pd_drop = tk.OptionMenu(calculus_frame, zoomR2pd_select, *zoomR2pd_list, command=update_2_form_zoom)
zoomR2pd_drop.grid(row=10, column=1)

# 2-form zoom insize
zoomR2_length_select = tk.DoubleVar()
zoomR2_length_list = [0.1, 0.2, 0.3, 0.4, 0.5]
zoomR2_length_select.set(zoomR2_length_list[2])
tk.Label(calculus_frame, text='Inset Fractional Size:').grid(row=11, column=0)
zoomR2_length_drop = tk.OptionMenu(calculus_frame, zoomR2_length_select, *zoomR2_length_list, command=update_2_form_zoom)
zoomR2_length_drop.grid(row=11, column=1)

def set_inset_target_calc():
    pass

# Set target
x_m_entry_calc = tk.Entry(calculus_frame, width=12)
y_m_entry_calc = tk.Entry(calculus_frame, width=12)
x_m_entry_calc.grid(row=12, column=0)
y_m_entry_calc.grid(row=12, column=1)
# and a button to submit these:
Set_target_btn_calc = tk.Button(calculus_frame, text='Set Target', command=set_inset_target_calc)
Set_target_btn_calc.grid(row=12, column=2, padx=20)

# 2-form dropdown
# select_form_2 = tk.StringVar()
# select_form_2.set(list_form_2_names[0])

# select_form_2_drop_label = tk.Label(calculus_frame, text='Select Pre-Defined 2-Form:')
# select_form_2_drop_label.grid(row=1, column=0, columnspan=3)

# select_form_2_drop = ttk.Combobox(calculus_frame, value=list_form_2_names, width=40)
# select_form_2_drop.current(0)
# select_form_2_drop.grid(row=2, column=0, columnspan=3)
# select_form_2_drop.bind("<<ComboboxSelected>>", selection_form_2_response)

stop = timeit.default_timer()
print('Time: ', stop - start)

tk.mainloop()


