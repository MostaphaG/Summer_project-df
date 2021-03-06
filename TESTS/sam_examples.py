# -*- coding: utf-8 -*-
"""
Created on Thu Dec 30 14:03:11 2021

@author: single user
"""
#%%
# Testing log scaling for 1-forms

import dformpy as fp
import numpy as np
import matplotlib.pyplot as plt
import timeit

r = np.linspace(-5, 5, 21)
xg, yg = np.meshgrid(r, r)

u = np.exp(xg*yg)
v = np.sin(yg)

fig1 = plt.figure(figsize=(12, 6))
ax1 = fig1.add_subplot(121)
ax2 = fig1.add_subplot(122)

f1 = fp.form_1(xg, yg, u, v)
f1.plot(ax1)

f1.log_scaling()
f1.plot(ax2)

# %%

r = np.linspace(-5, 5, 21)
xg, yg = np.meshgrid(r, r)

z = np.exp(xg**2 + yg**2)

fig1 = plt.figure(figsize=(12, 6))
ax1 = fig1.add_subplot(121)
ax2 = fig1.add_subplot(122)

f2 = fp.form_2(xg, yg, z)
f2.plot(ax1)

f2.log_scaling()
f2.plot(ax2)

#%%

# Test whether the curl is dodgy

import dformpy as fp
import numpy as np
import matplotlib.pyplot as plt
import timeit

r = np.linspace(-5, 5, 21)
x,y = np.meshgrid(r, r)
u = x*np.cos(y)
v = -y*np.sin(x)

vf1 = fp.vector_field(x,y,u,v)

fig1 = plt.figure(figsize=(12,6))
ax1 = fig1.add_subplot(121)
ax2 = fig1.add_subplot(122)


vf1.plot(ax1)
vf1.give_eqn('y*sin(x)', '-x*cos(y)')

vf1z = vf1.curl(target=(1,0), mag=5)
vf1z.plot(ax2)

U = vf1z.F_x
V = vf1z.F_y

dpd = len(U[:,0])


t = 0

if U[0,0] < 0:
    for i in range(dpd):
        if U[0,i] >= 0:
            t += 1 
            
elif U[0,0] > 0:
    for i in range(dpd):
        if U[0,i] <= 0:
            t += 1
            
elif U[0,0] == 0:
    t += 1

if V[0,0] < 0:
    for i in range(dpd):
        if V[i,0] >= 0:
            t += 1 
            
elif V[0,0] > 0:
    for i in range(dpd):
        if V[i,0] <= 0:
            t += 1
            
elif V[0,0] == 0:
    t += 1

# t = 1 means the curl is not accurate (zero or needs more zooming)

# %%

# Can we div and curl using even dpd?

import dformpy as fp
import numpy as np
import matplotlib.pyplot as plt
import timeit

r = np.linspace(-5, 5, 21)
x,y = np.meshgrid(r, r)
u = x*np.cos(y)
v = -y*np.sin(x)

vf1 = fp.vector_field(x,y,u,v)

fig1 = plt.figure(figsize=(18,6))
ax1 = fig1.add_subplot(131)
ax2 = fig1.add_subplot(132)
ax3 = fig1.add_subplot(133)


vf1.plot(ax1)
vf1.give_eqn('y*sin(x)', '-x*cos(y)')

vf1d = vf1.div((1,1), 20, 6)
vf1d.plot(ax2)

vf1c = vf1.curl((1,1), 20, 20)
vf1c.plot(ax3)

# Seems to work okay

# %%

# s_min testing

import dformpy as fp
import numpy as np
import matplotlib.pyplot as plt
import timeit

r = np.linspace(-5, 5, 21)
x,y = np.meshgrid(r, r)
u = y*np.sin(x)
v = -x*np.cos(y)
z = np.cos(x)*np.sin(1/y)

F1 = fp.form_1(x,y,u,v)
F2 = fp.form_2(x,y,z)

fig1 = plt.figure(figsize=(18,6))
ax1 = fig1.add_subplot(131)
ax2 = fig1.add_subplot(132)
ax3 = fig1.add_subplot(133)

F1.s_max = 10
F1.s_min = 1

F2.s_max = 10
F2.s_min = 4

F1.plot(ax1)
F2.plot(ax2)

# %%

#  1-form customisation testing

import dformpy as fp
import numpy as np
import matplotlib.pyplot as plt
import timeit

r = np.linspace(0, 5, 20)
x,y = np.meshgrid(r, r)
u = y*np.sin(x)
v = -x*np.cos(y)
# z = np.cos(x)*np.sin(1/y)

F1 = fp.form_1(x,y,u,v)
F1.give_eqn('y*sin(x)', '-x*cos(y)')

fig1 = plt.figure(figsize=(8,8))
ax1 = fig1.add_subplot(111)

# F1.colour('b')
# F1.arrow_heads()
# F1.head_width(0.5)
# F1.head_height(0.5)
# F1.log_scaling()
# F1.max_sheets(10)
# F1.sheet_size(0.02)
# F1.surround_space(20)
# F1.set_density(30)

F1.plot(ax1)

# %% Script to find the laplacian of a given function (0-form)

import dformpy as fp
import numpy as np
import matplotlib.pyplot as plt
import timeit

r = np.linspace(-2,2,20)
x,y = np.meshgrid(r, r)
# u = y*np.sin(x)
# v = -x*np.cos(y)
z = y*np.sin(x)

fig1 = plt.figure(figsize=(12,6))
ax1 = fig1.add_subplot(121)
ax2 = fig1.add_subplot(122)
ax1.set_aspect('equal')
ax2.set_aspect('equal')
# ax3 = fig1.add_subplot(133)

f = fp.form_0(x,y,z)
f.give_eqn('(x**2+y**2)**(-0.5)')
f.plot(ax1)

df = f.ext_d()
# df.plot(ax2)

hdf = df.hodge()
# hdf.plot(ax3)

dhdf = hdf.ext_d()

hdhdf = dhdf.hodge()

hdhdf.plot(ax2)

# %% Log scaling for VF

import dformpy as fp
import numpy as np
import matplotlib.pyplot as plt
import timeit

r = np.linspace(-2, 2, 25)
x,y = np.meshgrid(r, r)

u = x*(x**2 + y**2)**(-3)
v = y*(x**2 + y**2)**(-3)

fig1 = plt.figure(figsize=(12,12))
ax1 = fig1.add_subplot(121)
ax2 = fig1.add_subplot(122)
# ax3 = fig1.add_subplot(223)
# ax4 = fig1.add_subplot(224)

ax1.set_aspect('equal')
ax2.set_aspect('equal')
# ax3.set_aspect('equal')
# ax4.set_aspect('equal')

VF = fp.vector_field(x,y,u,v)
VF.plot(ax1)
VF.log_scaling()
VF.plot(ax2)

# VF.base = 5
# VF.plot(ax2)

# F1 = fp.form_1(x,y,u,v)
# F1.plot(ax1)
# F1.log_scaling()
# F1.plot(ax2)

# %% 0 -form log scaling

import dformpy as fp
import numpy as np
import matplotlib.pyplot as plt

r = np.linspace(-2, 2, 25)
x,y = np.meshgrid(r, r)

u = x*(x**2 + y**2)**(-1.5)

fig1 = plt.figure(figsize=(10,5))
ax1 = fig1.add_subplot(121)
ax2 = fig1.add_subplot(122)

F0 = fp.form_0(x,y,u)
F0.density_increase(5)
F0.give_eqn('x*(x**2 + y**2)**(-1.5)')
F0.plot(ax1)
F0.log_scaling()
F0.plot(ax2)

# %%

import dformpy as fp
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-5, 5, 25)
y = np.linspace(-3, 5, 18)
xg ,yg  = np.meshgrid(x , y)

u = yg*np.sin(xg)/((xg - 5)*(xg+2))
v = -xg*np.cos(yg)/(xg)

fig1 = plt.figure(figsize=(10,5))
ax1 = fig1.add_subplot(121)
ax2 = fig1.add_subplot(122)

F0 = fp.form_2(xg, yg, u)
F0.plot(ax1)
F0.log_scaling()
F0.plot(ax2)

# %%

import numpy as np

x = 5
y = 10

A = np.array([[1,2,3],[4,5,6],[7,8,9]])

R = np.zeros(shape=(x,y))



