# -*- coding: utf-8 -*-
"""
Created on Thu Dec 30 14:03:11 2021

@author: single user
"""
#%%
# Testing log scaling for 1-forms

import formpy as fp
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

import formpy as fp
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

vf1z = vf1.curl(target=(1,0), zoom=5)
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

import formpy as fp
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









