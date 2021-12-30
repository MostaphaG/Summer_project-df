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

fig1=plt.figure(figsize=(12,6))
ax1 = fig1.add_subplot(121)
ax2 = fig1.add_subplot(122)

f1 = fp.form_1(xg,yg,u,v)
f1.plot(ax1)

f1.log_scaling()
f1.plot(ax2)

#%%

r = np.linspace(-5, 5, 21)
xg, yg = np.meshgrid(r, r)

z = np.exp(xg**2+yg**2)

fig1=plt.figure(figsize=(12,6))
ax1 = fig1.add_subplot(121)
ax2 = fig1.add_subplot(122)

f2 = fp.form_2(xg,yg,z)
f2.plot(ax1)

f2.log_scaling()
f2.plot(ax2)
