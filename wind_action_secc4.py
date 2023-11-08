# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 14:21:00 2023

@author: cfcpc2
"""

import numpy as np

c_dir = 1.0
c_season = 1.0
v_b0 = 27.0  # km/h
p = 0.01
K = 0.2
n = 0.5
rho = 1.25  # kg/m3
z_max = 200.0  # m
z_0 = 0.003  # m
z_min = 1.0  # m
z_0II = 0.005  # m
k_I = 1.00
A_ref = 800  # m2
c_d = 1.00
c_f = 1.55
c_0 = 1.00

# Basic velocity pressure given in Expression (4.10):
v_b = c_dir * c_season * v_b0

# The 10 minutes mean wind velocity having prob p
c_prob = ((1 - K * np.log(-np.log(1 - p))) / (1 - K * np.log(-np.log(0.98))) ** n)

# Terrain factor depending on the roughness length z0
k_r = 0.19 * (z_0 / z_0II) ** 0.07

# Basic velocity pressure
q_b = 0.5 * rho * v_b ** 2  # Fixed operator and added 0.5

sigma_v = k_r * v_b * k_I

def c_r(z):
    condlist = [z >= z_min and z <= z_max, z <= z_min]
    funclist = [k_r * np.log(z / z_0), k_r * np.log(z_min / z_0)]
    exp = np.piecewise(z, condlist, funclist)
    return exp

def v_m(z):
    exp = c_r(z) * c_0 * v_b
    return exp

def I_v(z):
    condlist = [z >= z_min and z <= z_max, z <= z_min]
    funclist = [sigma_v / v_m(z), sigma_v / v_m(z_min)]
    exp = np.piecewise(z, condlist, funclist)
    return exp

def q_p(z):
    exp = (1 + 7 * I_v(z)) * 0.5 * rho * v_m(z) ** 2  # Fixed operator and added 0.5
    return exp

def c_e(z):
    exp = q_p(z) / q_b
    return exp

print (c_r(10.0))

print(f'c_prob={c_prob}')


