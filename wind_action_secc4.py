# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 14:21:00 2023

@author: cfcpc2
"""

import numpy as np
from sympy import symbols, Eq, Function,UnevaluatedExpr, Mul, Rational, sqrt
from sympy import Piecewise, nan, N, And, log
from sympy import *
from sympy import N
init_printing()
import matplotlib.pyplot as plt


# Section 4 Wind velocity and velocity presure


## Fubntions to show decimals
def round_expr(expr, num_digits=3):
	return expr.xreplace({n : round(n, num_digits) for n in expr.atoms(Number)})

def round_equation(eq, num_digits=3):
	lhs = eq.lhs
	rhs = eq.rhs
	rounded_rhs = round_expr(rhs, num_digits)
	return Eq(lhs, rounded_rhs)


# Symbolics
						 
ue=UnevaluatedExpr

c_dir = symbols('c_dir')#1.0
c_season = symbols('c_season')#1.0
v_b0 = symbols('v_b0')#27.0  # km/h
p = symbols('p')#0.01
K = symbols('K')#0.2
n = symbols('n')#0.5
rho = symbols('rho')#1.25  # kg/m3
z=symbols('z')
z_max = symbols('z_max')#200.0  # m
z_0 = symbols('z_0')#0.003  # m
z_min = symbols('z_min')#1.0  # m
z_0II = symbols('z_0II')#0.005  # m
k_I = symbols('k_I')#1.00
A_ref = symbols('A_ref')#800  # m2
c_d = symbols('c_d')#1.00
c_f = symbols('c_f')#1.55
c_0 = symbols('c_0')#1.00
v_b=symbols('v_b')
k_r=symbols('k_r')
c_r=symbols('c_r',cls=Function)(z)
z_e, c_pe = symbols('z_e c_pe')
z_i, c_p_i = symbols('z_i c_p_i')
c_d, c_s,c_f, A_ref = symbols('c_d c_s c_f A_ref')
A_fr, c_fr=symbols('A_fr c_fr')
F_fr=symbols('F_fr', cls=Function)(z_e) 
c_s, z_s, B, R, k_p=symbols('c_s z_s B R k_p')
c_d=symbols('c_d ')
c_sd=symbols('c_sd ')
c_prob=symbols('c_prob')
v_m=symbols('v_m',cls=Function)(z)
q_b=symbols('q_b')
sigma_v=symbols('sigma_v')
c_e=symbols('c_e', cls=Function)(z)
I_v=symbols('I_v',cls=Function)(z)
q_p=symbols('q_p', cls=Function)(z)




# 4.2 Basic values



def v_b_func(**kwargs):   #The basic wind velocity shall be calculated from Expression (4.1).
	kwargs = {eval(key): UnevaluatedExpr(value) for key, value in kwargs.items()}
	expr=c_dir * c_season * v_b0
	expr = expr.subs(kwargs)
	_eq=Eq(v_b,expr)
	return _eq


def c_prob_func(**kwargs): # The 10 minutes mean wind velocity having prob p
	kwargs = {eval(key): UnevaluatedExpr(value) for key, value in kwargs.items()}

	expr=((1 - Mul(K , log(-log(1 - p)),evaluate=False))/ (1 - Mul(K , log(-log(0.98,evaluate=False),evaluate=False),evaluate=False)))**n
	expr = expr.subs(kwargs)
	_eq=Eq(c_prob,expr)
	return _eq

# 4.3 Mean wind
#4.3.1 Variation with height

def v_m_func(z=z,**kwargs): # The mean wind velocity vm(z) at a height z 
	kwargs = {eval(key): UnevaluatedExpr(value) for key, value in kwargs.items()}
	expr = c_r * c_0 * v_b
	# Apply eval() to all keys
	expr = expr.subs(kwargs)
	_eq=Eq(v_m,expr)
	return _eq

#4.3.2 Terrain roughness

def c_r_func(z=z,**kwargs): #The roughness factor, cr(z)
	kwargs = {eval(key): UnevaluatedExpr(value) for key, value in kwargs.items()}
	condlist = [And(z >= z_min,z <= z_max), z <= z_min]
	funclist = [Mul(k_r, log(N(z / z_0,2), evaluate=False),evaluate=False), Mul(k_r, log(z_min / z_0, evaluate=False),evaluate=False)]
	expr = Piecewise(*zip(funclist,condlist))
	expr = expr.subs(kwargs)
	_eq=Eq(c_r,expr)
	return _eq


def k_r_func(**kwargs):
	kwargs.pop('k_r', None)
	kwargs = {eval(key): UnevaluatedExpr(value) for key, value in kwargs.items()}
	expr = 0.19 * (z_0 / z_0II) ** 0.07
	expr = expr.subs(kwargs)
	_eq=Eq(k_r,expr)
	return _eq

# Table 4.1  Terrain categories and terrain parameters
def Terrain_Category(terrain_type):

	if terrain_type == "0":
		z_0 = 0.003
		z_min = 1
	if terrain_type == "I":
		z_0 = 0.01
		z_min = 1
	if terrain_type == "II":
		z_0 = 0.05
		z_min = 2
	if terrain_type == "III":
		z_0 = 0.03
		z_min = 5
	if terrain_type == "IV":
		z_0 = 1.0
		z_min = 10.0
		
	return z_0, z_min


#4.4 Wind turbulence

def sigma_v_func(**kwargs):
	kwargs = {eval(key): UnevaluatedExpr(value) for key, value in kwargs.items()}
	expr = k_r * v_b * k_I
	expr = expr.subs(kwargs)
	_eq=Eq(sigma_v,expr)
	return _eq

def I_v_func(**kwargs):
	kwargs = {eval(key): UnevaluatedExpr(value) for key, value in kwargs.items()}

	condlist = [And(z >= z_min,z <= z_max), z <= z_min]
	#funclist = [sigma_v / v_m, sigma_v / v_m.subs(z, z_min)]
	funclist = [k_I/(c_0*log(z/z_0)),k_I/(c_0*log(z_min/z_0)) ]
	expr = Piecewise(*zip(funclist, condlist))
	expr = expr.subs(kwargs)
	_eq=Eq(I_v,expr)
	return _eq


#4.5 Peak velocity pressure


def q_b_func(**kwargs):
	kwargs = {eval(key): UnevaluatedExpr(value) for key, value in kwargs.items()}
	expr= 0.5 * rho * v_b ** 2
	expr = expr.subs(kwargs)
	_eq=Eq(q_b,expr)
	return _eq


def q_p_func(**kwargs):
	kwargs = {eval(key): UnevaluatedExpr(value) for key, value in kwargs.items()}
	expr=(1+7*I_v)*(1/2)*rho*v_m**2   # Fixed operator and added 0.5
	expr = expr.subs(kwargs)
	_eq=Eq(q_p,expr)
	return _eq


def c_e_func(**kwargs):
	kwargs = {eval(key): UnevaluatedExpr(value) for key, value in kwargs.items()}
	expr = q_p / q_b
	expr = expr.subs(kwargs)
	_eq=Eq(c_e,expr)
	return _eq


#%%
#Section 5 Wind actions

#5.2 Wind pressure on surfaces

z_e, c_pe = symbols('z_e c_pe')
W_e= symbols('W_e', cls=Function)(z_e)

def W_e_func(**kwargs): # The wind pressure acting on the external surfaces, we 
	kwargs = {eval(key): UnevaluatedExpr(value) for key, value in kwargs.items()}
	exp= q_p.subs(z,z_e)*c_pe
	expr = expr.subs(kwargs)
	_eq=Eq(W_e,expr)
	return _eq


z_i, c_p_i = symbols('z_i c_p_i')
W_i= symbols('W_i', cls=Function)(z_i)

def W_i_func(**kwargs): #  The wind pressure acting on the internal surfaces
	kwargs = {eval(key): UnevaluatedExpr(value) for key, value in kwargs.items()}
	exp= q_p.subs(z,z_i)*c_p_i
	expr = expr.subs(kwargs)
	_eq=Eq(W_e,expr)
	return _eq

#5.3 Wind forces

c_d, c_s,c_f, A_ref = symbols('c_d c_s c_f A_ref')
F_w=symbols('F_w', cls=Function)(z_e) 

def F_w_func(**kwargs):#The wind forces for the whole structure
	kwargs = {eval(key): UnevaluatedExpr(value) for key, value in kwargs.items()}
	exp=c_s*c_d*c_f*q_p.subs(z,z_e)*A_ref
	expr = expr.subs(kwargs)
	_eq=Eq(F_w,expr)
	return _eq
	


#%%
# Section 6 Structural factor cscd

#6.3 Detailed procedure 
#6.3.1 Structural factor cscd

A_fr, c_fr=symbols('A_fr c_fr')
c_sd=symbols('c_sd ')
c_s, z_s, B, R, k_p=symbols('c_s z_s B R k_p')
c_d=symbols('c_d ')
F_fr=symbols('F_fr', cls=Function)(z_e) 


def c_sd_func(**kwargs):
	kwargs = {eval(key): UnevaluatedExpr(value) for key, value in kwargs.items()}
	exp = (1 + 2 *k_p* I_v.subs(z, z_s)*sqrt(B**2+R**2))/(1+7*I_v.subs(z, z_s))
	expr = expr.subs(kwargs)
	_eq=Eq(c_sd,expr)
	return _eq
def c_s_func(**kwargs): #The size factor cs 
	kwargs = {eval(key): UnevaluatedExpr(value) for key, value in kwargs.items()}
	exp = (1 + 7 * I_v.subs(z, z_s)  * sqrt(B**2))/(1+7*I_v.subs(z, z_s))
	expr = expr.subs(kwargs)
	_eq=Eq(c_s,expr)
	return _eq
	

def c_d_func(**kwargs): #The dynamic factor cd
	kwargs = {eval(key): UnevaluatedExpr(value) for key, value in kwargs.items()}
	exp = (1 + 2 *k_p* I_v.subs(z, z_s)*sqrt(B**2+R**2))/(1+7*I_v.subs(z, z_s)*sqrt(B**2))
	expr = expr.subs(kwargs)
	_eq=Eq(c_d,expr)
	return _eq

def F_fr_func(**kwargs):
	kwargs = {eval(key): UnevaluatedExpr(value) for key, value in kwargs.items()}
	exp=c_fr*q_p.subs(z,z_e)*A_fr
	expr = expr.subs(kwargs)
	_eq=Eq(F_fr,expr)
	return _eq







