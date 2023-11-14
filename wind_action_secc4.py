# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 14:21:00 2023

@author: cfcpc2
"""

import numpy as np
from sympy import symbols, Eq, Function,UnevaluatedExpr, Mul, Rational
from sympy import Piecewise, nan
from sympy import *
from sympy import N
init_printing()
import matplotlib.pyplot as plt

def round_expr(expr, num_digits=3):
	return expr.xreplace({n : round(n, num_digits) for n in expr.atoms(Number)})

def round_equation(eq, num_digits=3):
	lhs = eq.lhs
	rhs = eq.rhs
	rounded_rhs = round_expr(rhs, num_digits)
	return Eq(lhs, rounded_rhs)

                         
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






# Basic velocity pressure given in Expression (4.10):
#v_b = c_dir * c_season * v_b0
def v_b_func(c_dir=c_dir, c_season=c_season, v_b0=v_b0):
	c_dir=UnevaluatedExpr(c_dir)
	c_season=UnevaluatedExpr(c_season)
	v_b0=UnevaluatedExpr(v_b0)
	return Eq(v_b, c_dir * c_season * v_b0, evaluate=False)
	
	


# The 10 minutes mean wind velocity having prob p
#c_prob = ((1 - K * np.log(-np.log(1 - p))) / (1 - K * np.log(-np.log(0.98))) )** n

c_prob=symbols('c_prob')

def c_prob_func(K=K,p=p,n=n):
	K=UnevaluatedExpr(K)
	p=UnevaluatedExpr(p)
	n=UnevaluatedExpr(n)
	
	val=((1 - Mul(K , log(-log(1 - p)),evaluate=False))/ (1 - Mul(K , log(-log(0.98,evaluate=False),evaluate=False),evaluate=False))) **n
	#val=((1 - K * np.log(p)) / 2) ** n
	#val=((1 - K ) / (1 - K * np.log(0.98))) ** n
	#val = ((1 - K * log(-log(1 - p)))**n, (1 - K * log(-log(0.98)))**n)
	
	return Eq(c_prob,val, evaluate=False)

# Terrain factor depending on the roughness length z0
k_r=symbols('k_r')

def k_r_func(z_0=z_0,z_0II=z_0II):
	z_0=UnevaluatedExpr(z_0)
	z_0II=UnevaluatedExpr(z_0II)
	val = 0.19 * (z_0 / z_0II) ** 0.07
	return Eq(k_r,val,evaluate=False)
	

# Basic velocity pressure
q_b=symbols('q_b')

def q_b_func(rho=rho, v_b=v_b):
	rho=UnevaluatedExpr(rho)
	v_b=UnevaluatedExpr(v_b)

	val= 0.5 * rho * v_b ** 2  # Fixed operator and added 0.5
	return Eq(q_b, val, evaluate=False)

# sigma_v
sigma_v=symbols('sigma_v')

def sigma_v_func(k_r=k_r, v_b=v_b, k_I=k_I):
	k_r=UnevaluatedExpr(k_r)
	v_b=UnevaluatedExpr(v_b)
	k_I=UnevaluatedExpr(k_I)
	
	val = k_r * v_b * k_I
	return Eq(sigma_v, val, evaluate=False)
	
c_r=symbols('c_r',cls=Function)(z)

def c_r_func(z=z, z_min=z_min, z_max=z_max, z_0=z_0, k_r=k_r, UE=False):
	if UE:
		z=UnevaluatedExpr(z)
		z_min=UnevaluatedExpr(z_min)
		z_max=UnevaluatedExpr(z_max)
		z_0=UnevaluatedExpr(z_0)
		k_r=UnevaluatedExpr(k_r)

	condlist = [And(z >= z_min,z <= z_max), z <= z_min]
	funclist = [Mul(k_r, log(N(z / z_0,2), evaluate=False),evaluate=False), Mul(k_r, log(z_min / z_0, evaluate=False),evaluate=False)]
	exp = Piecewise(*zip(funclist,condlist))
	return Eq(c_r, exp, evaluate=False)

v_m=symbols('v_m',cls=Function)(z)
def v_m_func(z=z,c_r=c_r, c_0=c_0,v_b=v_b):
	z=UnevaluatedExpr(z)
	c_r=UnevaluatedExpr(c_r)
	c_0=UnevaluatedExpr(c_0)
	v_b=UnevaluatedExpr(v_b)
	
	exp = c_r * c_0 * v_b
	return Eq(v_m, exp, evaluate=False)

I_v=symbols('I_v',cls=Function)(z)

def I_v_func(z=z,z_min=z_min, z_max=z_max, z_0=z_0, sigma_v=sigma_v, v_m=v_m, UE=False):
#def I_v_func(z=z, z_min=z_min, z_max=z_max, z_0=z_0, UE=False):
	if UE:
		z=UnevaluatedExpr(z)
		sigma_v=UnevaluatedExpr(sigma_v)
		v_m=UnevaluatedExpr(v_m)
		z_min=UnevaluatedExpr(z_min)
		z_max=UnevaluatedExpr(z_max)
		z_0=UnevaluatedExpr(z_0)

	condlist = [And(z >= z_min,z <= z_max), z <= z_min]

	funclist = [sigma_v / v_m,  sigma_v / v_m.subs(z, z_min)]
	#funclist = [k_I/(c_0*log(z/z_0)),k_I/(c_0*log(z_min/z_0)) ]
	exp = Piecewise(*zip(funclist, condlist))
	return Eq(I_v, exp, evaluate=False)


q_p=symbols('q_p', cls=Function)(z)
def q_p_func(z=z, I_v=I_v, rho=rho, v_m=v_m):
	I_v=UnevaluatedExpr(I_v)
	rho=UnevaluatedExpr(rho)
	v_m=UnevaluatedExpr(v_m)

	exp = (1 + 7 * I_v) * 0.5 * rho * v_m ** 2  # Fixed operator and added 0.5
	return Eq(q_p, exp, evaluate=False)

c_e=symbols('c_e', cls=Function)(z)
def c_e_func(z=z,q_p=q_p, q_b=q_b):
	q_p=UnevaluatedExpr(q_p)
	q_b=UnevaluatedExpr(q_b)

	exp = q_p / q_b
	return Eq(c_e, exp,evaluate=False)

#%%


# Define values
c_dir = 1.0
c_season = 1.0
v_b0 = 27.0  # km/h
p = 0.01
K = 0.2
n = 0.5
rho = 1.25
z_max = 200.0  # m
z_0 = 0.003  # m
z_min = 1.0  # m
z_0II = 0.005  # m
k_I = 1.00
A_ref = 800.0  # m2
c_d = 1.00
c_f = 1.55
c_0 = 1.00
#parameters=[c_dir,c_season,v_b0,p,K,n,rho,z_max,z_0,z_min,z_0II, k_I, A_ref, c_d, c_f,c_0]
#parameters=[c_dir=c_dir,c_season=c_season,v_b0=v_b0,p=p,K=K,n=n,rho=rho,z_max=z_max,z_0=z_0,z_min=z_min,z_0II=z_0II, k_I=k_I, A_ref=A_ref, c_d=c_d, c_f=c_f,c_0=c_0]
def c_ez(z,c_dir=c_dir,c_season=c_season,v_b0=v_b0,p=p,K=K,n=n,rho=rho,z_max=z_max,z_0=z_0,z_min=z_min,z_0II=z_0II, k_I=k_I, A_ref=A_ref, c_d=c_d, c_f=c_f,c_0=c_0):
	

# Basic velocity pressure given in Expression (4.10):

	def v_b():
		exp=c_dir * c_season * v_b0
		return exp


	def c_prob():
		exp = ((1 - K * np.log(-np.log(1 - p))) / (1 - K * np.log(-np.log(0.98))) )** n
		
		return exp

	def k_r():

		exp = 0.19 * (z_0 / z_0II) ** 0.07
		return exp

	def q_b():

		exp= 0.5 * rho * v_b() ** 2  # Fixed operator and added 0.5
		return exp

	def sigma_v():

		exp = k_r() * v_b() * k_I
		return exp

	def c_r(z):
		condlist = [np.logical_and(z >= z_min, z <= z_max), z <= z_min]
		funclist = [k_r() * np.log(z / z_0), k_r() * np.log(z_min / z_0)]
		if z != 0 and z_min != 0 and z_0 != 0:
			funclist = [k_r() * np.log(z / z_0), k_r() * np.log(z_min / z_0)]
		else:
			pass  
		#funclist = [k_r() * np.log(np.maximum(z, 1e-6) / z_0), k_r() * np.log(np.maximum(z_min, 1e-6) / z_0)]
		exp = np.piecewise(z, condlist, funclist)
		return exp

	def v_m(z):

		exp = c_r(z) * c_0 * v_b()
		return exp

	def I_v(z):
		condlist = [np.logical_and(z >= z_min, z <= z_max), z <= z_min]
		funclist = [sigma_v() / v_m(z), sigma_v() / v_m(z_min)]
		exp = np.piecewise(z, condlist, funclist)
		return exp

	def q_p(z):

		exp = (1 + 7 * I_v(z)) * 0.5 * rho * v_m(z) ** 2  # Fixed operator and added 0.5
		return exp

	val = q_p(z) / q_b()
	return val




# Generate 1000 points linearly spaced between 0 and 100
z_values = np.linspace(0, 100, 1000)

# Calculate c_e for each value of z
c_e_values = [c_ez(z) for z in z_values]

# =============================================================================
# # Plotting the results
# plt.plot(c_e_values,z_values)
# plt.xlabel('c_z')
# plt.ylabel('z')
# plt.title('Plot of c_e vs z')
# plt.grid(True)  # Add a grid for better readability
# plt.show()
# =============================================================================

