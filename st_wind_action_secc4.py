#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 09:59:43 2023

@author: namnguyen
"""

import streamlit as st
from wind_action_secc4 import *
import matplotlib.pyplot as plt


c_dir = st.sidebar.number_input('Direction coeffition $c_{dir}$ =', value= 1.0, min_value=0.0, step=0.1, format="%.3f")#1.0
c_season = st.sidebar.number_input('Season coeffition $c_{season}$ =', value= 1.0, min_value=0.0, step=0.1, format="%.3f")#1.0
v_b0 = st.sidebar.number_input('base velocity #$v_{b0}[km/h]$ =', value= 27.0, min_value=0.0, step=1.0, format="%.3f")#27.0  # km/h
p = st.sidebar.number_input('p =', value= 0.01, min_value=0.0, step=0.01, format="%.3f")#0.01
K = st.sidebar.number_input('K=', value= 0.2, min_value=0.0, step=0.01, format="%.3f")#0.2
n = st.sidebar.number_input('n =', value= 0.5, min_value=0.0, step=0.01, format="%.3f")#0.5
rho = st.sidebar.number_input('Density of  air $\\rho [kg/m^3]$ =', value= 1.25, min_value=0.0, step=0.01, format="%.3f")
z_max = st.sidebar.number_input('Maximum height $z_{max} [m]$ =', value= 200.0, min_value=10.0, step=1.0, format="%.3f")#200.0  # m
z_0 = st.sidebar.number_input('Base height $z_{0} [m]$ =', value= 0.003, min_value=0.0, step=0.001, format="%.3f")#0.003  # m
z_min = st.sidebar.number_input('Minimum height $z_{min} [m]$ =', value= 1.0, min_value=0.0, step=0.10, format="%.3f")#1.0  # m
z_0II = st.sidebar.number_input('$z_{0II} [m]$ =', value= 0.005, min_value=0.0, step=0.001, format="%.3f")#symbols('z_0II')#0.005  # m
k_I = st.sidebar.number_input('$k_I=$',value= 1.00, min_value=0.0, step=0.01, format="%.3f")#1.00
A_ref = st.sidebar.number_input('$A_{ref} [m^2]=$',value= 800.0, min_value=1.0, step=1., format="%.3f")#800  # m2
c_d = st.sidebar.number_input('$c_d=$',value= 1.00, min_value=0.0, step=0.01, format="%.3f")#1.00
c_f = st.sidebar.number_input('$c_f=$',value= 1.55, min_value=0.0, step=0.01, format="%.3f")#1.55
c_0 = st.sidebar.number_input('$c_0=$',value= 1.00, min_value=0.0, step=0.01, format="%.3f")#1.00
#v_b=symbols('v_b', cls=Function)(v_b0)

z=st.number_input('$z=$',value= 1.00, min_value=0.0, step=0.01, format="%.3f")

st.title('**Eurocode 1: Actions on structures - Part 1-4: General actions - Wind actions**')

st.header('**Section 4: Wind velocity and velocity pressure**')

st.markdown('---')

st.write('Basic velocity pressure given in Expression (4.10), $v_{b}$')

v_b1 = v_b_func()
st.latex(latex(v_b1))

# Calculate v_b using v_b_func() with specific values and display in LaTeX format
v_b2 = v_b_func(c_dir=c_dir, c_season=c_season, v_b0=v_b0)
st.latex(latex(v_b2))


v_b = N(v_b2.doit(),3)

st.latex(latex(v_b)+f"(m/s)")

st.markdown('---')

st.write('The 10 minutes mean wind velocity having the probability p:')
c_prob1=c_prob_func()
st.latex(latex(c_prob1))

c_prob2=c_prob_func(K=K,p=p,n=n)

st.latex(latex(c_prob2))

c_prob=N(c_prob2.doit(),3)
st.latex(latex(c_prob))

st.markdown('---')

st.write('Terrain factor depending on the roughness length $z_0$')

k_r1=k_r_func()
st.latex(latex(k_r1))

k_r2=k_r_func(z_0=z_0,z_0II=z_0II)

st.latex(latex(k_r2))

k_r=N(k_r2.doit(),3)
st.latex(latex(k_r))

st.markdown('---')

st.write('Basic velocity pressure $q_b$')

q_b1=q_b_func()
st.latex(latex(q_b1))

q_b2=q_b_func(rho=rho, v_b=v_b.rhs)

st.latex(latex(q_b2))

q_b=N(q_b2.doit(),3)
st.latex(latex(q_b))

st.markdown('---')

st.write('The turbulent component of wind velocity has a mean value of 0 and a standard deviation $\sigma_v$')

sigma_v1=sigma_v_func()
st.latex(latex(sigma_v1))

sigma_v2=sigma_v_func(k_r=k_r.rhs, v_b=v_b.rhs, k_I=k_I)

st.latex(latex(sigma_v2))

sigma_v=N(sigma_v2.doit(),3)
st.latex(latex(sigma_v))

st.markdown('---')

st.write('$c_r(z)$')
c_r1=c_r_func()
st.latex(latex(c_r1))

c_r2=c_r_func(z=z, z_min=z_min, z_max=z_max, z_0=z_0, k_r=k_r.rhs,UE=True)

st.latex(latex(c_r2))

c_r3=c_r_func(z=z, z_min=z_min, z_max=z_max, z_0=z_0, k_r=k_r.rhs,UE=False)

st.latex(latex(c_r3))

c_r=N(c_r2.doit(),3)

st.latex(latex(c_r))

st.markdown('---')

st.write("$v_m$")
v_m1=v_m_func()
st.latex(latex(v_m1))

v_m2=v_m_func(z=z,c_r=c_r.rhs, c_0=c_0,v_b=v_b.rhs)

st.latex(latex(v_m2))

v_m=N(v_m2.doit(),3)

st.latex(latex(v_m)+f'(m/s^2)')


st.markdown('---')
st.write('$I_v$')

I_v1=I_v_func()
st.latex(latex(I_v1))

I_v =I_v_func(z=z,z_min=z_min, z_max=z_max, z_0=z_0, sigma_v=sigma_v.rhs, v_m=v_m.rhs, UE=False)

#st.latex(latex(I_v2))
#I_v=I_v2.doit()
st.latex(latex(I_v))

st.markdown('---')

st.write("$q_p(z)$")
q_p1=q_p_func()
st.latex(latex(q_p1))

q_p2=q_p_func(z=z, I_v=I_v.rhs, rho=rho, v_m=v_m.rhs)

st.latex(latex(q_p2))

q_p=N(q_p2.doit(),3)

st.latex(latex(q_p))


st.markdown('---')

st.write("$c_e(z)$")
c_e1=c_e_func()
st.latex(latex(c_e1))

c_e2=c_e_func(z=z,q_p=q_p.rhs, q_b=q_b.rhs)

st.latex(latex(c_e2))

c_e=N(c_e2.doit(),3)

st.latex(latex(c_e))

# Generate 1000 points linearly spaced between 0 and 100
z = np.linspace(0, 100, 1000)

# Calculate e_e using the c_e_func
c_e = c_e_func(z,q_p=q_p.rhs, q_b=q_b.rhs)

st.write(c_e)

# =============================================================================
# # Plotting the results
# plt.plot(z, c_e)
# plt.xlabel('z')
# plt.ylabel('c_e')
# plt.title('Plot of c_e_func')
# 
# # Display the plot in Streamlit
# st.pyplot()
# 
# =============================================================================




