#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 09:59:43 2023

@author: namnguyen
"""
import pandas as pd
import numpy as np
import streamlit as st
from wind_action_secc4 import *
import matplotlib.pyplot as plt
init_printing()


c_dir = st.sidebar.number_input('Enter c_dir', value=1.0, min_value=0.0, step=0.1, format="%.3f")
c_dir = st.sidebar.number_input('Direction coefficient $c_{dir}=$', value= 1.0, min_value=0.0, step=0.1, format="%.3f")#1.0
c_season = st.sidebar.number_input('Season coeffition $c_{season} =$', value= 1.0, min_value=0.0, step=0.1, format="%.3f")#1.0
v_b0 = st.sidebar.number_input('base velocity $v_{b0}[km/h] =$', value= 27.0, min_value=0.0, step=1.0, format="%.3f")#27.0  # km/h
p = st.sidebar.number_input('p =', value= 0.01, min_value=0.0, step=0.01, format="%.3f")#0.01
K = st.sidebar.number_input('K=', value= 0.2, min_value=0.0, step=0.01, format="%.3f")#0.2
n = st.sidebar.number_input('n =', value= 0.5, min_value=0.0, step=0.01, format="%.3f")#0.5
rho = st.sidebar.number_input('Density of  air $\\rho [kg/m^3] =$', value= 1.25, min_value=0.0, step=0.01, format="%.3f")
z_max = st.sidebar.number_input('Maximum height $z_{max} [m]$ =', value= 200.0, min_value=10.0, step=1.0, format="%.3f")#200.0  # m
z_0 = st.sidebar.number_input('Base height $z_{0} [m]$ =', value= 0.3, min_value=0.0, step=0.001, format="%.3f")#0.003  # m
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

st.markdown('**4.1 Basis for calculation**')

st.markdown('**4.2 Basic values**')

st.write('The basic wind velocity shall be calculated from Expression (4.1). ')
v_b1 = v_b_func()
st.latex(latex(v_b1))

# Calculate v_b using v_b_func() with specific values and display in LaTeX format
v_b2 = v_b_func(c_dir=c_dir, c_season=c_season, v_b0=v_b0)
st.latex(latex(v_b2))


v_b = N(v_b2.doit(),3)

st.latex(latex(v_b)+f"(m/s)")

st.markdown(f"""
where: 
	
$v_b$ is the basic wind velocity, defined as a function of wind direction and time of year at $10$ m 
above ground of terrain category II 

$v_{{b0}}$ is the fundamental value of the basic wind velocity, see (1)P 

$c_{{dir}}$ is the directional factor, see Note 2. 

$c_{{season}}$ is the season factor, see Note 3.
""")

st.markdown('---')

st.write('The 10 minutes mean wind velocity having the probability p:')
c_prob1=c_prob_func()
st.latex(latex(c_prob1))

c_prob2=c_prob_func(K=K,p=p,n=n)

st.latex(latex(c_prob2))

c_prob=N(c_prob2.doit(),3)
st.latex(latex(c_prob))
st.markdown(f"""
where: 
	
$$K$$ is the shape parameter depending on the coefficient of variation of the extreme-value distribution. 
$$n$$ is the exponent. 
NOTE 5 The values for $K$ and $n$ may be given in the National Annex. The recommended values are $0,2$ 
for $K$ and $0,5$ for $n$.
""")

st.markdown('---')
st.markdown('**4.3 Mean wind** ')
st.markdown('**4.3.1 Variation with height**')

st.write('The mean wind velocity $v_m(z)$ at a height z above the terrain depends on the terrain roughness and orography and on the basic wind velocity, $v_b$, and should be determined using Expression (4.3)')

st.write("$v_m$")


k_r2=k_r_func(z_0=z_0,z_0II=z_0II)

k_r=N(k_r2.doit(),3)

c_r3=c_r_func(z=z, z_min=z_min, z_max=z_max, z_0=z_0, k_r=k_r.rhs,UE=False)


c_r=N(c_r3.doit(),3)

v_m1=v_m_func()
st.latex(latex(v_m1))

v_m2=v_m_func(z=z,c_r=c_r.rhs, c_0=c_0,v_b=v_b.rhs)

st.latex(latex(v_m2))

v_m=N(v_m2.doit(),3)

st.latex(latex(v_m)+f'(m/s^2)')

st.markdown(f"""

where: 
$$c_r(z)$$ is the roughness factor, given in 4.3.2 

$$c_o(z)$$ is the orography factor, taken as 1,0 unless otherwise specified in 4.3.3
""")

st.markdown('**4.3.2 Terrain roughness**')


st.write('The roughness factor $c_r(z)$')


k_r2=k_r_func(z_0=z_0,z_0II=z_0II)

k_r=N(k_r2.doit(),3)



c_r1=c_r_func()
st.latex(latex(c_r1))

c_r2=c_r_func(z=z, z_min=z_min, z_max=z_max, z_0=z_0, k_r=k_r.rhs,UE=True)

st.latex(latex(c_r2))

c_r3=c_r_func(z=z, z_min=z_min, z_max=z_max, z_0=z_0, k_r=k_r.rhs,UE=False)

st.latex(latex(c_r3))

c_r=N(c_r2.doit(),3)

st.latex(latex(c_r))

st.markdown(f"""

where:
	 
$z_0$ is the roughness length 

""")


st.write('Terrain factor depending on the roughness length $z_0$')

k_r1=k_r_func()
st.latex(latex(k_r1))

k_r2=k_r_func(z_0=z_0,z_0II=z_0II)

st.latex(latex(k_r2))

k_r=N(k_r2.doit(),3)
st.latex(latex(k_r))

st.markdown(f"""

where:
$$z_{{0,II}}$$ = 0,05 m (terrain category II, Table 4.1) 

$z_{{min}}$ is the minimum height defined in Table 4.1 

z_{{max}} is to be taken as $200$ m 

""")

st.markdown('**4.3.3 Terrain orography**')


st.markdown('**4.3.4 Large and considerably higher neighbouring structures**')



st.markdown('**4.3.5 Closely spaced buildings and obstacles**')



st.markdown('---')
st.markdown('**4.4 Wind turbulence**')
st.write('The turbulent component of wind velocity has a mean value of $0$ and a standard deviation $\sigma_v$')

sigma_v1=sigma_v_func()
st.latex(latex(sigma_v1))

sigma_v2=sigma_v_func(k_r=k_r.rhs, v_b=v_b.rhs, k_I=k_I)

st.latex(latex(sigma_v2))

sigma_v=N(sigma_v2.doit(),3)
st.latex(latex(sigma_v))


st.markdown('---')
st.write(' The turbulence intensity $I_v$ at height $z$ is defined as')
S
I_v1=I_v_func()
st.latex(latex(I_v1))

I_v =I_v_func(z=z,z_min=z_min, z_max=z_max, z_0=z_0, sigma_v=sigma_v.rhs, v_m=v_m.rhs, UE=False)

#st.latex(latex(I_v2))

st.latex(latex(I_v))


st.markdown(f"""

where: 
	
$k_I$ is the turbulence factor. The value of $k_I$ may be given in the National Annex. The recommended value 
for $k_I$ is $1,0$. 

$c_o$ is the orography factor as described in 4.3.3 

$z_0$ is the roughness length, given in Table 4.1 

""")

st.markdown('---')

st.markdown('**4.5 Peak velocity pressure**')

st.write("The peak velocity pressure q $q_p(z)$")

q_b2=q_b_func(rho=rho, v_b=v_b.rhs)

q_b=N(q_b2.doit(),3)

q_p1=q_p_func()
st.latex(latex(q_p1))

q_p2=q_p_func(z=z, I_v=I_v.rhs, rho=rho, v_m=v_m.rhs)

st.latex(latex(q_p2))

q_p=N(q_p2.doit(),3)

st.latex(latex(q_p))

st.markdown(f"""

where: 
	
$\\rho$ is the air density, which depends on the altitude, temperature and barometric pressure to be 
expected in the region during wind storms 

$c_e(z)$ is the exposure factor given in Expression (4.9) 

""")

c_e1=c_e_func()
st.latex(latex(c_e1))

c_e2=c_e_func(z=z,q_p=q_p.rhs, q_b=q_b.rhs)

st.latex(latex(c_e2))

c_e=N(c_e2.doit(),3)

st.latex(latex(c_e))

st.markdown('$q_b$ is the basic velocity pressure given in Expression (4.10)')

q_b1=q_b_func()
st.latex(latex(q_b1))

q_b2=q_b_func(rho=rho, v_b=v_b.rhs)

st.latex(latex(q_b2))

q_b=N(q_b2.doit(),3)
st.latex(latex(q_b))



# Generate 1000 points linearly spaced between 0 and 100
z_values = np.linspace(0, 100, 1000)
c_ez=c_ez(z,c_dir=c_dir,c_season=c_season,v_b0=v_b0,p=p,K=K,n=n,rho=rho,z_max=z_max,z_0=z_0,z_min=z_min,z_0II=z_0II, k_I=k_I, A_ref=A_ref, c_d=c_d, c_f=c_f,c_0=c_0)
# Calculate c_e for each value of z
c_ez_values = [c_ez for z in z_values]

# Plotting the results
plt.plot( c_e_values,z_values)
plt.xlabel('$c_e(z)$')
plt.ylabel('z')
plt.title(f'''Plot of  $z$ vs $c_e(z)$''')

# Display the plot in Streamlit
st.pyplot()
showPyplotGlobalUse = False
st.set_option('deprecation.showPyplotGlobalUse', False)



