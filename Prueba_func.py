# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 10:37:53 2023

@author: cfcpc2
"""

import streamlit as st
from sympy import symbols, Eq, Function,UnevaluatedExpr, Mul, Rational, sqrt, Pow
from sympy import *#Piecewise, nan, N, And, log
import matplotlib.pyplot as plt



a,b,c,d, x =symbols('a b c d x')


a=UnevaluatedExpr(a)
b=UnevaluatedExpr(b)
c=UnevaluatedExpr(c)

x=UnevaluatedExpr(x)

g= symbols('g', cls=Function)(x)
f= symbols('f', cls=Function)(x)

def g_func(x=x, **db):
	db = {eval(key): value for key, value in db.items()}
	expr= d*x**3
	#expr= expr.subs(db)
	return Eq(g, expr, evaluate= False)

def f_func(x=x, **db):
	db = {eval(key): value for key, value in db.items()}
	expr= g+a*x**2+b*x+c
	#expr= expr.subs(db)
	return Eq(f, expr, evaluate= False)




a_val= st.number_input("input value of $a$:", value=1.0)#, min_value=0.0,step=0.1, format="%.2f")
b_val= st.number_input("input value of $b$:", value=2.0)#, min_value=0.0,step=0.1, format="%.2f")
c_val= st.number_input("input value of $c$:", value=3.0)#, min_value=0.0, step=0.1,format="%.2f")
d_val= st.number_input("input value of $d$:", value=4.0)#, min_value=0.0, step=0.1,format="%.2f")
db={}
db={'a':a_val, 'b':b_val,'c':c_val, 'd':d_val,'g':g_func()}

#st.latex(latex(f_func()))

#db={'a':1, 'b':2,'c':3}
st.latex(latex(g_func()))
st.latex(latex(f_func()))
st.latex(latex(f_func().subs(db)))


x_val=st.number_input("input value of $x$:", value=10.0,step=1.0)#, min_value=0.0, step=0.1,format="%.2f")
st.latex(latex(f_func().subs(x,x_val).subs(db)))
st.latex(latex(f_func().subs(x,x_val).subs(db).doit()))






















