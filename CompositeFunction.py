# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 12:08:51 2023

@author: cfcpc2
"""
import numpy as np
from sympy import symbols, Eq, Function,UnevaluatedExpr, Mul, Rational, sqrt
class CompositeFunction:
	def __init__(self, F):
		self.F = F

	def apply(self, g, x):
		result_of_g = g(x)
		result_of_F = self.F(result_of_g)
		return result_of_F


def g(x):
	#x=UnevaluatedExpr(x)
	return  x + 2


def F(y):

	return y * 3


x_value=5

compositefunc= CompositeFunction()

result_g=compositefunc.g(x_value)

print(result_g )