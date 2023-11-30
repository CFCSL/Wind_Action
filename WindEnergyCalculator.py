# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 12:25:35 2023

@author: cfcpc2
"""

import numpy as np

class WindEnergyCalculator:
    def __init__(self, params):
        self.c_dir, self.c_season, self.v_b0, self.p, self.K, self.n, self.rho, \
        self.z_max, self.z_0, self.z_min, self.z_0II, self.k_I, self.A_ref, \
        self.c_d, self.c_f, self.c_0, self.c_pe, self.c_pi, self.c_fr, self.A_fr = params

    def v_b(self):
        exp = self.c_dir * self.c_season * self.v_b0
        return exp

    def c_prob(self):
        exp = ((1 - self.K * np.log(-np.log(1 - self.p))) / (1 - self.K * np.log(-np.log(0.98)))) ** self.n
        return exp

    def k_r(self):
        exp = 0.19 * (self.z_0 / self.z_0II) ** 0.07
        return exp

    def q_b(self):
        exp = 0.5 * self.rho * self.v_b() ** 2
        return exp

    def sigma_v(self):
        exp = self.k_r() * self.v_b() * self.k_I
        return exp

    def c_r(self, z):
        condlist = [np.logical_and(z >= self.z_min, z <= self.z_max), z <= self.z_min]
        funclist = [self.k_r() * np.log(z / self.z_0), self.k_r() * np.log(self.z_min / self.z_0)]
        if z != 0 and self.z_min != 0 and self.z_0 != 0:
            funclist = [self.k_r() * np.log(z / self.z_0), self.k_r() * np.log(self.z_min / self.z_0)]
        else:
            pass
        exp = np.piecewise(z, condlist, funclist)
        return exp

    def v_m(self, z):
        exp = self.c_r(z) * self.c_0 * self.v_b()
        return exp

    def I_v(self, z):
        condlist = [z >= self.z_min, z <= self.z_min]
        funclist = [self.sigma_v() / self.v_m(z), self.sigma_v() / self.v_m(self.z_min)]
        exp = np.piecewise(z, condlist, funclist)
        return exp

    def q_p(self, z):
        exp = (1 + 7 * self.I_v(z)) * 0.5 * self.rho * self.v_m(z) ** 2
        return exp

    def c_ez(self, z):
        val = self.q_p(z) / self.q_b()
        return val

    def W_e(self, z):
        exp = self.q_p(z) * self.c_pe
        return exp

    def W_i(self, z):
        exp = self.q_p(z) * self.c_pi
        return exp

    def F_w(self, z):
        exp = self.c_fr * self.c_d * self.c_f * self.q_p(z) * self.A_ref
        return exp

    def F_fr(self, z):
        exp = self.c_fr * self.q_p(z) * self.A_fr
        return exp

# Example usage
params = (0.5, 0.8, 10, 0.2, 1.5, 4, 1.225, 100, 0.01, 0.005, 0.005, 0.05, 1000, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 100)
wind_calculator = WindEnergyCalculator(params)

# Example usage of functions
z_value = 50
result_W_e = wind_calculator.W_e(z_value)
result_F_fr = wind_calculator.F_fr(z_value)

print(f"W_e({z_value}) = {result_W_e}")
print(f"F_fr({z_value}) = {result_F_fr}")