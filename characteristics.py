# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 10:15:19 2019

@author: Matt
"""
# This tracks some of the flight parameters. It was based on the parameters 
# that shouldn't change between engines and then it did.

class Chars(object):
    def __init__(self):

        self.M = 0.8
        self.h = 10566
        self.eta_m = 0.99
        self.CPR = 0.96
        self.eta_nozzle = 0.99
        self.eta_is_intake = 0.99 
        self.T_a = 220
        self.p_a = 23842
        self.R = 287
        self.LHV = 43e6
        self.Cp_air = 1000
        self.kappa_air = 1.4
        self.Cp_gas = 1150
        self.kappa_gas = 1.33
