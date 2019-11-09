# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 10:15:00 2019

@author: Matt
4"""
# Main file, run this to then do the things. Commented lines are for previous 
# assignmnents as they help to verify that changes made to engine.py have not
# ruined the thermodynamic cycle calculations.
import math
import matplotlib.pyplot as plt
import numpy as np
import engine
import augeng
import characteristics
import plotting
import numpy
from scipy.optimize import minimize

#def get_val(value):
#    return jtvals[value], cfmvals[value], leapvals[value]
        
# Assignment 1
#JT8D = engine.Engine(1.62,	81.2,	1.9,	3.5,	1303.15,	16.9,	0.83,	0.985,	0.88)
#CFM96 = Engine([5.2,	52.2,	1.6,	8,	1270,	32.5,	0.88,	0.91,	0.995, 1.55])
#Leap1B = Engine([8.6,	41.5,	1.5,	10,	1450,	39,	0.92,	0.92,	0.995, 1.75])
#GE90 = engine.Engine(8,	57.6,	1.6,	19,	1450,	42.56,	0.9,	0.99,	0.92)
#
#JT8D.set_flight_condition(M_inf=0, p_a=101325, T_a=288.15)
#jtans, jtvals = JT8D.do_math()
#cfmans, cfmvals = CFM96.do_math()
#leapans, leapvals = Leap1B.do_math()
#
#GE90.set_flight_condition(0.8, 22632, 216)
#geans, gevals = GE90.do_math()
#
#final_vals = ([jtans[0], cfmans[0], leapans[0]],
#              [jtans[1], cfmans[1], leapans[1]])
        
# Assignment 2+3
JT8D = engine.Engine(1.62, 90.2, 1.9, 3.5, 1150, 17, 0.85, 0.985, 0.88) 
Leap1B = engine.Engine(8.6, 50, 1.5, 10, 1450, 40, 0.92, 0.995, 0.92)
#
#JT8D.set_flight_condition(M_inf=0, p_a=101325, T_a=288.15)
jtans, jtvals = JT8D.do_math()
leapans, leapvals = Leap1B.do_math()
jteffs, leapeffs = jtans[2], leapans[2]

# Assignment 1 T-s diagrams
#JT8D.T_sDiagram('JT8D engine cycle', jtans[3], jtans[4])
#CFM96.T_sDiagram('CFM96 engine cycle')
#Leap1B.T_sDiagram('Leap-1B engine cycle')
#plotting.plot_all(JT8D, None, Leap1B)

#JT8D.T_sDiagram('JT8D engine cycle')
#Leap1B.T_sDiagram('Leap1B engine cycle')

turbojet = augeng.turbojet(75, 2.1, 1700, 21, 0.85, 0.995, 0.9, 1-0.05, 0.98, 1-0.03, 0.98)
dryans, dryvals = turbojet.thermocycle()
dry_thrust = dryans[0]
aug_thrust = dry_thrust * 1.4

# max number convergence exit: mdot_ab = 1.93; my F = 89.895 kN
# loop for ab. This uses scipy.optimize.minimize() (aka python's fmincon) to 
# look for the optimum mass flow to meet the augmented thrust requirement.
# abfun (and itbfun) are functions that take a mass flow and caclulate the 
# difference between the thermocycle with that much mass in the 2nd burner and
# the augmented thrust. This lets minimize() find the closest value.
abfun = lambda mdot_ab: numpy.abs(aug_thrust - turbojet.thermocycle(0, mdot_ab)[0][0])
info_ab = minimize(abfun, (1))
mdot_ab_opt = info_ab.x

# loop for itb
# This does the same for ab but now with itb mass flow. A mass flow of zero for 
# itb means that the itb changes nothing, and same for ab. 
itbfun = lambda mdot_itb: numpy.abs(aug_thrust - turbojet.thermocycle(mdot_itb, 0)[0][0])
info_itb = minimize(itbfun, (1))
mdot_itb_opt = info_itb.x

abans, abvals = turbojet.thermocycle(0, mdot_ab_opt)
itbans, itbvals = turbojet.thermocycle(mdot_itb_opt, 0)

# Assignment 3 T-s diagrams
#turbojet.T_sDiagram('dry', dryans[-2], dryans[-1])
#turbojet.T_sDiagram('ab', abans[-2], abans[-1])
#turbojet.T_sDiagram('itb', itbans[-2], itbans[-1])

slopes = plotting.plot_all(dryans[-2:]+['dry'], abans[-2:]+['ab'], itbans[-2:]+['itb'])