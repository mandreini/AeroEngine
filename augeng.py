# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 10:54:09 2019

@author: Matt
"""

import math
import characteristics
import engine

class turbojet(engine.Engine):
    # Initialise a turbojet. This is a child of the engine.Engine class.
    def __init__(self, mdot_core, boosterPR, TIT, PR_overall, compress_eta_is, eta_cc, expansion_eta_is, PR_itb, eta_itb, PR_ab, eta_ab):
        self.mdot_core = mdot_core
        self.boosterPR = boosterPR
        self.TIT = TIT
        self.PR_overall = PR_overall
        self.compress_eta_is = compress_eta_is
        self.eta_cc = eta_cc
        self.expansion_eta_is = expansion_eta_is
        self.PR_itb = PR_itb
        self.eta_itb = eta_itb
        self.PR_ab = PR_ab
        self.eta_ab = eta_ab
        
        self.M_inf = chars.M
        self.p_a = chars.p_a
        self.T_a = chars.T_a
        self.v_inf = chars.M * math.sqrt(chars.kappa_air*chars.R*chars.T_a)
        
        self.delta_s = lambda Cp, T2, T1, p2, p1: Cp*math.log(T2/T1, math.e) - chars.R*math.log(p2/p1, math.e)
        
    def thermocycle(self, mdot_itb=0, mdot_ab=0):
        # Perform the calculations of the thermodynamic cycle.
        # When mdot_itb or mdot_ab are equal to zero, then the combustion 
        # chamber 'runs,' but generates no change, so pressure and temperature
        # remain the same.
        T_a = self.T_a
        p_a = self.p_a
        M_inf = self.M_inf
        v_inf = self.v_inf
        T_ta = T_a*(1+(chars.kappa_air-1)/2*M_inf**2)
        p_ta = p_a*(1+(chars.kappa_air-1)/2*(M_inf**2))**(chars.kappa_air/(chars.kappa_air-1))
        T_t2 = T_ta
        p_t2 = p_a * (1 + chars.eta_is_intake*(chars.kappa_air-1)/2*M_inf**2)**(chars.kappa_air/(chars.kappa_air-1))
        
        # corrected mass flow
        mdot_core = self.mdot_core
        delta=p_ta/101325;
        theta=T_ta/288;
        mdot_core=mdot_core*delta/theta**.5
        deltas_2 = self.delta_s(chars.Cp_air, T_t2, T_a, p_t2, p_a)
        
        # LPC
        mdot_25 = mdot_core
        p_t25, T_t25, W_LPC, deltas_25 = self.compression(p_t2, T_t2, self.boosterPR, mdot_25)
        
        # HPC
        p_HPC = self.PR_overall / self.boosterPR
        mdot_3 = mdot_core
        p_t3, T_t3, W_HPC, deltas_3 = self.compression(p_t25, T_t25, p_HPC, mdot_3)
        
        # Combustion
        p_t4, T_t4, mdot_f, deltas_4 = self.combustion(mdot_3, T_t3, p_t3)
        mdot_4 = mdot_3 + mdot_f
        
        # HPT
        mdot_45 = mdot_4
        p_t45, T_t45, deltas_45 = self.expansion(p_t4, T_t4, W_HPC, mdot_45)
        
        # ITB - if mdot_itb = 0, then no changes to temperature or pressure
        T_t47 = T_t45 + self.eta_ab * mdot_itb * chars.LHV / mdot_45 / chars.Cp_gas
        p_t47 = p_t45 * self.PR_itb if mdot_itb != 0 else p_t45
        deltas_47 = self.delta_s(chars.Cp_gas, T_t47, T_t45, p_t47, p_t45)
        
        # LPT
        mdot_5 = mdot_45 + mdot_itb
        p_t5, T_t5, deltas_5 = self.expansion(p_t47, T_t47, W_LPC, mdot_5)
        
        # Afterburner - if mdot_ab = 0, then no changes to temperature or pressure
        T_t7 = T_t5 + self.eta_ab * mdot_ab * chars.LHV / mdot_5 / chars.Cp_gas
        p_t7 = p_t5 * self.PR_ab if mdot_ab != 0 else p_t5
        mdot_7 = mdot_5 + mdot_ab
        deltas_7 = self.delta_s(chars.Cp_gas, T_t7, T_t5, p_t7, p_t5)
        
        # Core Nozzle
        pt7ptcr = (1-1/chars.eta_nozzle*((chars.kappa_gas-1)/(chars.kappa_gas+1)))**(-chars.kappa_gas/(chars.kappa_gas-1))
        nozzlePressureRatio = p_t7/p_a
        nozzleStatus = pt7ptcr < nozzlePressureRatio  # True if choked, False if unchoked
        
        if nozzleStatus:
            T_8 = T_t7*(2/(chars.kappa_gas+1))
            p_8 = p_t7/pt7ptcr
            v_8 = math.sqrt(chars.kappa_gas*chars.R*T_8)
            
        else:
            p_8 = p_a
            T_8 = T_t7*(1-chars.eta_nozzle*(1-(p_8/p_t7)**((chars.kappa_gas-1)/chars.kappa_gas)))
            v_8 = math.sqrt(2*chars.Cp_gas*(T_t7-T_8))
        
        mdot_8 = mdot_7        
        rho_8 = p_8/T_8/chars.R
        A_8 = mdot_8/rho_8/v_8
        deltas_8 = self.delta_s(chars.Cp_gas, T_8, T_t7, p_8, p_t7)
        F_core = mdot_5*(v_8-self.v_inf) + A_8*(p_8-p_a) 
        TSFC = mdot_f / F_core

        v_eqj = A_8*(p_8-p_a)/mdot_8 + v_inf
        v_jeff_core = v_8 + A_8/mdot_8*(p_8-p_a)
        
        # P_gg
        T_tgg = T_t7
        p_tgg = p_t45*(((1-(1-(T_tgg/T_t45))/self.expansion_eta_is))**(chars.kappa_gas/(chars.kappa_gas-1)))
        W_gg = mdot_4*chars.Cp_gas*T_tgg*(1-(p_a/p_tgg)**((chars.kappa_gas-1)/chars.kappa_gas)) - 1/2*mdot_core*v_inf**2        
        
        # Efficiencies - these efficiencies replace all combustion related
        # values with Sum(values for each).
        eta_comb = (mdot_core * chars.Cp_gas * (T_t4-T_t3) + mdot_5*chars.Cp_gas*(T_t47-T_t45) + mdot_7*chars.Cp_gas*(T_t7-T_t5)) / ((mdot_f+mdot_itb+mdot_ab) * chars.LHV)
        eta_thdy = W_gg / chars.Cp_gas / (mdot_core * (T_t4-T_t3) + mdot_5*(T_t47-T_t45) + mdot_7*(T_t7-T_t5))
        eta_jetgen = (0.5*mdot_8*(v_jeff_core**2-v_inf**2))/W_gg
#        eta_thermal = eta_comb * eta_thdy * eta_jetgen
        eta_thermal = v_inf*(mdot_8*(v_jeff_core-v_inf)) / (0.5*mdot_8*(v_jeff_core**2-v_inf**2))
        
#        eta_prop = v_inf*(mdot_8*(v_jeff_core-v_inf)) / (0.5*mdot_8*(v_jeff_core**2-v_inf**2))
        eta_prop = eta_comb * eta_thdy * eta_jetgen
        eta_total = F_core*v_inf/(mdot_f+mdot_itb+mdot_ab)/chars.LHV
        
        effs = eta_comb, eta_thdy, eta_jetgen, eta_thermal, eta_prop, eta_total
        # The eta_thermal and eta_prop were switched in AET2, so I've switched 
        # them here to make them 'right,' but the formulas look a bit strange, 
        # so the commented out versions that were supposed to be the switched 
        # ones are left.

        deltases = [deltas_2, deltas_25, deltas_3, deltas_4, deltas_45, deltas_47, deltas_5, deltas_7, deltas_8, 0]
        Tes = [T_a, T_t2, T_t25, T_t3, T_t4, T_t45, T_t47, T_t5, T_t7, T_8]

        return ([F_core, TSFC, effs, Tes, deltases], locals())
    
chars = characteristics.Chars()
