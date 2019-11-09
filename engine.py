# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 12:38:47 2019

@author: Matt
"""
import math
import matplotlib.pyplot as plt
import numpy as np
import characteristics

class Engine(object):
    
    def __init__(self, BPR, mdot_core, FPR, HPCPR, TIT, PR_overall, compress_eta_is, eta_cc, expansion_eta_is):
        self.BPR = BPR
        self.mdot_core = mdot_core
        self.FPR = FPR
        self.HPCPR = HPCPR
        self.TIT = TIT
        self.PR_overall = PR_overall
        self.compress_eta_is = compress_eta_is
        self.expansion_eta_is = expansion_eta_is
        self.eta_cc = eta_cc
        
        self.M_inf = chars.M
        self.p_a = chars.p_a
        self.T_a = chars.T_a
        self.v_inf = chars.M * math.sqrt(chars.kappa_air*chars.R*chars.T_a)
        
        self.delta_s = lambda Cp, T2, T1, p2, p1: Cp*math.log(T2/T1, math.e) - chars.R*math.log(p2/p1, math.e)
        
    def set_flight_condition(self, M_inf, p_a, T_a):
        # This is partially obsolete, and ends up causing problems with 
        # characteristics.py, and should be fixed for AET5 if variable flight conditions
        self.M_inf = M_inf
        self.p_a = p_a
        self.T_a = T_a
        self.v_inf = M_inf * math.sqrt(chars.kappa_air*chars.R*self.T_a)
        
    def T_sDiagram(self, title, Tes, deltases):
        # Plot the T-s diagram. The commented lines are for using a fan and
        # for generating a logarithmic curve for the combustion chambers.
        # However, the curve was not working properly due to Python variable
        # handling requiring an extra step, so I commented it out to ignore
        # creating the curved lines. I will return and fix them if they 
        # comment on AET3 saying that they want the curves!
#        Tes = self.Tes#[:-2]
#        deltases = self.deltases#[:-2]
        s_vals = [sum(deltases[:i]) for i in range(len(deltases))]

#        isenfit = np.polyfit(s_vals[4:6], np.log(Tes[4:6]), 1)
#        xs = np.linspace(s_vals[4], s_vals[5])
#        ys = np.exp(isenfit[1]) * np.exp(isenfit[0] * xs)
#        xs, ys = list(xs), list(ys)        
        
#        s_vals = s_vals[:5] + xs + s_vals[6:]
#        Tes = Tes[:5] + ys + Tes[6:]

        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.plot(s_vals, Tes)
        plt.xlabel('s [J/K]')
        plt.ylabel('T [K]')
        plt.title(title)

        # This is to label the plot for when there is a logarithmic line, as 
        # the logarithmic line will then be the length of the original T-s 
        # stations + 50 as the curve takes 50 points to generate it.
#        stations = ['' , '', '2 [%i]' % round(Tes[2]), '21 [%i]' % round(Tes[3]), '25 [%i]' % round(Tes[4]), '3 [%i]' % round(Tes[5]), '4 [%i]' % round(Tes[-5]), '45 [%i]' % round(Tes[-4]), 'gg [%i]' % round(Tes[-3]), '5 [%i]' % round(Tes[-2]), '8 [%i]' % round(Tes[-1])]
#        i = 0
#        for xy in zip(s_vals[:5]+s_vals[5+len(xs)-1:], Tes[:5]+Tes[5+len(ys)-1:]):
#            i = i + 1
#            ax.annotate('%s' % stations[i], xy=xy, textcoords='data')
#        
#        plt.grid()
        plt.savefig('%s.png' % title)
#        plt.show()
        
    def compression(self, p_t1, T_t1, Pi, mdot):
        # Compute a compression (i.e. Fan, LPC, HPC) stage.
        # Pi = pressure ratio
        p_t2 = p_t1 * Pi
        T_t2 = T_t1*(1+(1/self.compress_eta_is)*(((Pi)**((chars.kappa_air-1)/chars.kappa_air))-1))
        delta_s = self.delta_s(chars.Cp_air, T_t2, T_t1, p_t2, p_t1)
        W = mdot*chars.Cp_air*(T_t2-T_t1)
        return p_t2, T_t2, W, delta_s
    
    def expansion(self, p_t1, T_t1, W_compress, mdot):
        # Compute an expansion (i.e. HPT, LPT) stage
        # W_compress is work done by related compression side, i.e. LPT = LPC + fan if necessary
        T_t2 = T_t1-W_compress/(chars.eta_m*mdot*chars.Cp_gas)
        p_t2 = p_t1*(((1-(1-(T_t2/T_t1))/self.expansion_eta_is))**(chars.kappa_gas/(chars.kappa_gas-1)))
        delta_s = self.delta_s(chars.Cp_gas, T_t2, T_t1, p_t2, p_t1)
        
        return p_t2, T_t2, delta_s
    
    def combustion(self, mdot, T_t1, p_t1):
        # Compute a combustion stage
        mdot_f = mdot*chars.Cp_gas*(self.TIT-T_t1)/(self.eta_cc*chars.LHV)
        p_t2 = p_t1 * chars.CPR
        delta_s = self.delta_s(chars.Cp_gas, self.TIT, T_t1, p_t2, p_t1)
        
        return p_t2, self.TIT, mdot_f, delta_s
        
    def do_math(self):        
        # inlet
        T_a = self.T_a
        p_a = self.p_a
        M_inf = self.M_inf
        v_inf = self.v_inf
        T_ta = T_a*(1+(chars.kappa_air-1)/2*M_inf**2)
        p_ta = p_a*(1+(chars.kappa_air-1)/2*(M_inf**2))**(chars.kappa_air/(chars.kappa_air-1))
        T_t2 = T_ta
        p_t2 = p_a * (1 + chars.eta_is_intake*(chars.kappa_air-1)/2*M_inf**2)**(chars.kappa_air/(chars.kappa_air-1))
        if M_inf == 0:
            p_t2 = p_ta * chars.eta_is_intake  # for the JT8D on ground
                
        # corrected mass flow
        mdot_core = self.mdot_core
        delta=p_ta/101325;
        theta=T_ta/288;
        mdot_core=mdot_core*delta/theta**.5
        mdot_bypass = mdot_core * self.BPR
        mdot_real = mdot_core + mdot_bypass
    
        mdot_2 = mdot_real
        deltas_2 = self.delta_s(chars.Cp_air, T_t2, T_a, p_t2, p_a)
            
        # Fan
        mdot_21 = mdot_core
        mdot_13 = mdot_bypass
        p_t21, T_t21, W_fan, deltas_21 = self.compression(p_t2, T_t2, self.FPR, mdot_21+mdot_13)
        
        # LPC
        p_LPC = self.PR_overall / self.FPR / self.HPCPR
        mdot_25 = mdot_core
        p_t25, T_t25, W_LPC, deltas_25 = self.compression(p_t21, T_t21, p_LPC, mdot_25)
        
        # HPC
        mdot_3 = mdot_core
        p_t3, T_t3, W_HPC, deltas_3 = self.compression(p_t25, T_t25, self.HPCPR, mdot_3)
        
        # Combustion
        p_t4, T_t4, mdot_f, deltas_4 =  self.combustion(mdot_3, T_t3, p_t3)
        mdot_4 = mdot_3 + mdot_f
        
        # HPT
        mdot_45 = mdot_4
        p_t45, T_t45, deltas_45 = self.expansion(p_t4, T_t4, W_HPC, mdot_45)
        
        # LPT
        mdot_5 = mdot_45
        p_t5, T_t5, deltas_5 = self.expansion(p_t45, T_t45, W_LPC+W_fan, mdot_5)
        
        # Core Nozzle
        T_t7 = T_t5
        p_t7 = p_t5
        deltas_7 = 0
        mdot_7 = mdot_5
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
            
        # Bypass Nozzle
        if self.BPR != 0:
            T_t16 = T_t21
            p_t16 = p_t21
            deltas_16 = 0
            mdot_16 = mdot_13
            pt16ptcr = (1-1/chars.eta_nozzle*((chars.kappa_air-1)/(chars.kappa_air+1)))**(-chars.kappa_air/(chars.kappa_air-1))
            bypassPressureRatio = p_t16/p_a
            bypassStatus = pt16ptcr < bypassPressureRatio  #  True if choked, False if unchoked
            
            if bypassStatus:
                T_18 = T_t16*(2/(chars.kappa_air+1))
                p_18 = p_t16/pt16ptcr
                v_18 = math.sqrt(chars.kappa_air*chars.R*T_18)
            else:
                p_18 = p_a
                T_18 = T_t16*(1-chars.eta_nozzle*(1-(p_18/p_t16)**((chars.kappa_air-1)/chars.kappa_air)))
                v_18 = math.sqrt(2*chars.Cp_air*(T_t16-T_18))
            
            mdot_18 = mdot_16        
            rho_18 = p_18/T_18/chars.R
            A_18 = mdot_18/rho_18/v_18
            deltas_18 = self.delta_s(chars.Cp_air, T_18, T_t16, p_18, p_t16)
            F_bypass = mdot_bypass*(v_18-v_inf) + A_18*(p_18-p_a)  
        else:
            F_bypass = 0
            
        # Thrust of both
        F_total = F_core + F_bypass
        TSFC = mdot_f / F_total
        
        v_eqj = A_8*(p_8-p_a)/mdot_8 + v_inf
        v_eqj_bp = A_18*(p_18-p_a)/mdot_bypass + v_inf
        v_jeff_core = v_8 + A_8/mdot_8*(p_8-p_a)
        v_jeff_bp = v_18+A_18/mdot_18*(p_18-p_a)
        
        # P_gg
        T_tgg = T_t45 - (W_LPC+W_fan*mdot_core/(mdot_core+mdot_bypass))/mdot_4/chars.Cp_gas/chars.eta_m
        p_tgg = p_t45*(((1-(1-(T_tgg/T_t45))/self.expansion_eta_is))**(chars.kappa_gas/(chars.kappa_gas-1)))
        W_gg = mdot_4*chars.Cp_gas*T_tgg*(1-(p_a/p_tgg)**((chars.kappa_gas-1)/chars.kappa_gas)) - 1/2*mdot_core*v_inf**2        
        deltas_gg = self.delta_s(chars.Cp_gas, T_tgg, T_t45, p_tgg, p_t45)        

        eta_comb = mdot_core * chars.Cp_gas * (T_t4-T_t3) / (mdot_f * chars.LHV)
        eta_thdy = W_gg / (mdot_core * chars.Cp_gas * (T_t4-T_t3))
        
        eta_jetgen = (0.5*mdot_8*(v_jeff_core**2-v_inf**2) + 0.5*mdot_18*(v_jeff_bp**2-v_inf**2))/W_gg
        eta_thermal = eta_comb * eta_thdy * eta_jetgen
        
        eta_prop = v_inf*(mdot_8*(v_jeff_core-v_inf) + mdot_18*(v_jeff_bp-v_inf)) / (0.5*mdot_8*(v_jeff_core**2-v_inf**2) + 0.5*mdot_18*(v_jeff_bp**2-v_inf**2))
        eta_total = F_total*v_inf/mdot_f/chars.LHV
        
        deltases = [deltas_2, deltas_21, deltas_25, deltas_3, deltas_4, deltas_45, deltas_gg, deltas_5, deltas_7, deltas_8, deltas_16, deltas_18]
        Tes = [T_a, T_t2, T_t21, T_t25, T_t3, T_t4, T_t45, T_tgg, T_t5, T_8, T_t16, T_18]
        dtes = [Tes[i+1]-Tes[i] for i in range(len(Tes[:-1]))]
        
        effs = eta_comb, eta_thdy, eta_jetgen, eta_thermal, eta_prop, eta_total
        return (F_total, TSFC, effs, Tes, deltases), locals()
        
chars = characteristics.Chars()