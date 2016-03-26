# -*- coding: utf-8 -*-
"""
Created on Sun Mar 13 09:53:57 2016

@author: Carlo
"""

from __future__ import division, print_function

import numpy as np
import math as m
import matplotlib.pyplot as plt

inf = float('inf')
T = 300 #K
k = 8.617E-5 # ev/K
kT = k * T  # eV  (0.02585 @ 300 K)
eps_0 = 8.854*10**-18 # q/cm3 / V/nm2
q = 1.6E-19

def V0(E_gp, Del_Ec, FE_p, EF_N):
    return E_gp + Del_Ec - FE_p - EF_N

def Depletion_Width(eps0, eps1, N_a, N_D, V0):
    x = (2*eps0*V0/(q*N_a*N_D*(N_D + (eps0/eps1)*N_a)))**0.5
    return x
  
'''
MAIN LOOP
'''

def calc_bands(layers, points):
    layer_list = [layer for layer in layers]
       
    ### Generate x_axis for each layer
    ### Generate eps for each point in layer
    x_list = []
    eps = []
    previous_pos = 0
    tot_thickness = 0
    for layer in layer_list:
        layer_x = []
        eps_x = []
        tot_thickness += layer.thickness
        for x in np.linspace(previous_pos, tot_thickness, points):
            layer_x.append(x)
            eps_x.append(layer.matl.eps * eps_0)
        previous_pos = layer.thickness
        x_list.append(layer_x)
        eps.append(eps_x)
    
    ### Band Gap and discontinuities
    E_gap = [layer.matl.EG for layer in layer_list]
    Del_E_gap = [E_gap[index + 1] - E_gap[index] for index in range(len(layer_list)-1)]
    Del_EC = [Del_E_gap[index] * 0.67 for index in range(len(Del_E_gap))]
    Del_EV = [Del_E_gap[index] * 0.33 for index in range(len(Del_E_gap))]
    
    ### Doping 
    N = [layer.doping for layer in layer_list]
    N_C = [layer.matl.NC for layer in layer_list]
    N_V = [layer.matl.NV for layer in layer_list]
    
    ### Fermi level to Conduction/Valence Band
    FE = [-k*T*np.log(N[i]/N_V[i]) for i in range(len(N) - 1)]
    EF = [-k*T*np.log(N[2*i + 1]/N_C[2*i + 1]) for i in range(len(N) - 1)]
    
    ### Potential Difference
    V_0 = [V0(E_gap[i], Del_EC[i], FE[i], EF[i]) for i in range(len(layer_list) - 1)]
    
    ### Depletion Widths
    Dep_width_loc = []
    for index in range(len(layer_list) - 1):
        Dep_width_loc.append(Depletion_Width(eps[index][0], eps[index + 1][0], N[index], N[index + 1]
        , V_0[index])*N[index + 1])
        Dep_width_loc.append(Depletion_Width(eps[index][0], eps[index + 1][0], N[index], N[index + 1]
        , V_0[index])*N[index])
    
    ### Energy values (eV) of Valence band
    y_Ev = []
    y_Ec = []
    y_Evac = []
    previous_thickness = 0
    total_thickness = 0
    for i, layer in enumerate(x_list):
        y_Ev_temp = []
        y_Ec_temp = []
        y_Evac_temp = []
        total_thickness += layer_list[i].thickness
        if layer_list[i].n_or_p == 'p':
            for j, x in enumerate(layer):
                if x < -Dep_width_loc[i] + total_thickness and x >= previous_thickness:
                    temp = -FE[i]
                    y_Ev_temp.append(temp)
                    y_Ec_temp.append(temp + E_gap[i])
                    y_Evac_temp.append(temp + E_gap[i] + layer_list[i].matl.chi)
                    
                elif x >= -Dep_width_loc[i] + total_thickness and x < total_thickness:
                    temp = -((N[i]/(2*eps[i][j]))*(x + Dep_width_loc[i] - total_thickness + previous_thickness)**2)*q - FE[i]
                    y_Ev_temp.append(temp)
                    y_Ec_temp.append(temp + E_gap[i])
                    y_Evac_temp.append(temp + E_gap[i] + layer_list[i].matl.chi)
                                                                                
                elif x == total_thickness - previous_thickness:
                    temp = (-(Dep_width_loc[i] + x - total_thickness + previous_thickness)**2 * N[i]/(2*eps[i][j]))*q - FE[i]
                    y_Ev_temp.append(temp)
                    y_Ec_temp.append(temp+ + E_gap[i])
                    y_Evac_temp.append(temp + E_gap[i] + layer_list[i].matl.chi)                                       
        
        elif layer_list[i].n_or_p == 'n': 
            for j, x in enumerate(layer):
                if x == previous_thickness:
                    temp = -Del_EV[i-1] - (((Dep_width_loc[i-1] + x - previous_thickness)**2)*N[i-1]/(2*eps[i-1][j]))*q - FE[i-1]
                    y_Ev_temp.append(temp)
                    y_Ec_temp.append(temp + E_gap[i])
                    y_Evac_temp.append(temp + E_gap[i] + layer_list[i-1].matl.chi - Del_EC[i-1])
                    
                elif x > previous_thickness and x <= previous_thickness + Dep_width_loc[i] :
                    y_Ev_temp.append(-Del_EV[i-1] - (((Dep_width_loc[i-1])**2)*N[i-1]/(2*eps[i-1][j]))*q - FE[i-1]
                    - (2*(x - previous_thickness)*Dep_width_loc[i] - (x - previous_thickness)**2)*q*N[i]/(2*eps[i][j]))
                    
                    y_Ec_temp.append(-Del_EV[i-1] - (((Dep_width_loc[i-1])**2)*N[i-1]/(2*eps[i-1][j]))*q - FE[i-1]
                    - (2*(x - previous_thickness)*Dep_width_loc[i] - (x - previous_thickness)**2)*q*N[i]/(2*eps[i][j]) + E_gap[i])
                    
                    y_Evac_temp.append(-Del_EV[i-1] - (((Dep_width_loc[i-1])**2)*N[i-1]/(2*eps[i-1][j]))*q - FE[i-1]
                    - (2*(x - previous_thickness)*Dep_width_loc[i] - (x - previous_thickness)**2)*q*N[i]/(2*eps[i][j]) + E_gap[i] + 
                    layer_list[i-1].matl.chi - Del_EC[i-1])
                    
                elif x > previous_thickness + Dep_width_loc[i] and x <= total_thickness:
                    temp = -Del_EV[i-1] - V_0[i-1] - FE[i-1]
                    y_Ev_temp.append(temp)
                    y_Ec_temp.append(temp + E_gap[i])
                    y_Evac_temp.append(temp + E_gap[i] + layer_list[i-1].matl.chi - Del_EC[i-1])
                    
        y_Evac.append(y_Evac_temp)            
        y_Ec.append(y_Ec_temp)            
        y_Ev.append(y_Ev_temp) 
        previous_thickness = total_thickness
    print (Del_E_gap)
    return x_list, y_Ev, y_Ec, y_Evac, Dep_width_loc

def plot_bands(calc, layers):
    x_list, y_Ev, y_Ec, y_Evac, Dep_width = calc    
    
    # Convert to nm
    x_points = []
    for i in x_list:
        x_temp = []
        for x in i:
            x_temp.append(x*1E9)
        x_points.append(x_temp)
    
    prev_thickness = 0 
    total_thickness = 0
    fig = plt.figure(1)
    ax1 = fig.add_subplot(111)
    Avg_EG = np.mean([layer.matl.EG for layer in layers])    
    
    ax1.set_ylim([40*np.max(y_Ev), 1.5*np.max(y_Ec)]) 
    ax1.set_ylim([-2*Avg_EG, 1.25*np.max(y_Evac)])
    ax1.set_xlim([-50, x_points[-1][-1] + 50])
    ax1.hlines(0, 0, x_points[-1][-1], "red")
    
    ### Main plot: Energy Band Diagram
    for index in range(len(x_list)):
        ax1.plot(x_points[index], y_Ev[index], 'black')     # Valence Band
        ax1.plot(x_points[index], y_Ec[index], 'black')     # Conduction Band
        ax1.plot(x_points[index], y_Evac[index], 'black')
        ax1.vlines(x_points[index][0], -10, 10, linestyle = '--')   # interface
        ax1.vlines(x_points[index][-1], -10, 10, linestyle = '--')  # interface
        total_thickness += layers[index].thickness*1E9
        if layers[index].n_or_p == 'p':
            temp = -Dep_width[index]*1E9 + total_thickness
            ax1.vlines(temp, -10, 10, linestyle = '--', color = 'blue') # x_p
        elif layers[index].n_or_p == 'n':
            temp = Dep_width[index]*1E9 + prev_thickness
            ax1.vlines(temp, -10, 10, linestyle = '--', color ='red')   #x_N
        prev_thickness = total_thickness
    
    ax1.set_xlabel("Position, nm")
    ax1.set_ylabel("Energy, eV")
    
    title = ""
    for layer in layers:
        title += (layer.n_or_p + " " + layer.matl.name + " " + str(layer.doping) + " ")
    ax1.set_title(title)
    
    ### Electric Field vs Position
    fig2 = plt.figure(2)
    ax2 = fig2.add_subplot(111)
    E_field = Elec_field(x_list, layers, Dep_width)
    ax2.set_ylim(1.5*np.min(E_field), 1E-8)
    ax2.set_xlim(-50, total_thickness + 50)
    for index in range(len(x_list)):
        ax2.plot(x_points[index], E_field[index], 'black')
#        ax2.vlines(x_points[index][0], -10, 10, linestyle = '--')   # interface
#        ax2.vlines(x_points[index][-1], -10, 10, linestyle = '--')  # interface
    
    
    
def Elec_field(x_list, layers, Dep_width):
    prev_thickness = 0
    tot_thickness = 0
    E_field = []
    for ind, x in enumerate(x_list):
        E_temp = []
        tot_thickness += layers[ind].thickness
        if layers[ind].n_or_p == 'p':
            for i in x:
                if i < -Dep_width[ind] + tot_thickness - prev_thickness and i >= 0:
                    E_temp.append(0)
                elif i >= -Dep_width[ind] + tot_thickness and i <= tot_thickness:
                    E_temp.append(-layers[ind].doping*(i + Dep_width[ind] - tot_thickness + prev_thickness)/layers[ind].matl.eps*eps_0)
        elif layers[ind].n_or_p == 'n':
            for i in x:
                if i >= prev_thickness and i <= Dep_width[ind] + prev_thickness:
                    E_temp.append(layers[ind].doping*(i - Dep_width[ind] - prev_thickness)/layers[ind].matl.eps*eps_0)
                elif i > Dep_width[ind] + prev_thickness:
                    E_temp.append(0)
        prev_thickness = tot_thickness        
        E_field.append(E_temp)
    return E_field
        
        
    
    
def Gen_plot(): 
    layer0 = Layer(matl=GaAs, n_or_p= 'p', doping=1e18, thickness=350)
    layer1 = Layer(matl=AlGaAs, n_or_p= 'n', doping=2e17, thickness=200)
    
    layers = [layer0, layer1]

    plot_bands(calc_bands(layers, 100), layers)


class Material:
    """
    Semiconductor material with the following properties...
    
    NC = conduction-band effective density of states in cm^-3
    
    NV = valence-band effective density of states in cm^-3
    
    EG = Band gap in eV
    
    chi = electron affinity in eV (i.e. difference between conduction
          band and vacuum level)
    
    eps = static dielectric constant (epsilon / epsilon0)
    
    ni = intrinsic electron concentration in cm^-3 (defined by p = n = ni
    when the undoped material is in thermal equilibrium)
    
    Evac_minus_Ei is the [positive] energy difference (in eV) between the
        vacuum level and the "intrinsic" fermi level, i.e. the fermi level
        at which p=n.
    
    name = a string describing the material (for plot labels etc.)
    """

    def __init__(self, NC, NV, EG, me, mh, chi, eps, name=''):
        if NC is None:
            self.NC = 2.51E-19 * (me*T/300)**(1.5)
        elif NC is not None:
            self.NC = NC
        if NV is None:
            self.NV = 2.51E-19 * (mh*T/300)**(1.5)
        elif NC is not None:
            self.NV = NV
        self.EG = EG
        self.me = me
        self.mh = mh
        self.chi = chi
        self.eps = eps
        self.name = name
        
        # Sze equation (29), p21...
        self.ni = m.sqrt(self.NC * self.NV * m.exp(-self.EG / kT))
        # Sze equation (27), p20...
        self.Evac_minus_Ei = (self.chi + 0.5 * self.EG
                              + 0.5 * kT * m.log(self.NC / self.NV))

#Sze Appendix G

GaAs = Material(NC=4.3E17,
                NV=8.87E18,
                EG=1.424,
                me = 0.067,
                mh = 0.087,
                chi=4.07,
                eps=13.1,
                name='GaAs')
                
AlGaAs = Material(NC = 6.94E17,
                  NV = 1.13E19,
                  EG = 1.7981,
                  me = 0.585,
                  mh = 0.1024,
                  chi = 3.74,   #x = 0.3
                  eps = 12.2,
                  name= "AlGaAs") 

Si = Material(NC=3.2e19,
              NV=1.8e19,
              EG=1.12,
              me = 0.9163,
              mh = 0.153,
              chi=4.05,
              eps=11.9,
              name='Si')

Ge = Material(NC=1e19,
              NV=5e18,
              EG=0.661,
              me = 1.59,
              mh = 0.043,
              chi=4,
              eps=16.2,
              name='Ge')

AlN = Material(NC=6.3e18,
              NV=4.8e20,
              EG=6.2,
              me = 0.4,
              mh = 0.24,
              chi= 0.6,
              eps=9.14,
              name='AlN')

ZnO = Material(NC=6.3e18,
              NV=4.8e20,
              EG=3.37,
              me = 0.24,
              mh = 0.59,
              chi= 0.6,
              eps=7.926, 
              name='ZnO)

# dieletric constant = 7.926 = sqrt(e_|| * e_norm) Retrieved from: http://www.springer.com/cda/content/document/cda_downloaddocument/9783540888468-c1.pdf?SGWID=0-0-45-692199-p173851917

class Layer:
    """
    Layer of semiconductor with the following properties...
    
    matl = a material (an object with Material class)
    
    n_or_p = a string, either 'n' or 'p', for the doping polarity
    
    doping = density of dopants in cm^-3
    
    thickness = thickness of the layer in nm
    """
    def __init__(self, matl, n_or_p, doping, thickness):
        self.matl = matl
        self.n_or_p = n_or_p
        self.doping = doping
        self.thickness = thickness*1E-9
