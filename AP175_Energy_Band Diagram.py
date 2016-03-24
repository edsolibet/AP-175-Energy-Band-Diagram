# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 10:12:08 2016

@author: Carlo
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

X = 0.3
T = 300 #K
eps_0 = 8.854*10**-18 #cm^-3 kg^-1 s^4 A^2
m_0 = 9.11*10**-31
q = 1.6E-19


def Eg(x):
    return 1.424 + 1.247*x #eV

def eps(x):
    return (13.1 - 3.0*x)*eps_0

def mh(x):
    return (0.50 + 0.29*x)*m_0
    
def me(x):
    return (0.0665 + 0.083*x)*m_0

''''------------------Doping Conentrations---------------------------------'''
N_a = 1*10**18 #cm^-3  p side
N_D = 2*10**17 #cm^-3  N side

'''------------------ Band Discontinuities--------------------------------- '''
E_gp = Eg(0)
E_gN = Eg(X)
Del_Eg = E_gN - E_gp
Del_Ec = 0.67 * Del_Eg #GaAs - AlGaAs
Del_Ev = 0.33 * Del_Eg #GaAs - AlGaAs

'''------------------ p-GaAs Region-----------------------------------------'''
N_cp = 2.51*10**19 * ((me(0)/m_0)*(T/300))**1.5 #cm^-3
N_vp = 2.51*10**19 * ((mh(0)/m_0)*(T/300))**1.5 #cm^-3
p = N_a
FE_p = -8.617*10**-5 * T * np.log(N_a/N_vp)

'''-------------------N AlGaAs region---------------------------------------'''

N_cN = 2.51*10**19 * ((me(X)/m_0)*(T/300))**1.5
N_vN = 2.51*10**19 * ((mh(X)/m_0)*(T/300))**1.5
N = N_D
EF_N = -8.617*10**-5 * T * np.log(N_D/N_cN)

'''-------------------Contact Potential-------------------------------------'''
V0 = (E_gp + Del_Ec - FE_p - EF_N)
    
'''-------------------Depletion Width---------------------------------------'''
x_p = (2*eps(0)*V0/(1.6*10**-19 * N_a*N_D*(N_D + (eps(0)/eps(X))*N_a)))**0.5 * N_D
x_N = x_p*N_a/N_D
x_w = x_p + x_N

'''-----------------------x a-xis division ---------------------------------'''
x = np.linspace(-2.0E-7, 2.0E-7, 151)
z_index = np.where(x == 0)[0][0]
xp = list(x[:z_index+1])
xN = list(x[z_index:])

'''------------------------Potential Energy---------------------------------'''
def Potential_p(xp, x_p, x_N):
    pass

def Potential_N(xN, x_p, x_N):
    pass



'''-------------------Valence Band Energy-----------------------------------'''
        

        
def Ev_p(xp, x_p, x_N):
    y_Ev_p  = []
    for index, i in enumerate(xp):
        if i <= -x_p:
            y_Ev_p.append(-FE_p)
        elif i > -x_p and i < 0:
            y_Ev_p.append((-(N_a/(2*eps(0)))*(i + x_p)**2)*q - FE_p)
        elif i == 0:
            y_Ev_p.append((-((x_p**2)*N_a/(2*eps(0))))*q - FE_p)
    return y_Ev_p
    
def Ev_N(xN, x_p, x_N):
    y_Ev_N  = []
    for index, i in enumerate(xN):
        if i == 0:
            y_Ev_N.append(-Del_Ev - ((x_p**2)*N_a/(2*eps(0)))*q - FE_p )
        elif i > 0 and i <= x_N:
            y_Ev_N.append(-Del_Ev + (-((x_p**2)*N_a/(2*eps(0))) - (N_D/(2*eps(X)))*(2*i*x_N - i**2))*q - FE_p)
        elif i > x_N:
            y_Ev_N.append(-Del_Ev - V0 - FE_p)
    return y_Ev_N

'''--------------------Conduction Band Energy-------------------------------'''

def Ec_p(y_Ev_p):
    y_Ec_p = []
    for i in y_Ev_p:
        y_Ec_p.append(i + E_gp )
    return y_Ec_p

def Ec_N(y_Ev_N):
    y_Ec_N = []
    for i in y_Ev_N:
        y_Ec_N.append(i + E_gN )
    return y_Ec_N

'''--------------------Vacuum Level-----------------------------------------'''
def VacLev(x, x_p, x_N):
    if x < 0    :
        return ((N_a/(2*eps(0)))*(x + x_p)**2)*q
    else:
        return (((x_p**2)*N_a/(2*eps(0))) + (N_D/(2*eps(X)))*(2*x*x_N - x**2))*q

'''--------------------Plot-------------------------------------------------'''

y_Ev_p = Ev_p(xp, x_p, x_N) 
y_Ev_N = Ev_N(xN, x_p, x_N) 
y_Ec_p = Ec_p(Ev_p(xp, x_p, x_N))
y_Ec_N = Ec_N(Ev_N(xN, x_p, x_N))
y_lim = [y_Ec_p[0]*(-1.4), y_Ec_p[0]*(1.3)]
x_lim = [x[0], x[-1]]

vacuum_level = []

fig = plt.figure(1)
ax1 = fig.add_subplot(111)

ax1.hlines(0, x[0], x[-1]) # E = 0

ax1.annotate(r"$x_p$", xy = (-x_p, y_lim[0]), xytext = (-x_p, y_lim[0]*0.70))
ax1.annotate(r"$x_N$", xy = (x_N, y_lim[0]), xytext = (x_N, y_lim[0]*0.70))
ax1.vlines(-x_p, y_lim[0], y_lim[0]*0.75)
ax1.vlines(0, y_lim[0], y_lim[0]*0.75)
ax1.vlines(x_N, y_lim[0], y_lim[0]*0.75)

ax1.text(x[0], 0.1, r"$F_p$")
ax1.text(x[-1]*0.9, -0.15, r"$F_N$")
#ax1.text(x[0], y_Ev_p[0]*3, r"$E_{vp} = %0.1f$ meV" %(-FE_p*1000), fontsize = 10)
ax1.annotate(r"$E_{vp} = %0.1f$ meV" %(-FE_p*1000), xy = (x[0], y_Ev_p[0]), xytext = (x[0], y_Ev_p[0]*3), 
             fontsize =10, arrowprops=dict(facecolor='black', shrink=0.05))
#ax1.text(x[0], y_Ec_p[0]*1.1, r"$E_{cp} = %0.2f$ meV" %(y_Ec_p[0]*1000), fontsize = 10)
ax1.annotate(r"$E_{cp} = %0.2f$ meV" %(y_Ec_p[0]*1000), xy = (x[0], y_Ec_p[0]), xytext = (x[0], y_Ec_p[0]*1.1), 
             fontsize = 10, arrowprops=dict(facecolor='black', shrink=0.05))
#ax1.text(x[-1]*0.8, y_Ec_N[-1]*4, r"$E_{cN} = %0.1f$ meV" %(EF_N*1000), fontsize = 10)
ax1.annotate(r"$E_{cN} = %0.1f$ meV" %(EF_N*1000), xy = (x[-1], y_Ec_N[-1]), xytext = (x[-1]*0.8, y_Ec_N[-1]*4), 
             fontsize = 10, arrowprops=dict(facecolor='black', shrink=0.05))
#ax1.text(x[-1]*0.8, y_Ev_N[-1]*0.95, r"$E_{vN} = %0.2f$ meV" %(-y_Ec_N[0]*1000), fontsize = 10)
ax1.annotate( r"$E_{vN} = %0.2f$ meV" %(-y_Ec_N[0]*1000), xy = (x[-1], y_Ev_N[-1]), xytext = (x[-1]*0.8, y_Ev_N[-1]*0.95), 
             fontsize = 10, arrowprops=dict(facecolor='black', shrink=0.05))
ax1.plot(xp, y_Ev_p)    # Valence Band Energy in p region
ax1.hlines(y_Ev_p[-1], xp[-1], xN[-1]/5) # Valence Band discontinuity
ax1.text(xN[-1]/5, (y_Ev_p[-1] + y_Ev_N[0])/2, r"$\Delta E_v = %0.1f $ meV" %(Del_Ev*1000), fontsize = 10)
ax1.hlines(y_Ev_N[0], xp[-1], xN[-1]/5)  # Valence Band discontinuity
ax1.vlines(0, y_Ev_p[-1], y_Ev_N[0]) #Valence Band Discontinuity
ax1.plot(xN, y_Ev_N)    # Valence Band Energy in N region
ax1.plot(xp, y_Ec_p)    # Conduction Band Energy in p region
ax1.hlines(y_Ec_p[-1], xp[-1], xN[-1]/5) # Conduction Band Discontinuity
ax1.vlines(0, y_Ec_p[-1], y_Ec_N[0]) #Conduction Band Discontinuity
ax1.text(xN[-1]/5, (y_Ec_p[-1] + y_Ec_N[0])/2, r"$\Delta E_c = %0.1f $ meV" %(Del_Ec*1000), fontsize = 10)
ax1.hlines(y_Ec_N[0], xp[-1], xN[-1]/5) # Conduction Band Discontinuity
ax1.plot(xN, y_Ec_N)    # Conduction Band Energy in N region


ax1.set_ylim(y_lim[0], y_lim[1])
ax1.set_xlim(x_lim[0], x_lim[1])
ax1.set_xticks(np.linspace(1.1*x_lim[0], 1.1*x_lim[-1], 5))
ax1.set_yticks(np.linspace(y_lim[0], y_lim[1], 9))



