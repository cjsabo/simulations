import matplotlib.pyplot as pl
import numpy as np
import math as m
import pandas as pd

#Fundamental constants

#electron charge (C)
q = 1.602e-19

#boltzmann constant
k = 1.38e-23

#temperature
T = 300

#thermal voltage
Vt = (k*T)/q

#permittivity of vacuum (F/cm)
eps0 = 8.85e-14

#Constants for problem 1

#flatband voltage (V)
V_FB = .5

#oxide layer width (nm)
x_ox = 100

#resitivity of silicon (ohm*cm)
rho = 1

#relative dielectric constant of silicon
eps_Si = 11.7

#relative dielectric constant of silicon dioxide
eps_SiO2 = 3.9

#calculated values for problem 1

#capacitance per unit area across the oxide
C_ox = (eps0*eps_SiO2)/(100e-7)

#doping density (cm^-3, from graph of resistivity vs. doping density)
N_a = 2e16

#intrinsic carrier density (cm^-3)
ni = 1e10

#phi_f (V)
phi_F = Vt*np.log(N_a/ni)

#threshold voltage
V_T = V_FB + 2*phi_F + np.sqrt(4*eps0*eps_Si*q*N_a*phi_F)/C_ox

#plotting

#voltage applied at the gate
V_G_acc = np.linspace(-6, V_FB, 1000)
C_acc = [C_ox]*len(V_G_acc)

phi_s = np.linspace(0, 2*phi_F, 1000)

V_G_dep = [0]*len(phi_s)
C_dep = [0]*len(phi_s)

for i in range(0, len(phi_s)):
    V_G_dep[i] = V_FB + phi_s[i] + np.sqrt(2*eps0*eps_Si*q*N_a*phi_s[i])/C_ox
    x_d = np.sqrt((2*eps0*eps_Si*phi_s[i])/(q*N_a))
    C_dep[i] = 1/(1/C_ox + x_d/(eps0*eps_Si))
    
V_G_inv = np.linspace(V_G_dep[len(V_G_dep)-1], 6, 1000)
C_inv = [C_ox]*len(V_G_inv)

V_G = np.concatenate([V_G_acc, V_G_dep, V_G_inv])
C_approx = np.concatenate([C_acc, C_dep, C_inv])

#Low frequency model
phi_s_LF = np.linspace(-0.15, 2.35*phi_F, 1000)
E_sq = [0]*len(phi_s_LF)
V_G_LF = [0]*len(phi_s_LF)
C_sLF = [0]*len(phi_s_LF)
C_LF = [0]*len(phi_s_LF)

for i in range(0, len(phi_s_LF)):
    E_sq[i] = 2*(phi_s_LF[i]/abs(phi_s_LF[i]))*np.sqrt(((q*ni*Vt)/(eps0*eps_Si))*(m.cosh((phi_s_LF[i] - phi_F)/Vt) + (phi_s_LF[i]/Vt)*m.sinh(phi_F/Vt) - m.cosh(phi_F/Vt)))
    V_G_LF[i] = V_FB + phi_s_LF[i] + (x_ox*1e-7)*E_sq[i]*(eps_Si/eps_SiO2)
    C_sLF[i] = 2*((q*ni)/E_sq[i])*(m.sinh((phi_s_LF[i] - phi_F)/Vt) + m.sinh(phi_F/Vt))
    C_LF[i] = 1/(1/C_ox + 1/C_sLF[i])
    
pl.plot(V_G, C_approx, 'b')
pl.plot(V_G_LF, C_LF, 'r')
pl.xlabel("Gate Voltage (V)")
pl.ylabel("Capacitance (F)")
pl.legend(["Approx. Capacitance", "Low Freq. Model"])
pl.title("Capacitance vs. Gate Voltage")
pl.rcParams["figure.dpi"] = 300
pl.grid()
pl.show()