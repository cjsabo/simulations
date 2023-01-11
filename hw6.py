import matplotlib.pyplot as pl
import numpy as np
import math as m
import pandas as pd
import scipy as sp

def qwe(h,m,l,q):
    E = (h**2)/(8*m*l**2)
    return E/q

#constants

#fundamental electron mass (kg)
m0 = 9.11e-31

#electron charge (C)
q = 1.602e-19

#speed of light (cm/s)
c = 3e10

#planck's constant (weird ass units)
h = 6.626e-34

#boltzmann's constant (also weird ass units)
k = 1.38e-23

#temperature (K)
T = 300

#Big K for GaAs (cm^-1 eV^-1)
K = 11700

#GaAs bangap (eV)
Eg_GaAs = 1.42

#Thermal voltage
Vt = (k*T)/q


#Parameters for edge emitter ONLY

#cavity length (microns) 
L_ee = 300

#laser ridge width/diameter (m)
W_ee = 3e-6

#spontaneous emission factor 
beta_ee = 1e-5

#mirror reflectivities 
R1_ee = .1
R2_ee = .9

#Area (cm^2)
A_ee = W_ee*1e2*L_ee*1e-4


#Parameters for surface emitter ONLY

#laser ridge width/diameter (m)
W_se = 5e-6

#spontaneous emission factor 
beta_se = 1e-2

#mirror reflectivities 
R1_se = .995
R2_se = .999

#Area (cm^2) (area of a circle)
A_se = (((W_se/2)*1e2)**2)*np.pi

#Shared parameters

#non-radiative recombination time (s)
tau_nr = 100e-9

#Bimolecular recombination constant (cm^2/s)
B = 5e-5

#Auger recombination constant (cm^4/s)
C = 6.25e-18

#Group velocity (cm/s)
v_gr = 9e9

#waveguide losses (cm^-1)
alpha = 3

#confinement factor
gamma = .02

#effective electron mass (kg)
m_e = .067*m0

#effective hole mass (kg)
m_h = .45*m0

#depletion layer width (m)
x_d = .4e-6

#quantum well width (m)
L_x = 100e-10

#effective refractive index
n_eff = c/v_gr

print("\r\n")

#Calculated values

#1 transparent density (from Google Sheet goal seek, cm^-2)
N_tr = 1.13e12


#effective density of states in conduction band of quantum well (cm^-2)
Nc_qw = ((4*np.pi*m_e*k*T)/(h**2))/(1e4)

#effective density of states in valence band of quantum well (cm^-2)
Nv_qw = ((4*np.pi*m_h*k*T)/(h**2))/(1e4)

#First energy level of CONDUCTION band (eV)
E1n = qwe(h, m_e, L_x, q)

#First energy level of VALENCE band (eV)
E1p = qwe(h, m_h, L_x, q)

#Bandgap of entire quantum well (eV), also 5. Photon Energy
Eg_qw = Eg_GaAs + E1n + E1p

#max gain (cm^-1)
gmax = K*np.sqrt(E1n+E1p)
print("Max gain (cm^-1): ", gmax, "\r\n")

#wavelength of light (m)
lamb = (c/(Eg_qw*q/h))*(1e-2)

#length of surface emitter cavity (microns) (1/neff or neff)
L_se = (3/2)*(lamb/1e-6)*(1/n_eff)


#2 peak gain (is he just looking for expression?)


#3 differential gain coefficient at trasnparency density (cm)
curly_l = gmax*(np.exp(-1*N_tr/Nc_qw)/Nc_qw + np.exp(-1*N_tr/Nv_qw)/Nv_qw)


#4 photon lifetime (s)

#edge emitter photon lifetime (s)
tau_ph_ee = 1/(v_gr*(alpha+(1/(L_ee*1e-4))*np.log(1/np.sqrt(R1_ee*R2_ee))))

#surface emitter photon lifetime (s)
tau_ph_se = 1/(v_gr*(alpha+(1/(L_se*1e-4))*np.log(1/np.sqrt(R1_se*R2_se))))


#6 modal gain at threshold (when lasing starts)

#modal gain edge emitter (cm^-1)
mg_ee = alpha + (1/(2*L_ee*1e-4))*np.log(1/(R1_ee*R2_ee))
#print(mg_ee/gamma)

#modal gain surface emitter (cm^-1)
mg_se = alpha + (1/(2*L_se*1e-4))*np.log(1/(R1_se*R2_se))
#print(mg_se/gamma)


#7 threshold carrier density (might do in Google Sheets) (CHANGE THESE, NEED EXACT)

#edge emitter (from Google Sheet, cm^-2)
N0_ee = 6.289e12

#surface emitter (from Google Sheet, cm^-2) (had to increase number of wells from 1 to 3) (double check with Dr. VanZeghbroeck to find out how to solve)
N0_se = 3.157e12

#need to recalculate differential gain coefficient at threshold now
curly_l_ee = gmax*(np.exp(-1*N0_ee/Nc_qw)/Nc_qw + np.exp(-1*N0_ee/Nv_qw)/Nv_qw)

curly_l_se = gmax*(np.exp(-1*N0_se/Nc_qw)/Nc_qw + np.exp(-1*N0_se/Nv_qw)/Nv_qw)

#8 threshold current density

#current density edge emitter (A/cm^2)
Jth_ee = q*(B*(N0_ee**2) + N0_ee/(tau_nr))
Ith_ee = Jth_ee*A_ee

#current density surface emitter (A/cm^2)
Jth_se = q*(B*(N0_se**2) + N0_se/(tau_nr))
Ith_se = Jth_se*A_se


#9 laser power at twice threshold current

#photon density for edge emitter
S0_ee = ((2*Jth_ee - Jth_ee)/q)*(1/(v_gr*curly_l_ee*gamma*(N0_ee - N_tr)))
#Power out surface emitter
P1_ee = Eg_qw*q*S0_ee*W_ee*1e2*v_gr*np.log(1/np.sqrt(R1_ee))

#photon density for surface emitter
S0_se = ((2*Jth_se - Jth_se)/q)*(1/(v_gr*curly_l_se*gamma*(N0_se - N_tr)))
#Power out surface emitter
P1_se = Eg_qw*q*S0_se*W_se*1e2*v_gr*np.log(1/np.sqrt(R1_se))

#10 differential efficiency

#edge emitter
DE_ee = Eg_qw*(np.log(1/np.sqrt(R1_ee)))/(np.log(1/np.sqrt(R1_ee*R2_ee))+alpha*L_ee*1e-4)
QE_ee = DE_ee/Eg_qw

#surface emitter
DE_se = Eg_qw*np.log(1/np.sqrt(R1_se))/(np.log(1/np.sqrt(R1_se*R2_se))+alpha*L_se*1e-4)
QE_se = DE_se/Eg_qw


#11 small signal equivalent circuit values of C, L and corresponding resonant raidal frequency w_p

#edge emitter
#Area of edge emitter (cm^2)
A_ee = W_ee*1e2*L_ee*1e-4
#m parameter
m_ee = (N0_ee*np.exp(N0_ee/Nc_qw))/(Nc_qw*(np.exp(N0_ee/Nc_qw) - 1)) + (N0_ee*np.exp(N0_ee/Nv_qw))/(Nc_qw*(np.exp(N0_ee/Nv_qw) - 1))
#Capacitance
C_ee = (q*N0_ee*A_ee)/(m_ee*Vt)
#Inductance
Ind_ee = (1/C_ee)*(tau_ph_ee/(gamma*curly_l_ee*S0_ee*v_gr))
#resonant frequency
w_p_ee = 1/np.sqrt(C_ee*Ind_ee)

#surface emitter
#Area of Surface emitter (cm^2) (area of a circle)
A_se = (((W_se/2)*1e2)**2)*np.pi
#m parameter
m_se = (N0_se*np.exp(N0_se/Nc_qw))/(Nc_qw*(np.exp(N0_se/Nc_qw) - 1)) + (N0_se*np.exp(N0_se/Nv_qw))/(Nc_qw*(np.exp(N0_se/Nv_qw) - 1))
#Capacitance
C_se = (q*N0_se*A_se)/(m_se*Vt)
#Inductance
Ind_se = (1/C_se)*(tau_ph_se/(gamma*curly_l_se*S0_se*v_gr))
#resonant frequency
w_p_se = 1/np.sqrt(C_se*Ind_se)


#print out answer to all the questions
print("#1 Transparency Density (cm^-2): ", N_tr)
print("\r")
print("#2 Peak Gain (expression in notes) \n")
print("#3 Differential gain coefficient (cm): ", curly_l)
print("\r")
print("#4 Photon Lifetime (ps)")
print("Edge Emitter: ", tau_ph_ee*1e12)
print("Surface Emitter: ", tau_ph_se*1e12)
print("\r")
print("#5 Photon Energy (eV): ", Eg_qw)
print("\r")
print("#6 Modal Gain at threshold (cm^-1)")
print("Edge Emitter: ", mg_ee)
print("Surface Emitter: ", mg_se)
print("\r")
print("#7 Carrier Density at threshold (cm^-2)")
print("Edge Emitter: ", N0_ee)
print("Surface Emitter: ", N0_se)
print("\r")
print("#8 Threshold Current Density (A/cm^2)")
print("Edge Emitter: ", Jth_ee)
print("Surface Emitter: ", Jth_se)
print("\r")
print("#9 Laser power at twice threshold current (mW)")
print("Edge Emitter: ", P1_ee*1e3)
print("Surface Emitter: ", P1_se*1e3)
print("\r")
print("#10 Differential efficiency and Quantum efficiency: ")
print("Edge Emitter DE: ", DE_ee)
print("Surface Emitter DE: ", DE_se)
print("\r")
print("Edge Emitter QE: ", QE_ee)
print("Surface Emitter QE: ", QE_se)
print("\r")
print("#11 Small signal C,L and w_p (pF, nH, Grad/s and GHz)")
print("Edge Emitter Capacitance: ", C_ee*1e12)
print("Edge Emitter Inductance: ", Ind_ee*1e9)
print("Edge Emitter frequency: ", w_p_ee*1e-9, " ", w_p_ee*1e-9/(2*np.pi))
print("\r")
print("Surface Emitter Capacitance: ", C_se*1e12)
print("Surface Emitter Inductance: ", Ind_se*1e9)
print("Surface Emitter frequency: ", w_p_se*1e-9, " ", w_p_se*1e-9/(2*np.pi))


print("\r")