import matplotlib.pyplot as pl
import numpy as np
import math as m
import pandas as pd
import scipy as sp

def sign(phi):
    if(phi == 0):
        return 1
    else:
        return phi/abs(phi)
def E(Vt,LD,phi):
    field = sign(phi)*(Vt/LD)*np.sqrt(2*(np.exp(phi/Vt) - (phi/Vt) -1))
    return field

def seek(phi_s, Va, phi_i, phi_AlGaAs, phi_sp, phi_GaAs, low, high):
    if((phi_GaAs + phi_sp + phi_AlGaAs) > (phi_i - Va)):
        return (phi_s + low)/2, low, phi_s
    elif((phi_GaAs + phi_sp + phi_AlGaAs) < (phi_i - Va)):
        return (phi_s+high)/2, phi_s, high
    else:
        return phi_s, low, high
    

#Linespace of plots (nm)
x = np.linspace(0,430e-9,1000)

xa = np.linspace(0, 100e-9, 1000)
xb = np.linspace(100e-9, 110e-9, 1000)
xc = np.linspace(110e-9, 430e-9, 10000)

#Constants

#Applied voltage
Va = 0;

#electronic charge
q = 1.602e-19

#Boltzmann constant
k = 1.38e-23

#temperature (K)
T = 300

#permittivity of vacuum (F/cm)
eps_0 = 8.85e-14

#Free electron mass (kg)
m_0 = 9.11e-31

#thermal voltage
Vt = 2.58e-2

#Al(.4)Ga(.6)As parameters

#layer thickness
d_AlGaAs = 100e-9

#doping density (cm^-3)
Nd_AlGaAs = 1e17

#bandgap
Eg_AlGaAs = 1.924

#DOS effective electron mass
me_AlGaAs = .794*m_0

#electron affinity
chi_AlGaAs = 3.77

#effective density of states (cm^-3)
Nc_AlGaAs = 7.46e17

#relative dielectric constant
eps_r_AlGaAs = 11.7

#intrinsic carrier concentration
ni_AlGaAs = 0

#GaAs parameters

#doping density (cm^-3)
Nd_GaAs = 1e14

#bandgap
Eg_GaAs = 1.424

#DOS electron effective mass
me_GaAs = .067*m_0

#electron affinity
chi_GaAs = 4.07

#effective density of states
Nc_GaAs = 4.7e17

#relative dielectric constant
eps_r_GaAs = 12.9

#Spacer layer thickness
d_sp = 10e-9

#Chrome workfunction
phi_m = 4.6

print("\r\n")

#Calculated values

#change in valence band and conduction band energy

dEc = abs(.79*.4)
print("delta Ec (eV): ", dEc)

dEv = abs(-.46*.4)
print("delta Ev (eV): ", dEv)

#EC-EF

#AlGaAs
EC_EF_AlGaAs = -1*Vt*np.log(Nd_AlGaAs/Nc_AlGaAs)
print("EC-EF (AlGaAs, meV): ", EC_EF_AlGaAs*1000)

#GaAs
EC_EF_GaAs = -1*Vt*np.log(Nd_GaAs/Nc_GaAs)
print("EC-EF (GaAs, meV): ", EC_EF_GaAs*1000)

#Fermi Energies and built in potential

EF_AlGaAs = (phi_m - chi_AlGaAs) - Vt*np.log(Nc_AlGaAs/Nd_AlGaAs)
print("Fermi energy of AlGaAs layer (eV): ", EF_AlGaAs)

EF_GaAs = (phi_m - chi_GaAs) - Vt*np.log(Nc_GaAs/Nd_GaAs)
print("Fermi energy of GaAs layer (eV): ", EF_GaAs)

phi_i = (phi_m - chi_GaAs - EF_GaAs)
print("Built in potential (V): ", phi_i)

Ei_AlGaAs = ((phi_m - chi_AlGaAs) + (phi_m - chi_AlGaAs - Eg_AlGaAs))/2 + Vt*np.log(Nd_AlGaAs/Nc_AlGaAs)/2
print("Intrinsic Energy AlGaAs (eV): ", Ei_AlGaAs)

Ei_GaAs = ((phi_m-chi_GaAs) + (phi_m-chi_GaAs-Eg_GaAs))/2 + Vt*np.log(Nd_GaAs/Nc_GaAs)/2
print("Intrinsic Energy GaAs (eV): ", Ei_GaAs)

#Dielectric Constant

#AlGaAs
eps_AlGaAs = eps_r_AlGaAs*eps_0
print("Dielectric Constant AlGaAs (F/cm): ", eps_AlGaAs)

#GaAs
eps_GaAs = eps_r_GaAs*eps_0
print("Dielectric Constant GaAs (F/cm): ", eps_GaAs)

#Debye length GaAs
LD = np.sqrt((eps_GaAs*Vt)/(q*Nd_GaAs)) #in cm
print("Debye length (nm): ", LD/(1e-7)) #print out in nm

print("\r\n")

#Flatband diagram plot calculations
cond = []
val = []
EF_a = []
EF_c = []
Ei = []
Va_plot = [Va]*len(x)
for i in range(len(x)):
    if(x[i] >= 0 and x[i] < 110e-9):
        cond.append(phi_m - chi_AlGaAs)
        val.append(phi_m - chi_AlGaAs - Eg_AlGaAs)
        Ei.append(Ei_AlGaAs)
    else:
        cond.append(phi_m - chi_GaAs)
        val.append(phi_m - chi_GaAs - Eg_GaAs)
        Ei.append(Ei_GaAs)
        
for i in range(len(xa)):
    EF_a.append(EF_AlGaAs)

for i in range(len(xc)):
    EF_c.append(EF_GaAs)
        
#Plotting the potential and field
phi_start = 0
tolerance = .000001
low = Va - 3
high = Va + 3
phi_AlGaAs = 0
phi_sp = 0
phi_GaAs = 0

#GaAs layer
pot_GaAs = [0]*len(xc)
efield_GaAs = [0]*len(xc)

while (True):
    for i in range(len(xc)):
        if(i == 0):
            pot_GaAs[i] = phi_start
        else:
            efield = E(Vt, LD, pot_GaAs[i-1])
            pot_GaAs[i] = (-1*efield*100*(xc[i] - xc[i-1]) + pot_GaAs[i-1])

    for i in range(len(xc)):
        efield_GaAs[i] = (E(Vt, LD, pot_GaAs[i]))

    #Spacer layer
    efield_sp = [0]*len(xb)
    pot_sp = [0]*len(xb)
    for i in range(len(xb)-1,-1,-1):
        if(i == len(xb)-1):
            pot_sp[i] = pot_GaAs[0]
            efield_sp[i] = E(Vt, LD, pot_GaAs[0])*(eps_GaAs/eps_AlGaAs)
        else:
            pot_sp[i] = pot_sp[i+1] - efield_sp[i+1]*100*(xb[i] - xb[i+1])
            efield_sp[i] = efield_sp[i+1]

    #AlGaAs layer
    efield_AlGaAs = [0]*len(xa)
    pot_AlGaAs = [0]*len(xa)

    #Charge density in AlGaAs
    rho_AlGaAs = q*Nd_AlGaAs

    #total charge on doped layer
    Q_AlGaAs = rho_AlGaAs*100*d_AlGaAs

    for i in range(len(xa)-1, -1, -1):
        if(i == len(xa)-1):
            efield_AlGaAs[i] = efield_sp[0]
        else: #check this again
            efield_AlGaAs[i] = (rho_AlGaAs/eps_AlGaAs)*(100)*(xa[i] - xa[i+1]) + efield_AlGaAs[i+1]

    for i in range(len(xa)-1, -1, -1):
        if(i == len(xa)-1):
            pot_AlGaAs[i] = pot_sp[0]
        else:
            pot_AlGaAs[i] = pot_AlGaAs[i+1] - efield_AlGaAs[i+1]*100*(xa[i] - xa[i+1])

    #plot for energy band diagram
    cond_bd_AlGaAs = [0]*len(xa)
    val_bd_AlGaAs =[0]*len(xa)
    cond_bd_sp = [0]*len(xb)
    val_bd_sp =[0]*len(xb)
    cond_bd_GaAs = [0]*len(xc)
    val_bd_GaAs =[0]*len(xc)
    Ei_bd_AlGaAs = [0]*len(xa)
    Ei_bd_sp = [0]*len(xb)
    Ei_bd_GaAs = [0]*len(xc)

    for i in range(len(xa)):
        cond_bd_AlGaAs[i] = (cond[10] - pot_AlGaAs[i])
        val_bd_AlGaAs[i] = (val[10] - pot_AlGaAs[i])
        Ei_bd_AlGaAs[i] = Ei_AlGaAs - pot_AlGaAs[i]

    for i in range(len(xb)):
        cond_bd_sp[i] = (cond[10] - pot_sp[i])
        val_bd_sp[i] = (val[10] - pot_sp[i])
        Ei_bd_sp[i] = Ei_AlGaAs - pot_sp[i]

    for i in range(len(xc)):
        cond_bd_GaAs[i] = (cond[900] - pot_GaAs[i])
        val_bd_GaAs[i] = (val[900] - pot_GaAs[i])
        Ei_bd_GaAs[i] = Ei_GaAs - pot_GaAs[i]
    
    phi_GaAs = phi_start
    phi_sp = pot_sp[0] - pot_sp[len(xb)-1]
    phi_AlGaAs = pot_AlGaAs[0] - pot_AlGaAs[len(xa)-1]
    
    if(abs(-1*phi_i-Va-(phi_GaAs + phi_sp + phi_AlGaAs)) <= tolerance):
        break
    else:
        phi_start, low, high = seek(phi_start, Va, -1*phi_i, phi_AlGaAs, phi_sp, phi_GaAs, low, high)

#accumulated charge in GaAs (this needs work)
totcharge = 0
for i in range(len(xc)):
    if(xc[i] <= LD):
        totcharge = totcharge - (eps_GaAs*efield_GaAs[i])
    
print("Potential across AlGaAs layer: ", phi_AlGaAs)
print("Potential across spacer layer: ", phi_sp)
print("Potential across GaAs layer: ", phi_GaAs)
print("Field just to left of interface: ", efield_sp[len(xb)-1])
print("Field just to right of interface: ", efield_GaAs[0])
print("Field at metal-AlGaAs interface: ", efield_AlGaAs[0])
print("Charge in GaAs layer (C/cm^2): ", totcharge)
print("\r\n")

#Flatband diagram plot
pl.figure(1)
pl.plot(x/1e-9, cond, "b")
pl.plot(x/1e-9, val, "m")
pl.plot(xa/1e-9, EF_a, "c--")
pl.plot(xc/1e-9, EF_c, "r--")
pl.plot(x/1e-9, Ei, "g--")
pl.title("Flatband diagram of HEMT structure")
pl.legend(["Conduction band edge", "Valence band edge", "Fermi energy (AlGaAs)",  "Fermi energy (GaAs)", "Intrinsic Energy"])
pl.ylabel("Energy (eV)")
pl.xlabel("Position (nm)")
pl.grid()

#Plotting the potential
pl.figure(2)
pl.plot(xc/1e-9, pot_GaAs, "r")
pl.plot(xb/1e-9, pot_sp, "g")
pl.plot(xa/1e-9, pot_AlGaAs, "b")
pl.title("Potential vs. Position")
pl.ylabel("Potential (V)")
pl.xlabel("Position (nm)")
pl.legend(["GaAs layer", "Spacer Layer", "AlGaAs layer"])
pl.grid()

#Plotting the field
pl.figure(3)
pl.plot(xa/1e-9, efield_AlGaAs, "b")
pl.plot(xb/1e-9, efield_sp, "g")
pl.plot(xc/1e-9, efield_GaAs, "r")
pl.title("Field vs. Position")
pl.ylabel("Electric field (V/cm)")
pl.xlabel("Position (nm)")
pl.legend(["AlGaAs layer", "Spacer Layer", "GaAs layer"])
pl.grid()

#Plotting band diagram
pl.figure(4)
pl.plot(xa/1e-9, cond_bd_AlGaAs, "b")
pl.plot(xb/1e-9, cond_bd_sp, "b")
pl.plot(xc/1e-9, cond_bd_GaAs, "b")
pl.plot(xa/1e-9, val_bd_AlGaAs, "r")
pl.plot(xb/1e-9, val_bd_sp, "r")
pl.plot(xc/1e-9, val_bd_GaAs, "r")
pl.plot(x/1e-9, Va_plot, "m--")
pl.plot(xa/1e-9, Ei_bd_AlGaAs, "c--")
pl.plot(xb/1e-9, Ei_bd_sp, "c--")
pl.plot(xc/1e-9, Ei_bd_GaAs, "c--")
pl.title("Band diagram of HEMT structure")
pl.ylabel("Energy (eV)")
pl.xlabel("Position (nm)")
pl.grid()

pl.show()