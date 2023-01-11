import matplotlib.pyplot as pl
import numpy as np
import math as m
import pandas as pd

#Constants
#Charge of an electron (C)
q = 1.602e-19

#Dielectric constant of air (C^2/(N*m^2))
e_0 = 8.85e-12

#Epsilon infinity and zero
e_inf = 6.52
e_zero = 16

#Fundamental mass of electron (kg)
m_o = 9.11e-31

#Transverse effective mass
m_e_t = .42*m_o

#longitudinal effective mass
m_e_l = .29*m_o

#Effective mass in material
m_e = 3/((1/m_e_l) + (2/m_e_t))

#Electron mobility of substrate layer
mew_sub = 40.0/10000

#Electron mobility of epi layer
mew_epi = 1000.0/10000

#density of electrons for epi layer (m^-3)
n_e_epi = 4e20

#density of electrons for substrate layer (m^-3)
n_e_sub = 6.5e24

#Plasma frequency epi layer
w_p_epi = np.sqrt(((q**2)*n_e_epi)/(e_inf*e_0*m_e))

#Plasma frequency substrate layer 
w_p_sub = np.sqrt(((q**2)*n_e_sub)/(e_inf*e_0*m_e))

#Transverese phonon frequency (rad/s given cm^-1)
w_t = (797)*((3e10)*2*np.pi)

#Longitudinal phonon frequency (rad/s given cm^-1)
w_l = (970)*((3e10)*2*np.pi)

#Thickness of epi layer
d = 7e-6*0

#Plasma Damping term for substrate will go here
gamma_sub = (q/(m_e * mew_sub))

#Plasma damping term for epi layer will go here
gamma_epi = q/(m_e * mew_epi)

#Phonon damping term
bgamma = (2)*(3e10)*(2*np.pi)

#Functions
#Range of frequencies
w = np.linspace(.1*(3e10*2*np.pi), 2000*(3e10*2*np.pi), 1000)

#Imported Data
impdata = pd.read_csv("/Users/cobisabo/Documents/FTIRdata.csv",header = None, names = ['Wavenumber', 'Sub', 'EpiSub'], skiprows = [0])

#Separated columns out
wn_imp = impdata['Wavenumber']
sub_imp = impdata['Sub']
episub_imp = impdata['EpiSub']

#dielectric constant of epi layer
e_epi = e_inf*(1 + (w_l**2 - w_t**2)/(w_t**2 - w**2 - w_t*bgamma*1j) - (w_p_epi**2)/(w**2 + w*gamma_epi*1j))

#dielectric constant of substrate layer
e_sub = e_inf*(1 + (w_l**2 - w_t**2)/(w_t**2 - w**2 - w_t*bgamma*1j) - (w_p_sub**2)/(w**2 + w*gamma_sub*1j))

#Refractive index epi layer
n_epi = np.sqrt(e_epi)

#Refractive index substrate layer
n_sub = np.sqrt(e_sub)

#Refractive index air
n_air = 1

#Angle in air
theta0 = 10*(np.pi/180)

#Incident angle of light epi layer(radians)
theta1 = np.arcsin((n_air*np.sin(theta0))/n_epi)

#Angle in substrate layer
theta2 = np.arcsin((n_air*np.sin(theta0))/n_sub)

#Fabry Perot equations

#Incident Wavelength
iwl = (2*np.pi*3e8)/w

#Just some random constant for Fabry Perot
delta = (4.0*np.pi*n_epi*d*np.cos(theta1))/iwl

#Reflection coefficients TE mode
r_01_TE = (n_air*np.cos(theta0) - n_epi*np.cos(theta1))/(n_air*np.cos(theta0) + n_epi*np.cos(theta1))

r_12_TE = (n_epi*np.cos(theta1) - n_sub*np.cos(theta2))/(n_epi*np.cos(theta1) + n_sub*np.cos(theta2))

#Reflection coefficients TM mode
r_01_TM = (n_epi*np.cos(theta0) - n_air*np.cos(theta1))/(n_epi*np.cos(theta0) + n_air*np.cos(theta1))

r_12_TM = (n_sub*np.cos(theta1) - n_epi*np.cos(theta2))/(n_sub*np.cos(theta1) + n_epi*np.cos(theta2))

#Reflected fraction of incident field for TE
A_RI_TE = (r_01_TE + r_12_TE*np.exp(delta*1j))/(1 + r_01_TE*r_12_TE*np.exp(delta*1j))

#Reflected fraction of incident field for TM
A_RI_TM = (r_01_TM + r_12_TM*np.exp(delta*1j))/(1 + r_01_TM*r_12_TM*np.exp(delta*1j))

#Reflectance for TE mode
R_TE = abs(A_RI_TE)**2

#Reflectance for TE mode
R_TM = abs(A_RI_TM)**2

#Plots
#Calculated Data plots
pl.subplot(3,2,1)
#pl.plot(w/(2*np.pi*3e10), e_sub, "r-")
#pl.plot(w/(2*np.pi*3e10), e_epi, "g-")
pl.title("Dielectric Constant vs. Wavenumber")
pl.ylabel("Dielectric Constant")
pl.xlim([.001, 2200])
pl.ylim([-5,15])
pl.grid()

#TE mode plot
pl.subplot(3,2,3)
pl.plot(w/(2*np.pi*3e10), R_TE, "r-")
#pl.title("Reflectivity vs. Wavenumber")
pl.ylabel("Reflectivity (TE mode)")
pl.xlim([600, 2200])
pl.ylim([0,1])
pl.yticks(np.arange(0,1.1,.1))
pl.grid()

#TM mode plot
pl.subplot(3,2,4)
pl.plot(w/(2*np.pi*3e10), R_TM, "b-")
pl.title("Reflectivity vs. Wavenumber")
pl.xlabel(r'Wavenumber $cm^{-1}$')
pl.ylabel("Reflectivity (TM mode)")
pl.xlim([600, 2200])
pl.ylim([0,1])
pl.yticks(np.arange(0,1.1,.1))
pl.grid()

#Imported Data plots

#Imported data WITH Epi layer
pl.subplot(3,2,5)
pl.plot(wn_imp, episub_imp, "g-")
#pl.plot(w/(2*np.pi*3e10), R_TM, "r-")
pl.xlabel(r'Wavenumber $cm^{-1}$')
pl.ylabel("Reflectivity (EpiSub imp.)")
pl.xlim([.001, 2000])
pl.ylim([0,1])
pl.yticks(np.arange(0,1.1,.1))
pl.grid()

#Imorted data WITHOUT Epi Layer
pl.subplot(3,2,6)
pl.plot(wn_imp, sub_imp, "c-")
pl.plot(w/(2*np.pi*3e10), R_TM, "r-")
pl.xlabel(r'Wavenumber $cm^{-1}$')
pl.ylabel("Reflectivity (Sub imp.)")
pl.xlim([.001, 2000])
pl.ylim([0,1])
pl.yticks(np.arange(0,1.1,.1))
pl.grid()

pl.show()