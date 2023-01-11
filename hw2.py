import matplotlib.pyplot as pl
import numpy as np
import math as m

#problem 2
T = 600
m0 = 9.11e-31
me = 1.08*m0
mh = .81*m0
k = 1.38e-23
h = 6.626e-34
#eV
Eg0 = 1.166
a = .437e-3
b = 636
q = 1.602e-19
Nd = 1e17
Na = 5e16
Tn = T/300
N = Na + Nd

print("Problem 2")
#Effective density of states conduction band, cm^-3
Nc = (2*((2*np.pi*me*k*T)/h**2)**(3/2))*(1/(100**3))
print("Nc = ", Nc)

#Effective density of states valence band, cm^-3
Nv = (2*((2*np.pi*mh*k*T)/h**2)**(3/2))*(1/(100**3))
print("Nv = ", Nv)

#Temperature dependence of the bandgap
Eg = Eg0 - (a*T**2)/(b+T)
print("Eg = ", Eg)

#intrinsic carrier density
ni = np.sqrt(Nc*Nv)*np.exp(-1*(Eg*q)/(2*k*T))
print("ni = ", ni)

#electron density
n = Nd-Na
print("n = ", n)

#hole density
p = (ni**2)/n
print("p = ", p)

#intrinsic Fermi energy
Ei = Eg/2 + .5*k*T*q*np.log(Nv/Nc)
print("Ei = ", Ei)

#Fermi energy
Ef = Ei + k*T*q*np.log(n/ni)
print("Ef = ", Ef)

#electron mobility
mun = 88*((Tn)**(-.57)) + (1250*(Tn**-2.33))/(1+(N/(1.26e17*Tn**2.4))*.88*Tn**-.146)
print("mun = ", mun)

#hole mobility
mup = 54.3*((Tn)**(-.57)) + (407*(Tn**-2.33))/(1+(N/(2.35e17*Tn**2.4))*.88*Tn**-.146)
print("mup = ", mup)

#resistivity
rho = 1/(q*(mun*n + mup*p))
print("rho = ", rho)

print("\r\n")

#problem 3
print("Problem 3")

