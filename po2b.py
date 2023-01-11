import matplotlib.pyplot as pl
import numpy as np
import math as m
import cmath as cm
import pandas as pd
import qutip as qt
import sympy as sp

#only works for normal incidence
def getMat(l, n1, n2, d):
    V = [[0.5*(1+(n1/n2))*np.exp(-1j*(2*np.pi/l)*(n2-n1)*d), 0.5*(1+(n1/n2))*((n2-n1)/(n2+n1))*np.exp(-1j*(2*np.pi/l)*(n2+n1)*d)], [0.5*(1+(n1/n2))*((n2-n1)/(n2+n1))*np.exp(1j*(2*np.pi/l)*(n2+n1)*d), 0.5*(1+(n1/n2))*np.exp(1j*(2*np.pi/l)*(n2-n1)*d)]]
    return V

#problem 6

#refractive indices
n0 = 1
nh = 2.3
nl = 1.38
ns = 1.52
lam0 = 1300e-9
dh = lam0/(4*nh)
dl = lam0/(4*nl)
k0 = 2*np.pi*n0/lam0
Rt = .9934

#for pairs of layers
N = np.linspace(1,15,15)

#wavelength range (nm)
lam = np.linspace(1000e-9, 1600e-9, 1500)

#position (nm)
z = np.linspace(-500e-9, dh+dl, 2000)

E = np.zeros(len(z))
Eang = np.zeros(len(z))
b1 = 0
b2 = 0
E0pl = 0
E0ml = 0
E1pl = 0
E1ml = 0
V1 = getMat(lam0, n0, nl, 0)
V2 = np.matmul(getMat(lam0, nl, nh, 0), V1)

for i in range(len(z)):
    if z[i] < 0:
        E0p = np.exp(-1j*k0*z[i])
        E0m = Rt*np.exp(1j*k0*z[i])
        E[i] = abs(E0p+E0m)
        Eang[i] = np.angle(E0p+E0m)
        E0pl = E0p
        E0ml = E0m
    elif z[i] >= 0 and z[i] < dl:
        E1p = (V1[0][0]*E0pl + V1[0][1]*E0ml)*np.exp(-1j*k0*nh*z[i])
        E1m = (V1[1][0]*E0pl + V1[1][1]*E0ml)*Rt*np.exp(1j*k0*nh*z[i])
        E1pl = E1p
        E1ml = E1m
        E[i] = abs(E1p + E1m)
        Eang[i] = np.angle(E1p + E1m)
    elif z[i] >= dl and z[i] <= dl+dh:
        E2p = (V2[0][0]*E1pl + V2[0][1]*E1ml)*np.exp(-1j*k0*nl*z[i])
        E2m = (V2[1][0]*E1pl + V2[1][1]*E1ml)*Rt*np.exp(1j*k0*nl*z[i])
        E[i] = abs(E2p + E2m)
        Eang[i] = np.angle(E2p + E2m)
    else:
        E[i] = 0
        Eang[i] = 0

R_total = np.zeros((len(N), len(lam)))
R_squared = np.zeros((len(N), len(lam)))
R_phase = np.zeros((len(N), len(lam)))

for i in range(len(N)):
    for u in range(len(lam)):
        d_tot = (N[i])*dh + N[i]*dl
        V = getMat(lam[u], nh, ns, d_tot)
        for j in range(int(N[i]-1)):
            d_tot = d_tot-dh
            V = np.matmul(V, getMat(lam[u], nl, nh, d_tot))
            d_tot = d_tot-dl
            V = np.matmul(V, getMat(lam[u], nh, nl, d_tot))
        d_tot = d_tot - dh
        V = np.matmul(V, getMat(lam[u], nl, nh, d_tot))
        V = np.matmul(V, getMat(lam[u], n0, nl, 0))
        R_squared[i,u] = abs(-1*V[1,0]/V[1,1])**2
        R_phase[i,u] = np.angle(-1*V[1,0]/V[1,1])

pl.figure(1)
pl.subplot(2,1,1)
pl.plot(lam*1e9, R_squared[0,:], 'r')
pl.plot(lam*1e9, R_squared[1,:], 'b')
pl.plot(lam*1e9, R_squared[6,:], 'g')
pl.plot(lam*1e9, R_squared[14,:], 'c')
pl.ylabel("Reflectivity")
pl.title("Refelctivity vs. wavelength")
pl.legend(['N = 1', 'N = 2', 'N = 7', 'N = 15'])

pl.subplot(2,1,2)
pl.plot(lam*1e9, R_phase[0,:], 'r')
pl.plot(lam*1e9, R_phase[1,:], 'b')
pl.plot(lam*1e9, R_phase[6,:], 'g')
pl.plot(lam*1e9, R_phase[14,:], 'c')
pl.xlabel("Wavelength (nm)")
pl.ylabel("Angle (R)")
pl.title("Angle vs. wavelength")

pl.figure(2)
pl.subplot(2,1,1)
pl.plot(lam*1e9, R_squared[5,:])
pl.ylabel("Reflectivity")
pl.title("Reflectivity vs. wavelength (N = 6)")
pl.subplot(2,1,2)
pl.plot(lam*1e9, R_phase[5,:])
pl.xlabel("Wavelength (nm)")
pl.ylabel("Angle(R)")
pl.title("Angle vs. wavelength (N = 6)")

pl.figure(3)
pl.subplot(2,1,1)
pl.plot(z*1e6, E)
pl.title("E vs. position")
pl.ylabel("E")
pl.xlim([-.5, .21])
pl.subplot(2,1,2)
pl.plot(z*1e6, Eang)
pl.title("Angle(E) vs. position")
pl.ylabel("Angle(E)")
pl.xlabel("Position (um)")
pl.xlim([-.5, .21])

pl.figure(4)
pl.subplot(2,1,1)
pl.plot(lam*1e9, R_squared[5,:])
#pl.xlabel("Wavelength (nm)")
pl.ylabel("Reflectivity")
pl.title("Reflectivity vs. Wavelength (N = 6, missing layer)")
pl.subplot(2,1,2)
pl.plot(lam*1e9, R_phase[5,:])
pl.xlabel("Wavelength (nm)")
pl.ylabel("Angle(R)")
pl.title("Angle(R) vs. Wavelength (N = 6, missing layer)")

pl.show()



