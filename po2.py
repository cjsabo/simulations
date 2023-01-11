import matplotlib.pyplot as pl
import numpy as np
import math as m
import cmath as cm
import pandas as pd
import qutip as qt
import sympy as sp

#Physical optics homework 2

#problem 2

#part a
nair = 1
nglass = 1.58
n1 = 1.38
n2 = 2.30
lam = 1;

l = np.linspace(.6*lam,1.6*lam, 1000)
t = np.linspace(0, np.pi/2, 1000)

lambda0, theta0 = np.meshgrid(l,t)

d0 = 0
t1 = -1*lam/(4*n1)
t2 = -1*lam/(4*n2)

theta1 = np.arcsin((nair/n1)*np.sin(theta0))
theta2 = np.arcsin((n1/n2)*np.sin(theta1))
theta3 = np.arcsin((n2/nglass)*np.sin(theta2))

k0z = (2*np.pi*nair/lambda0)*np.cos(theta0)
k1z = (2*np.pi*n1/lambda0)*np.cos(theta1)
k2z = (2*np.pi*n2/lambda0)*np.cos(theta2)
k3z = (2*np.pi*nglass/lambda0)*np.cos(theta3)

p32_TE = k2z/k3z
p32_TM = (nglass**3 * k2z)/(n2**3 * k3z)

R32_TE = (1-p32_TE)/(1+p32_TE)
R32_TM = (1-p32_TM)/(1+p32_TE)

p21_TE = k1z/k2z
p21_TM = (n2**3 * k1z)/(n1**3 * k2z)

R21_TE = (1-p21_TE)/(1+p21_TE)
R21_TM = (1-p21_TM)/(1+p21_TE)

p10_TE = k0z/k1z
p10_TM = (n1**3 * k0z)/(nair**3 * k1z)

R10_TE = (1-p10_TE)/(1+p10_TE)
R10_TM = (1-p10_TM)/(1+p10_TE)

V32_TE = [[np.exp(1j*k3z*(t1+t2)),R32_TE*np.exp(1j*k3z*(t1+t2))],[R32_TE*np.exp(-1j*k3z*(t1+t2)), np.exp(-1j*k3z*(t1+t2))]]

V32_TM = [[np.exp(1j*k3z*(t1+t2)),R32_TM*np.exp(1j*k3z*(t1+t2))],[R32_TM*np.exp(-1j*k3z*(t1+t2)), np.exp(-1j*k3z*(t1+t2))]]

V21_TE = [[np.exp(-1j*k2z*t2), R21_TE*np.exp(-1j*k2z*t2)],[R21_TE*np.exp(1j*k2z*t2), np.exp(1j*k2z*t2)]]

V21_TM = [[np.exp(-1j*k2z*t2), R21_TM*np.exp(-1j*k2z*t2)],[R21_TM*np.exp(1j*k2z*t2), np.exp(1j*k2z*t2)]]

V10_TE = [[np.exp(-1j*k1z*t1), R10_TE*np.exp(-1j*k1z*t1)],[R10_TE*np.exp(1j*k1z*t1), np.exp(1j*k1z*t1)]]

V10_TM = [[np.exp(-1j*k1z*t1), R10_TM*np.exp(-1j*k1z*t1)],[R10_TM*np.exp(1j*k1z*t1), np.exp(1j*k1z*t1)]]

matrixconstant_TE = (-1/8)*(1+p32_TE)*(1+p21_TE)*(1+p10_TE)
matrixconstant_TM = (-1/8)*(1+p32_TM)*(1+p21_TM)*(1+p10_TM)

V30_TE = np.matmul(V32_TE,np.matmul(V21_TE,V10_TE))
V30_TM = np.matmul(V32_TM,np.matmul(V21_TM,V10_TM))

R_TE = np.ones((1000,1000))
R_TM = np.ones((1000,1000))

for i in range(1000):
    for j in range(1000):
        R_TE[i,j] = abs(-1*V30_TE[0,1,i,j]/V30_TE[0,0,i,j])**2
        R_TM[i,j] = abs(-1*V30_TM[0,1,i,j]/V30_TM[0,0,i,j])**2
        
pl.figure(1)
ax = pl.axes(projection = '3d')
ax.contour3D(lambda0, theta0, R_TE, 50, cmap=pl.get_cmap('hot'))
ax.set_xlabel('lambda/lambda0')
ax.set_ylabel('incident angle (radians)')
ax.set_title('Reflectivity (TE)')

pl.figure(2)
ax = pl.axes(projection = '3d')
ax.contour3D(lambda0, theta0, R_TM, 50, cmap=pl.get_cmap('hot'))
ax.set_xlabel('lambda/lambda0')
ax.set_ylabel('incident angle (radians)')
ax.set_title('Reflectivity (TM)')


#print("Amplitude reflectivity (TE): ", -1*V30_TE[0,1]/V30_TE[0,0])
#print("Reflectivity squared (TE): ", abs(-1*V30_TE[0,1]/V30_TE[0,0])**2)


#problem 4

#part d
d = np.linspace(0,1e-6,1000)
lambda1 = 600e-9
n1_4 = 1.55
n2_4 = 1
n3_4 = 1.55
theta1 = 55*(cm.pi/180)
k_mag = (2*cm.pi*n1_4/lambda1)
k_x = k_mag*cm.cos(theta1)
p_x = k_mag*cm.cos(theta1)
q_x = k_mag*cm.sqrt((n2_4/n1_4)**2 - cm.sin(theta1)**2)

B = ((q_x-p_x)/(q_x+p_x))*np.exp(2*1j*q_x*d)

rs = ((1+B)*k_x - (1-B)*q_x)/((1-B)*q_x + (1+B)*k_x)

rs_squared = abs(rs)**2

pl.figure(3)
pl.plot(d*1e6,rs_squared,'b')
pl.xlim([0, 1])
pl.ylim([0, 1])
pl.xlabel('d (microns)')
pl.ylabel('Reflectivity')
pl.title('Reflectivity vs. d')
pl.grid()
pl.show()


