import matplotlib.pyplot as pl
import numpy as np
import math as m
import pandas as pd
import qutip as qt
from mpl_toolkits.mplot3d import Axes3D
import scipy as sp

#Physical Optics final project
#Polarimeter, rotating quarter-wave plate and fixed polarizer

#rotation matrix
def Rotation(theta):
    rot = np.matrix([[1,0,0,0], [0,np.cos(2*theta),np.sin(2*theta),0], [0,-1*np.sin(2*theta),np.cos(2*theta),0],[0,0,0,1]])
    return rot

#quarter-wave plate matrix
def qwp():
    wp = np.matrix([[1,0,0,0],[0,1,0,0],[0,0,0,-1],[0,0,1,0]])
    return wp

#quarter waveplate with change in lambda
def qwp1(l,dl):
    wp = np.matrix([[1,0,0,0],[0,1,0,0],[0,0,np.cos((np.pi/2) + (dl/l)*(np.pi/2)),-1*np.sin((np.pi/2) + (dl/l)*(np.pi/2))],[0,0,np.sin((np.pi/2) + (dl/l)*(np.pi/2)),np.cos((np.pi/2) + (dl/l)*(np.pi/2))]])
    return wp

#horizontal polarizer matrix
def hp():
    p = .5*np.matrix([[1,1,0,0],[1,1,0,0],[0,0,0,0],[0,0,0,0]])
    return p

#full matrix
def full(theta):
    f = .5*np.matrix([[1,(np.cos(2*theta))**2,np.sin(2*theta)*np.cos(2*theta), np.sin(theta)],[1,(np.cos(2*theta))**2,np.sin(2*theta)*np.cos(2*theta), np.sin(theta)], [0,0,0,0],[0,0,0,0]])
    return f
 
#initial state of polarization
isop = np.matrix([[1],[.59],[-.15],[.25]])

#angular velocity in which to rotate waveplate
w_0 = 2*np.pi*1e6

#1 period 
T = 2*np.pi/(w_0)

#time
t = np.linspace(0,T,1000)

#hold intensities
int_out = np.array(np.zeros(len(t)))

s1 = np.array(np.zeros(len(t)))
s2 = np.array(np.zeros(len(t)))
s3 = np.array(np.zeros(len(t)))

#wavelength designed for
l0 = 600e-9

#deviation from that wavelength
dl = 0

for i in range(len(t)):
    x = np.matmul(hp(),np.matmul(Rotation(-1*w_0*t[i]),np.matmul(qwp1(l0,dl),np.matmul(Rotation(w_0*t[i]),isop))))
    #x1 = np.matmul(full(w_0*t[i]),isop)
    s1[i] = x[1]
    s2[i] = x[2]
    s3[i] = x[3]
    int_out[i] = x[0]   
    
#Determine initial state of polarization
points = [s1,s2,s3]

#parameter A
A = (2/T)*sp.integrate.trapz(int_out, t)
#print(A)

#Parameter B
B = (4/(T))*sp.integrate.trapz(int_out*np.array(np.sin(2*w_0*t)), t)
#print(B)

#Parameter C
C = (4/(T))*sp.integrate.trapz((int_out*np.array(np.cos(4*w_0*t))), t)
#print(C)

#parameter D
D = (4/(T))*sp.integrate.trapz((int_out*np.array(np.sin(4*w_0*t))), t)

#Original Stokes Parameters
print('S0: ', A-C)
print('S1: ', 2*C)
print('S2: ', 2*D)
print('S3: ', B)

pl.figure(1)
pl.plot(t*1e6,int_out)
pl.xlabel('Time (us)')
pl.ylabel('Intensity I(t)')
pl.title('Intensity vs. Time')

fig, axs = pl.subplots(figsize=(4, 4), subplot_kw=dict(projection='3d'))
bloch = qt.Bloch(fig=fig, axes=axs)
bloch.point_marker = ['o']
bloch.point_size = [10]
bloch.view = [40,35]
bloch.add_points(points)
bloch.xlabel = ['$S1$','']
bloch.ylabel = ['$S2$','']
bloch.zlabel = ['$S3$','']
bloch.render()


pl.show()
