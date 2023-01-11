import matplotlib.pyplot as pl
import numpy as np
import math as m
import pandas as pd
import qutip as qt
import scipy as sp

def rotMat(theta):
    rot = np.matrix([[np.cos(theta), np.sin(theta)],[-1*np.sin(theta), np.cos(theta)]])
    return rot

def func(x):
    p = x[0]
    c = x[1]
    F1 = np.zeros(2)
    F1[0] = (.177-1j*.354)*(np.cos(p)*np.cos(c) - 1j*np.sin(p)*np.sin(c)) + .306*(np.sin(p)*np.cos(c) + 1j*np.cos(p)*np.sin(c))
    F1[1] = (.306-1j*.612)*(np.cos(p)*np.cos(c) - 1j*np.sin(p)*np.sin(c)) + .530*(np.sin(p)*np.cos(c) + 1j*np.cos(p)*np.sin(c))
    return F1

#Physical optics hw 4

#problem 1

#matrix for polarizer
tp1 = 60*(np.pi/180)
tph = np.matrix([[1,0],[0,0]])
p1 = np.matmul(rotMat(-1*tp1), np.matmul(tph,rotMat(tp1)))

#matrix for quarter wave plate
tqwp1 = 30*(np.pi/180)
qwph = np.exp(-1j*np.pi/4)*np.matrix([[1,0],[0,1j]])
qwp1 = np.matmul(rotMat(-1*tqwp1), np.matmul(qwph, rotMat(tqwp1)))

#total polarization matrix
pol1 = np.matmul(p1, qwp1)

#solve for Ex and Ey
Ex = 1/np.sqrt(1+abs((-.177 + 1j*.354)/.306)**2)
Ey = ((-.177+1j*.354)/.306)*Ex
alpha1 = np.arctan(-1*abs(Ey/Ex))
phiyx1 = np.arctan(.708/-.354)
psi1 = .5*np.arctan(np.cos(phiyx1)*np.tan(2*alpha1))
chi1 = np.arcsin(np.sin(2*alpha1)*np.sin(phiyx1))/2
print("Ex = ", Ex)
print("Ey = ", Ey)
print("alpha = ", alpha1*(180/np.pi))
print("Phiy - Phix  = ", phiyx1*(180/np.pi))
print("Psi = ", psi1*(180/np.pi))
print("Chi = ", chi1*(180/np.pi))
print("\r\n")

#problem 2
m1 = np.matrix([[1,0],[0,0]])
m2 = (1/np.sqrt(2))*np.matrix([[1,1j], [1j,1]])
m3 = np.matrix([[1,0],[0,-1]])

m2b = (1/np.sqrt(2))*np.matrix([[1,-1j], [-1j,1]])

mf = np.matmul(m1, np.matmul(m2b,np.matmul(m3,np.matmul(m2,m1))))
print("Final matrix: ", mf)

#problem 4


