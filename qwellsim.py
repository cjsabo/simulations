import matplotlib.pyplot as pl
import numpy as np
import math as m
import pandas as pd
import scipy.special as sp
import scipy.integrate as intl
import sympy as sym

#Function to perform Euler method
def euler(q, m, hbar, x, V0, E, wf, wfp, wfdp, En, L):
    for k in range(0, 1000):
        if(x[k] < -1*L/2):
            if(k == 0):
                wf[k] = 0.01e-25
                wfp[k] = 0.01e-25
                wfdp[k] = (2*m0/hbar**2)*(V0 + q*E*x[k] - En)*wf[k]
            else:
                wf[k] = (wf[k-1] + (x[k] - x[k-1])*wfp[k-1])
                wfp[k] = (wfp[k-1] + (x[k] - x[k-1])*wfdp[k-1])
                wfdp[k] = ((2*m0/hbar**2)*(V0 + q*E*x[k] - En)*wf[k])
        elif(x[k] >= -1*L/2 and x[k] <= L/2):
            wf[k] = (wf[k-1] + (x[k] - x[k-1])*wfp[k-1])
            wfp[k] = (wfp[k-1] + (x[k] - x[k-1])*wfdp[k-1])
            wfdp[k] = ((2*m0/hbar**2)*(q*E*x[k] - En)*wf[k])
        else:
            wf[k] = (wf[k-1] + (x[k] - x[k-1])*wfp[k-1])
            wfp[k] = (wfp[k-1] + (x[k] - x[k-1])*wfdp[k-1])
            wfdp[k] = ((2*m0/hbar**2)*(V0 + q*E*x[k] - En)*wf[k])
    return wf, wfp, wfdp

#Find value for n that will make right side equal to 0, or close to it
def opt_n(wf, wfdp, E, x, m, hbar, q):
    En = -1*((wfdp[0]/wf[0])*((hbar**2)/(2*m)) - (V0 + q*E*x[len(wf)-1]))
    an = (-1*En)/((((q**2) * (E**2) * (hbar**2))/(2*m))**(1/3))
    n = ((-1*an)**3)/(2/(3*np.pi)) + (1/4)
    return n
    
#hbar
hbar = (6.626e-34)/(2*np.pi)

#charge of an electron
q = 1.602e-19

#mass of a free electron
m0 = 9.11e-31

#Plotting function from x = -5 to 5
x = np.linspace(-5e-9, 5e-9, 1000)

#height at edge of well
V0 = 1.602e-19

#width of well
L = 3e-9

#Electric Field
E = 1e7

#Potential (for now, a harmonic oscillator)
V = x**2
for i in range(0, 1000):
    if(x[i] < -1*L/2):
        V[i] = V0 + q*E*x[i]
    elif(x[i] >= -1*L/2 and x[i] <= L/2):
        V[i] = q*E*x[i] + (q*E*L/2)
    else:
        V[i] = V0 + q*E*x[i]
        
#Energy level to plot
n = .878

#Constant an used for energy levels
an = -1*((3*np.pi/2)*(n - 1/4))**(2/3)

#Energy for given n
if(E != 0):
    En = -1*an*(((q**2) * (E**2) * (hbar**2))/(2*m0))**(1/3)
else:
    En = (hbar**2/(2*m0))*((n*np.pi/L)**2)

#Wavefunction array with initial condition already set
wf = [0.01e-25]*1000

#Derivative of wavefunction with initial condition set
wfp = [0.01e-25]*1000

#2nd derivative of wavefunction with initial condition set
wfdp = [.01]*1000

#Calculating values for wavefunction, derivative and second derivative
wf, wfp, wfdp = euler(q, m0, hbar, x, V0, E, wf, wfp, wfdp, En, L)

#Solve for n that gives you best curve
#n = opt_n(wf, wfdp, E, x, m0, hbar, q)
#print(n)

an = -1*((3*np.pi/2)*(n - 1/4))**(2/3)

#Energy for given n
if(E != 0):
    En = -1*an*(((q**2) * (E**2) * (hbar**2))/(2*m0))**(1/3)
else:
    En = (hbar**2/(2*m0))*((n*np.pi/L)**2)

wf, wfp, wfdp = euler(q, m0, hbar, x, V0, E, wf, wfp, wfdp, En, L)
        
#Airy functions
ai, aip, bi, bip = sp.airy(((2*m0/(hbar**2 * q**2 * E**2))**(1/3))*(q*E*x - En))

pl.plot(x, V/V0, "r--")
#pl.plot(x, wf, "b-")
#pl.plot(x, ai, "g-")
#pl.ylim([-.5e-19, 2e-19])
#pl.xlim([-5e-9, 6e-9])

pl.grid()
pl.show()