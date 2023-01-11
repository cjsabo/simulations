import matplotlib.pyplot as pl
import numpy as np
import math as m
import pandas as pd
import qutip as qt

#Homework 5 code

#problem 1

#part a

t = np.linspace(-2,2,1000)
o_01 = 1
lo_01 = 1
tgate = 1

V = o_01*np.cos(lo_01*t)*np.exp(-1*(t**2)/(2*tgate**2))
pl.figure(1)
pl.plot(t,V)
pl.xlabel('Time (s)')
pl.ylabel('V(t)')


tg = 5
w = np.linspace(-10, 10, 1000)
O_w = (o_01/2)*np.exp((-1*((w-lo_01)**2)*tg**2)/2) + (o_01/2)*np.exp((-1*((w+lo_01)**2)*tg**2)/2)
pl.figure(2)
pl.plot(w, O_w)
pl.xlabel("Frequency (w)")
pl.ylabel("F(V(t))")



#part d

#gate time
tg = np.linspace(0,50e-9,1000)

#dephasing time
T1 = 1e-6

#lambda (can change this value)
lam = .01

#Transition frequency
O_01 = 2*np.pi*(.2e9)

#function
P1 = (1/(1+lam**2))*((np.sin(.5*O_01*np.sqrt(1+lam**2)*tg))**2)*np.exp(-1*tg/T1)

#plot
pl.figure(3)
pl.plot(tg*1e9,P1)
pl.xlabel('Gate Time (ns)')
pl.ylabel('Excited State population')

pl.show()

