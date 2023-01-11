import matplotlib.pyplot as pl
import numpy as np
import math as m
import pandas as pd
import qutip as qt
from mpl_toolkits.mplot3d import Axes3D

#code for HW4

#problem 3

#part a
R = 1
B0 = 1
#h = 6.626e-34
h = 1
e = 1.602e-19
hbar = h/(2*np.pi)
Ec = h*.2e9
Ej = h*5e9
phi0 = h/(2*e)

#Hamiltonian
n = []
for i in range(-10,11):
    n.append(i)


ng = np.linspace(-2,2, 200)
ev = []
ef = []
H0 = np.zeros([len(ng), len(n), len(n)])
for i in range(len(ng)):
    for j in range(len(n)):
        if n[j] == -10:
            H0[i] = H0[i] + 4*Ec*((n[j] - ng[i])**2)*qt.basis(len(n),n[j]+10)*qt.basis(len(n), n[j]+10).dag()
        elif n[j] == 10:
            H0[i] = H0[i] + 4*Ec*((n[j] - ng[i])**2)*qt.basis(len(n),n[j]+10)*qt.basis(len(n), n[j]+10).dag()
        else:
            H0[i] = H0[i] + 4*Ec*((n[j] - ng[i])**2)*qt.basis(len(n),n[j]+10)*qt.basis(len(n), n[j]+10).dag()-(Ej/2)*(qt.basis(len(n), n[j]+11)*qt.basis(len(n), n[j]+10).dag() + qt.basis(len(n), n[j]+9)*qt.basis(len(n), n[j]+10).dag())
    v,f = np.linalg.eig(H0[i])
    idx = np.argsort(v)
    v = v[idx]
    f = f[:][idx]
    ev.append(v)
    ef.append(f)

H1 = np.zeros([len(n), len(n)])
ng1 = 0
for j in range(len(n)):
        if n[j] == -10:
            H1 = H1 + 4*Ec*((n[j]-ng1)**2)*qt.fock(len(n),n[j], n[0])*qt.basis(len(n), n[j]+10).dag()
        elif n[j] == 10:
            H1 = H1 + 4*Ec*((n[j]-ng1)**2)*qt.fock(len(n),n[j], n[0])*qt.fock(len(n), n[j], n[0]).dag()
        else:
            H1 = H1 + 4*Ec*((n[j]-ng1)**2)*qt.fock(len(n),n[j], n[0])*qt.fock(len(n), n[j], n[0]).dag()-(Ej/2)*(qt.fock(len(n), n[j], n[0])*qt.fock(len(n), n[j]+1, n[0]).dag() + qt.fock(len(n), n[j]+1, n[0])*qt.fock(len(n), n[j], n[0]).dag())
            
v1, f1 = np.linalg.eig(H1)
idx1 =  np.argsort(v1)
v1 = v1[idx1]
f1 = f1[:,idx1]

ev0 = np.zeros(len(ng))
ev1 = np.zeros(len(ng))
ev2 = np.zeros(len(ng))
ev3 = np.zeros(len(ng))
ev4 = np.zeros(len(ng))
for i in range(len(ng)):
    ev0[i] = ev[i][0]
    ev1[i] = ev[i][1]
    ev2[i] = ev[i][2]
    ev3[i] = ev[i][3]
    ev4[i] = ev[i][4]
floor = min(ev0)

#problem 1
n1 = [0,1,2,3]
theta = np.linspace(0, 2*np.pi, 1000)

pl.figure(1)
pl.plot(ng, (ev0-floor)*1e-9)
pl.plot(ng, (ev1-floor)*1e-9)
pl.plot(ng, (ev2-floor)*1e-9)
pl.plot(ng, (ev3-floor)*1e-9)
pl.plot(ng, (ev4-floor)*1e-9)
pl.xlabel("ng")
pl.ylabel("E/E0")
pl.title("Energy Dispersion of Cooper Pair Box")

pl.figure(2)
pl.subplot(5,1,1)
pl.bar(n,f1[:,0])
pl.ylabel("psi0")
pl.title("Cooper Pair Box Wavefunctions")
pl.subplot(5,1,2)
pl.bar(n,f1[:,1])
pl.ylabel("psi1")
pl.subplot(5,1,3)
pl.bar(n,f1[:,2])
pl.ylabel("psi2")
pl.subplot(5,1,4)
pl.bar(n,f1[:,3])
pl.ylabel("psi3")
pl.subplot(5,1,5)
pl.bar(n,f1[:,4])
pl.xlabel("n")
pl.ylabel("psi4")

pl.figure(3)
for i in range(len(n1)):
    pl.plot(theta, n1[i]*theta)
pl.xlabel("Theta (radians)")
pl.ylabel("Phase (radians)")
pl.title("Phase vs. Theta")
pl.legend(["n = 0", "n = 1", "n = 2", "n = 3"])

pl.show()
    