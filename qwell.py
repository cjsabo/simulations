import numpy as np
import math as m
import matplotlib.pyplot as plt

def pot(V, E, x, q, L):
    p = 0.0
    if(x < -1*L/2):
        p = V + q*E*x
    elif(x >= -1*L/2 and x <= L/2):
        p = q*E*x
    else:
        p = V + q*E*x
    return p

#hbar
hbar = 1.055e-34

#charge of an electron
q = 1.602e-19

#Electric field
E = 1e7

#mass of a free electron 
m0 = 9.11e-31

#effective electron mass
eem = .067*m0

#Effective hole mass
ehm = .45*m0

#height of well (conduction band) (J)
V0 = .3*q

#Height of barrier (valence band)
Vh = .2*q

#Length of well
L = 10e-9

#Band gap
Eg = 1.9*q

#number of points
N = 500

#range of x
x = np.linspace(-20e-9, 20e-9, N)

#range of photon energy
Eph = np.linspace(0, 3.0e-19, N)

#delta x
dx = x[1] - x[0]

#kinetic energy matrix
T = np.zeros((N-2,N-2))

#filling kinetic energy matrix
for i in range(0,N-2):
    for j in range(0,N-2):
        if(i == j):
            T[i,j] = -2
        elif(np.abs(i-j) == 1):
            T[i,j] = 1
        else:
            T[i,j] = 0

#Potential energy matrices

#With Field
Vcond = np.zeros((N-2,N-2))
Vval = np.zeros((N-2, N-2))

#Without Field
Vcnf = np.zeros((N-2,N-2))
Vhnf = np.zeros((N-2,N-2))

#Fill in potential energy matrix
for i in range(0,N-2):
    for j in range(0,N-2):
        if(i == j):
            Vcond[i,j] = pot(V0, E, x[i+1], q, L)
            Vval[i,j] = pot(Vh, E, x[i+1], q, L)
            Vcnf[i,j] = pot(V0, 0, x[i+1], q, L)
            Vhnf[i,j] = pot(Vh, 0, x[i+1], q, L)
        else:
            Vcond[i,j] = 0
            Vval[i,j] = 0
            Vcnf[i,j] = 0
            Vhnf[i,j] = 0

#Graph for potential function
Vgraphc = [0]*N
Vgraphv = [0]*N
Vgraphcnf = [0]*N
Vgraphhnf = [0]*N

for p in range(0,N):
    if(x[p] < -1*L/2):
        Vgraphc[p] = (V0 + q*E*x[p])/q
        Vgraphv[p] = (V0 - Eg + q*E*x[p])/q
        Vgraphcnf[p] = V0/q
        Vgraphhnf[p] = (V0 - Eg)/q
    elif(x[p] >= -1*L/2 and x[p] <= L/2):
        Vgraphc[p] = (q*E*x[p])/q
        Vgraphv[p] = (q*E*x[p] + V0 + Vh - Eg )/q
        Vgraphcnf[p] = 0
        Vgraphhnf[p] = (V0 - Eg + Vh)/q
    else:
        Vgraphc[p] = (V0 + q*E*x[p])/q
        Vgraphv[p] = (V0 - Eg + q*E*x[p])/q
        Vgraphcnf[p] = V0/q
        Vgraphhnf[p] = (V0 - Eg)/q

#Create Hamiltonian Matrix
Hcond = -1*((hbar**2)/(2*eem*dx**2))*T + Vcond
Hval = -1*((hbar**2)/(2*ehm*dx**2))*T + Vval
Hcnf = -1*((hbar**2)/(2*eem*dx**2))*T + Vcnf
Hvnf = -1*((hbar**2)/(2*ehm*dx**2))*T + Vhnf

#Find eigenvalues and and eigenvectors

#Conduction band with field
valc, vecc = np.linalg.eig(Hcond)
z = np.argsort(valc)
z = z[0:2]
energiesc = (valc[z])

#Valence band with field
valv, vecv = np.linalg.eig(Hval)
v = np.argsort(valv)
v = v[0:2]
energiesv = (valv[v])

#Conduction band no field
valcnf, veccnf = np.linalg.eig(Hcnf)
cnf = np.argsort(valcnf)
cnf = cnf[0:2]
energiescnf = (valcnf[cnf])

#Valence band with field
valvnf, vecvnf = np.linalg.eig(Hvnf)
vnf = np.argsort(valvnf)
vnf = vnf[0:2]
energiesvnf = (valvnf[vnf])
                
#Plotting first 2 wave functions for both electrons and holes

#With Field
psic1 = [d + abs(energiesc[0])/q for d in vecc[:,z[0]]]
psic1 = np.append(psic1, abs(energiesc[0]/q))
psic1 = np.insert(psic1, 0, abs(energiesc[0]/q))

psic2 = [d + abs(energiesc[1])/q for d in vecc[:,z[1]]]
psic2 = np.append(psic2, abs(energiesc[1]/q))
psic2 = np.insert(psic2, 0, abs(energiesc[1]/q))

psiv1 = [d - (Eg-V0-Vh+abs(energiesv[0]))/q for d in vecv[:,v[0]]]
psiv1.reverse()
psiv1 = np.append(psiv1, -1*(Eg-V0-Vh-abs(energiesv[0]))/q)
psiv1 = np.insert(psiv1, 0, -1*(Eg-V0-Vh-abs(energiesv[0]))/q)

psiv2 = [d - (Eg-V0-Vh+abs(energiesv[1]))/q for d in vecv[:,v[1]]]
psiv2.reverse()
psiv2 = np.append(psiv2, -1*(Eg-V0-Vh-abs(energiesv[1]))/q)
psiv2 = np.insert(psiv2, 0, -1*(Eg-V0-Vh-abs(energiesv[1]))/q)

#Without field
psicnf1 = [d + abs(energiescnf[0])/q for d in veccnf[:,cnf[0]]]
psicnf1 = np.append(psicnf1, abs(energiescnf[0]/q))
psicnf1 = np.insert(psicnf1, 0, abs(energiescnf[0]/q))

psicnf2 = [d + abs(energiescnf[1])/q for d in veccnf[:,cnf[1]]]
psicnf2 = np.append(psicnf2, abs(energiescnf[1]/q))
psicnf2 = np.insert(psicnf2, 0, abs(energiescnf[1]/q))

psivnf1 = [d - (Eg-V0-Vh+abs(energiesvnf[0]))/q for d in vecvnf[:,vnf[0]]]
psivnf1.reverse()
psivnf1 = np.append(psivnf1, -1*(Eg-V0-Vh-abs(energiesvnf[0]))/q)
psivnf1 = np.insert(psivnf1, 0, -1*(Eg-V0-Vh-abs(energiesvnf[0]))/q)

psivnf2 = [d - (Eg-V0-Vh+abs(energiesv[1]))/q for d in vecvnf[:,vnf[1]]]
psivnf2.reverse()
psivnf2 = np.append(psivnf2, -1*(Eg-V0-Vh-abs(energiesvnf[1]))/q)
psivnf2 = np.insert(psivnf2, 0, -1*(Eg-V0-Vh-abs(energiesvnf[1]))/q)

#Matrix element without field
Mnf1 = np.dot(psicnf1, psivnf1)
Mnf2 = np.dot(psicnf2, psivnf2)

#Matrix element with field
M1 = np.dot(psic1, psiv1)
M2 = np.dot(psic2, psiv2)

#Absorption coefficient with field found by ratio of matrix elements
alpha1 = M1/Mnf1
alpha2 = M2/Mnf2

#Alpha function to be plotted
alpha = [0]*N

for i in range(0,N):
    if(Eph[i] < (Eg - V0 - Vh + energiesc[0] + energiesv[0])):
        alpha[i] = 0
    elif(Eph[i] >= (Eg - V0 - Vh + energiesc[0] + energiesv[0]) and Eph[i] <= (Eg - V0 - Vh + energiesc[1] + energiesv[1])):
        alpha[i] = alpha1
    elif(Eph[i] > (Eg - V0 - Vh + energiesc[1] + energiesv[1])):
        alpha[i] = alpha2

        
print(energiescnf[0]/q)
print(energiescnf[1]/q)
print(energiesvnf[0]/q)
print(energiesvnf[1]/q)

print(energiesc[0]/q)
print(energiesc[1]/q)
print(energiesv[0]/q)
print(energiesv[1]/q)
#Plotting potential well and wavefunctions together

#With Field
plt.subplot(1,3,1)
plt.plot(x, psic1, "b-")
plt.plot(x, psiv1, "b-")
plt.plot(x, psic2, "g-")
plt.plot(x, psiv2, "g-")
plt.plot(x, Vgraphc, "r--")
plt.plot(x, Vgraphv, "r--")
plt.xlim([-10e-9, 10e-9])
plt.title("Psi's vs. Pos.")
plt.xlabel("Position")
plt.ylabel("Energy (eV)")
plt.grid()

#Without Field
plt.subplot(1,3,2)
plt.plot(x, psicnf1, "b-")
plt.plot(x, psivnf1, "b-")
plt.plot(x, psicnf2, "g-")
plt.plot(x, psivnf2, "g-")
plt.plot(x, Vgraphcnf, "r--")
plt.plot(x, Vgraphhnf, "r--")
plt.title("Psi's vs. Pos.")
plt.xlabel("Position")
plt.xlim([-10e-9, 10e-9])
plt.grid()

#Absorption Coefficient
plt.subplot(1,3,3)
plt.plot(Eph, alpha, "c-")
plt.title("Abs. vs. Eph")
plt.xlabel("Photon Energy (eV)")
plt.ylabel("Alpha")

plt.show()
