import matplotlib.pyplot as plt
import numpy as np
import math as m
import pandas as pd
import qutip as qt
from scipy.sparse import diags, linalg, kron

L     = 300e-9 # m length the box 
N_pts = 101 # number of points for the phase space

e = 1.602e-19 #charge of an electron
N = np.linspace(0,25,26)
r = 98e-9 #radius of quantum dot (calculated)
eps0 = 8.85e-12 #permittivity of vaccuum
epsr = 12.9 #relative dielectric constant for GaAs
dist = 10e-9 #distance between "plates"
C = (eps0*epsr*np.pi*(r**2))/(10e-9) #capacitance of gate
Ec = (e**2)/(2*C) #charging energy (capacitive energy)
Ec_ev = Ec/e
#print(C)
#print(Ec)
#print(Ec_ev)

omega = 2*np.pi*(200e9) #angular frequency

phi_pts_x= np.linspace(-L/2, L/2, N_pts) # phi vector for x

dp        = phi_pts_x[2] - phi_pts_x[1] # distance between points
d1_coeff  = (1. / (2. * dp)) # coefficient for the first derivative
d2_coeff  = (1. / (dp ** 2)) # coefficient for the second derivative

# defining operators for open boundary conditions
id_op = diags(np.ones(phi_pts_x.size), 0, shape=(N_pts,N_pts))
x_op  = diags(phi_pts_x, 0, shape=(N_pts,N_pts))
x2_op = diags(phi_pts_x**2, 0, shape=(N_pts,N_pts))
d1x   = diags([-d1_coeff, d1_coeff], [-1,1], shape=(N_pts,N_pts))
d2x   = diags([d2_coeff, -2.0 * d2_coeff, d2_coeff], [-1,0,1], shape=(N_pts,N_pts))

# first derivative operator
#print(d1x.toarray() / d1_coeff)

# second derivative operator
#print(d2x.toarray() / d2_coeff)

# x^2 operator
#print(x2_op.toarray())

hbar = 6.582e-16 # eV s
m_e  = 5.685e-12 # eV s2/m2 free electron mass
#hbar = (6.626e-34)/(2*np.pi)
#m_e = 9.11e-31

m    = m_e * 0.063 # eV s2/m2 electron mass in GaAs

E1 = hbar ** 2 / 2 / m

# the Hamiltonian of the box
H = - E1 * kron(d2x, id_op) - E1 * kron(id_op,d2x) + .5*m*(omega**2)*(kron(x2_op,id_op) + kron(id_op,x2_op))

# solve numerically the Hamiltonian for the first k=11 energy levels
evals, ekets = linalg.eigsh(H, k=26, which='SA')

# sort the eigenenergies and eigenvalues
sort_idxs    = evals.argsort()
evals        = np.sort(evals)
zero_energy  = evals[0]
evals        = evals - zero_energy # Define the 0 energy with respect to the lowest lying state.
ekets        = [ekets[:,idx].reshape([N_pts, N_pts]) for idx in sort_idxs] # reshaping the wavefunctions for 2D

mu = np.zeros(len(N)) #chemical potential
for i in range(len(mu)):
    mu[i] = 2*Ec_ev*N[i]+(evals[i])
    
#for i in range(len(mu)-1):
#    print(mu[i+1]-mu[i])
    
sum1, sum2 = 0,0
for i in range((N_pts)):
    for j in range((N_pts)):
        sum2 += abs(ekets[0][i,j])**2

count = -1
while((sum1/sum2) <= .999):
    sum1 = 0
    count += 1
    for i in range(-1*count,count,1):
        for j in range(-1*count, count, 1):
            sum1 += abs(ekets[0][50+i,50+j])**2
print(count)

fig, axs = plt.subplots(figsize=(4,3))

im1 = axs.imshow(np.real(ekets[0][:,:]),
           extent = [phi_pts_x[0] * 1e9, phi_pts_x[-1] * 1e9, phi_pts_x[0] * 1e9, phi_pts_x[-1] * 1e9], 
           cmap='RdBu', vmin=-0.05, vmax=0.05)

u = 0
for i in range(N_pts):
    for j in range(N_pts):
        if (np.real(ekets[0][i,j]))**2 > u:
            u = (np.real(ekets[0][i,j]))**2

#print(u)
#for i in range(N_pts):
#    for j in range(N_pts):
#        if (np.real(ekets[0][i,j]))**2 >= .001*(u):
#            print(i, " ", j)
            
fig.colorbar(im1, ax=axs, label = r'$\Re\Psi$')
axs.set_xlabel(r'$x$ (nm)')
axs.set_ylabel(r'$y$ (nm)')

plt.tight_layout()
#plt.show()

# plot the first 10 orbital energies

fig, axs = plt.subplots(figsize=(6,4))

axs.plot(np.arange(0,26), evals[0:26] * 1e3, '-o')

axs.set_xlabel(r'idx')
axs.set_ylabel(r'$E_\mathrm{orb}$ (meV)')

plt.tight_layout()

fig, axs = plt.subplots(figsize=(6,4))
axs.plot(N,mu*1e3,'-o')
axs.set_xlabel('Number of electrons (N)')
axs.set_ylabel('$\mu (N)$ (meV)')

plt.show()




