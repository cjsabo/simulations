import matplotlib.pyplot as pl
import numpy as np
import math as m
import pandas as pd
import qutip as qt
from mpl_toolkits.mplot3d import Axes3D

#Code for HW3

#problem 3

#constants
w0 = 2*np.pi*5e9                                      #resonant frequency
wL = w0                                               #drive frequency
bigO = 2*np.pi*5e6                                    #Rabi drive frequency
T1 = 5e-6                                             #seconds, relaxation time
T_phi = 1e-6                                          #seconds, dephasing time
t_gate = np.linspace(0,100e-9, 1001)                  #seconds, gate time

#Pauli matrices
sx = qt.sigmax()
sy = qt.sigmay()
sz = qt.sigmaz()
sm = qt.Qobj([[0,1],[0,0]])
sp = qt.Qobj([[0,0],[1,0]])

#jump operators
Lm = np.sqrt(1/T1)*sm
L_phi = np.sqrt(1/(2*T_phi))*sz

#make initial state the |1> state
psi_init = qt.basis(2,1)

#Rabi Hamiltonian
H_rabi = -.5*(w0-wL)*sz + .5*bigO*sx

#master equation solver
result = qt.mesolve(H_rabi, psi_init, t_gate, c_ops=[Lm, L_phi])

#extract point for bloch sphere
r_x = []
r_y = []
r_z = []

for i in range(0,len(t_gate), 10):
    state = result.states[i]
    r_x.append(qt.expect(sx,state))
    r_y.append(qt.expect(sy,state))
    r_z.append(qt.expect(sz,state))
    
points = [r_x, r_y, r_z]

fig, axs = pl.subplots(figsize=(4, 4), subplot_kw=dict(projection='3d'))

bloch = qt.Bloch(fig=fig, axes=axs)
bloch.point_marker = ['o']
bloch.point_size = [10]
bloch.view = [40,35]
bloch.add_points(points)
bloch.render()


#calculate fidelity of gate
psi_fin_ideal = qt.basis(2,0)
psi_fin = psi_fin_ideal*psi_fin_ideal.dag()

psi_fin_actual = result.states[len(t_gate)-1]

arg = (psi_fin.sqrtm()*psi_fin_actual*psi_fin.sqrtm()).sqrtm()

#final fidelity
fid = arg.tr()**2
print(fid)

pl.show()
