import matplotlib.pyplot as pl
import numpy as np
import math as m
from qutip import *


hbar = (6.626e-34)/(2*m.pi)
w0 = 2*m.pi*1
#Problem 1

#part a
#identity
i = Qobj([[1,0], [0,1]])
print('Eigenenergies of identity: ', i.eigenenergies())
print('Eigenvectors of identity: ', i.eigenstates())

#sigmax
sx = Qobj([[0,1], [1,0]])
print('Eigenenergies of sigmax: ', sx.eigenenergies())
print('Eigenvectors of sigmax: ', sx.eigenstates())

#sigmay
sy = Qobj([[0,-1j], [1j,0]])
print('Eigenenergies of sigmay: ', sy.eigenenergies())
print('Eigenvectors of sigmay: ', sy.eigenstates())

#sigmaz
sz = Qobj([[1,0], [0,-1]])
print('Eigenenergies of sigmaz: ', sz.eigenenergies())
print('Eigenvectors of sigmaz: ', sz.eigenstates())

#part b
#Commutation xy
print('Commutation of x and y:')
sxy = sx*sy;
syx = sy*sx;
print(sxy, ' ', syx, ' ', sxy-syx)

#Commutation yz
print('Commutation of y and z:')
syz = sy*sz;
szy = sz*sy;
print(syz, ' ', szy, ' ', syz-szy)

#Commutation zx
print('Commutation of z and x:')
szx = sz*sx;
sxz = sx*sz;
print(szx, ' ', sxz, ' ', szx-sxz)

#part c
commxy = sxy - syx
zk = Qobj([[1],[0]])
zb = Qobj([[1,0]])
unc = abs(commxy.matrix_element(zb,zk)/2)
print('uncertainty: ', unc)



#problem 2
H2 = (-1/2)*hbar*sz
print(H2.eigenenergies(), H2.eigenstates())

#problem 5
syy = tensor(sy,sy)
print(syy)
