import matplotlib.pyplot as pl
import numpy as np
import math as m
import pandas as pd
import qutip as qt
import sympy as sym

#Physical Optics Homework 3

#problem 4

#part a
sigma = 1
D = 1
z = np.linspace(0,.5,100)
z2 = np.linspace(0,5,1000)

sigz = np.sqrt(sigma**2 + 1j*2*D*z)

pl.figure(1)
pl.plot(z, abs(sigz))
pl.xlabel('z')
pl.ylabel('width')
pl.title('Width vs. position')

#part b
pl.figure(2)
bc = .5
a = (1/(1/(sigma**2) + 1j*bc)) + 1j*2*D*z2
a1 = 1/a
sigz2 = np.sqrt(1/(-1j*bc + a1))


pl.plot(z2, abs(sigz2))
pl.xlabel('z')
pl.ylabel('width')
pl.title('Width vs. position (chirped)')


#extra credit problem
l = np.linspace(.6, 1.6, 1000)      #wavelength in microns
dl = .01
n = np.sqrt((.6961663*l**2)/(l**2 - .0684043**2) + (.4079426*l**2)/(l**2 - .1162414**2) + (.8974794*l**2)/(l**2 - 9.896161**2) + 1)
nplus = np.sqrt((.6961663*(l+dl)**2)/((l+dl)**2 - .0684043**2) + (.4079426*(l+dl)**2)/((l+dl)**2 - .1162414**2) + (.8974794*(l+dl)**2)/((l+dl)**2 - 9.896161**2) + 1)
nmin = np.sqrt((.6961663*(l-dl)**2)/((l-dl)**2 - .0684043**2) + (.4079426*(l-dl)**2)/((l-dl)**2 - .1162414**2) + (.8974794*(l-dl)**2)/((l-dl)**2 - 9.896161**2) + 1)

#part a
pl.figure(3)
pl.plot(l,n,'r')
pl.xlabel("Wavelength (microns)")
pl.ylabel("Refractive index")
pl.title("Refractive index vs. wavelength")

#part b
pl.figure(4)
T = (1/3e8)*(n - l*(nplus - nmin)/(2*dl))
pl.plot(l,T)
pl.xlabel('Wavelength (microns)')
pl.ylabel('time (s/km)')
pl.title('Propagation time vs. wavelength')

#part c
pl.figure(5)
T1 = (3e8)*T
T1plus = np.gradient(T1)/np.gradient(l)
pl.plot(l,(1/3e8)*T1plus)
pl.xlabel('wavelength')
pl.ylabel('dispersion')
pl.title('dispersion vs. wavelength')

pl.figure(6)
x = np.linspace(0,5,1000)
y = x**4
dx = np.gradient(x)
dy = np.gradient(y)
pl.plot(x, dy/dx)


pl.show()
