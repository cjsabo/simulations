import numpy as np
import math as m
import matplotlib.pyplot as plt

def eigenener(m, Ej, Ec):
    return -1*Ej + np.sqrt(8*Ej*Ec)*(m + 1/2) - (Ec/12)*(6*(m**2) + 6*m + 3)

#range of phi
phi = np.linspace(-1*np.pi, np.pi, 1000)

#charge of an electron
q = 1.602e-19

#Planck's constant
h = 6.626e-34

#epsilon naught
ep0 = 8.85e-12

#epsilon of substrate material. material being used:Silicon dioxide
epox = 3.9

#effective relative dielectric constant
eeff = 2

#normalized plancks constant
hbar = h/(2*np.pi)

#external magnetic field
B = 1e-3

#flux quantum
phi0 = h/(2*q)

#SQUID loop area
A = phi0/B

#SQUID loop side length and perimeter
s = np.sqrt(A)
p = 4*s

#Josephson junction current density (A/m^2)
Jc = (500)*(10000)

#Tunnel barrier
tb = 2e-9

#Ratio of EJ to Ec
Ejc = 50.0

#frequency
f = 8e9

#loop variable
i = 1
while i == 1:
    f = int(input("Please enter a frequency to operate at: "))
    print("\r")
    #Charging energy (Ec)
    Ec = (h*f)/(np.sqrt(Ejc*8)) 
    #Josephson Junction energy (EJ)
    Ej = Ejc * Ec 
    #Solving for critical current
    Ic = (2*q*Ej)/h 
    print("Charging Energy (meV): ", (Ec/q)*1000)
    print("Josephson Junciton Energy (meV): ", (Ej/q)*1000)
    print("Critical Current/Super current (nA): ", Ic/1e-9, "\r\n")
    #solving for total capacitance
    C = (q**2)/(2*Ec) 
    #solving for length of the capacitor
    l = (C)/(ep0*(eeff))
    print("Total Capacitance (fF):", C/1e-15)
    print("Length of interdigitated capacitor (mm):", l/1e-3, "\r\n")
    #length of the resonator
    lr = (3e8)/(np.sqrt(eeff)*f)
    print("Length of the resonator (mm): ", lr/1e-3, "\r\n")
    #value for ground state energy, first excited and second excited
    E0 = eigenener(0, Ej, Ec)
    E1 = eigenener(1, Ej, Ec)
    E2 = eigenener(2, Ej, Ec)
    print("Ground state energy (meV): ", (E0/q)*1000)
    print("First excited state energy (meV): ", (E1/q)*1000)
    print("Second excited state energy (meV): ", (E2/q)*1000)
    print("Anharmonicity factor: ", (((E2-E1)-(E1-E0))/(E1-E0)), "\r\n")
    #Josephson Junction dimensions
    print("Josephson Junction current density (A/cm^2): ", Jc/10000)
    print("Josephson junction size (um^2): ", (Ic/Jc)*(1e12), "\r\n")
    #SQUID loop size
    print("SQUID loop size (um^2): ", (s/(1e-6))**2, "\r\n")
    i = int(input("Would you like to try another frequency? (1 for yes, 0 for no) "))
    print("\r")
