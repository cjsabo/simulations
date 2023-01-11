import matplotlib.pyplot as pl
import numpy as np
import math as m
import pandas as pd
import qutip as qt
import sympy as sp

#function for problem 2
def getRho(g,t):
    Rho = qt.Qobj([[1/4, (-1/4)*np.exp(-1j*2*g*t), (1/4)*np.exp(-1j*2*g*t), -1/4], [(-1/4)*np.exp(1j*2*g*t), 1/4, -1/4, (1/4)*np.exp(1j*2*g*t)], [(1/4)*np.exp(1j*2*g*t), -1/4, 1/4, (-1/4)*np.exp(1j*2*g*t)], [-1/4, (1/4)*np.exp(1j*2*g*t), (-1/4)*np.exp(-1j*2*g*t), 1/4]])
    return Rho

def getPT(g,t):
    pt = qt.Qobj([[(1/2)*np.cos(2*g*t), 1/2], [1/2, (1/2)*np.cos(2*g*t)]])
    return pt

def H_rabit(t,w,s):
    H_op = .5*w*s
    return [H_op, np.ones(len(t))]

def H_coupled(t, g, sxa, sxb, sya, syb, sm_A, sp_A, sm_B, sp_B):
    H_1 = [-0.5 * g * sp_A * sm_B, np.ones(len(t))] 
    H_2 = [-0.5 * g * sm_A * sp_B, np.ones(len(t))]
    H_int = [.5*g*(sxa*sxb + sya*syb), np.ones(len(t))]
    return H_int
    

#problem 1

#part b 
f = 1
w0 = 2*np.pi*f
t = np.linspace(0, 1/f, 1000)

pl.figure(1)
pl.plot(t, np.sin(w0*t))
pl.xlabel("time (s)")
pl.ylabel("Tr")
pl.title("Tr vs. time")

#part d and e
p = np.linspace(0,1,1000)
tr_pw = np.zeros(1000)
doe = np.zeros(len(p))

for i in range(len(p)):
    rho_w = qt.Qobj([[(1-p[i])/4,0,0,0], [0, (p[i]+1)/4, (1-p[i])/4, 0], [0, (1-p[i])/4, (p[i]+1)/4, 0], [0, 0, 0, (1-p[i])/4]])
    rho_w2 = rho_w*rho_w
    tr_pw[i] = rho_w2.tr()
    doe[i] = qt.entropy_linear(rho_w)
    

pl.figure(2)
pl.plot(p, tr_pw)
pl.xlabel("p")
pl.ylabel("Tr(p_w^2)")
pl.title("Tr vs. p")

pl.figure(3)
pl.plot(p, doe)
pl.xlabel("p")
pl.ylabel("Degree of entanglement")
pl.title("Degree of entanglement vs. p")

#problem 2

t1 = np.linspace(0,5,1000)
g1 = 2
purity = np.zeros(len(t))
concurr = np.zeros(len(t))

for i in range(len(t1)):
    state = qt.Qobj([[.5*np.exp(-1j*g1*t1[i])], [-.5*np.exp(1j*g1*t1[i])], [.5*np.exp(1j*g1*t1[i])], [-.5*np.exp(-1j*g1*t1[i])]])
    rab = state*state.dag()
    rab2 = rab*rab
    purity[i] = rab2.tr()
    concurr[i] = np.sin(2*g1*t1[i])
    

pl.figure(4)
pl.plot(t1, purity)
pl.plot(t1, concurr)
pl.xlabel("Time (s)")
pl.ylabel("Purity and Concurrence")
pl.legend(["Purity", "Concurrence"])

#problem 3

psi0 = qt.basis(2,0)
psi1 = qt.basis(2,1)
psip = (1/np.sqrt(2))*(psi0 + psi1)
psim = (1/np.sqrt(2))*(psi0 - psi1)
pp = psip*psip.dag()
pm = psim*psim.dag()

M0 = psi0*psi0.dag()
M1 = psi1*psi1.dag()

s0 = qt.qeye(2)
sx = qt.sigmax()
sy = qt.sigmay()
sz = qt.sigmaz()

wq = 1
t3 = np.linspace(0,10, 1000)

init_s = psi0

#solve rabi hamiltonian numerically
result = qt.mesolve(H_rabit(t3,wq,sy), init_s, t3, e_ops = [M0, M1])

dynamics = np.array(result.expect)

p0 = dynamics[0]
p1 = dynamics[1]




pl.figure(5)
pl.plot(t3, p0)
pl.plot(t3, p1)
pl.legend(["Ground state", "Excited state"])
pl.xlabel("Time (s)")
pl.ylabel("Population probability")

count = 0
ind1 = 0
ind2 = 0
for i in range(len(t3)):
    if count < 1:
        if abs(p0[i]-.5) <= .0005:
            ind1 = i
            count = count+1

ind2 = ind1*3

res = qt.mesolve(H_rabit(t3,wq,sy), init_s, t3)
print(res.states[ind1], "\r\n ", res.states[ind2])
print(t3[ind1], "\r\n", t3[ind2])

sig = res.states[ind1]*res.states[ind1].dag()

inter = pp.sqrtm()*sig*pp.sqrtm()
inter2 = inter.sqrtm()
fid = inter2.tr()**2
print(fid)


#part b
delta = 0
g = 3
psi_initial = qt.tensor(psip, psim)

sx_A = qt.tensor(sx, qt.qeye(2))
sx_B = qt.tensor(qt.qeye(2), sx)

sy_A = qt.tensor(sy, qt.qeye(2))
sy_B = qt.tensor(qt.qeye(2), sy)

MA1 = qt.tensor(M1, qt.qeye(2))
MB1 = qt.tensor(qt.qeye(2), M1)

sp = sx + 1j * sy
sm = sx - 1j * sy

sp_A = qt.tensor(sp, s0)
sm_A = qt.tensor(sm, s0)
sp_B = qt.tensor(s0, sp)
sm_B = qt.tensor(s0, sm)

mat = ((1/2)*g*(qt.tensor(sx,sx) + qt.tensor(sy,sy)))




t4 = np.linspace(0,3,1000)

#res2 = qt.mesolve(H_coupled(t4, g, sx_A, sx_B, sy_A, sy_B), psi_initial, t4, e_ops=[MA1, MB1])

res2 = qt.mesolve(H_coupled(t4, g, sx_A, sx_B, sy_A, sy_B, sm_A, sp_A, sm_B, sp_B), psi_initial, t4)

ws = qt.Qobj([[1/2], [1j/2], [-1j/2], [-1/2]])
#print(ws)

wi = 0
count = 0
for i in range(len(t4)):
    if abs(res2.states[i][0] - ws[0]) <= .001 and abs(res2.states[i][1] - ws[1]) <= .001 and abs(res2.states[i][2] -ws[2]) <= .001 and abs(res2.states[i][3]-ws[3]) <= .001 and count < 1:
        wi = i
        count = count+1
    #print(res2.states[i][1])

print(t4[wi])
print(res2.states[wi])

sig1 = res2.states[wi]*res2.states[wi].dag()
ppA = sp_A*sp_A.dag()
pmA = sm_A*sm_A.dag()
ppB = sp_B*sp_B.dag()
pmB = sm_B*sm_B.dag()

in1a = ppA.sqrtm()*sig1*ppA.sqrtm()
in2a = in1a.sqrtm()

fidA = in2a.tr()**2

in1b = ppB.sqrtm()*sig1*ppB.sqrtm()
in2b = in1b.sqrtm()
fidB = in2b.tr()**2

print(fidA, " ", fidB)
print(qt.concurrence(res2.states[wi]))

pl.show()
