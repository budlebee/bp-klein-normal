import numpy as np
import matplotlib.pyplot as plt
from copy import copy

# for armchair junctiong

# declare parameters
# dimension default: energy is eV, length is 10^-10
me = 0.51099895000*(10**6) / (299792458**2)  # 전자의 질량 MeV/c**2
hbar = 6.582119569 * 10**(-16)
vy = 5.6*10**5
m = 1.42*me

mx = 0.15*me
my = 1.18*me

eg1 = 0.3  # bandgap energy
eg2 = 0.3
potential = 1
energy = 0.5

D = 0  # potential depth D

kxMax = np.sqrt((energy - eg1/2)*2*m/(hbar**2))

kstep = 10000
kx = np.linspace(kxMax, -1, kstep)  # np.linspace(0, 1, kstep)
ky = (1/(hbar*vy)) * np.sqrt(energy**2 -
                             (eg1/2 + ((hbar*kx)**2)/(2*m))**2)  # np.linspace(0, 1, kstep)

qx = kx
qy = -(1/(hbar*vy))*np.sqrt((energy-potential)**2 -
                            (eg2/2 + (((hbar*kx)**2)/(2*m)))**2)  # np.linspace(0, 1, kstep)

inciAngle = []
results = []
thetaS = 0
thetaSR = 0
phiS = 0
phiSR = 0
lamb = 1
lambdot = -1

qxr = 0
qyr = 0
vx = 1
vy = 1
# vy/vx = sqrt(mx/my)
gamma = np.sqrt(mx/my)
phiV = np.arctan((gamma**2)*ky/kx)
phiS = np.arctan((1/gamma)*np.tan(phiV))

A = np.exp(-1j*qx*D)*(np.exp(1j*thetaS+1j*thetaSR) + np.exp(1j*phiS+1j*phiSR) -
                      lamb*lambdot*np.exp(1j*thetaSR + 1j*phiSR) - lamb*lambdot*np.exp(1j*thetaS+1j*phiS)) - np.exp(-1j*qxr*D)*(np.exp(1j*thetaS+1j*thetaSR) + np.exp(1j*phiS+1j*phiSR) - lamb*lambdot*np.exp(1j*thetaSR+1j*phiS) - lamb*lambdot*np.exp(1j*thetaS + 1j*phiSR))

t = (lamb*lambdot*np.exp(-1j*kx*D)*(np.exp(1j*thetaSR) -
     np.exp(1j*thetaS))*(np.exp(1j*phiS)-np.exp(1j*phiSR)))/A

T = t * np.conjugate(t)
# alpha 값에 의해 kx ky 값이 좀 바뀌고, V와 D를 특정값으로 한뒤 phiV 에 따른 변화를 plot 하기.


#thetaK = np.angle(complex((eg1/2 + hbar**2 * kx**2), (hbar * vy * ky)))
#thetaQ = np.angle(complex((eg2/2 + hbar**2 * qx**2), (hbar * vy * qy)))
# transmission = - ((2*np.sin(thetaK) * np.sin(thetaQ)) /
#                  (1+np.cos(thetaK + thetaQ)))

for b in range(0, kstep):

    thetaK = np.angle(
        complex(((eg2-eg1)/2 + ((hbar * kx[b])**2)/(2*m)), (hbar * vy * ky[b])))
    # print(thetaK)
    inciAngle.append(np.angle(
        complex(((eg2-eg1)/2 + ((hbar * kx[b])**2)/(2*m)), (hbar * vy * ky[b])), deg=True))
    thetaQ = np.angle(
        complex(((eg2-eg1)/2 + ((hbar * qx[b])**2)/(2*m)), (hbar * vy * qy[b])))
    if b == 0:
        print(kx[b])
        print(np.angle(
            complex((eg1/2 + hbar**2 * kx[b]**2), (hbar * vy * ky[b])), deg=True))
    transmission = - ((2*np.sin(thetaK) * np.sin(thetaQ)) /
                      (1+np.cos(thetaK + thetaQ)))
    results.append(copy(transmission))

# print(results)

# plot

#plt.plot(inciAngle, results)

# plt.show()
