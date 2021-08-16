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

eg1 = 0.3  # bandgap energy
eg2 = 0.3
potential = 1
energy = 0.5


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

plt.plot(inciAngle, results)

plt.show()
