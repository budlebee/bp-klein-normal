import numpy as np
import matplotlib.pyplot as plt
from copy import copy

# for armchair junctiong

# declare parameters
# dimension default: energy is eV, length is 10^-10
me = 0.51099895000*(10**6) / (299792458**2)  # 전자의 질량 MeV/c**2
hbar = 6.582119569 * 10**(-16)
K0 = 2.35*10**(-9)
vy = 5.6*10**5
m = 1.42*me
mx = 0.15*me
my = 1.18*me
v0 = vy  # hbar*K0 / mx


eg1 = 0.3  # bandgap energy
eg2 = 0.3
potential = 1
energy = 0.5

D = 200*10**(-9)  # potential depth D. set 200nm.
V0 = 25*10**(-3)  # 25meV

E0 = 10*10**(-3)  # 10meV

kxMax = np.sqrt((energy - eg1/2)*2*m/(hbar**2))

kstep = 10000
# kx = np.linspace(kxMax, -1, kstep)  # np.linspace(0, 1, kstep)
# ky = (1/(hbar*vy)) * np.sqrt(energy**2 -
#                             (eg1/2 + ((hbar*kx)**2)/(2*m))**2)  # np.linspace(0, 1, kstep)

#qx = kx
# qy = -(1/(hbar*vy))*np.sqrt((energy-potential)**2 -
#                            (eg2/2 + (((hbar*kx)**2)/(2*m)))**2)  # np.linspace(0, 1, kstep)

inciAngle = []
results = []


lamb = 1  # upper band is +1
lambdot = -1  # lower band is -1


vx = 1
vy = 1
# vy/vx = sqrt(mx/my)
gamma = np.sqrt(mx/my)
alpha = np.arctan(1/gamma)  # potential 장벽이 얼마나 기울어졌느냐. 변수 취급해야하지만 일단 편의상 0으로.
# phiV = np.arctan((gamma**2)*ky/kx)  # 하지만 kx, ky 변화시키지 않고 그냥 phiV 만 변화시킬거임.
phiStep = 1000
phiV = np.linspace(-np.pi/2, np.pi/2, phiStep)
kx = (2*E0*gamma/(lamb*hbar*v0))*np.sqrt(1/(gamma**2 + (np.tan(phiV))**2))
ky = (1/gamma)*np.sqrt(((2*E0)/(lamb*hbar*v0))**2-kx**2)
phiK = np.arctan(ky/kx)
k = np.sqrt(kx**2 + ky**2)
phiKdot = phiK - alpha
phiKRdot = np.arctan((kx*np.sin(alpha) + ky*np.cos(alpha)) /
                     (-kx*np.cos(alpha) + ky*np.sin(alpha)))
kxdot = k*np.cos(phiKdot)
kydot = k*np.sin(phiKdot)
phiS = np.arctan((1/gamma)*np.tan(phiV))
# alpah 값이 달라지면 둘이 다르지만, 일단 alpha 가 0인 상태로 테스트만.
phiSR = np.arctan((gamma)*np.tan(phiKRdot))

# thetaV는? theta 가 포텐셜 장벽 내부의 phi다.

qy = ky
qx = np.sqrt((2*(V0-E0)/(hbar*v0))**2 - (gamma*qy)**2)
qyr = qy
qxr = -qx
q = np.sqrt(qx**2 + qy**2)
thetaK = np.arctan(qy/qx)
thetaKdot = thetaK - alpha
qxdot = q*np.cos(thetaKdot)
qydot = q*np.sin(thetaKdot)
qxrdot = qxr*np.cos(alpha) + qyr*np.sin(alpha)
thetaKRdot = np.arctan((qx*np.sin(alpha) + qy*np.cos(alpha)) /
                       (-qx*np.cos(alpha) + qy*np.sin(alpha)))
thetaV = np.arctan((gamma**2)*qy/qx)  # qy와 qx 를 kx, ky 로 표현하자.
thetaS = np.arctan((1/gamma)*np.tan(thetaV))
thetaSR = np.arctan((gamma)*np.tan(thetaKRdot))


A = np.exp(-1j*qxdot*D)*(np.exp(1j*thetaS + 1j*thetaSR) + np.exp(1j*phiS + 1j*phiSR) -
                         lamb*lambdot*np.exp(1j*thetaSR + 1j*phiSR) - lamb*lambdot*np.exp(1j*thetaS + 1j*phiS)) - np.exp(-1j*qxrdot*D)*(np.exp(1j*thetaS+1j*thetaSR) + np.exp(1j*phiS+1j*phiSR) - lamb*lambdot*np.exp(1j*thetaSR+1j*phiS) - lamb*lambdot*np.exp(1j*thetaS + 1j*phiSR))

t = (lamb*lambdot*np.exp(-1j*kxdot*D)*(np.exp(1j*thetaSR) -
     np.exp(1j*thetaS))*(np.exp(1j*phiS)-np.exp(1j*phiSR)))/A

T = np.absolute(t)**2  # np.absolute(t)**2  # t * np.conjugate(t)
# alpha 값에 의해 kx ky 값이 좀 바뀌고, V와 D를 특정값으로 한뒤 phiV 에 따른 변화를 plot 하기.


#thetaK = np.angle(complex((eg1/2 + hbar**2 * kx**2), (hbar * vy * ky)))
#thetaQ = np.angle(complex((eg2/2 + hbar**2 * qx**2), (hbar * vy * qy)))
# transmission = - ((2*np.sin(thetaK) * np.sin(thetaQ)) /
#                  (1+np.cos(thetaK + thetaQ)))


fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
ax.plot(phiV, T)

plt.show()


# for b in range(0, kstep):
#
#    thetaK = np.angle(
#        complex(((eg2-eg1)/2 + ((hbar * kx[b])**2)/(2*m)), (hbar * vy * ky[b])))
#    # print(thetaK)
#    inciAngle.append(np.angle(
#        complex(((eg2-eg1)/2 + ((hbar * kx[b])**2)/(2*m)), (hbar * vy * ky[b])), deg=True))
#    thetaQ = np.angle(
#        complex(((eg2-eg1)/2 + ((hbar * qx[b])**2)/(2*m)), (hbar * vy * qy[b])))
#    if b == 0:
#        print(kx[b])
#        print(np.angle(
#            complex((eg1/2 + hbar**2 * kx[b]**2), (hbar * vy * ky[b])), deg=True))
#    transmission = - ((2*np.sin(thetaK) * np.sin(thetaQ)) /
#                      (1+np.cos(thetaK + thetaQ)))
#    results.append(copy(transmission))

# print(results)

# plot

#plt.plot(inciAngle, results)

# plt.show()
