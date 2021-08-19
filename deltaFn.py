import numpy as np
import matplotlib.pyplot as plt

# define variables
#kx = 1
#ky = 1
hbar = 4.135667696 * (10**(-15))  # plank const eV*s.
vy = 5.6*10**5  # armchair direction velocity factor. in bp, vx != vy.
V0 = 0.5*10**5  # 1  # potential factor of delta Fn.
Eg = 0.3  # energy gap in bp.
me = 0.51099895000*(10**6) / (299792458**2)  # 전자의 질량 MeV/c**2
m = 1.42*me
mx = 0.15*me
my = 1.18*me
gamma = np.sqrt(mx/my)  # 실제 공간상의 속도가 xy 가 서로 다르기에 factor 로 보정해줘야됨.
# my = 1  # effective mass of particle
#phiK = np.arctan(ky/kx)
#E = np.sqrt((Eg/2 + hbar**2*kx**2/(2*m))**2 + (hbar*vy*ky)**2)
E = 0.5  # 0.5eV # energy of electron
# phiV = np.arctan(gamma**2 * ky/kx)  # real space incident angle.

#phiStep = 100
#phiV = np.linspace(-np.pi/2, np.pi/2, phiStep)
#a = (hbar**4)/(4*m**2)
#c = (Eg**2)/4 - E**2
#b = (Eg*hbar**2)/(2*m)+(hbar*vy*np.tan(phiV)/gamma)**2

# kx = np.sqrt((-b + np.sqrt(b**2 - 4*a*c))/(2*a))  # kx representation of phiV

kxMax = np.sqrt((2*m/hbar**2) * (E-Eg/2))
print(kxMax)
kx = np.linspace(1/1000, kxMax, 1000)
ky = (1/(hbar*vy)) * np.sqrt(E**2 -
                             (Eg/2 + ((hbar*kx)**2)/(2*m))**2)
phiV = np.arctan(gamma**2 * ky/kx)
Trx = kx**2 / (kx**2 + (m*V0/hbar)**2)
Try = ky**2 / (ky**2 + (m*V0/hbar)**2)

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
plt.plot(phiV, np.sqrt(Trx**2 + Try**2)/np.sqrt(2))
plt.show()
