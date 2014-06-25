import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import math

################################################################################################
# This code was written by Azadeth and Juan N. Garavito-Camargo as part
# of the VOSS14 student project. Advisor: Michelle Trenti
################################################################################################

# Defining some constants ------------------------------------------------------

Dh = 3000
thetax = 10.0
thetay = 12.0
h = 0.73
mp = 8.6E8/h 
zi = 8.0
zf = 10.0
omega_m = 0.25
L = 500
MW_Mstar = 5E10 

# Reading mass to light information ------------------------------------------------

Lum = np.loadtxt("../data/mass_to_light_z9.txt")

Hmass = Lum[:,0]
Luminosity = Lum[:,1]
Mstar = Lum[:,2]

#---------------------------------------------------------------------------------

# Function that tranforms mass halo to luminosity --------------------------

def Luminos(nP, Hmass, Lumi, Mstar, mp):
    M_h = nP*mp
    minimum = abs(Hmass - M_h)
    index = np.where(minimum==np.amin(minimum))
    L = Lumi[index]
    M = Mstar[index]
    return L, M

# Reading data --------------------------------------------------------------------------

snap13 = np.loadtxt("../data/Millenium13")
z13 =snap13[0,12]
np13 =snap13[:,13]
index = np.where(np13>6000)
np13 = np13[index]
"""
Lum13 = []
Mstar13 = []
for i in range(len(np13)):
    Lu, Ms = Luminos(np13[i] ,Hmass, Luminosity, Mstar, mp)
    Lum13.append(Lu[0])
    Mstar13.append(Ms[0])

index = np.where(Mstar13>MW_Mstar)
"""
X13 = snap13[index,17]
Y13 = snap13[index,18]
Z13 = snap13[index,19]
#print np.amin(X13), np.amax(X13)
#print np.amin(Y13), np.amax(Y13)
#print np.amin(Z13), np.amax(Z13)



snap14 = np.loadtxt("../data/Millenium14")
z14 =snap14[0,12]
np14 =snap14[:,13]
index = np.where(np14>6000)
np14 = np14[index]
"""
Lum14 = []
Mstar14 = []

for i in range(len(np14)):
    Lu, Ms = Luminos(np14[i] ,Hmass, Luminosity, Mstar, mp)
    Lum14.append(Lu[0])
    Mstar14.append(Ms[0])

index = np.where(Mstar14>MW_Mstar)
"""
X14 = snap14[index,17]
Y14 = snap14[index,18]
Z14 = snap14[index,19]

snap15 = np.loadtxt("../data/Millenium15")
z15 =snap15[0,12]
np15 =snap15[:,13]
index = np.where(np15>6000)
np15 = np15[index]
"""
Lum15 = []
Mstar15 = []
for i in range(len(np15)):
    Lu, Ms = Luminos(np15[i] ,Hmass, Luminosity, Mstar, mp)
    Lum15.append(Lu[0])
    Mstar15.append(Ms[0])

index = np.where(Mstar15 > MW_Mstar)
"""
X15 = snap15[index,17]
Y15 = snap15[index,18]
Z15 = snap15[index,19]


#Luminos(np13[0] ,Hmass, Luminosity, Mstar, mp)

#lt.scatter(np13*mp, Mstar13, c='red')
#lt.xscale("log")
#plt.yscale("log")
#plt.show()

# This is the integral we have to do to compute de comovil distance
def my_func(z, omega_m):
    return 1 / np.sqrt(omega_m*(1+z)**3 + (1-omega_m))

Dc,err = quad(my_func, 0, np.mean(zi + zf) ,omega_m)
Dc = Dc*Dh

# Computing Diameter angular distance

def DAD(theta, z, Dc):
    theta =(theta/60.0)*np.pi/180.0
    L = theta * Dc
    return L

Lx = DAD(thetax, (zi + zf)*0.5 , Dc)
Ly = DAD(thetay, (zi + zf)*0.5 , Dc)
print Lx
#zsnap = [Z13, Z14, Z15]
zsnap = [9.27, 8.54, 7.88]

# This function order the snapshots

def snapZ(zsnap, zi, zf):
    x = []
    zsnap = np.sort(zsnap)
    for i in range(len(zsnap)-1):
        x.append((zsnap[i+1]+zsnap[i])/2.0)
    x.insert(0, zi)
    x.append(zf)
    return x

zavg = snapZ(zsnap, zi, zf)
Lz = []

# Compute the comovil between snapshots

for i in range(len(zavg)-1):
    Dc,err = quad(my_func, zavg[i], zavg[i+1], omega_m)
    Dc = Dc*Dh
    Lz.append(Dc)
