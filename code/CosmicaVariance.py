import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import math
from scipy.stats import poisson

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
np_lim = 7.258733E10/mp


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
print max(np13)*mp
index = np.where(np13>np_lim)
np13 = np13[index]

Lum13 = []
Mstar13 = []
for i in range(len(np13)):
    Lu, Ms = Luminos(np13[i] ,Hmass, Luminosity, Mstar, mp)
    Lum13.append(Lu[0])
    Mstar13.append(Ms[0])

index2 =  np.where(Mstar13>MW_Mstar)
print len(index2[0])
X13 = snap13[index,17]
Y13 = snap13[index,18]
Z13 = snap13[index,19]
#print np.amin(X13), np.amax(X13)
#print np.amin(Y13), np.amax(Y13)
#print np.amin(Z13), np.amax(Z13)



snap14 = np.loadtxt("../data/Millenium14")
z14 =snap14[0,12]
np14 =snap14[:,13]
print np.amax(np14)*mp
index = np.where(np14>np_lim)
np14 = np14[index]

Lum14 = []
Mstar14 = []
for i in range(len(np14)):
    Lu, Ms = Luminos(np14[i] ,Hmass, Luminosity, Mstar, mp)
    Lum14.append(Lu[0])
    Mstar14.append(Ms[0])

index2 = np.where(Mstar14>MW_Mstar)
print len(index2[0])
X14 = snap14[index,17]
Y14 = snap14[index,18]
Z14 = snap14[index,19]

snap15 = np.loadtxt("../data/Millenium15")
z15 =snap15[0,12]
np15 =snap15[:,13]
print max(np15)*mp
index = np.where(np15>np_lim)
np15 = np15[index]
Lum15 = []
Mstar15 = []
for i in range(len(np15)):
    Lu, Ms = Luminos(np15[i] ,Hmass, Luminosity, Mstar, mp)
    Lum15.append(Lu[0])
    Mstar15.append(Ms[0])

index2 = np.where(Mstar15 > MW_Mstar)
print len(index2[0])
X15 = snap15[index,17]
Y15 = snap15[index,18]
Z15 = snap15[index,19]


#print len(X13[0])
#print len(X14[0])
#print len(X15[0])



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
#rint Lz, zavg    

########################################################################################
# PENCIL BEAM IMPLEMENTATION ---------------------------------------------------------
########################################################################################

def origin(L):
    a1 = np.random.random(1)*L
    b1 = np.random.random(1)*L
    c1 = np.random.random(1)*L
    return a1, b1, c1


def random_cube(L, Lx, Ly, Lz, a1, b1, c1):
    x = []
    y = []
    z = []
    #---------Origin 
    #a1 = np.random.random(1)*L
    #b1 = np.random.random(1)*L
    #c1 = np.random.random(1)*L
    #a1 = 1.0
    #b1 = 1.0
    #c1 = 1.0
    # -------- bottom xedge
    a2 = (a1 + Lx)
    b2 = b1 #
    c2 = c1
    # ---------- right yedge
    a3 = a2 # 
    b3 = (b1 + Ly)
    c3 = c1
    #----------- top xedge
    a4 = a1 #*cos(theta) 
    b4 = (b1 + Ly) #*sin(theta)
    c4 = c1
    # --------- # left yedge
    a5 = a1
    b5 = b1
    c5 = c1 + Lz
    #------------ second face of the triangle
    a6 = (a1 + Lx )#*cos(theta) 
    b6 = b1 #*sin(theta)
    c6 = c1 + Lz
    # ---------- right yedge
    a7 = (a1 + Lx)# *cos(theta) 
    b7 = (b1 + Ly) # *sin(theta)
    c7 = c1 + Lz
    #----------- top xedge
    a8 = a1 #*cos(theta) 
    b8 = (b1 + Ly) #*sin(theta)
    c8 = c1 + Lz

    x.append(a1)
    y.append(b1)
    z.append(c1)
    x.append(a2)
    y.append(b2)
    z.append(c2)
    x.append(a3)
    y.append(b3)
    z.append(c3)
    x.append(a4)
    y.append(b4)
    z.append(c4)
    x.append(a5)
    y.append(b5)
    z.append(c5)
    x.append(a6)
    y.append(b6)
    z.append(c6)
    x.append(a7)
    y.append(b7)
    z.append(c7)
    x.append(a8)
    y.append(b8)
    z.append(c8)

    #fig = plt.figure()
    #ax = fig.gca(projection='3d')
    #ax.set_xlim([0, L])
    #ax.set_ylim([0, L])
    #ax.set_zlim([0, L])
    #ax.scatter(x, y, z)
    #X, Y, Z = data(L)
    #ax.scatter(X, Y, Z, c='r', alpha=0.5)
    #plt.show()
    #ax.set_xlabel("x")
    #ax.set_ylabel("y")
    #ax.set_zlabel("z")
    return a1, a2, a3, a4, a5, a6, a7, a8, b1, b2, b3, b4, b5, b6, b7, b8, c1, c2, c3, c4, c5, c6, c7 ,c8

def pbc(L, Lx, Ly, Lz, ax, ay, az, X13, Y13, Z13, nP):
    a1, a2, a3, a4, a5, a6, a7, a8, b1, b2, b3, b4, b5, b6, b7, b8, c1, c2, c3, c4, c5, c6, c7 ,c8 = random_cube(L, Lx, Ly, Lz, ax, ay, az)
    #a1, a2, a3, a4, a5, a6, a7, a8, b1, b2, b3, b4, b5, b6, b7, b8, c1, c2, c3, c4, c5, c6, c7 ,c8 = random_cube_general(L, theta, phi)
    r = []
    x = [a1, a2, a3, a4, a5, a6, a7, a8, b1, b2, b3, b4, b5, b6, b7, b8, c1, c2, c3, c4, c5, c6, c7 ,c8]
    for i in range(len(x)):
        if x[i] > L:
            y = x[i]-L
            r.append(y)
        else:
            r.append(x[i])
            
     
    X = X13
    Y = Y13
    Z = Z13
    
    a1 = r[0]
    a2 = r[2]
    b1 = r[8]
    b3 = r[10]
    c1 = r[16]
    c6 = r[21]
    
    if ((a1<=a2) & (b1 <=b3) & (c1<=c6)):
        index = np.where((X<a2) & (X>a1) & (Y<b3) & (Y>b1) & (Z<c6) & (Z>c1)) 
    
    elif((a1>a2) & (b1 <=b3) & (c1<=c6)):
        index = np.where(((X<a2) | (X>a1)) & (Y<b3) & (Y>b1) & (Z<c6) & (Z>c1))
        
    elif((a1<=a2) & (b1 >b3) & (c1<=c6)):
        index = np.where((X<a2) & (X>a1) & ((Y<b3) | (Y>b1)) & (Z<c6) & (Z>c1))
    
    elif((a1<=a2) & (b1<=b3) & (c1>c6)):
        index = np.where((X<a2) & (X>a1) & (Y<b3) & (Y>b1) & ((Z<c6) | (Z>c1)))
        
    elif((a1>a2) & (b1>b3) & (c1<=c6)):
        index = np.where(((X<a2) | (X>a1)) & ((Y<b3) | (Y>b1)) & ((Z<c6) & (Z>c1)))
        
    elif((a1>a2) & (b1<=b3) & (c1>c6)):
        index = np.where(((X<a2) | (X>a1)) & ((Y<b3) & (Y>b1)) & ((Z<c6) | (Z>c1)))
    
    elif((a1<=a2) & (b1>b3) & (c1>c6)):
        index = np.where(((X<a2) & (X>a1)) & ((Y<b3) | (Y>b1)) & ((Z<c6) | (Z>c1)))
        
    elif((a1>a2) & (b1>b3) & (c1>c6)):
        index = np.where(((X<a2) | (X>a1)) & ((Y<b3) | (Y>b1)) & ((Z<c6) | (Z>c1)))
            
    X1 = X[index]
    Y1 = Y[index]
    Z1 = Z[index]
    #mass = nP[index]
    #print X1, Y1, Z1, mass
    return len(X1)
    #print len(X13)
    #print X, Y, Z, mass

Nhalos = []


while i<100000:
    x, y, z = origin(L)
    #print "origins", x, y, z


    #rint "--------- Halos-------------"
    n1 = pbc(L, Lx, Ly, Lz[0], x, y, z, X13, Y13, Z13, np13)
    n2 = pbc(L, Lx,  Ly, Lz[1], x, y, z+Lz[0], X14, Y14, Z14, np14)
    n3 = pbc(L, Lx, Ly, Lz[2], x, y, z+Lz[0] + Lz[1], X15, Y15, Z15, np15)
    Nhalos.append(n1 + n2 +n3)
    i+=1
#print Nhalos


nhalos = np.linspace(0, np.amax(Nhalos), np.amax(Nhalos)+1)

print "Nhalos mean", np.mean(Nhalos)
print "Standard deviation", np.std(Nhalos)
Poisson = []
print "Len nhalos", len(nhalos)
l = np.mean(Nhalos)
for i in range(len(nhalos)):
	Poisson.append(poisson.pmf(nhalos[i] , l))
f = open("Nhalo.dat", "w")
for i in range(len(Nhalos)):
    f.write(str(Nhalos[i]))
f.close()


plt.hist(Nhalos, range=(0,np.amax(Nhalos)), bins=(np.amax(Nhalos)+3.0)/9.0, normed=True)
#plt.hist(Nhalos,  bins=(np.amax(Nhalos)+3.0)/9.0, normed=True)
plt.plot(nhalos,Poisson, c='r', linewidth='2.5')
plt.xlabel(r"$\mathrm{N}$", fontsize=25)
plt.ylabel(r"$\mathrm{P(N)}$", fontsize=25)
#plt.text(20, 0.02, r"$\mathrm{L = 1.58^8L_{\odot}}$", fontsize=20)
#plt.text(20, 0.017, r"$\mathrm{M_{halo} = 7.28^{10}M_{\odot}}$", fontsize=20)
plt.legend((r"$Poisson\ distribution$",r"$Millennium\ Simulation$"),  "upper right")
plt.show()


