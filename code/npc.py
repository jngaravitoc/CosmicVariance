# This code is made by Juan Nicolas Garavito and Azadeh Fattahi under 
# The supervision of Michelle Trenti.

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys


def data(filename):
	data = np.loadtxt(filename)
	x = data[:,0]
	y = data[:,1]
	z = data[:,2]
	return x, y, z

def random_cube(L, Lz, Lx):
    x = []
    y = []
    z = []
    #---------Origin 
    a1 = np.random.random(1)*L
    b1 = np.random.random(1)*L
    c1 = np.random.random(1)*L
    #a1 = 1.0
    #b1 = 1.0
    #c1 = 1.0
    # -------- bottom xedge
    a2 = (a1 + Lx)
    b2 = b1 #*sin(theta)
    c2 = c1
    # ---------- right yedge
    a3 = a2 #*cos(theta) 
    b3 = (b1 + Lx)
    c3 = c1
    #----------- top xedge
    a4 = a1 #*cos(theta) 
    b4 = (b1 + Lx) #*sin(theta)
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
    b7 = (b1 + Lx) # *sin(theta)
    c7 = c1 + Lz
    #----------- top xedge
    a8 = a1 #*cos(theta) 
    b8 = (b1 + Lx) #*sin(theta)
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
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.set_xlim([0, L])
    ax.set_ylim([0, L])
    ax.set_zlim([0, L])
    ax.scatter(x, y, z)
    X, Y, Z = data(L)
    #ax.scatter(X, Y, Z, c='r', alpha=0.5)
    plt.show()
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    return a1, a2, a3, a4, a5, a6, a7, a8, b1, b2, b3, b4, b5, b6, b7, b8, c1, c2, c3, c4, c5, c6, c7 ,c8

# Our pencil beam has the same distance in the x and y axys.

def pbc(L, Lz,Lx):
    a1, a2, a3, a4, a5, a6, a7, a8, b1, b2, b3, b4, b5, b6, b7, b8, c1, c2, c3, c4, c5, c6, c7 ,c8 = random_cube(L, Lz, Lx)
    #a1, a2, a3, a4, a5, a6, a7, a8, b1, b2, b3, b4, b5, b6, b7, b8, c1, c2, c3, c4, c5, c6, c7 ,c8 = random_cube_general(L, theta, phi)
    r = []
    x = [a1, a2, a3, a4, a5, a6, a7, a8, b1, b2, b3, b4, b5, b6, b7, b8, c1, c2, c3, c4, c5, c6, c7 ,c8]
    for i in range(len(x)):
        if x[i] > L:
            y = x[i]-L
            r.append(y)
        else:
            r.append(x[i])
            
     
    X, Y, Z = data(L)
    
    if ((a1<=a2) & (b1 <=b3) & (c1<=c6)):
        index =np.where((X<a2) & (X>a1) & (Y<b3) & (Y>b1) & (Z<c6) & (Z>c1)) 
    
    elif((a1>a2) & (b1 <=b3) & (c1<=c6)):
        index = np.where(((X>a2) | (X<a1)) & (Y<b3) & (Y>b1) & (Z<c6) & (Z>c1))
      
    elif((a1<=a2) & (b1 >b3) & (c1<=c6)):
        index = np.where((X<a2) & (X>a1) & ((Y<b3) | (Y>b1)) & (Z<c6) & (Z>c1))
        
    elif((a1<=a2) & (b1<=b3) & (c1>c6)):
        index = np.where((X<a2) & (X>a1) & (Y>b3) & (Y<b1) & ((Z>c6) | (Z<c1)))
        
    elif((a1>a2) & (b1>b3) & (c1<=c6)):
        index = np.where(((X>a2) | (X<a1)) & ((Y<b3) | (Y>b1)) & ((Z<c6) & (Z<c1)))
        
    elif((a1>a2) & (b1>b3) & (c1>c6)):
        index = np.where(((X>a2) | (X<a1)) & ((Y<b3) | (Y>b1)) & ((Z>c6) | (Z<c1)))
                
          
    X1 = X[index]
    Y1 = Y[index]
    Z1 = Z[index]
    
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.scatter(X, Y, Z, c='r', alpha=0.5)
    ax.scatter(X1, Y1, Z1, c='green')
    ax.scatter(r[0:8], r[8:16], r[16:24])
    #ax.scatter(X, Y, c='r',alpha=0.5)
    #ax.scatter(X1, Y1, c='green')
    #ax.scatter(r[0:8],r[8:16])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    #ax.plot(a1, a2)
    ax.set_xlim([0, L])
    ax.set_ylim([0, L])
    ax.set_zlim([0, L])
    return X1, Y1, Z1

filename = sys.argv[1]
pbc(10.0, 8, 8)



