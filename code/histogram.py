import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data.dat")

x = data[:,0]

print x

plt.hist(x, bins=40)
plt.show()
