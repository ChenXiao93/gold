import numpy as np
import matplotlib.pyplot as plt

f1 = np.loadtxt("../A10/gauge_field")
f2 = np.loadtxt("../A10/total_current")
plt.subplot(2,1,1)
plt.plot(f1[:,1], f1[:,2])


plt.subplot(2,1,2)
plt.plot(f2[:,1], f2[:,2])

plt.show()