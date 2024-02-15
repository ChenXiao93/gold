import numpy as np
import matplotlib.pyplot as plt

f = np.loadtxt("gauge_field")

plt.plot(f[:,1],f[:,2])

plt.show()