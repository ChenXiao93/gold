import numpy as np
import matplotlib.pyplot as plt
import MyPackage.Basic_Function  as my

t,n =my.ReadData("./n_ex",1,2)

plt.plot(t*0.02419,n)
fs = 14


plt.xlabel("Propagation time [fs]", fontsize = fs)
plt.ylabel("Gauge field A [Ha]", fontsize = fs)
plt.title("Incident laser pulse. 3500nm, 0.354eV", fontsize = fs+2)


plt.show()
