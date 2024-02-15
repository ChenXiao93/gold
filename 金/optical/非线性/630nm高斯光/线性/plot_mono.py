import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op

fs = 12
plt.figure()

Wavelength = 630 # nm
Omega = cx.Wavelength2Energy(Wavelength)  # eV
eta = 0.1# eV
t_ind = cx.ReadOne("./pbe", 1)
A_ind= cx.ReadOne("./pbe", 2)
t_ext = cx.ReadOne("./A1", 0)
A_ext = cx.ReadOne("./A1", 1)

eps1,eps2 = op.Get_Epsilon_Gauge(t_ext, A_ext, t_ind, A_ind, Omega,eta)

print(eps1, eps2)

plt.plot(t_ext, A_ext, "r-")
plt.plot(t_ind, A_ind, "b--")

plt.show()
