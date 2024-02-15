import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op
import OctopusTools.Laser as la
import OctopusTools.chi3 as cc


amplitude = [40]
wavelength = 800
# Omega = cx.Wavelength2Energy(wavelength)  # eV
Omega = np.arange(1.4, 1.7, 0.01)  # eV
eta = 0.0 # eV
t = 1500

for a in amplitude:
    t_ext, A_ext, A_env = la.GaussianPulse(PropTime=t + 0.05,
                                           PulseDuration=int(t * 2 / 5),
                                           Wavelength=800,
                                           Amplitude=a,
                                           TimeStep=0.2,
                                           Phase=0)
    file = "./A%s/gauge_field" % (a)
    t_ind = cx.ReadOne(file, 1)
    A_ind = cx.ReadOne(file, 2)
    # 规范场方法
    eps1,eps2 = op.Get_Epsilon_Gauge(t_ext, A_ext, t_ind, A_ind, Omega, eta=0.05)
    # 电流方法
    current = cx.ReadOne("./A%s/total_current" % a, 2)
    epsj1,epsj2 = op.Get_Epsilon_Current_damping(t_ext, A_ext, current, Omega, eta=0.05)


    plt.subplot(1,2,1)
    plt.plot(Omega, eps1)
    plt.plot(Omega, epsj1)

    plt.subplot(1,2,2)
    plt.plot(Omega, eps2)
    plt.plot(Omega, epsj2)

plt.show()