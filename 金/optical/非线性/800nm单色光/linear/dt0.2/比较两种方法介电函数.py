import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op
import OctopusTools.Laser as la

dt = 0.2

plt.figure()
fs = 12
Wavelength = 800 # nm
# Omega = cx.Wavelength2Energy(Wavelength)  # eV
Omega = np.arange(0.5,3,0.01)  # eV
t_ext, A_ext, A_env = la.GaussianPulse(PropTime=500 + 0.05,
                                PulseDuration=int(500 * 2 / 5),
                                Wavelength=800,
                                Amplitude=1,
                                TimeStep=dt,
                                Phase=0)

#  方法一： 规范场
eta = 0.00 # [Ha/hbar]
t_ind = cx.ReadOne("./gauge_field_dt%s" % (dt), 1)
A_ind= cx.ReadOne("./gauge_field_dt%s" % (dt), 2)
eps1,eps2 = op.Get_Epsilon_Gauge(t_ext, A_ext, t_ind, A_ind, Omega,eta)
# 方法二： 电流法
t_cur = cx.ReadOne("./total_current_dt%s" % (dt), 1)
j_cur= cx.ReadOne("./total_current_dt%s" % (dt), 2)

eps1_j,eps2_j = op.Get_Epsilon_Current_damping(t_ext, A_ext, j_cur, Omega, eta)



#  比较画图
plt.subplot(1,2,1)
plt.plot(Omega, eps1, label = "A")
plt.plot(Omega, eps1_j, label = "J")
plt.legend()
plt.subplot(1,2,2)
plt.plot(Omega,eps2, label = "A")
plt.plot(Omega, eps2_j, label = "J")
plt.legend()


plt.show()
