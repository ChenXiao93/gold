import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op
import OctopusTools.Laser as la

dt = 0.2
# dt = 0.05


plt.figure()
fs = 12
Wavelength = 800 # nm
Omega = cx.Wavelength2Energy(Wavelength)  # eV
eta = 0.0015 # [Ha/hbar]
t_ind = cx.ReadOne("./dt%s/gauge_field_dt%s" % (dt, dt), 1)
A_ind= cx.ReadOne("./dt%s/gauge_field_dt%s" % (dt, dt), 2)
t_ext, A_ext, A_env = la.GaussianPulse(PropTime=500 + 0.05,
                                PulseDuration=int(500 * 2 / 5),
                                Wavelength=800,
                                Amplitude=1,
                                TimeStep=dt,
                                Phase=0)

eps1,eps2 = op.Get_Epsilon_Gauge(t_ext, A_ext, t_ind, A_ind, Omega,eta)

print("%.4f    %.4f"   % (eps1, eps2) )

# plt.plot(t_ext, A_ext, "r-")
plt.plot(t_ind, A_ind, "-", label = "$A_{ind}$")
plt.plot(t_ind, np.exp(-eta * t_ind), "-",label = "$e^{-\eta t}$")
plt.plot(t_ind, A_ind * np.exp(-eta * t_ind), "--", label = "$A_{ind}e^{-\eta t} $")
# plt.plot(t_ext,A_env, "--")
plt.legend(fontsize=  fs)

plt.xlabel("t [$\hbar/Ha$]", fontsize = fs)
plt.ylabel("$A_{ind}$ [Ha]", fontsize = fs)

plt.show()
