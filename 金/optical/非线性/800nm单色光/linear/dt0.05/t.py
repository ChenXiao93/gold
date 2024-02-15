import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op
import OctopusTools.Laser as la


t_ext, A_ext, A_env = la.GaussianPulse(PropTime=500 + 0.05,
                                PulseDuration=int(500 * 2 / 5),
                                Wavelength=800,
                                Amplitude=1,
                                TimeStep=0.2,
                                Phase=0)

plt.plot(t_ext, A_ext, label = "$A_{ext}$")
plt.plot(t_ext, A_env, "--",label = "Envelope")
fs = 12
plt.legend(fontsize = fs)
plt.xlabel("t [$\hbar/Ha$]", fontsize = fs)
plt.ylabel("$A_{ext}$ [Ha]", fontsize = fs)

plt.show()