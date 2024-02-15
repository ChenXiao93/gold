import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op
import OctopusTools.Laser as la

Amp = [10,20,30,40]

wavelength = 300
energy = cx.Wavelength2Energy(wavelength)
omega_au = energy / 27.2114

for a in Amp:
    time, envy = la.env_sin2_full(head=500, mid=9000, tail=500, dt=0.2)
    aa = envy * np.cos( time * omega_au) * a
    # tt, aa, envy = la.GaussianPulse(PropTime=t+0.05, PulseDuration=int(t * 2 / 5), Wavelength=300, Amplitude=a, TimeStep=0.2,
    #                                     Phase=0)
    plt.plot(time,aa)


    cx.save("./A%s" % a, time, aa, "### t [hbar/Ha]    A [Ha/e]")
plt.show()

