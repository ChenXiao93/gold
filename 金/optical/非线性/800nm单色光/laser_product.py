import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op
import OctopusTools.Laser as la

Amp = [10,22,33,40]
Time = [15000, 40000]

for t in Time:
    for a in Amp:

        tt, aa, envy = la.GaussianPulse(PropTime=t+0.05, PulseDuration=int(t * 2 / 5), Wavelength=800, Amplitude=a, TimeStep=0.2,
                                            Phase=0)
        cx.save("./LASER/t%s/a%s" % ( t,a), tt, aa, "### t [hbar/Ha]    A [Ha/e]")


