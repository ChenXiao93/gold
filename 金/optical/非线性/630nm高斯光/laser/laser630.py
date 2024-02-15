import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
##############  Generate monochromatic pulse #################################################################
def GaussianPulse(PropTime, PulseDuration, Wavelength, Amplitude, TimeStep, tail):
    ################# Input parameters###########################
    # PropTime is the total propagation time [hbar/Ha]
    # PulseDuration is the time range when the intensity decrease to 1/e [habr/Ha]
    PeakPosition = PropTime / 2   # Pulse peak position tau [fs]
    omega_si =  cx.Wavelength2Omega(Wavelength)  # Wavelength [nm]  to  Omega[1/s]
    omega_au = omega_si  * cx.timeau2s  #

    time = np.arange(0, PropTime + tail, TimeStep)
    envelope = Amplitude * np.exp(-np.square((time - PeakPosition) / (PulseDuration/2)))
    gauge_field = envelope * np.cos(omega_au * (time-PeakPosition) )
    return (time, gauge_field, envelope)   # t[hbar/Ha] A [Hartree/e]
#########################################################################
Amp = [1,4,7,10]
wavelength = 630
energy = cx.Wavelength2Energy(wavelength)
omega_au = energy / 27.2114
t = 1000
dt = 0.05
for a in Amp:
    tt, aa, envy = GaussianPulse(PropTime=t+dt, PulseDuration=int(t * 2/5), Wavelength=wavelength,
                                    Amplitude=a, TimeStep=dt,tail = 0)
    plt.plot(tt,aa)
    # plt.plot(tt, envy)
    # cx.save("./t%s/A%s" % (t, a), tt, aa,"### t [hbar/Ha]    A [Ha/e]")

plt.show()

