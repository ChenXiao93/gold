import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op
import OctopusTools.Laser as la

##############  Generate monochromatic pulse #################################################################
def GaussianPulse(PropTime, PulseDuration, Wavelength, Amplitude, TimeStep, tail):
    ################# Input parameters###########################
    # PropTime is the total propagation time [hbar/Ha]
    # PulseDuration is the time range when the intensity decrease to 1/e [habr/Ha]
    PeakPosition = PropTime / 2   # Pulse peak position tau [fs]
    omega_si =  cx.Wavelength2Omega(Wavelength)  # Wavelength [nm]  to  Omega[1/s]
    omega_au = omega_si  * cx.timeau2s  #
    energy = cx.Wavelength2Energy(Wavelength) # wavelength [nm] to  Energy [eV]

    time = np.arange(0, PropTime + tail, TimeStep)
    envelope = Amplitude * np.exp(-np.square((time - PeakPosition) / (PulseDuration/2) ))
    gauge_field = envelope * np.cos(omega_au * (time-PeakPosition) )
    return (time, gauge_field, envelope)   # t[hbar/Ha] A [Hartree/e]
#########################################################################
Amp = [1,4,7,10]
wavelength = 300
energy = cx.Wavelength2Energy(wavelength)
omega_au = energy / 27.211386
t = 1500
dt = 0.2
###################
for a in Amp:
    tt, aa, envy = GaussianPulse(PropTime=t+0.2, PulseDuration=int(t * 2 / 5), Wavelength=wavelength,
                                    Amplitude=a, TimeStep=0.2,tail = 100)
    # aa = envy * np.cos(tt * omega_au) * a

    plt.plot(tt,aa)
    # plt.plot(tt, envy)
    # cx.save("./t10000/A%s" % a, tt, aa, "### t [hbar/Ha]    A [Ha/e]")

plt.show()

