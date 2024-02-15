import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx

def pump_pulse_800(a):

    wavelength = 800  # nm
    omega = 2 * np.pi * cx.Light_speed / (wavelength * 1E-9)  # 1/s
    dt= 0.015    # time step magnitude  [hbar/Hartree]
    T = 500
    T1 = 200
    t = np.arange(0,T,dt)   #  [hbar/Hartree] normal propagation time
    t1 = np.arange(T, T+T1, dt)  # add a flat tail
    tt = np.append(t,t1)  # total propagation time

    A = a * np.sin( omega * t * cx.timeau2s)  * np.square(np.sin(np.pi *  t /T)) #  * epsilon_800 # omega(1/s) t [hbar/Hartree]
    A1 = np.zeros_like(t1)
    AA = np.append(A,A1)
    return(tt,AA)

##############################
wavelength = 800 # nm
energy = cx.Wavelength2Energy(wavelength) # eV
omega_au = energy / cx.hartree2ev  # a.u.

for i in Intenisty:
    Iau = np.sqrt(i) / np.sqrt(3.509470 * 1E16)
    a = Iau/omega_au*cx.c_au
    # a_modified = a * eps1_3500
    t,A = pump_pulse_3500(a)

    plt.plot(t,A)
    plt.xlabel("Propagation time [$\hbar/Ha$]", fontsize=14)
    plt.ylabel("Gauge field A [Ha]", fontsize=14)
    # cx.save("./A%i" % a, t, A, "t [hbar/Ha]  A [Ha]")
plt.show()

