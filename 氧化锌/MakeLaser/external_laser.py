import numpy as np
import OctopusTools.Basic as my
import matplotlib.pyplot as plt
import os

def pump_pulse_800(a):

    wavelength = 800  # nm
    omega = 2 * np.pi * my.Light_speed / (wavelength * 1E-9)  # 1/s
    dt= 0.015    # time step magnitude  [hbar/Hartree]
    T = 500
    T1 = 200
    t = np.arange(0,T,dt)   #  [hbar/Hartree] normal propagation time
    t1 = np.arange(T, T+T1, dt)  # add a flat tail
    tt = np.append(t,t1)  # total propagation time

    A = a * np.sin( omega * t * my.timeau2s)  * np.square(np.sin(np.pi *  t /T)) #  * epsilon_800 # omega(1/s) t [hbar/Hartree]
    A1 = np.zeros_like(t1)
    AA = np.append(A,A1)
    my.save("../newdata3.3ev/laser_800nm/A%s" % a, tt,AA, "t [hbar/Ha]  A [Ha]")

    return(tt,AA)

def pump_pulse_3500(a):

    wavelength = 3500  # nm
    omega = 2 * np.pi * my.Light_speed / (wavelength * 1E-9)  # 1/s
    dt= 0.015   # time step magnitude  [hbar/Hartree]
    T = 1200
    T1 = 300
    t = np.arange(0,T,dt)   #  [hbar/Hartree] normal propagation time
    t1 = np.arange(T, T+T1, dt)  # add a flat tail
    tt = np.append(t,t1)  # total propagation time

    A = a * np.sin( omega * t * my.timeau2s)  * np.square(np.sin(np.pi *  t /T)) #  * epsilon_800 # omega(1/s) t [hbar/Hartree]
    A1 = np.zeros_like(t1)
    AA = np.append(A,A1)
    # my.save("../newdata3.3ev/laser_3500nm/A%s" % a, tt,AA, "t [hbar/Ha]  A [Ha]")

    return(tt,AA)
##############################
AMPLITUDE = [10,20,25,30,35,40,50,70,100]


for a in AMPLITUDE:
    # t, A = pump_pulse_800(a)
    t,A = pump_pulse_3500(a)
    plt.plot(t,A)
    plt.xlabel("Propagation time [$\hbar/Ha$]", fontsize=14)
    plt.ylabel("Gauge field A [Ha]", fontsize=14)

plt.show()
