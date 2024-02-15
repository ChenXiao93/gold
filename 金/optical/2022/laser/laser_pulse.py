import numpy as np
import matplotlib.pyplot as plt
import chenxiao.Basic as cx

def ASineSquarePulse(A0):
    ################# Input parameters###########################
    dt = 0.1  # time step of pulse dt [au]
    t = 800  # Total propagation time t [au]
    T = 800  # Pulse duration T [au], Better T = t
    tau = 400  # Pulse peak position tau [au]
    wavelength = 630  # [nm]
    energy = cx.Wavelength2Energy(wavelength)  # photon energy \hbar \omega[eV]
    phi = 0  # phase \phi
    ################# used constants  ############################
    h = 4.1356676969 * 1e-15  # planck constant eV·s
    hbar = h / (2 * np.pi)  # reduced planck constant
    c = 299792458  # m/s   light speed
    epsilon0 = 8.854187817E-12  # [F/m]  vacuum permitivity
    ##########  Other parameters  ###########################
    time = np.arange(0, t, dt)  # [a.u.]
    omega = cx.Wavelength2Omega(wavelength)  # [1/s]
    At = A0 * np.sin(omega * time * cx.timeau2s + phi) * np.square(np.sin(np.pi * (time) / T))  # T == len(t)
    A_envelope = A0 * np.square(np.sin(np.pi * (time) / T))
    ##############  Add straight line ###################
    time_add = np.arange(t, t + 100, dt)
    At_add = np.zeros_like(time_add)
    time = np.concatenate((time, time_add))
    At = np.concatenate((At, At_add))
    A_envelope = np.concatenate((A_envelope, At_add))
    Et = cx.A2E(time, At)
    ################## Plotting! ##########################################
    fs = 14  # fontsize
    plt.figure(figsize=(12, 4))
    ############ Electric field ########################################
    plt.subplot(1, 2, 1)
    plt.plot(time, At, label="Real Pulse")
    plt.plot(time, A_envelope, label="$sin^2$ Envelope")
    plt.legend(fontsize=fs, loc=3)  # loc = 3 left bottom
    plt.xlabel("Propagation time t [$\hbar/Ha$]", fontsize=fs)
    plt.ylabel("Gauge field A [$Ha$]", fontsize=fs)
    plt.title("dt = %s [au], %i steps" % (dt, len(time)), fontsize=fs + 2)
    plt.subplot(1, 2, 2)
    plt.plot(time, Et)
    plt.xlabel("Propagation time t [$\hbar/Ha$]", fontsize=fs)
    plt.ylabel("Electric field field E [$V/\mathring{A}$]", fontsize=fs)

    ####################   Save  ######################################
    # f = open("A%s" % A0, "w+")
    # f.write("##Time [a.u.]      A [Ha]   \n")
    # for i in range(len(time)):
    #     f.write("%.2f      %.12f   \n" % (time[i], At[i]))
    # f.close()
    ####################################################################
    plt.show()
    return ()  # t[fs], E[V/Angstrom]


ASineSquarePulse(1)