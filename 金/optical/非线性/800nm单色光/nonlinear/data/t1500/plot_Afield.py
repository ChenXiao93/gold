import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op
import OctopusTools.Laser as la

amplitude = [10,22,33,40]
fs = 12
for a in amplitude:
    file = np.loadtxt("./A%s/gauge_field" % (a))
    t = file[:, 1]  #  [hbar/Hartree]
    t_fs = file[:,1] * 0.02419 # fs
    A = file[:,2] # A [Ha]
    E = cx.A2E(t,A) # E [V/Angstrom]

    plt.subplot(1,2,1)
    plt.plot(t_fs, A, label = "A=%s" % a)
    plt.xlabel("t [fs]", fontsize=fs)
    plt.ylabel("A [Ha]", fontsize=fs)
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.plot(t_fs, E, label="E=%.2f" % max(E))
    plt.xlabel("t [fs]", fontsize=fs)
    plt.ylabel("E [$V/\mathring{A}$]", fontsize=fs)
    plt.legend()


plt.show()
