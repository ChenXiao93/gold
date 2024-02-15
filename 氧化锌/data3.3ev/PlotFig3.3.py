import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.Keldysh as k

def GetIntensity(wavelength):
    Intensity = []
    for a in names.get('A' + str(wavelength)):
        file = np.loadtxt("laser%s/A%s" % (wavelength, a))
        E = max(cx.A2E(file[:, 0], file[:, 1]))
        epsilon  = names.get('epsilon' + str(wavelength))
        n = np.sqrt(epsilon)
        I = cx.E2I(E=E, Refractive_index=n)  # intensity in air!
        Intensity.append(I)
    return (Intensity)

def PlotNex(wavelength):
    plt.figure()
    start = 500
    meannex = []
    for a in names.get('A' + str(wavelength)):
        file = np.loadtxt("%snm/A%s/n_ex" % (wavelength, a))
        t = file[:, 1]
        nex = file[:, 2]
        mean_n = np.mean(nex[start:])
        V = 335.5089 * np.power((5.292e-11), 3)  # [b^3] to [m^3]   grep "Cell volume"  slurm*
        time = names.get('FWHM' + str(wavelength)) * cx.timeau2s  # s
        mean_n = mean_n / (V * time)
        meannex.append(mean_n)
        plt.plot(t, nex, label="A = %s" % (a))
    plt.legend()
    return (meannex)


def Keldysh(wavelength):
    epsilon = names.get('epsilon' + str(wavelength))
    n = np.sqrt(epsilon)
    energy = cx.Wavelength2Energy(wavelength)
    Intensity = np.array(GetIntensity(wavelength))
    TDDFT_nex = PlotNex(wavelength)

    keldysh_parameter, Keldysh_nex = [], []

    print("Intensity [W/cm^2], KelydshParameter, Keldysh_nex [m^-3 s^-1],TDDFT_nex [m^-3 s^-1]")
    for i in range(len(Intensity)):
        gamma, W = k.Kelydsh_W(Intensity=Intensity[i] , photon_energy=energy, band_gap=3.28, RefractiveIndex=1)
        # np.sqrt(names.get('epsilon' + str(wavelength)))
        keldysh_parameter.append(gamma)
        Keldysh_nex.append(W)

        print("%.2e,          %.2f,            %.1e ,           %.1e " % (Intensity[i], gamma, W, TDDFT_nex[i]))

    ################  ax1  ##################################
    # n_ex_el follows Keldysh model
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    l1 = ax1.scatter(Intensity, Keldysh_nex, facecolors="none", edgecolors="red", linewidths=1.5, label="Keldysh")
    l2 = ax1.scatter(Intensity, TDDFT_nex, c="b", linewidths=1.5, label="TDDFT")
    # l3 = ax1.loglog(fit_I,fitted_ww, "--", c = "k")
    l4 = ax1.scatter(0.5e12, 1.6e39, label = "exp laser threshold", color = "green")
    ax1.set_xscale("log")
    ax1.set_yscale("log")

    fs = 16
    ax1.legend(handles=[l1, l2, l4], fontsize=fs - 3, loc='upper left')
    ax1.set_ylabel("Excited electrons $[m^{-3} \cdot s^{-1}]$", fontsize=fs)
    # ax1.set_xticks(fontsize = fs)
    ax1.set_xlabel("Intensity $[W/cm^2]$", fontsize=fs)
    # ax1.set_xlim(min(Intensity), max(Intensity))
    #######   ax2  ##########################################
    # upper x axis shows the Keldysh parameter
    ax2 = ax1.twiny()
    ax2.set_xscale("log")
    ax2.invert_xaxis()
    ax2.set_xlim(max(keldysh_parameter), min(keldysh_parameter))
    ax2.set_xlabel("Keldysh parameter $\gamma$", fontsize=fs)
    format_g = np.around(keldysh_parameter, decimals=2)
    ax2.set_xticks(format_g)
    ax2.set_xticklabels(format_g)

    plt.tight_layout()
    return ()


###############################################################################
## Define the amplitude list for two cases.
names = locals()
names['A' + str(800)] = [2, 3, 4, 5, 6, 8, 10, 13, 20, 27, 40, 55, 80]
names['A' + str(3500)] = [10, 15, 20, 25, 30, 35, 45, 60, 85, 120]
names['FWHM' + str(800)] = 500  # atom unit [hbar/Ha]
names['FWHM' + str(3500)] = 1200
names['epsilon' + str(800)] = 3.286 # https://refractiveindex.info/?shelf=main&book=ZnO&page=Bond-o
names['epsilon' + str(3500)] = 10  # https://refractiveindex.info/?shelf=main&book=ZnO&page=Bond-o

WaveLen = 800
Keldysh(wavelength=WaveLen)

plt.show()