import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.Keldysh as k

def GetIntensity(wavelength):
    Intensity = []
    for a in names.get('A' + str(wavelength)):
        file = np.loadtxt("./%s/laser" % (a))
        print("A = %.2f" % max(file[:,2]))
        E = max(cx.A2E(file[:, 1], file[:,2]))
        E = E / 1  #  GaugeFieldPropagate = no, so no epsilon modification!
        I = cx.E2I(E=E, Refractive_index=1)  # intensity in air!
        Intensity.append(I)
    return (Intensity)

def PlotNex(wavelength):
    plt.figure()
    global start, end
    meannex = []
    for a in names.get('A' + str(wavelength)):
        file = np.loadtxt("./%s/n_ex" % ( a))
        t = file[:, 1]
        nex = file[:, 2]
        mean_n = np.max(nex[start:end])  # use minimum
        t1 = t[start:end]
        nex1 = np.ones_like(t1) * mean_n
        V = 335.5089 * np.power((5.292e-11), 3)  # [b^3] to [m^3]   grep "Cell volume"  slurm*
        time = names.get('FWHM' + str(wavelength)) * cx.timeau2s  # s
        mean_n = mean_n / (V * time)
        meannex.append(mean_n)
        plt.plot(t, nex, label="A = %s" % (a))
        plt.plot(t1, nex1, "k--")

    plt.legend()
    return (meannex)


def Keldysh(wavelength):
    f = np.loadtxt("../zno_data.csv", delimiter=",")

    epsilon = names.get('epsilon' + str(wavelength))
    n = np.sqrt(epsilon)
    energy = cx.Wavelength2Energy(wavelength)
    Intensity = np.array(GetIntensity(wavelength))
    TDDFT_nex = PlotNex(wavelength)

    keldysh_parameter, Keldysh_nex = [], []

    print("Intensity [W/cm^2], KelydshParameter, Keldysh_nex [m^-3 s^-1],TDDFT_nex [m^-3 s^-1]")
    for i in range(len(Intensity)):
        gamma, W = k.Kelydsh_W(Intensity=Intensity[i] , photon_energy=energy, band_gap=2.79, RefractiveIndex=1)
        # np.sqrt(names.get('epsilon' + str(wavelength)))
        keldysh_parameter.append(gamma)
        Keldysh_nex.append(W)

        print("%.2e,          %.2f,            %.1e ,           %.1e " % (Intensity[i], gamma, W, TDDFT_nex[i]))

    #####################  Linear fitting ####################
    # USE THE FIRST 4 DATA OF OCTOPUS DO THE LINEAR FITTING.
    def linearfitting(I, n, I1):
        linear_fitting_real = np.polyfit(np.log10(I), np.log10(n), 1)
        parameter_real = np.poly1d(linear_fitting_real)
        a2 = parameter_real[1]  # k1
        a1 = parameter_real[0]  # b
        fit_I = np.array(I1)
        fitted_n = np.power(fit_I, a2) * np.power(10, a1)

        return (fit_I, fitted_n)
    fitted_I, fitted_n = linearfitting(I=Intensity[0:3], n=Keldysh_nex[0:3], I1=Intensity[0:5])
    ################  ax1  ##################################
    # n_ex_el follows Keldysh model
    fig = plt.figure(figsize=(6, 5))
    ax1 = fig.add_subplot(111)
    l1 = ax1.scatter(Intensity, Keldysh_nex, facecolors="none", edgecolors="red", linewidths=1.5, label="Keldysh")
    l2 = ax1.scatter(Intensity, TDDFT_nex, c="b", linewidths=1.5, label="TDDFT")
    l3 = ax1.loglog(fitted_I,fitted_n, "--", c = "k")
    l4 = ax1.scatter(0.6e12, 1.6e39, label="exp laser threshold", color="green",  marker="*", linewidths=3)
    l5 = ax1.scatter(f[:, 2], f[:, 3], facecolors="none",edgecolors="orange", linewidths=1.5, label="SBE")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    fs = 16
    ax1.legend(handles=[l1, l2, l4, l5], fontsize=fs - 3, loc='lower right')
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
    plt.savefig("nex3500.eps", dpi=900, bbox_inches='tight')
    return ()

###############################################################################
names = locals()
names['A' + str(3500)] = ["1e11","3e11","6e11","1e12","3e12","6e12","1e13"]
names['FWHM' + str(3500)] = 500
names['epsilon' + str(3500)] = 3.1 # https://refractiveindex.info/?shelf=main&book=ZnO&page=Bond-o
start = int(1200/0.80)
end = int(1400/0.80)
WaveLen = 3500
Keldysh(wavelength=WaveLen)

plt.show()