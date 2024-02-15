import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.Keldysh as k

def GetEpsilon(A):

    time_laser, A_laser = cx.ReadData("./laser_800nm/A%s" % A, 0, 1)
    time_gauge, A_gauge = cx.ReadData("./A%s/gauge_field" % A, 1, 4)
    Amp = max(A_laser)
    omega = cx.Wavelength2Omega(800) * cx.timeau2s  # [1/s] --> 1 / [hbar/Ha]
    Eta = 0.02  # [Ha/hbar]
    dAdt_laser = cx.dA_over_dt(time_laser, A_laser)
    dAdt_gauge = cx.dA_over_dt(time_gauge, A_gauge)
    total_real, total_imag = cx.FT_eta(time_gauge, dAdt_gauge, Eta, omega)
    laser_real, laser_imag = cx.FT_eta(time_laser, dAdt_laser, Eta, omega)
    Denominator = np.square(total_real + laser_real) + np.square(total_imag + laser_imag)
    eps1 = (laser_real * (total_real + laser_real) + laser_imag * (total_imag + laser_imag)) / Denominator
    eps2 = (total_real * laser_imag - laser_real * total_imag) / Denominator
    print("A = %.2f,  eps1 = %.2f, eps2 = %.2f"  % (Amp, eps1, eps2))

    return(eps1,eps2)

def GetIntensity(wavelength):

    Intensity = []
    for a in names.get('A' + str(wavelength)):
        eps1,eps2 = GetEpsilon(a)
        file = np.loadtxt("./laser_800nm/A%s" % (a))
        E = max(cx.A2E(file[:, 0], file[:, 1]))
        E = E / np.sqrt(eps1**2+eps2**2)
        # E = E / eps1
        I = cx.E2I(E=E, Refractive_index=1)  # intensity in air!
        Intensity.append(I)

    return (Intensity)

def PlotNex(wavelength):

    global start, end

    plt.figure()
    meannex = []
    for a in names.get('A' + str(wavelength)):
        file = np.loadtxt("./A%s/n_ex" % ( a))
        t = file[:, 1]
        nex = file[:, 2]
        mean_n = np.min(nex[start:end])  # use minimum
        t1 = t[start:end]
        nex1 = np.ones_like(t1) * mean_n
        V = 335.5089 * np.power((5.292e-11), 3)  # [b^3] to [m^3]   grep "Cell volume"  slurm*
        time = names.get('FWHM' + str(wavelength)) * cx.timeau2s  # s
        mean_n = mean_n / (V * time)
        meannex.append(mean_n)
        plt.plot(t, nex, label="A = %s" % (a))
        plt.plot(t1, nex1, "k--")
        print("A = %s, n = %.2e" % (a, mean_n))

    plt.legend()

    return (meannex)


#####################  Linear fitting ####################
# USE THE FIRST 4 DATA OF OCTOPUS DO THE LINEAR FITTING.
def linearfitting(I,n,I1):

    linear_fitting_real = np.polyfit(np.log10(I), np.log10(n), 1)
    parameter_real = np.poly1d(linear_fitting_real)
    a2 = parameter_real[1]  # k1
    a1 = parameter_real[0]  # b
    fit_I = np.array(I1)
    fitted_n = np.power(fit_I, a2) * np.power(10, a1)

    return(fit_I, fitted_n)




def Keldysh(wavelength):

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

    ########  I, n is used to do linear fitting,  I1 is the intensity range that you want to get!!!
    fitted_I, fitted_n = linearfitting(I = Intensity[0:6], n =  Keldysh_nex[0:6], I1 = Intensity[0:8])
    ################  ax1  ##################################
    # n_ex_el follows Keldysh model
    fig = plt.figure(figsize=(6, 5))
    ax1 = fig.add_subplot(111)
    l1 = ax1.scatter(Intensity, Keldysh_nex, facecolors="none", edgecolors="red", linewidths=1.5, label="Keldysh")
    l2 = ax1.scatter(Intensity, TDDFT_nex, c="b", linewidths=1.5, label="TDDFT")
    l3 = ax1.loglog(fitted_I,fitted_n, "--", c = "k")
    l4 = ax1.scatter(0.5e12, 1.6e39, label="exp laser threshold", color="green",  marker="*", linewidths=3)
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

    plt.savefig("nex.png", dpi = 900, bbox_inches='tight')
    return ()

#################Define the amplitude list for two cases. ####################
names = locals()
names['A' + str(800)] = [10,20,25,30,35,50,70,100]
names['FWHM' + str(800)] = 300# atom unit [hbar/Ha]   1/e max
start = int(400/0.75)
end = int(450/0.75)
# ##################   Main Program    ######################################
WaveLen = 800
# PlotNex(wavelength=WaveLen)
Keldysh(wavelength=WaveLen)

plt.show()