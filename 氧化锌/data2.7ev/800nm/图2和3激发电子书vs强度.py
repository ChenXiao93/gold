import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.Keldysh as k
from matplotlib.ticker import FuncFormatter
fs = 16
fig3 = plt.figure(figsize=(6, 4))  ## 不同n_ex
ax3 = fig3.add_subplot(111)
fig1 = plt.figure(figsize=(6, 5))  ## 主要图
ax1 = fig1.add_subplot(111)
fig4 = plt.figure(figsize=(6, 5))  ## 主要图放大
ax4 = fig4.add_subplot(111)

def GetEpsilon(A):

    time_laser, A_laser = cx.ReadData("./laser/A%s" % A, 0, 1)
    time_gauge, A_gauge = cx.ReadData("./gauge_field/gaugefield_A%s" % A, 1, 4)
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
        file = np.loadtxt("./laser/A%s" % (a))
        E = max(cx.A2E(file[:, 0], file[:, 1]))
        E = E / np.sqrt(eps1**2)
        I = cx.E2I(E=E, Refractive_index=1)  # intensity in air!
        Intensity.append(I)
    return (Intensity)

def PlotNex(wavelength):

    start = names.get('start' + str(wavelength))
    meannex = []
    for a in names.get('A' + str(wavelength)):
        file = np.loadtxt("./n_ex/n_A%s" % a)
        t = file[:, 1]
        nex = file[:, 2]
        mean_n = np.min(nex[start:-1])  # use minimum
        t1 = t[start:-1]
        nex1 = np.ones_like(t1) * mean_n
        V = 335.5089 * np.power((5.292e-11), 3)  # [b^3] to [m^3]   grep "Cell volume"  slurm*
        time = names.get('FWHM' + str(wavelength)) * cx.timeau2s  # s
        mean_n = mean_n / (V * time)  ### 平均时间用的是半高全宽！！！！
        meannex.append(mean_n)

        eps1,eps2 = GetEpsilon(a)
        file = np.loadtxt("./laser/A%s" % (a))
        E = max(cx.A2E(file[:, 0], file[:, 1]))
        E = E / np.sqrt(eps1**2)
        I = cx.E2I(E=E, Refractive_index=1) *1E-12 # [TW/cm^2]]
        ## 不同A的n_ex
        ax3.plot(t* 0.02419, nex, label="I = %.2f" % I)
        ax3.plot(t1* 0.02419, nex1, "k--")
        ax3.set_xlabel("t [fs]", fontsize = fs)
        ax3.set_ylabel("$n_{ex} [1/cell]$", fontsize = fs)

    return (meannex)


#####################  Linear fitting ####################
# USE THE FIRST 6 DATA OF OCTOPUS DO THE LINEAR FITTING.
def linearfitting(I,n,I1):

    linear_fitting_real = np.polyfit(np.log10(I), np.log10(n), 1)
    parameter_real = np.poly1d(linear_fitting_real)
    a2 = parameter_real[1]  # k1
    a1 = parameter_real[0]  # b
    fit_I = np.array(I1)
    fitted_n = np.power(fit_I, a2) * np.power(10, a1)

    return(fit_I, fitted_n)

def Keldysh(wavelength):

    f = np.loadtxt("../zno_data.csv", delimiter=",")
    sbe_I = f[:,0]
    sbe_n = f[:,1]

    energy = cx.Wavelength2Energy(wavelength)
    Intensity = np.array(GetIntensity(wavelength))
    TDDFT_nex = PlotNex(wavelength)
    keldysh_parameter, Keldysh_nex = [], []
    print("Intensity [W/cm^2], KelydshParameter, Keldysh_nex [m^-3 s^-1],TDDFT_nex [m^-3 s^-1]")
    for i in range(len(Intensity)):
        gamma, W = k.Kelydsh_W(Intensity=Intensity[i] , photon_energy=energy, band_gap=2.79, RefractiveIndex=1)
        keldysh_parameter.append(gamma)
        Keldysh_nex.append(W)
        print("%.2e,          %.2f,            %.1e ,           %.1e " % (Intensity[i], gamma, W, TDDFT_nex[i]))
    ########  I, n is used to do linear fitting,  I1 is the intensity range that you want to get!!!####
    fitted_Ikm, fitted_km = linearfitting(I = Intensity[0:6], n =  Keldysh_nex[0:6], I1 = Intensity[0:6])
    fitted_Itddft, fitted_tddft = linearfitting(I=Intensity[0:6], n=TDDFT_nex[0:6], I1=Intensity[0:6])
    fitted_Isbe, fitted_nsbe = linearfitting(I=f[0:6,0], n=f[0:6,1], I1=f[0:6,0])
    ################  ax1 主图 ##################################
    l1 = ax1.scatter(Intensity, Keldysh_nex, facecolors="none", edgecolors="red", linewidths=1.5, label="Keldysh")
    l2 = ax1.scatter(Intensity, TDDFT_nex, c="b", linewidths=1.5, label="TDDFT")
    l4 = ax1.scatter(0.5e12, 1.6e39, label="exp laser threshold", color="green",  marker="*", linewidths=5)
    l5 = ax1.scatter(sbe_I, sbe_n ,  facecolors="none",edgecolors="orange", linewidths=1.5, label = "SBE")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.legend(handles=[l1, l2, l4, l5], fontsize=fs - 3, loc='upper left')
    ax1.set_ylabel("Excited electrons $[m^{-3} \cdot s^{-1}]$", fontsize=fs)
    ax1.set_xlabel("Intensity $[W/cm^2]$", fontsize=fs)
    ax2 = ax1.twiny()
    ax2.set_xscale("log")
    ax2.invert_xaxis()
    ax2.set_xlim(max(keldysh_parameter), min(keldysh_parameter))
    ax2.set_xlabel("Keldysh parameter $\gamma$", fontsize=fs)
    format_g = np.around(keldysh_parameter, decimals=2)
    ax2.set_xticks(format_g)
    ax2.set_xticklabels(format_g)
    plt.tight_layout()
    # plt.savefig("nex.eps", dpi = 900, bbox_inches='tight')
    #######  zoom in 放大图  ax4 #################
    start = 0
    end = 3
    l41 = ax4.scatter(Intensity[start:end], Keldysh_nex[start:end], facecolors="none", edgecolors="red", linewidths=1.5, label="Keldysh")
    l42 = ax4.scatter(Intensity[start:end], TDDFT_nex[start:end], c="b", linewidths=1.5, label="TDDFT")
    ax4.loglog(fitted_Ikm[start:end],fitted_km[start:end], "-", c = "r")
    ax4.loglog(fitted_Ikm[start:end], fitted_tddft[start:end], "-", c="b",)
    ax4.loglog(fitted_Isbe[start:6], fitted_nsbe[start:6], "-", c="orange")
    l44 = ax4.scatter(0.5e12, 1.6e39, label="exp laser threshold", color="green",  marker="*", linewidths=5)
    l45 = ax4.scatter(sbe_I[start:6], sbe_n[start:6] ,  facecolors="none",edgecolors="orange", linewidths=1.5, label = "SBE")
    ax4.set_xscale("log")
    ax4.set_yscale("log")
    ax4.legend(handles=[l41, l42, l44, l45], fontsize=fs - 3, loc='upper left')
    ax4.set_ylabel("Excited electrons $[m^{-3} \cdot s^{-1}]$", fontsize=fs)
    ax4.set_xlabel("Intensity $[W/cm^2]$", fontsize=fs)
    # ax4.set_xlim([1.9E11, 2E12])
    # ax4.set_ylim([1E38, 1.6E40 ])

    ax42 = ax4.twiny()
    ax42.set_xscale("log")
    ax42.invert_xaxis()  ## 反转上x轴
    ax42.set_xlabel("Keldysh parameter $\gamma$", fontsize=fs)
    format_g4 = np.around(keldysh_parameter, decimals=1)
    ax42.set_xticks(keldysh_parameter[start:end])
    ax42.set_xticklabels(format_g4[start:end])

    return ()
#################Define the amplitude list for two cases. ####################
names = locals()
names['A' + str(800)] = [5,10,15,20,25,30,55,90]
names['FWHM' + str(800)] = 0.585*600# atom unit [hbar/Ha]   1/e max
names['start' + str(800)] = int(500/0.75)
##################   Main Program    ######################################
WaveLen = 800
Keldysh(wavelength=WaveLen)
ax3.legend()
plt.show()