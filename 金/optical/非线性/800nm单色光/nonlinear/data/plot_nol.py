import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op
import OctopusTools.Laser as la
import OctopusTools.chi3 as cc

time = [500,1500,4000]
amplitude = [10,22,33,40]
wavelength = 800
Omega = cx.Wavelength2Energy(wavelength)  # eV
eta = 0.01 # eV

N0,N2,K0,K2,CHI3R,CHI3I = [],[],[],[],[],[]
for t in time:
    Intensity, epsilon1, epsilon2, nn, kk = [], [], [], [], []
    for a in amplitude:
        t_ext, A_ext, A_env = la.GaussianPulse(PropTime=t + 0.05,
                                               PulseDuration=int(t * 2 / 5),
                                               Wavelength=800,
                                               Amplitude=a,
                                               TimeStep=0.2,
                                               Phase=0)
        file = "./t%s/A%s/gauge_field" % (t,a)
        t_ind = cx.ReadOne(file, 1)
        A_ind = cx.ReadOne(file, 2)
        eps1,eps2 = op.Get_Epsilon_Gauge(t_ext, A_ext, t_ind, A_ind, Omega, eta)
        n,k = op.simple_refractive_index(eps1, eps2)
        E0 = cx.A2E(t_ext, A_ext) # V/Angstrom
        I = cx.E2I(E = np.abs(max(E0) / eps1), Refractive_index = 1) # [W/cm^2]

        epsilon1.append(eps1)
        epsilon2.append(eps2)
        nn.append(n)
        kk.append(k)
        Intensity.append(I)

        print("T = %s,A = %s,eps1 = %.4f, eps2 = %.4f, n = %.4f, k = %.4f, I = %.2e[W/cm^2] "
              % (t,a,eps1, eps2,n,k,I))

    fitted_nn, para_nn = cc.nfit(Intensity, nn, 1)
    print(nn)
    print(fitted_nn)
    fitted_kk, para_kk = cc.nfit(Intensity, kk, 1)

    n0 = para_nn[0]
    n2 = para_nn[1] * 1E-4  # [m^2/W]
    k0 = para_kk[0]
    k2 = para_kk[1] * 1E-4  # [m^2/W]
    #  计算三阶电极化率
    chi3r = (1 * n0 * cx.Vaccum_permitivity * cx.Light_speed / 3) * (n0 * n2 - k0 * k2)  # [m^2/V^2]
    chi3i = (1 * n0 * cx.Vaccum_permitivity * cx.Light_speed / 3) * (n0 * k2 + k0 * n2)  # [m^2/V^2]

    N0.append(n0)
    N2.append(n2)
    K0.append(k0)
    K2.append(k2)
    CHI3R.append(chi3r)
    CHI3I.append(chi3i)

    print("n0=%.4f, n2 = %.4e [m^2/W]" % (n0, n2))
    print("k0=%.4f, k2 = %.4e [m^2/W]" % (k0, k2))
    print("chi3r=%.4e [m^2/V^2], chi3i = %.4e [m^2/V^2]" % (chi3r, chi3i))
    #########################
    # plt.figure()
    # plt.subplot(1,2,1)
    # plt.plot(Intensity,epsilon1,"--o")
    # plt.subplot(1,2,2)
    # plt.plot(Intensity,epsilon2,"--o")

    plt.figure()
    fs = 12
    plt.subplot(1,2,1)
    plt.plot(Intensity,nn,"o")
    plt.plot(Intensity, fitted_nn, "--")
    plt.ylabel("n", fontsize = fs)
    plt.subplot(1,2,2)
    plt.plot(Intensity, kk, "o")
    plt.plot(Intensity,fitted_kk,"--")
    plt.title("T=%s au = %.1f fs" % (t, t / 27.2114) )
    plt.ylabel("k", fontsize=fs)

plt.figure()
plt.subplot(1,2,1)
plt.loglog(np.array(time) * 0.02419 * 1E-15, np.abs(CHI3R), "--o" )
plt.scatter(200*1e-15, 1e-19, color = "red")
print(time)
print(CHI3R)
print(CHI3I)

plt.subplot(1,2,2)
plt.loglog(np.array(time) * 0.02419 * 1E-15  * 2 / 5, np.abs(CHI3I), "--o" )


plt.show()