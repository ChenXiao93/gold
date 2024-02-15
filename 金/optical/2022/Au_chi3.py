e0 = cx.Vaccum_permitivity  # [F/m]
c = cx.Light_speed  # [m/s]
path = "/data117/chen/jupyter/2nonlinearity/Ir/data/"
Eta = 0  # [Ha/hbar]
wavelength = 700  # [nm]
omega = cx.Wavelength2Omega(wavelength) * cx.timeau2s  # [1/s] --> 1 / [hbar/Ha]
e0 = cx.Vaccum_permitivity  # [F/m]
n0 = 5.26  #####  WAITED TO BE CHANGED!!!!  ###################
dt = 0.1
fs = 14


def dA_over_dt(t, A):
    # t [hbar/Hartree]  A [Hartree]
    dt = t[2] - t[1]
    dAdt = [0]
    for i in range(len(A) - 1):
        dA = A[i + 1] - A[i]
        e = dA / dt
        dAdt.append(e)
    dAdt = np.array(dAdt)  # Ha^2 / hbar
    return (dAdt)


def FT(t, dAdt, Eta, omega):
    # omega can be one value, or even list.   [Ha/hbar]
    # Eta is for the imaginary part of response
    # t [hbar/Ha]; dA [Hartree]
    # return the real and imagnary part of FT result.
    dt = t[2] - t[1]  # dt [hbar / Ha]
    a = 0  # initialize the real of FT
    b = 0  # initialize the imag of FT
    for i in range(len(t)):
        a += dAdt[i] * dt * np.exp(-Eta * t[i]) * np.cos(omega * t[i])
        b += dAdt[i] * dt * np.exp(-Eta * t[i]) * np.sin(omega * t[i])
    return (a, b)


def nIcurve(duration):
    # Pulse_duration  < 2300 [a.u.]

    I = []
    n = []
    n_im = []

    for a in [1, 1.2, 1.4, 1.6, 1.8, 2, 2.1, 2.3, 2.6, 3]:
        #     for a in [ 1,1.2,1.4,1.6,1.8]:

        T_duration, A_ext, A_ind = CutDuration(Amplitude=a, time=5000, Pulse_duration=duration, Plotind=False,
                                               Plotext=False)

        dAdt_ext = dA_over_dt(T_duration, A_ext)
        dAdt_ind = dA_over_dt(T_duration, A_ind)

        ext_real, ext_imag = FT(T_duration, dAdt_ext, Eta, omega)
        ind_real, ind_imag = FT(T_duration, dAdt_ind, Eta, omega)

        Denominator = np.square(ind_real + ext_real) + np.square(ind_imag + ext_imag)
        epsilon_real = (ext_real * (ind_real + ext_real) + ext_imag * (ind_imag + ext_imag)) / Denominator
        epsilon_imag = (ind_real * ext_imag - ext_real * ind_imag) / Denominator

        n_real, n_imag = cx.refractive_index(epsilon_real, epsilon_imag)

        modification = epsilon_real
        #         modification  = 1

        Emax = max(cx.A2E(T_duration, A_ext)) / modification  # epsilon modification
        I_mod = cx.E2I(Emax, n_real) * 4  # [W/cm^2]  # cited from the 'n0n2k0k2' paper

        I.append(I_mod)
        n.append(ind_real)
        n_im.append(ind_imag)

    return (I, n, n_im)


def plotnI():
    plt.figure(figsize=(12, 6))

    for Tpulse in [100, 500, 1000, 1500, 2000]:
        I, n, n_im = nIcurve(Tpulse)

        plt.subplot(1, 2, 1)

        plt.plot(I, n, "--o", label="Tpulse = %s" % Tpulse)
        plt.ylabel("$FT{E_{Ind}}_{real}$", fontsize=fs)
        plt.xlabel("Intensity  [$W/cm^2$]", fontsize=fs)
        plt.legend(fontsize=fs)

        plt.subplot(1, 2, 2)
        plt.plot(I, n_im, "--o")
        plt.ylabel("$FT{E_{Ind}}_{imag}$", fontsize=fs)
        plt.xlabel("Intensity  [$W/cm^2$]", fontsize=fs)

    return ()


plotnI()

plt.show()
