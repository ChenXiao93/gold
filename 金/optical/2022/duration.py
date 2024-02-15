import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx

def env(x,y, pulse_duration):
    #  sin2 attenuation
    envy =   max(y) * (np.heaviside(x - pulse_duration, 1)   * \
             np.square( np.sin( (x - pulse_duration)  *  np.pi / tail /2    + np.pi/2)) + \
             np.heaviside(pulse_duration - x - 0.01, 1) ) *  np.heaviside( pulse_duration + tail - x, 1)
    newy = envy * y / max(envy)
    return(newy)

def FT(t, dAdt, omega):
# omega can be one value, or even list.   [Ha/hbar]
# t [hbar/Ha]; dA [Hartree]
# return the real and imagnary part of FT result.
    dt = t[2] - t[1]  # dt [hbar / Ha]
    a = 0  # initialize the real of FT
    b = 0  # initialize the imag of FT
    for i in range(len(t)):
        a += dAdt[i] * dt  * np.cos(omega * t[i])
        b += dAdt[i] * dt  * np.sin(omega * t[i])
    return(a,b)

def FT1(t, dAdt, omega):
# omega can be one value, or even list.   [Ha/hbar]
# t [hbar/Ha]; dA [Hartree]
# return the real and imagnary part of FT result.
    dt = t[2] - t[1]  # dt [hbar / Ha]
    a = 0  # initialize the real of FT
    b = 0  # initialize the imag of FT
    for i in range(len(t)):
        a += dAdt[i] * dt  * np.cos(omega * t[i])
        b += dAdt[i] * dt  * np.sin(omega * t[i])
    return(a,b)
#######################################################################################
def Field_Plot(Amplitude, pulse_duration):
    ######   A_ext and A_ind ###############
    time_ext, A_ext = cx.ReadData(path+"A%s/laser" % Amplitude, 1, 2)
    time_ind, A_ind = cx.ReadData(path+"A%s/gauge_field" % Amplitude, 1, 2)
    start = 0
    end = len(A_ind)
    time_ext, A_ext = cx.shortten(time_ext, A_ext, start, end)
    time_ind, A_ind = cx.shortten(time_ind, A_ind, start, end)

    ######## Add the envelope to cut it ######################
    A_ext = env(time_ext, A_ext, pulse_duration+head)
    A_ind = env(time_ind, A_ind, pulse_duration+head)

    ########    dA/dt for ext and ind  ######################
    dAdt_ext = cx.dA_over_dt(time_ext, A_ext)
    dAdt_ind = cx.dA_over_dt(time_ind, A_ind)

    #####   FT for full frequency ############################
    ext_real, ext_imag = FT(time_ext, dAdt_ext, omega)
    ind_real, ind_imag = FT(time_ind, dAdt_ind, omega)

    ########   Dielectruic function ##########################
    Denominator = np.square(ind_real + ext_real) + np.square(ind_imag + ext_imag)
    epsilon_real = (ext_real * (ind_real + ext_real) + ext_imag * (ind_imag + ext_imag)) / Denominator
    epsilon_imag = (ind_real * ext_imag - ext_real * ind_imag) / Denominator

    n_real, n_imag = cx.refractive_index(epsilon_real, epsilon_imag)

    return(n_real, n_imag )
##########################################################################################
def N_Cal(Amplitude, pulse_duration):
    ######   A_ext and A_ind ###############
    time_ext, A_ext = cx.ReadData(path+"A%s/laser" % Amplitude, 1, 2)
    time_ind, A_ind = cx.ReadData(path+"A%s/gauge_field" % Amplitude, 1, 2)
    start = 0
    end = len(A_ind)
    time_ext, A_ext = cx.shortten(time_ext, A_ext, start, end)
    time_ind, A_ind = cx.shortten(time_ind, A_ind, start, end)

    ######## Add the envelope to cut it ######################
    A_ext = env(time_ext, A_ext, pulse_duration+head)
    A_ind = env(time_ind, A_ind, pulse_duration+head)

    ########    dA/dt for ext and ind  ######################
    dAdt_ext = cx.dA_over_dt(time_ext, A_ext)
    dAdt_ind = cx.dA_over_dt(time_ind, A_ind)

    #####   FT for full frequency ############################
    ext_real, ext_imag = FT(time_ext, dAdt_ext, omega)
    ind_real, ind_imag = FT(time_ind, dAdt_ind, omega)

    ########   Dielectruic function ##########################
    Denominator = np.square(ind_real + ext_real) + np.square(ind_imag + ext_imag)
    epsilon_real = (ext_real * (ind_real + ext_real) + ext_imag * (ind_imag + ext_imag)) / Denominator
    epsilon_imag = (ind_real * ext_imag - ext_real * ind_imag) / Denominator
    ###### Refractive Index ##################################
    modification = epsilon_real
    n_real, n_imag = cx.refractive_index(epsilon_real, epsilon_imag)
    Emax = max(cx.A2E(time_ext, A_ext))  / modification  # epsilon modification
    I = cx.E2I(Emax, n_real) * 4  # [W/cm^2]  # cited from the 'n0n2k0k2' paper
    # I = cx.E2I(Emax, 1)
    # I = cx.E2I(Emax, n_real)
    # I = cx.E2I(Emax, n_real) * 4
    # I = cx.E2I(Emax, 1) * 4
    # print("%.2e   \n  %.2f   \n  %.2f" % (I, n_real, n_imag))
    return(I, n_real, n_imag)
#########################################################################
def N_plot(pulse_duration):
    I, n_real, n_imag = [], [], []
    for a in [1.2,1.6,2]:
        I_a, n_real_a, n_imag_a = N_Cal(a,pulse_duration)
        I.append(I_a)
        n_real.append(n_real_a)
        n_imag.append(n_imag_a)

    plt.subplot(1,2,1)
    plt.plot(I, n_real, "--o")

    plt.subplot(1,2,2)
    plt.plot(I, n_imag, "--o")
#### N-th order polynomial fitting  ##################################
def nfit(x, y, n):
    x = np.array(x)
    fitting = np.polyfit(x, y, n)
    parameters = np.poly1d(fitting)
    fitted_y = 0
    for i in range(n + 1):
        fitted_y += np.power(x, i) * parameters[i]

    return (fitted_y, parameters)
######################################################################################################
def chi3_cal(pulse_duration):
    I, n_real, n_imag = [], [], []
    for a in [1.2,1.6,2]:
        I_a, n_real_a, n_imag_a = N_Cal(a,pulse_duration)
        I.append(I_a)
        n_real.append(n_real_a)
        n_imag.append(n_imag_a)
    #######   fitting order = 1 #######
    fitted_yr, pr = nfit(I, n_real, 1)
    fitted_yi, pi = nfit(I, n_imag, 1)
    n0 = pr[0]
    n2 = pr[1]*1E-4
    k0 =  pi[0]
    k2 = pi[1]*1E-4

    chi3r = (4 * n0 * e0 * c / 3) * (n0*n2 - k0*k2)
    chi3i = (4 * n0 * e0 * c / 3) * (n0*k2 + k0*n2)

    return(chi3r,chi3i)
#########################   Parameters  #########################################
head  =1000
tail = 1000
fs = 14
path = "./envelop_data/"
wavelength_0 = 630  # [nm]
omega = cx.Wavelength2Omega(wavelength_0) * cx.timeau2s  # [1/s] --> 1 / [hbar/Ha]
e0 = cx.Vaccum_permitivity   #[F/m]
c = cx.Light_speed   # [m/s]
# omega = np.arange(omega_0 * 0.96, 1.04 * omega_0, omega_0*0.001) # To do the FT for more omega points
# Time = np.arange(0,3300,300)
Time= 1500
CHI3R, CHI3I = [], []
Amplitude = [0.1, 0.7, 1.2, 1.6, 2]
ep1, ep2 = [], []
#### Fig1  Enveloped cutted external field ##################
for a in Amplitude:
    epsilon_real, epsilon_imag = Field_Plot(a, pulse_duration = Time)
    ep1.append( epsilon_real)
    ep2.append(epsilon_imag)

plt.subplot(1,2,1)
plt.plot(Amplitude, ep1, "--o")
# plt.legend(fontsize = fs)
plt.xlabel("Amplitude A [Ha]", fontsize =  fs)
plt.ylabel("$n$", fontsize = fs)

plt.subplot(1,2,2)
plt.plot(Amplitude, ep2, "--o")
# plt.legend(fontsize = fs)
plt.xlabel("Amplitude A [Ha]", fontsize =  fs)
plt.ylabel("$k$", fontsize = fs)

plt.show()

#### Fig2  Enveloped cutted induced field ##################
# for a in Amplitude:
#     epsilon_real, epsilon_imag = Field_Plot(a, pulse_duration=Time)
#
#     plt.plot(Amplitude, epsilon_imag)
#
# # plt.legend(fontsize=fs)
# plt.xlabel("Amplitude A [Ha]", fontsize=fs)
# plt.ylabel("\epsilon_{imag} [a.u]$", fontsize=fs)
# plt.show()
############# Fig 3 ###############################################

# for pulse_duration in Time:
#     chi3r,chi3i = chi3_cal(pulse_duration)
#     CHI3R.append(chi3r)
#     CHI3I.append(chi3i)
#     print(chi3r, chi3i)

#####  Ref #############
# rotenberg = np.loadtxt("./exp/Rotenberg630real.txt", delimiter=',')
# rotenbergimag = np.loadtxt("./exp/Rotenberg630imag.txt", delimiter=',')
# #############################
# plt.subplot(1,2,1)
# plt.loglog(Time * 0.02419 * 1E-15, np.abs(CHI3R), "--o")
# plt.loglog(rotenberg[:,0],  rotenberg[:,1],"--o", label = "exp", markersize = 5)
# plt.subplot(1,2,2)
# plt.loglog(Time * 0.02419 * 1E-15, np.abs(CHI3I), "--o")
# plt.loglog(rotenbergimag[:,0],  rotenbergimag[:,1],"--o", label = "exp", markersize = 5)
#
# plt.show()












