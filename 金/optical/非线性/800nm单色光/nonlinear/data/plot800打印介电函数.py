import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op
import OctopusTools.Laser as la
import OctopusTools.chi3 as ch


def N_Cal(Amplitude, pulse_duration, index):
    ######   A_ext 外场 and A_ind 感应场 ###############
    time_ext, A_ext, envy = la.GaussianPulse(PropTime=pulse_duration + 0.05, PulseDuration=int(pulse_duration * 2 / 5),
                                    Wavelength=800, Amplitude=Amplitude,
                                    TimeStep=0.2,Phase=0)
    time_ind, A_ind = cx.ReadData("./t%s/A%s/gauge_field" % (pulse_duration,Amplitude), 1, 2)
    ax1[index].plot(time_ind * 0.02419, A_ind, label="ind")
    #### 加上自动化调控衰减！ ####
    eta_au = np.log( np.abs(A_ind[-1] / A_standard))  / time_ind[-1]  # atom unit
    time, envy1 = la.env_exp(mid=(pulse_duration*(1-ratio)), tail=(pulse_duration*ratio+0.05), dt=0.2, gamma=eta_au)
    A_ind = A_ind * envy1
    A_ext = A_ext * envy1
    ax1[index].plot(time_ind * 0.02419, envy1 * max(A_ind), label = "env")
    ax1[index].plot(time_ind * 0.02419, A_ind, label = "ind_env")
    ######## 求介电函数 ####################
    wavelength = 800  # nm
    Omega = cx.Wavelength2Energy(wavelength) #1123
    epsilon_real, epsilon_imag = op.Get_Epsilon_Gauge(time_ext, A_ext, time_ind, A_ind, Omega, eta = 0.00)
    print("eps1=%.4f  eps2=%.4f (A=%s)" % (epsilon_real,epsilon_imag,Amplitude ))
    ###### Refractive Index  求折射率和光强 ##################################
    n_real, n_imag = op.refractive_index_positive(epsilon_real, epsilon_imag)
    Emax = max(cx.A2E(time_ext, A_ext+A_ind))
    I = cx.E2I(E = Emax, Refractive_index = 1) * 4  # [W/cm^2]  # cited from the 'n0n2k0k2'
    return(I, n_real, n_imag)

print("eps1=-25.402  eps2=1.3365 (exp)")

ratio = 0.6
A_standard = 1E-12  # Ha
amp = [10,22, 33, 40]
time = 15000
I, n_real, n_imag = [], [], []
fig1 = plt.figure()
ax11 = fig1.add_subplot(221)
ax12 = fig1.add_subplot(222)
ax13 = fig1.add_subplot(223)
ax14 = fig1.add_subplot(224)
ax1 = [ax11,ax12,ax13,ax14]

fig2 = plt.figure()
ax21 = fig2.add_subplot(221)
ax22 = fig2.add_subplot(222)
ax23 = fig2.add_subplot(223)
ax24 = fig2.add_subplot(224)
ax2 = [ax21,ax22,ax23,ax24]

for i in range(len(amp)):
    I_a, n_real_a, n_imag_a = N_Cal(Amplitude = amp[i], pulse_duration = time, index = i)
    # print("I=%.2e [W/cm^2], n=%.4f, k=%.4f" %(I_a,n_real_a,n_imag_a))
    I.append(I_a)
    n_real.append(n_real_a)
    n_imag.append(n_imag_a)
    #######   线性拟合 #######







fitted_yr, pr = ch.nfit(I, n_real, 1)
fitted_yi, pi = ch.nfit(I, n_imag, 1)

fig3 = plt.figure()
ax31 = fig3.add_subplot(121)
ax32 = fig3.add_subplot(122)
ax31.plot(I, n_real,"o")
ax31.plot(I, fitted_yr,"--")
ax31.set_ylabel("n")
ax31.set_xlabel("I [$W/cm^2$]")
ax32.plot(I, n_imag,"o")
ax32.plot(I, fitted_yi,"--")
ax32.set_ylabel("k")
ax32.set_xlabel("I [$W/cm^2$]")

n0 = pr[0]
n2 = pr[1]*1E-4  # [m^2/W]
k0 =  pi[0]
k2 = pi[1]*1E-4 # [m^2/W]
beta = k2 / 300E-9  # [m/W]
chi3r = ( 4 * n0 * cx.Vaccum_permitivity  * cx.Light_speed / 3) * (n0*n2 - k0*k2)  # [m^2/V^2]
chi3i = ( 4 * n0 * cx.Vaccum_permitivity  * cx.Light_speed / 3) * (n0*k2 + k0*n2)  # [m^2/V^2]
n0exp, k0exp = 1.5258, 1.8878
print("n0exp=%.4f  k0exp=%.4f" % (n0exp,k0exp))
print("n0=%.4f  k0=%.4f \nn2=%.4em^2/W \nk2=%.4em^2/W \nbeta=%.2em/W \nchi3r=%.4em^2/V^2 \nchi3i=%.4em^2/V^2"
      % (n0,k0, n2,k2,beta,chi3r,chi3i) )


ax11.set_title("$A_{ind}$")
ax21.set_title("$A_{tot}$")
ax14.legend()
ax24.legend()
plt.show()