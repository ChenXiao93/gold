import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op
import OctopusTools.Laser as la
import OctopusTools.chi3 as ch

def N_Cal(Amplitude, pulse_duration):
    ######   A_ext 外场 and A_ind 感应场 ###############
    time_ext, A_ext = cx.ReadData("./laser/A%s" % Amplitude, 0, 1)
    time_ind, A_ind = cx.ReadData("./gauge_field/A%s" % Amplitude, 1, 2)
    ###### 先截断入射场和感应场，总体传播时间长度=脉冲持续时间+500头+500尾
    start = 0
    end = int((500+500+pulse_duration)/0.2)  #
    time_ext,  A_ext= cx.shortten(time_ext, A_ext, start, end)
    time_ind, A_ind = cx.shortten(time_ind, A_ind, start, end)
    if Amplitude == 10:
        ax1.plot(time_ind * 0.02419, A_ind, label="ind") # 画出截断后的原始感应场
    ######## Add the envelope to cut it 加包络 切断时间 ######################
    #### 加上自动化调控衰减！ ####
    ratio = 0.5
    A_standard = 1E-8  # Ha
    eta_au = np.log( np.abs(A_ind[-1] / A_standard))  / time_ind[-1]  # atom unit
    # print("eta = %.4f eV" % (eta_au*27.2114))
    time, envy = la.env_exp(mid=int(pulse_duration*(1-ratio)) + 500,
                            tail=int(pulse_duration*ratio) + 500, dt=0.2, gamma=eta_au)
    A_ind = A_ind * envy
    A_ext = A_ext * envy
    if Amplitude == 10:
        ax1.plot(time_ind * 0.02419, A_ind,label = "ind_env")  # 画出截断后的经过exp衰减修正的感应场
        ax1.plot(time_ind * 0.02419, envy * max(A_ind), label="env") # 画出exp衰减包络
    ######## 求介电函数 ####################
    wavelength = 300  # nm
    Omega = cx.Wavelength2Energy(wavelength) # eV
    epsilon_real, epsilon_imag = op.Get_Epsilon_Gauge(time_ext, A_ext, time_ind, A_ind, Omega, eta = 0.0)
    print("eps1=%.4f  eps2=%.4f" % (epsilon_real,epsilon_imag))
    ###### Refractive Index  求折射率和光强 ##################################
    n_real, n_imag = op.refractive_index_positive(epsilon_real, epsilon_imag)
    Emax = max(cx.A2E(time_ext, A_ext+A_ind))
    I = cx.E2I(E = Emax, Refractive_index = 1) * 4  # [W/cm^2]  # cited from the 'n0n2k0k2'
    return(I, n_real, n_imag)

def chi3_cal(pulse_duration):
    I, n_real, n_imag = [], [], []
    for a in [10,20,30,40]:
        I_a, n_real_a, n_imag_a = N_Cal(a,pulse_duration)
        I.append(I_a)
        n_real.append(n_real_a)
        n_imag.append(n_imag_a)
    #######   线性拟合 #######
    fitted_yr, pr = ch.nfit(I, n_real, 1)
    fitted_yi, pi = ch.nfit(I, n_imag, 1)
    n0 = pr[0]
    n2 = pr[1]*1E-4  # [m^2/W]
    k0 =  pi[0]
    k2 = pi[1]*1E-4 # [m^2/W]

    chi3r = ( 4 * n0 * cx.Vaccum_permitivity  * cx.Light_speed / 3) * (n0*n2 - k0*k2)  # [m^2/V^2]
    chi3i = ( 4 * n0 * cx.Vaccum_permitivity  * cx.Light_speed / 3) * (n0*k2 + k0*n2)  # [m^2/V^2]
    return(chi3r,chi3i)
##############################################################
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
print("eps1=-1.2360  eps2=5.7608  Johnson and Christy 1972")
Time= np.arange(0,5000,1000)
# Time = np.array([3000])
CHI3R, CHI3I = [],[]
for pulse_duration in Time:
    chi3r,chi3i = chi3_cal(pulse_duration)
    CHI3R.append(chi3r)
    CHI3I.append(chi3i)

ax1.legend()

# print(CHI3R,CHI3I)
plt.figure()
plt.subplot(1,2,1)
plt.loglog(Time * 0.02419 * 1E-15, np.abs(CHI3R), "--o")

plt.subplot(1,2,2)
plt.loglog(Time * 0.02419 * 1E-15, np.abs(CHI3I), "--o")

plt.show()