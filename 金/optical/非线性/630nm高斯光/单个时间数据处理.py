import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op
import OctopusTools.Laser as la
import OctopusTools.chi3 as ch


def N_Cal(Amplitude, pulse_duration, index, dt):
    ######   A_ext 外场 and A_ind 感应场 ###############
    time_ext, A_ext = cx.ReadData("./laser/t%s/A%s" % (pulse_duration,Amplitude), 0, 1)
    time_ind, A_ind = cx.ReadData("./induced/t%s/A%s/gauge_field" % (pulse_duration,Amplitude), 1, 2)
    start = 0
    end = int(pulse_duration/dt)  # 把ext场和ind场调整到长度一致. 100是额外添加的尾部，用于平复感应场的滞后波动
    time_ext,  A_ext= cx.shortten(time_ext, A_ext, start, end)
    time_ind, A_ind = cx.shortten(time_ind, A_ind, start, end)
    A_tot1 = A_ext + A_ind

    # E_ind = cx.A2E(time_ind, A_ind)
    # E_tot = cx.A2E(time_ind, A_tot1)
    # bili = max(E_tot) / max(E_ext)
    # print("bili = %.2f" % bili)

    # plt.figure()
    # plt.subplot(2,2,i+1)
    # plt.plot(time_ext * 0.02419,  A_ext, "blue", label = "ext")
    # plt.plot(time_ind*0.02419 , A_ind, label = "ind")
    ax1[index].plot(time_ind * 0.02419, A_ind, label="Ori")
    #### 加上自动化调控衰减！ ####
    eta_au = np.log( np.abs(A_ind[-1] / A_standard))  / time_ind[-1]  # atom unit
    # eta_au= 0.05
    eta_ev = eta_au * 27.2114  # eV
    etafull.append(eta_ev)
    # print("eta = %.4feV" % eta_ev)
    time, envy = la.env_exp(mid=int(pulse_duration*(1-ratio)), tail=int(pulse_duration*ratio) , dt=dt, gamma=eta_au)

    A_ind = A_ind * envy
    A_ext = A_ext * envy
    A_tot2 = A_ext + A_ind

    # Emax_tot = max(cx.A2E(time_ext, A_tot2)[start:end])
    ax1[index].plot(time_ind * 0.02419, envy * max(A_ind),label = "Env")
    ax1[index].plot(time_ind * 0.02419, A_ind,label = "Mod")

    ax2[index].plot(time_ind * 0.02419, A_tot1, label = "Ori")
    ax2[index].plot(time_ind * 0.02419, A_tot2,label = "Mod")
    ######## 求介电函数 ####################
    wavelength = 630  # nm
    Omega = cx.Wavelength2Energy(wavelength) #1123
    epsilon_real, epsilon_imag = op.Get_Epsilon_Gauge(time_ext, A_ext, time_ind, A_ind, Omega, eta = 0.00)
    print("eps1=%.4f  eps2=%.4f (A=%s)" % (epsilon_real,epsilon_imag,Amplitude ))
    ###### Refractive Index  求折射率和光强 ##################################
    modification = np.sqrt(epsilon_real**2 + epsilon_imag**2)   # 修正1  ｜\epsilon｜
    # modification = np.abs(epsilon_real)  # 修正2  ｜\epsilon_{real}｜
    n_real, n_imag = op.refractive_index_positive(epsilon_real, epsilon_imag)
    # Emax_tot = max(cx.A2E(time_ext, A_ext+A_ind)[start:end])
    # print(Emax_ext,Emax_ind,Emax_tot,Emax_ext/Emax_tot)
    # 真空中的光强！才可以用epsilon修正 来近似 介质内部的光强！  cx.E2I 里面是I= 1/2 xxx， 文献里是 I= 2 xxx， 所以有个4的倍数关系
    # 用总电场来计算，就无需修正了！！！ 总电场就是晶体内部电子所感受到的真实的电场！
    Emax = max(cx.A2E(time_ext, A_ext)[start:end])/modification
    I = cx.E2I(E = Emax, Refractive_index = 1)  # [W/cm^2]  # cited from the 'n0n2k0k2'
    return(I, n_real, n_imag)

print("eps1 = -11.980  eps2 = 1.1524 (exp)")
ratio = 0.001
A_standard = 1E-4 # Ha
amp = [1,4,7,10]
time = 10000
dt = 0.05
I, n_real, n_imag,etafull  = [], [], [],[]
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
    I_a, n_real_a, n_imag_a = N_Cal(Amplitude = amp[i], pulse_duration = time, index = i,dt=dt)
    I_a = I_a * 4   #  cx.E2I 里面是I= 1/2 xxx， 文献里是 I= 2 xxx， 所以有个4的倍数关系
    # print("I=%.2e [W/cm^2], n=%.4f, k=%.4f" %(I_a,n_real_a,n_imag_a))
    I.append(I_a)
    n_real.append(n_real_a)
    n_imag.append(n_imag_a)
    #######   线性拟合 #######







fitted_yr, pr = ch.nfit(I, n_real, 1)
fitted_yi, pi = ch.nfit(I, n_imag, 1)

fig3 = plt.figure()
ax31 = fig3.add_subplot(211)
ax32 = fig3.add_subplot(212)
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



for i in range(4):
    print("eta%s=%.2feV" %(i+1,etafull[i]) )

fs = 12


ax14.legend()
ax24.legend()

ax24.set_xlabel("t [fs]", fontsize = fs)
ax23.set_xlabel("t [fs]", fontsize = fs)
ax21.set_ylabel("$A_{tot} [Ha]$", fontsize = fs)
ax11.set_ylabel("$A_{ind} [Ha]$", fontsize = fs)

plt.show()