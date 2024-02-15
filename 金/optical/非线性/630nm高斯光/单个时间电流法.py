import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op
import OctopusTools.Laser as la
import OctopusTools.chi3 as ch
#########################################
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
######################################
def Get_Epsilon_Current(time, A_ext, current, omega,index, ratio,standard):

    ###########################
    # 第二种： 从电流出发
    # calculate complex dielecreic function epsilon by using current method!
    ## Input parameters:
    #  time  # [hbar/Ha]  propagation time
    #  A_ext # [Ha]  external field
    #  current # [au]    induced current  j(t)
    #  Omega  # eV,   can be a value or  a list
    ## Output parameters
    #epsilon_real  # dimensionless
    #epsilon_imag  # dimensionless
    ###########################
    dt = time[2] - time[1]  # dt [hbar / Ha]
    # octopus输出的电流J 需要乘 - 1 !
    current = current * (-1)
    # omega_au = omega * cx.timeau2s / cx.hbar
    ##### 画图 外场及其修正
    # ax1[index].plot(time * 0.02419, A_ext, label="A-ext")
    #### 画图 感应电流及其修正
    ax2[index].plot(time * 0.02419, current, label="$J_{ori}$")
    #### 加上自动化调控衰减！ ####
    t0 = time[-1] * (1 - ratio)
    # eta_au =  np.abs( np.log( np.abs(standard / current[-1])) / (time_ind[-1]-t0)) # 指数包络
    eta_au =  np.sqrt(np.abs(np.log( np.abs(standard / current[-1]))))  / (2 * (time_ind[-1]-t0))  # 高斯包络
    eta_ev = eta_au * 27.2114
    # eta_ev = 0.2
    # time, envy = la.env_exp(mid=(pulse_duration*(1-ratio)), tail=(pulse_duration*ratio), dt=dt, gamma=eta_au)
    envy = la.Env_Gaussian(time_ind, ratio, eta_ev)  # 高斯包络
    # envy = la.Env_exponential(time_ind, ratio, eta_ev)  # 指数包络
    current = current * envy
    ax2[index].plot(time * 0.02419, envy * max(current), label="$Env$")
    ax2[index].plot(time * 0.02419, current, label="$J_{mod}$")
    ######
    Ew, Jw = cx.At2Ew(omega, time, A_ext, current, mod = "True")
    sigmaw = Jw / Ew
    # sigmaw = cx.positivecomplex(sigmaw)
    epsilonw = 1 + 1j * sigmaw / omega_au
    return(epsilonw,eta_ev)

####################
wavelength = 630
omega = cx.Wavelength2Energy(wavelength) # eV
omega_au = omega/27.2114 # a.u.
amp = [1,4,7,10] # Ha
T= [3000] # a.u.
RATIO = [0.001]
STD = [1E-8]# Ha
print("eps1 = -11.980  eps2 = 1.1524 (exp)")

for t in range(len(T)):
    eps1_current, eps2_current, Intensity, N,K = [], [], [],[],[]
    for i in range(len(amp)):
        # 比较电流，完美一致， 除了第一个点之外，但是第一个点我们设的是0，所以对计算无影响！
        compare_current(pulse_duration = T[t],
                        Amplitude = amp[i],
                        dt = 0.05, index = i)
        # 通过电流计算介电函数  读取数据
        time_ind, A_ind = cx.ReadData("./induced/t%s/A%s/gauge_field" % (T[t],amp[i]), 1, 2)
        time_ext, A_ext, envy = la.GaussianPulse(PropTime=T[t] + 0.05,
                                                PulseDuration=int(T[t] * 2 / 5),
                                                Wavelength=630,
                                                Amplitude=amp[i], TimeStep=0.05, Phase=0)
        J_cal = cx.current_cal(t_ind=time_ind, A_ind=A_ind, V=114.5822)
        # 开始计算
        epsilonw,eta_ev = Get_Epsilon_Current(time_ext, A_ext, J_cal, omega,index = i,
                                            ratio=RATIO[t], standard=STD[t] * amp[i])
        print("eps1 = %.4f, eps2 = %.4f ， A= %s eta=%.2feV.  (电流法)" %(epsilonw.real, epsilonw.imag,amp[i],eta_ev ))
        eps1_current.append(epsilonw.real)
        eps2_current.append(epsilonw.imag)
        ######  求光强 ##################################
        Emax = max(cx.A2E(time_ext, A_ext + A_ind))
        I = cx.E2I(E=Emax, Refractive_index=1)  # [W/cm^2]  # cited from the 'n0n2k0k2'
        Intensity.append(I)
        ###### Refractive Index  折射率
        n_real, n_imag = op.refractive_index_positive(epsilonw.real, epsilonw.imag)
        N.append(n_real)
        K.append(n_imag)

    Intensity = np.array(Intensity)

    n0exp, k0exp = 0.16629, 3.4652
    print("n = %.4f, k = %.4f (exp)" % (n0exp,k0exp) )
    ##################################
    #######   线性拟合 #######
    order = 1
    fitted_yr, pr = ch.nfit(Intensity * 4, N, order)  # cx.E2I 里面是I= 1/2 xxx， 文献里是 I= 2 xxx， 所以有个4的倍数关系
    fitted_yi, pi = ch.nfit(Intensity * 4, K, order)
    n0 = pr[0]
    # n0=n0exp
    n2 = pr[1] * 1E-4  # [m^2/W]
    # n2 = np.abs(n2) * (-1)  # 确保n2<0
    k0 =  pi[0]
    # k0 = k0exp
    k2 = pi[1] * 1E-4  # [m^2/W]
    # k2 = np.abs(k2)  # 确保k2>0
    beta = k2 / (wavelength * 1E-9)  # [m/W]
    chi3r = ( 4 * n0 * cx.Vaccum_permitivity  * cx.Light_speed / 3) * (n0*n2 - k0*k2)  # [m^2/V^2]
    chi3i = ( 4 * n0 * cx.Vaccum_permitivity  * cx.Light_speed / 3) * (n0*k2 + k0*n2)  # [m^2/V^2]
    print("n0 = %.4f, k0=%.4f" % (pr[0],pi[0]))
    print("n2 = %.2e, k2=%.2e" % (n2,k2))
    print("chi3r = %.2e, chi3i=%.2e" % (chi3r, chi3i))





    ##################################
    # plt.figure()  # 4个点的介电函数
    # plt.subplot(1,2,1)
    # plt.plot(Intensity, eps1_current, "-o", label = "$\epsilon_1$")
    # plt.legend()
    # plt.subplot(1,2,2)
    # plt.plot(Intensity, eps2_current, "-o", label = "$\epsilon_2$")
    # plt.legend()

    plt.figure() # 4个点的 refractive index
    plt.subplot(1,2,1)
    plt.plot(Intensity, N, "o",label = "$n$")
    plt.plot(Intensity, fitted_yr, "--")
    plt.legend()
    plt.subplot(1,2,2)
    plt.plot(Intensity, K, "o",label = "$k$")
    plt.plot(Intensity, fitted_yi, "--")
    plt.legend()

    ax14.legend()
    ax24.legend()
plt.show()
