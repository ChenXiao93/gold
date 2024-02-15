import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op
import OctopusTools.Laser as la
import OctopusTools.chi3 as ch

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
    omega_au = omega * cx.timeau2s / cx.hbar
    ##### 画图 外场及其修正
    # ax1[index].plot(time * 0.02419, A_ext, label="A-ext")
    #### 画图 感应电流及其修正
    # ax2[index].plot(time * 0.02419, current, label="J")
    #### 加上自动化调控衰减！ ####
    t0 = time[-1] * (1 - ratio)
    # eta_au =  np.abs( np.log( np.abs(standard / current[-1])) / (time_ind[-1]-t0)) # 指数包络
    eta_au =  np.sqrt(np.abs(np.log( np.abs(standard / current[-1]))))  / (2 * (time_ind[-1]-t0))  # 高斯包络
    # eta_ev = eta_au * 27.2114
    eta_ev = 0.3
    # time, envy = la.env_exp(mid=(pulse_duration*(1-ratio)), tail=(pulse_duration*ratio), dt=dt, gamma=eta_au)
    envy = la.Env_Gaussian(time_ind, ratio, eta_ev)  # 高斯包络
    # envy = la.Env_exponential(time_ind, ratio, eta_ev)  # 指数包络
    current = current * envy
    # ax2[index].plot(time * 0.02419, envy * max(current), label="$Env$")
    # ax2[index].plot(time * 0.02419, current, label="$J_{mod}$")
    ######
    # J(t) FT J(w)
    Jw_real,Jw_imag = 0, 0
    for i in range(len(current)):
        Jw_real += current[i] * dt * np.cos(omega_au * time[i])  # J_w [au]  # 使用修正后的current进行计算！
        Jw_imag += current[i] * dt * np.sin(omega_au * time[i])
    ## 电流不可以为负？ 傅立叶变换之后均为正
    Jw_real = np.abs(Jw_real)
    Jw_imag = np.abs(Jw_imag)
    # print("J1 = %.2f, J2 = %.2f" % (Jw_real,Jw_imag))
    # A(t) FT A(w)
    Aw_real, Aw_imag = 0, 0  # initialize for FT
    for i in range(len(A_ext)):
        Aw_real += A_ext[i] * dt * np.cos(omega_au * time[i] )
        Aw_imag += A_ext[i] * dt * np.sin(omega_au * time[i] )
    # A(w) ->  E(w)
    E_imag = omega_au * Aw_real / 137
    E_real = -omega_au * Aw_imag / 137
    E_real = np.abs(E_real)
    E_imag = np.abs(E_imag)
    # print("E1 = %.2f, E2 = %.2f" % (E_real,E_imag))
    # E0 = max(A_ext) # 公式里所谓的E0，其实是A0
    # Conductivity  电导  sigma
    denominator = np.square(E_real)+np.square(E_imag)
    sigma_real = ( Jw_real* E_real + Jw_imag * E_imag ) / denominator
    sigma_imag = ( Jw_imag* E_real - Jw_real * E_imag ) / denominator
    # 电导必须为正，负数值没有意义
    sigma_real = np.abs(sigma_real)
    sigma_imag = np.abs(sigma_imag)
    # print("sigma1 = %.2f, sigma2 = %.2f" % (sigma_real, sigma_imag))
    #从电导求介电函数eps, 这里使用了原子单位没有4pi, 但cgs单位制多个4pi, V 在之前求电流时候已经引入了，这里不需要V了？
    # print(sigma_real, sigma_imag)
    epsilon_real = 1 -  sigma_imag  / (omega_au)
    epsilon_imag =   sigma_real / (omega_au)
    return(epsilon_real,epsilon_imag,eta_ev)

#########################################
amp = [1,4,7,10] # Ha
T= [1000,3000,5000,7000,10000] # a.u.
Ts= np.array(T) * cx.timeau2s * 1e12 * 2/5 # ps  # 半高全宽
RATIO = [0.2,0.3,0.3,0.3,0.3]
STD = [1E-8,1E-12,1E-8,1E-8,1E-8]# Ha
CHI3R,CHI3I = [],[]

for t in range(len(T)):
    print("eps1 = -11.980  eps2 = 1.1524 (exp)")
    eps1_current, eps2_current, Intensity, N,K = [], [], [],[],[]
    for i in range(len(amp)):
        # 通过电流计算介电函数
        time_ind, A_ind = cx.ReadData("./induced/t%s/A%s/gauge_field" % (T[t], amp[i]), 1, 2)
        time_ext, A_ext, envy = la.GaussianPulse(PropTime=T[t] + 0.05,
                                                PulseDuration=int(T[t] * 2 / 5),
                                                Wavelength=630,
                                                Amplitude=amp[i], TimeStep=0.05, Phase=0)
        J_cal = cx.current_cal(t_ind=time_ind, A_ind=A_ind, V=114.5822)
        omega = cx.Wavelength2Energy(630)  # eV
        eps1_j, eps2_j, eta_ev = Get_Epsilon_Current(time_ext, A_ext, J_cal, omega, index=i,
                                                    ratio=RATIO[t], standard=STD[t] * amp[i])
        print("eps1 = %.4f, eps2 = %.4f ， A= %s eta=%.2feV.  (电流法)" % (eps1_j, eps2_j, amp[i], eta_ev))
        eps1_current.append(eps1_j)
        eps2_current.append(eps2_j)
        ###### Refractive Index  求折射率和光强 ##################################
        n_real, n_imag = op.refractive_index_positive(eps1_j, eps2_j)
        # n_real, n_imag = op.eps2n_fomo_630nmgold(epsilon_real,epsilon_imag)
        Emax = max(cx.A2E(time_ext, A_ext + A_ind))
        I = cx.E2I(E=Emax, Refractive_index=1)  # [W/cm^2]  # cited from the 'n0n2k0k2'
        Intensity.append(I)
        n_real, n_imag = op.refractive_index_positive(eps1_j, eps2_j)
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
    # n0 = pr[0]
    n0=n0exp
    n2 = pr[1] * 1E-4  # [m^2/W]
    # n2 = np.abs(n2) * (-1)  # 确保n2<0
    # k0 =  pi[0]
    k0 = k0exp
    k2 = pi[1] * 1E-4  # [m^2/W]
    # k2 = np.abs(k2)  # 确保k2>0
    beta = k2 / 300E-9  # [m/W]
    chi3r = ( 4 * n0 * cx.Vaccum_permitivity  * cx.Light_speed / 3) * (n0*n2 - k0*k2)  # [m^2/V^2]
    chi3i = ( 4 * n0 * cx.Vaccum_permitivity  * cx.Light_speed / 3) * (n0*k2 + k0*n2)  # [m^2/V^2]
    print("n0 = %.4f, k0=%.4f" % (pr[0],pi[0]))
    print("n2 = %.2e, k2=%.2e" % (n2,k2))
    print("chi3r = %.2e, chi3i=%.2e" % (chi3r, chi3i))

    CHI3R.append(chi3r)
    CHI3I.append(chi3i)


########  画总图， 读取实验数据 ##############
f630r = np.loadtxt("./exp/Rotenberg630real.txt", delimiter=",")
f630i = np.loadtxt("./exp/Rotenberg630imag.txt", delimiter=",")
f532 = np.loadtxt("./exp/532nm.txt")
f600 = np.loadtxt("./exp/600nm.txt")
f800 = np.loadtxt("./exp/800nm.txt")
f1064= np.loadtxt("./exp/1064nm.txt")
f300tdft = np.loadtxt("./exp/300nmTDDFT.txt")


fig1 = plt.figure()
fs = 12
ax11 = fig1.add_subplot(121)
ax11.axes.loglog(Ts, np.abs(CHI3R),"bo", markerfacecolor = "white", label ="TDDFT630nm")
ax11.axes.loglog(f300tdft[:,0], np.abs(f300tdft[:,1]),"ro", markerfacecolor = "white", label= "TDDFT300nm")
ax11.axes.loglog(f630r[:,0]*1e12 , np.abs(f630r[:,1]),"*", label= "Rotenberg630nm")
ax11.axes.loglog(f532[:,0], np.abs(f532[:,1]),"d", label= "Exp532nm")
# ax11.axes.loglog(f600[0], np.abs(f600[1]),"p", label= "Exp600nm")
ax11.axes.loglog(f800[0], np.abs(f800[1]),"p", label= "Exp800nm")
ax11.axes.loglog(f1064[:,0], np.abs(f1064[:,1]),"p", label= "Exp1064nm")
ax11.set_xlabel("t [ps]", fontsize = fs)
ax11.set_ylabel("$|\chi^3_{real}| [m^2/V^2]$", fontsize = fs)
ax11.legend()


ax12 = fig1.add_subplot(122)
ax12.axes.loglog(Ts, np.abs(CHI3I),"bo",markerfacecolor = "white", label ="TDDFT630nm")
ax12.axes.loglog(f300tdft[:,0], np.abs(f300tdft[:,2]),"ro", markerfacecolor = "white", label= "TDDFT300nm")
ax12.axes.loglog(f630i[:,0]*1e12 , np.abs(f630i[:,1]),"*", label= "Rotenberg630nm")
ax12.axes.loglog(f532[:,0], np.abs(f532[:,2]),"d", label= "Exp532nm")
# ax12.axes.loglog(f600[0], np.abs(f600[2]),"p", label= "Exp600nm")
ax12.set_xlabel("t [ps]", fontsize = fs)
ax12.set_ylabel("$|\chi^3_{imag}| [m^2/V^2]$", fontsize = fs)
ax12.legend()

plt.show()
