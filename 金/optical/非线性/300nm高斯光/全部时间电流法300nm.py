import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.Laser as la
import OctopusTools.chi3 as ch

def Get_Epsilon_Current(t, A, J, w_ev, V, index, ratio, standard):
    ## ratio, standard 为了个性化调控 damping 修正！
    dt = t[2] - t[1]  # dt [hbar / Ha]
    w_au = w_ev /cx.hartree2ev
    #### 个性化调控 J damping 修正！ ####
    t0 = t[-1] * (1 - ratio)
    # eta_au =  np.abs( np.log( np.abs(standard / J[-1])) / (time_ind[-1]-t0)) # 指数包络
    eta_au =  np.sqrt(np.abs(np.log( np.abs(standard / J[-1]))))  / (2 * (time[-1]-t0))  # 高斯包络
    eta_ev = eta_au * cx.hartree2ev
    # time, envy = la.env_exp(mid=(pulse_duration*(1-ratio)), tail=(pulse_duration*ratio), dt=dt, gamma=eta_au)
    envy = la.Env_Gaussian(time, ratio, eta_ev)  # 高斯包络
    # envy = la.Env_exponential(time_ind, ratio, eta_ev)  # 指数包络
    J_damp = J * envy
    ################
    Et = cx.A2Eau(t,A) # a.u.
    Ew, Jw = 0,0 # initialize for Fourier Transform
    for i in range(len(A)):
        Ew += Et[i] * dt * np.exp(w_au * t[i] * 1j)  # 外场无需加修正
        Jw += J_damp[i] * dt * np.exp(w_au * t[i] * 1j)  # 感应场杂峰多，加修正
    sigma = Jw / (Ew * V)
    epsilon = 1 + 4 * np.pi * 1j * sigma / w_au
    return(epsilon,eta_ev)

def GaussianPulse(PropTime, PulseDuration, Wavelength, Amplitude, TimeStep, tail):
    PeakPosition = PropTime / 2   # Pulse peak position tau [fs]
    omega_si =  cx.Wavelength2Omega(Wavelength)  # Wavelength [nm]  to  Omega[1/s]
    omega_au = omega_si  * cx.timeau2s  #
    time = np.arange(0, PropTime + tail, TimeStep)
    envelope = Amplitude * np.exp(-np.square((time - PeakPosition) / (PulseDuration/2) ))
    gauge_field = envelope * np.cos(omega_au * (time-PeakPosition) )
    return (time, gauge_field, envelope)   # t[hbar/Ha] A [Hartree/e]
#########################################
amp = [1,4,7,10] # Ha
T= [500, 1500, 3000, 5000, 10000] # a.u.
Ts= np.array(T) * cx.timeau2s * 1e12 * 2/5 # ps  # FWHM, 因为高斯光产生的时候设置PulseDuration=int(t*2/5)
CHI3R,CHI3I = [],[]
dt = [0.2, 0.2, 0.2, 0.2, 0.05] # a.u.
wavelegth = 300 # nm
omega = cx.Wavelength2Energy(wavelegth)  # eV
Volume = 114.5822 # bohr^3

RATIO = [0.2, 0.2, 0.2, 0.2, 0.2]
STD = [1E-8,1E-8,1E-8,1E-8,1E-8]# Ha
#########################################
for t in range(len(T)):
    print("eps1 = -1.236  eps2 = 5.761 (ex)")  # 实验值
    eps1_current, eps2_current, Intensity, N,K = [], [], [],[],[] # 初始化
    for i in range(len(amp)):
        # 通过电流计算介电函数
        time, A_ext, envy = GaussianPulse(PropTime=T[t]+ dt[i], PulseDuration=int(T[t] * 2 / 5),
                                            Wavelength= wavelegth, Amplitude=amp[i],
                                            TimeStep=dt[i],tail = 100)
        A_ind = cx.ReadOne("./data/t%s/A%s/gauge_field" % (T[t], amp[i]), 2)
        ### 截断 对齐
        cutoff = min(len(A_ext), len(A_ind))
        A_ext = A_ext[0: cutoff]
        A_ind = A_ind[0: cutoff]
        time = time[0: cutoff]
        ###  计算电流
        J_cal = cx.current_cal(t_ind=time, A_ind=A_ind, V=Volume) # calculated 感应电流

        epsilonw, eta_ev = Get_Epsilon_Current(t = time, A = A_ext,
                                            J = J_cal, w_ev = omega,
                                            V = Volume , index = i,
                                            ratio = RATIO[t], standard = STD[t] * amp[i] )
        print("eps1 = %.4f, eps2 = %.4f ， A= %s eta=%.2feV.  (电流法)" % (epsilonw.real,epsilonw.imag, amp[i], eta_ev))
        eps1_current.append(epsilonw.real)
        eps2_current.append(epsilonw.imag)
        ######  求光强 ##################################
        Emax = max(cx.A2E(time, A_ext + A_ind))  # Etot
        I = cx.E2I(E=Emax, Refractive_index=1)  # [W/cm^2]  # cited from the 'n0n2k0k2'
        Intensity.append(I)
        ###### Refractive Index  折射率
        n_real, n_imag = cx.eps2n(epsilonw)
        N.append(n_real)
        K.append(n_imag)


    Intensity = np.array(Intensity)
    n0exp, k0exp = 1.5258, 1.8878
    print("n = %.4f, k = %.4f (exp)" % (n0exp,k0exp))
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
    beta = k2 / (wavelegth * 1E-9)  # [m/W]
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

fig1 = plt.figure()
fs = 12
ax11 = fig1.add_subplot(121)
ax11.axes.loglog(Ts, np.abs(CHI3R),"bo", markerfacecolor = "white", label ="TDDFT%snm" % wavelegth)
ax11.axes.loglog(f630r[:,0]*1e12 , np.abs(f630r[:,1]),"*", label= "Rotenberg630nm")
ax11.axes.loglog(f532[:,0], np.abs(f532[:,1]),"d", label= "Exp532nm")
ax11.axes.loglog(f600[0], np.abs(f600[1]),"p", label= "Exp600nm")
ax11.axes.loglog(f800[0], np.abs(f800[1]),"p", label= "Exp800nm")
ax11.axes.loglog(f1064[:,0], np.abs(f1064[:,1]),"p", label= "Exp1064nm")
ax11.set_xlabel("t [ps]", fontsize = fs)
ax11.set_ylabel("$|\chi^3_{real}| [m^2/V^2]$", fontsize = fs)
ax11.legend()


ax12 = fig1.add_subplot(122)
ax12.axes.loglog(Ts, np.abs(CHI3I),"bo",markerfacecolor = "white", label ="TDDFT%snm" % wavelegth)
ax12.axes.loglog(f630i[:,0]*1e12 , np.abs(f630i[:,1]),"*", label= "Rotenberg630nm")
ax12.axes.loglog(f532[:,0], np.abs(f532[:,2]),"d", label= "Exp532nm")
ax12.axes.loglog(f600[0], np.abs(f600[2]),"p", label= "Exp600nm")
ax12.set_xlabel("t [ps]", fontsize = fs)
ax12.set_ylabel("$|\chi^3_{imag}| [m^2/V^2]$", fontsize = fs)
ax12.legend()

plt.show()
