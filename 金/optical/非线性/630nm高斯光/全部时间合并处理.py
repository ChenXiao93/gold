import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op
import OctopusTools.Laser as la
import OctopusTools.chi3 as ch

def N_Cal(Amplitude, pulse_duration, index, dt, ratio,A_standard):
    ######   A_ext 外场 and A_ind 感应场 ###############
    time_ext, A_ext = cx.ReadData("./laser/t%s/A%s" % (pulse_duration,Amplitude), 0, 1)
    time_ind, A_ind = cx.ReadData("./induced/t%s/A%s/gauge_field" % (pulse_duration,Amplitude), 1, 2)
    start = 0
    end = int((pulse_duration)/dt)  # 把ext场和ind场调整到长度一致. 100是额外添加的尾部，用于平复感应场的滞后波动
    time_ext,  A_ext= cx.shortten(time_ext, A_ext, start, end)
    time_ind, A_ind = cx.shortten(time_ind, A_ind, start, end)
    #### 加上自动化调控衰减！ ####
    eta_au = np.log( np.abs(A_ind[-1] / A_standard))  / time_ind[-1]  # atom unit
    eta_ev = eta_au * 27.2114
    # time, envy = la.env_exp(mid=(pulse_duration*(1-ratio)), tail=(pulse_duration*ratio), dt=dt, gamma=eta_au)
    envy = la.Env_Gaussian(time_ind, ratio, eta_ev)
    A_ind = A_ind * envy
    # A_ext = A_ext * envy  # 无需对外场进行修正！
    ######## 求介电函数 ####################
    wavelength = 630  # nm
    Omega = cx.Wavelength2Energy(wavelength)
    epsilon_real, epsilon_imag = op.Get_Epsilon_Gauge(time_ext, A_ext, time_ind, A_ind, Omega, 0)
    epsilon_real = -1 * np.abs(epsilon_real)
    epsilon_imag = +1 * np.abs(epsilon_imag)
    ###### Refractive Index  求折射率和光强 ##################################
    n_real, n_imag = op.refractive_index_positive(epsilon_real, epsilon_imag)
    # n_real, n_imag = op.eps2n_fomo_630nmgold(epsilon_real,epsilon_imag)
    Emax = max(cx.A2E(time_ext, A_ext+A_ind))
    I = cx.E2I(E = Emax, Refractive_index = 1)  # [W/cm^2]  # cited from the 'n0n2k0k2'
    return(I, n_real, n_imag)

# RATIO = [0.2, 0.2, 0.2, 0.025, 0.3]
# AST = [1E-6, 1E-6,1E-6,1E-6,1E-6]# Ha
RATIO = [0., 0., 0., 0., 0.]
AST = [1E-8, 1E-8,1E-8,1E-8,1E-8]# Ha
amp = [1,4,7,10]  # Ha
T= [1000,3000,5000,7000,10000] # a.u.

CHI3R,CHI3I = [],[]
dt = 0.05
for t in range(len(T)):
    I, n_real, n_imag = [], [], []
    for i in range(len(amp)):
        I_a, n_real_a, n_imag_a = N_Cal(Amplitude = amp[i], pulse_duration = T[t], index = i, dt=dt,
                                        ratio=RATIO[t], A_standard = AST[t] * amp[i])
        I.append(I_a)
        n_real.append(n_real_a)
        n_imag.append(n_imag_a)
    I = np.array(I)
        #######   线性拟合 #######
    order =1
    fitted_yr, pr = ch.nfit(I*4, n_real, order) #cx.E2I 里面是I= 1/2 xxx， 文献里是 I= 2 xxx， 所以有个4的倍数关系
    fitted_yi, pi = ch.nfit(I*4, n_imag, order)
    # n0 = pr[0]
    n0 = 0.16629
    n2 = -1 * np.abs(pr[1]*1E-4)  # [m^2/W]
    # n2 = pr[1] * 1E-4  # [m^2/W]
    # k0 =  pi[0]
    k0 = 3.4652
    k2 = +1 * np.abs(pi[1]*1E-4) # [m^2/W]
    # k2 = pi[1] * 1E-4  # [m^2/W]
    beta = k2 / 300E-9  # [m/W]
    chi3r = ( 4 * n0 * cx.Vaccum_permitivity  * cx.Light_speed / 3) * (n0*n2 - k0*k2)  # [m^2/V^2]
    chi3i = ( 4 * n0 * cx.Vaccum_permitivity  * cx.Light_speed / 3) * (n0*k2 + k0*n2)  # [m^2/V^2]

    CHI3R.append(chi3r)
    CHI3I.append(chi3i)

CHI3R = np.array(CHI3R)
CHI3I = np.array(CHI3I)
Ts= np.array(T) * cx.timeau2s * 1e12 *  2/5 # ps  # 半高全宽


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

# for i in range(len(CHI3I)):
#     print("%.2e" % CHI3I[i])