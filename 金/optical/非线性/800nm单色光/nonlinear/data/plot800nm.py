import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op
import OctopusTools.Laser as la
import OctopusTools.chi3 as ch

def N_Cal(Amplitude, pulse_duration, index):
    time_ext, A_ext, envy = la.GaussianPulse(PropTime=pulse_duration + 0.05, PulseDuration=int(pulse_duration * 2 / 5),
                                    Wavelength=800, Amplitude=Amplitude,
                                    TimeStep=0.2,Phase=0)
    ######   A_ext 外场 and A_ind 感应场 ###############
    time_ind, A_ind = cx.ReadData("./t%s/A%s/gauge_field" % (pulse_duration,Amplitude), 1, 2)
    #### 加上自动化调控衰减！ ####
    eta_au = np.log( np.abs(A_ind[-1] / A_standard))  / time_ind[-1]  # atom unit
    time, envy1 = la.env_exp(mid=(pulse_duration*(1-ratio)), tail=(pulse_duration*ratio+0.05), dt=0.2, gamma=eta_au)
    A_ind = A_ind * envy1
    A_ext = A_ext * envy1
    ######## 求介电函数 ####################
    wavelength = 800  # nm
    Omega = cx.Wavelength2Energy(wavelength)
    epsilon_real, epsilon_imag = op.Get_Epsilon_Gauge(time_ext, A_ext, time_ind, A_ind, Omega, 0)
    ###### Refractive Index  求折射率和光强 ##################################
    n_real, n_imag = op.refractive_index_positive(epsilon_real, epsilon_imag)
    Emax = max(cx.A2E(time_ext, A_ext+A_ind))
    I = cx.E2I(E = Emax, Refractive_index = 1) * 4  # [W/cm^2]  # cited from the 'n0n2k0k2'
    return(I, n_real, n_imag)

ratio = 0.6
A_standard = 1E-8 # Ha
amp = [10,22,33,40]
T= [500,1500,4000,15000] # a.u.
CHI3R,CHI3I = [],[]
for time in T:
    I, n_real, n_imag = [], [], []
    for i in range(len(amp)):
        I_a, n_real_a, n_imag_a = N_Cal(Amplitude = amp[i], pulse_duration = time, index = i)
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
    beta = k2 / 300E-9  # [m/W]
    chi3r = ( 4 * n0 * cx.Vaccum_permitivity  * cx.Light_speed / 3) * (n0*n2 - k0*k2)  # [m^2/V^2]
    chi3i = ( 4 * n0 * cx.Vaccum_permitivity  * cx.Light_speed / 3) * (n0*k2 + k0*n2)  # [m^2/V^2]

    CHI3R.append(chi3r)
    CHI3I.append(chi3i)
CHI3R = np.array(CHI3R)
CHI3I = np.array(CHI3I)
Ts= np.array(T) * cx.timeau2s * 1e12  # ps

path="/Users/chenxiao/Project/paper/2023.5.gold/optical/linearity/300nm高斯光"
f630r = np.loadtxt(path+"/exp/Rotenberg630real.txt", delimiter=",")
f630i = np.loadtxt(path+"/exp/Rotenberg630imag.txt", delimiter=",")
f532 = np.loadtxt(path+"/exp/532nm.txt")
f600 = np.loadtxt(path+"/exp/600nm.txt")
f800 = np.loadtxt(path+"/exp/800nm.txt")
f1064= np.loadtxt(path+"/exp/1064nm.txt")

fig1 = plt.figure()
fs = 12
ax11 = fig1.add_subplot(121)
ax11.axes.loglog(Ts, np.abs(CHI3R),"o", label ="TDDFT300nm")
ax11.axes.loglog(f630r[:,0]*1e12 , np.abs(f630r[:,1]),"*", label= "Rotenberg630nm")
ax11.axes.loglog(f532[:,0], np.abs(f532[:,1]),"d", label= "Exp532nm")
ax11.axes.loglog(f600[0], np.abs(f600[1]),"p", label= "Exp600nm")
ax11.axes.loglog(f800[0], np.abs(f800[1]),"p", label= "Exp800nm")
ax11.axes.loglog(f1064[:,0], np.abs(f1064[:,1]),"p", label= "Exp1064nm")

ax11.set_xlabel("t [ps]", fontsize = fs)
ax11.set_ylabel("$|\chi^3_{real}| [m^2/V^2]$", fontsize = fs)
ax11.legend()


ax12 = fig1.add_subplot(122)
ax12.axes.loglog(Ts, np.abs(CHI3I),"o", label ="TDDFT300nm")
ax12.axes.loglog(f630r[:,0]*1e12 , np.abs(f630r[:,1]),"*", label= "Rotenberg630nm")
ax12.axes.loglog(f532[:,0], np.abs(f532[:,2]),"d", label= "Exp532nm")
ax12.axes.loglog(f600[0], np.abs(f600[2]),"p", label= "Exp600nm")

ax12.set_xlabel("t [ps]", fontsize = fs)
ax12.set_ylabel("$|\chi^3_{imag}| [m^2/V^2]$", fontsize = fs)
ax12.legend()

plt.show()