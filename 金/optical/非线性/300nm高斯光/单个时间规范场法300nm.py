import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.Laser as la


def gauge_mono(t, A_ext, A_ind, w_ev,damp):  # 电场法 激光单色场
    dt = t[2] - t[1]
    w_au = w_ev / cx.hartree2ev
    dAdt_ext = cx.dA_over_dt(t, A_ext)
    dAdt_ind = cx.dA_over_dt(t, A_ind)
    Extw, Indw = 0, 0  # initialize for Fourier Transform
    for i in range(len(A_ind)):
        Extw += dAdt_ext[i] * dt * np.exp(w_au * t[i] * 1j) * damp[i]
        Indw += dAdt_ind[i] * dt * np.exp(w_au * t[i] * 1j) * damp[i]
    inv_epsilon = 1 +  Indw / Extw
    epsilon = 1 / inv_epsilon
    return(epsilon)
####################
wavelength = 300 # nm
omega = cx.Wavelength2Energy(wavelength) # eV
amp = [1,4,7,10] # Ha
T= [500, 1500,3000,5000,10000] # a.u.
RATIO = [0.2, 0.5, 0.5, 0.5, 0.58]
STD = [1E-8,1E-8,1E-8,1E-8,1E-8] # Ha


for t in range(len(T)):
    print("#####################")
    print("eps1 = -1.237  eps2 = 6.031 (exp)")
    print("T = %s a.u." % T[t])
    for i in range(len(amp)):
        ######  A_ext 外场 and A_ind 感应场 ###############
        time, A_ext = cx.ReadData("./laser/t%s/A%s" % (T[t], amp[i]), 0, 1)
        A_ind = cx.ReadOne("./data/t%s/A%s/gauge_field" % (T[t], amp[i]), 2)
        # 把ext场和ind场调整到长度一致
        dt = time[2] - time[1]
        start = 0
        end = int(T[t] / dt)
        time = time[start:end]
        A_ext = A_ext[start:end]
        A_ind = A_ind[start:end]
        ## 衰减函数配置
        t0 = time[-1] * (1 - RATIO[i])
        standard = STD[t] * amp[i]
        # eta_au = np.abs(np.log(np.abs( standard / A_ind[-1])) / (time[-1] - t0))  # 指数包络
        eta_au =  np.sqrt(np.abs(np.log( np.abs( standard / A_ind[-1]))))  / (2 * (time[-1]-t0))
        eta_ev =  eta_au * cx.hartree2ev
        envy = la.Env_Gaussian(time, RATIO[t], eta_ev)  # 高斯包络
        # envy = la.Env_exponential(time, RATIO[t], eta_ev)  # 指数包络
        plt.plot(time,envy)
        ###### 电场法
        eps =  gauge_mono(time, A_ext, A_ind, omega,damp = envy)
        print("eps1_a=%.4f,eps2_a=%.4f, A =%s Ha, eta = %.2feV" %(eps.real, eps.imag, amp[i], eta_ev))

plt.show()
