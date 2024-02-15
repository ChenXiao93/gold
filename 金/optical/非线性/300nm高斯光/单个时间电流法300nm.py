import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.Laser as la
####################
def current_mono(t, E, J, w_ev, V, damp): # 电流法 激光单色场
    dt = t[2] - t[1]
    w_au = w_ev / cx.hartree2ev  # a.u.
    Ew, Jw = 0,0 # initialize for Fourier Transform
    for i in range(len(E)):
        Ew += E[i] * dt * np.exp(w_au * t[i] * 1j)  # 外场无需加修正
        Jw += J[i] * dt * np.exp(w_au * t[i] * 1j) * damp[i] # 感应场杂峰多，加修正
    sigma = Jw / (Ew * V)
    return(sigma)
##########################
fig1 = plt.figure()
fs = 12
ax1 = fig1.add_subplot(211)
ax2 = fig1.add_subplot(212)
####################
wavelength = 300 # nm
omega = cx.Wavelength2Energy(wavelength) # eV
omega_au = omega / 27.2114 # a.u.
amp = [10] # Ha
T= [500] # a.u.

print("eps1 = -1.237  eps2 = 6.031 (exp)")
print("T = %s a.u." % (T[0]))

for t in range(len(T)):
    for i in range(len(amp)):
        ######  A_ext 外场 and A_ind 感应场 ###############
        time, A_ext = cx.ReadData("./laser/t%s/A%s" % (T[t], amp[i]), 0, 1)
        E_ext = cx.A2Eau(time, A_ext) # *(-1)  # a.u.
        J_ind = cx.ReadOne("./data/t%s/A%s/total_current" % (T[t], amp[i]), 2)
        # 把ext场和ind场调整到长度一致
        dt = time[2] - time[1]
        start = 0
        end = int(T[t] / dt)
        time = time[start:end]
        E_ext = E_ext[start:end]
        J_ind = J_ind[start:end]
        ######
        ax1.plot(time, E_ext, label = "$E_{ext}$")
        ax2.plot(time, J_ind, label = "$J_{ind}$")
        #### env
        t0 = int(T[t] / 2)
        PulseDuration = int(T[t] * 2 / 5)
        env = np.exp(-np.square((time - t0) / (PulseDuration / 2)))
        ax1.plot(time, env  * max(E_ext))

        ## 衰减函数配置
        damp = env
        # damp = np.ones_like(time)   #  无衰减
        # damp = la.Env_exponential(time_ext, ratio = 0.2, eta_ev = 0.3)
        # ax2.plot(time, damp * max(J_ind), label="$envelope$")
        ax2.plot(time, env * max(J_ind))
        ax2.plot(time, J_ind * damp, label="$J_{mod}$")
        ##### 电流法
        sigma = current_mono(time, E_ext, J_ind, omega, V = 114.5822, damp = damp)
        sigma = cx.positivecomplex(sigma) # 强令sigma为正！
        epsilon = 1 +  4 * np.pi * 1j * sigma / (omega_au)
        print("eps1=%.3f, eps2=%.3f, sigma1=%.3f, sigma2=%.3f.  A= %s  (current)"
            %(epsilon.real, epsilon.imag, sigma.real, sigma.imag, amp[i]))




ax1.legend(fontsize = fs)
ax2.legend(fontsize = fs)
plt.show()
