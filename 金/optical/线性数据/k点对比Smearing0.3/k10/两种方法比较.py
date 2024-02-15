# Smearing = 0.3*eV
# SmearingFunction = fermi_dirac
# RestartFixedOccupations = yes
# ExtraStates = 5
# a = 4.08*angstrom
# nk=10
# prop = 500
import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op

def current_kick(time, current,omega, kick_amplitude):
    # 输入参数:
    #  time  # [hbar/Ha]  传播时间
    #  current # [au]    octopus直接输出的感应电流  j(t)
    #  omega  # eV,   光子能量 \hbar \omega
    dt = time[2] - time[1]  # dt [hbar / Ha]
    omega_au = omega * cx.timeau2s / cx.hbar  # 原子单位换算 cx.timeau2s / cx.hbar = 1 / 27.2114
    epsilon1, epsilon2 = [],[]
    current = -1 * current   # octopus输出的电流J 需要乘 - 1 !
    for w in range(len(omega_au)):
        Jw1,Jw2 = 0,0
        for i in range(len(current)):
            Jw1 += current[i] * dt * np.cos(omega_au[w] * time[i])  # J_w [au]
            Jw2 += current[i] * dt * np.sin(omega_au[w] * time[i])
        sigma1 = Jw1  * cx.c_au / kick_amplitude # \sigma = J / E
        sigma2 = Jw2  * cx.c_au / kick_amplitude
        #从电导求介电函数eps, 这里使用了原子单位没有4pi, 但cgs单位制多个4pi
        eps1 = 1 -  sigma1 / omega_au[w]
        eps2 =   sigma2 / omega_au[w]
        epsilon1.append(eps1)
        epsilon2.append(eps2)
    epsilon1 = np.array(epsilon1)
    epsilon2 = np.array(epsilon2)

    return(epsilon1, epsilon2)




start, end = 1, 6  # eV
Omega = np.arange(start, end, 0.02)  # eV
Wavelength = cx.Energy2Wavelength(Omega) #  Energy [eV]  to Wavelength [nm]
A0 = 1   # Ha
## 规范场kick方法
damp = 0.0 # eV
time = cx.ReadOne("./gauge_field", 1)
epsilon = op.gauge_field_kick(prop_time = time, kick_delay = 0, kick_amplitude=A0,
        dAprobe_dt = cx.ReadOne("./gauge_field", 5), damping_factor = damp, frequency_range = Omega)
## 电流法kick方法
current = cx.ReadOne("./total_current", 2)
epsj1, epsj2 = current_kick(time, current,Omega, A0)


#### 画图
fs = 12
plt.figure()
plt.subplot(2, 1, 1)
plt.plot(Wavelength, epsilon.imag, "-", label="Gauge")
plt.plot(Wavelength, epsj2, "-", label="Current")
plt.subplot(2,1,2)
plt.plot(Wavelength, epsilon.real, "-", label="Gauge")
plt.plot(Wavelength, epsj1, "-", label="Current")

# 画出实验参考值
f_exp = np.loadtxt("../../exp/JC.csv", delimiter=",", skiprows=1, encoding="UTF-8")
eps1_exp, eps2_exp = op.n2eps(f_exp[:,1],f_exp[:,2])
Omega_exp = f_exp[:,0] # eV
Omega_exp_nm = cx.Energy2Wavelength(Omega_exp) # [nm]

plt.subplot(2,1, 1)
# plt.plot(Omega_exp, eps2_exp, "--b", label="exp JC")
plt.plot(Omega_exp_nm, eps2_exp, "--b", label="exp JC")
plt.xlim([200,1200])
plt.ylim([0,20])
# plt.xlim([1,6])
# plt.ylim([0,10])
plt.ylabel("$Im(\epsilon)$", fontsize=fs)

plt.subplot(2,1, 2)
# plt.plot(Omega_exp, eps1_exp, "--r", label="exp JC")
plt.plot(Omega_exp_nm, eps1_exp, "--b", label="exp JC")
plt.legend(fontsize=fs)
plt.xlim([200,1200])
plt.ylim([-70,20])
# plt.xlim([1,6])
# plt.ylim([-60,0])
plt.ylabel("$Re(\epsilon)$", fontsize = fs)
plt.xlabel("$Wavelength$ [nm]", fontsize=fs)

plt.show()
