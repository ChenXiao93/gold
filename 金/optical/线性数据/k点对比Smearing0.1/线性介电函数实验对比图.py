import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op

def epsilon_kick_gaussian(prop_time, kick_amplitude, induced_efield, damping_factor, frequency_range):
    frequency_range = frequency_range / 27.2114  # 从 eV 改到 原子单位
    dt = prop_time[2] - prop_time[1]
    f_damp = np.exp(- np.square( prop_time / damping_factor))
    dAdt = induced_efield * (-cx.c_au) * f_damp
    epsilon = np.zeros_like(frequency_range, dtype=complex)
    for w in range(len(frequency_range)):
        real, imag = 0, 0  # initialize the real and imag of FT
        for i in range(len(prop_time)):
            real += dAdt[i] * dt * np.cos(frequency_range[w] * prop_time[i])
            imag += dAdt[i] * dt * np.sin(frequency_range[w] * prop_time[i])
        eps_inv = 1 + (real + 1j * imag) / kick_amplitude
        epsilon[w] = 1 / eps_inv
    return(epsilon)

fig1 = plt.figure()
ax11 = fig1.add_subplot(211)
ax12 = fig1.add_subplot(212)

fig2 = plt.figure()
ax21 = fig2.add_subplot(211)
ax22 = fig2.add_subplot(212)


fs = 12
start, end = 1, 6
Omega = np.arange(start, end, 0.01)  # eV
Wavelength = cx.Energy2Wavelength(Omega) #  Energy [eV]  to Wavelength [nm]
A0 = 1
f_exp = np.loadtxt("../exp/JC.csv", delimiter=",", skiprows=1, encoding="UTF-8")
eps1_exp, eps2_exp = op.n2eps(f_exp[:,1],f_exp[:,2])
Omega_exp = f_exp[:,0] # eV
Omega_exp_nm = cx.Energy2Wavelength(Omega_exp) # [nm]

time = cx.ReadOne("./k20", 1)
Aprobe= cx.ReadOne("./k20", 2)
Eprobe = -1*cx.ReadOne("./k20", 5) / cx.c_au # E= -1/c  dA/dt  原子单位c=137
damp1 = 0.15 # eV
damp2ev = 0.15 # eV
damp2 = 27.2114 / damp2ev   # au time


env1 =  np.exp(- damp1 / 27.2114 * time) * max(Aprobe)
env2 =  np.exp(- np.square( time / damp2)) * max(Aprobe)
Aprobe_damp1 = Aprobe * np.exp(- damp1 / 27.2114 * time)
Aprobe_damp2 = Aprobe * np.exp(- np.square( time / (damp2)))

ax11.plot(time,Aprobe, label = "Origin")
ax11.plot(time,Aprobe_damp1, label = "Modified: $\eta = %seV$" % damp1)
ax11.plot(time,env1, label = "Envelope: $exp(-\eta t)$")
ax11.legend(fontsize=fs)
ax11.set_ylabel("$A$ [Ha]", fontsize=fs)

ax12.plot(time,Aprobe, label = "ori")
ax12.plot(time,Aprobe_damp2, label = "$T = %i au = %.2feV$" % (damp2, damp2ev))
ax12.plot(time,env2, label = "Envelope: $exp(-(\\frac{t}{T})^2)$")
ax12.legend(fontsize=fs)
ax12.set_xlabel("$t$ [au]", fontsize=fs)


epsilon1 = op.epsilon_kick(prop_time = time, kick_delay = 0, kick_amplitude=A0,
        induced_efield = Eprobe, damping_factor = damp1, frequency_range = Omega)
epsilon2 = epsilon_kick_gaussian(prop_time = time, kick_amplitude=A0,
        induced_efield = Eprobe, damping_factor = damp2, frequency_range = Omega)






ax21.plot(Wavelength, epsilon1.imag, "b-", label= "Exponential" )
ax21.plot(Wavelength, epsilon2.imag, "r-", label="Gaussian")
ax21.plot(Omega_exp_nm, eps2_exp, "--",color = "orange", label="exp JC")
ax21.legend(fontsize=fs)
ax21.set_ylabel("$Im(\epsilon)$", fontsize=fs)
ax21.set_xlim([200,1000])
ax21.set_ylim([0,20])

ax22.plot(Wavelength, epsilon1.real, "b-", label= "Exponential" )
ax22.plot(Wavelength, epsilon2.real, "r-", label="Gaussian" )
ax22.plot(Omega_exp_nm, eps1_exp, "--",color = "orange", label="exp JC")
ax22.legend(fontsize=fs)
ax22.set_ylabel("$Re(\epsilon)$", fontsize = fs)
ax22.set_xlabel("$Wavelength$ [nm]", fontsize=fs)
ax22.set_xlim([200,1000])
ax22.set_ylim([-45,5])


plt.show()
