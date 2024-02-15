# Smearing = 0.1*eV
# SmearingFunction = fermi_dirac
# RestartFixedOccupations = yes
# ExtraStates = 4
# a = 4.08*angstrom
# nk=10，15，20  对比不同k值
# prop = 500
import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op

exp_path = ["../exp/JC.csv"]  # 实验路径
# path = ["k10", "k15","k20"]   # 模拟路径
path = ["k20"]
fs = 12
plt.figure()
# 画出tddft结果
start, end = 1, 6
Omega = np.arange(start, end, 0.01)  # eV
Wavelength = cx.Energy2Wavelength(Omega) #  Energy [eV]  to Wavelength [nm]
A0 = 1
damp = 0.05 # eV

for i in path:
        time = cx.ReadOne("./%s" % i, 1)
        Aprobe= cx.ReadOne("./%s" % i, 2)
        Eprobe = -1*cx.ReadOne("./%s" % i, 5) / cx.c_au # E= -1/c  dA/dt  原子单位c=137
        epsilon = op.epsilon_kick(prop_time = time, kick_delay = 0, kick_amplitude=A0,
                induced_efield = Eprobe, damping_factor = damp, frequency_range = Omega)

        plt.subplot(2,1, 1)
        plt.plot(Wavelength, epsilon.imag, "-", label= i )
        plt.ylabel("$Im(\epsilon)$", fontsize=fs)

        plt.subplot(2,1,2)
        plt.plot(Wavelength, epsilon.real, "-", label= i )
        plt.ylabel("$Re(\epsilon)$", fontsize = fs)
        plt.xlabel("$Wavelength$ [nm]", fontsize=fs)

# 画出实验参考值
f_exp = np.loadtxt(exp_path[0], delimiter=",", skiprows=1, encoding="UTF-8")
eps1_exp, eps2_exp = op.n2eps(f_exp[:,1],f_exp[:,2])
Omega_exp = f_exp[:,0] # eV
Omega_exp_nm = cx.Energy2Wavelength(Omega_exp) # [nm]

plt.subplot(2,1, 1)
plt.plot(Omega_exp_nm, eps2_exp, "--b", label="exp JC")
plt.ylim([0,10])
plt.xlim([200,1200])

plt.subplot(2,1, 2)
plt.plot(Omega_exp_nm, eps1_exp, "--r", label="exp JC")
plt.legend(fontsize=fs)
plt.ylim([-60,5])
plt.xlim([200,1200])


plt.show()
