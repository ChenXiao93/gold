# Smearing = 0.1*eV
# SmearingFunction = fermi_dirac
# RestartFixedOccupations = yes
# ExtraStates = 4
# a = 4.08*angstrom
# nk= 20
# prop = 500
import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op

exp_path = ["../exp/JC.csv"]  # 实验路径
path = np.arange(0.05,0.55,0.05)   # 模拟路径

fs = 12
plt.figure()
# 画出tddft结果
start, end = 1.2, 3.1
Omega = np.arange(start, end, 0.01)  # eV
Wavelength = cx.Energy2Wavelength(Omega) #  Energy [eV]  to Wavelength [nm]
Omega_au = Omega / 27.2114  # atom unit
A0 = 1
for i in path:
        time = cx.ReadOne("{:.2f}".format(i), 1)
        Aprobe= cx.ReadOne("{:.2f}".format(i), 2)
        Eprobe = -1*cx.ReadOne("{:.2f}".format(i), 5) / cx.c_au # E= -1/c  dA/dt  原子单位c=137
        Eind_si = Eprobe * 51.44 # 1 [Ha] = 27.2114 [eV], 1 [Bohr] = 0.529 [Angstrom]   27.2114/0.529 =51.44
        fdamp = op.damping(prop_time = time,  kick_delay = 0, damping_factor = 0.3)
        epsilon = op.epsilon_kick(prop_time = time, kick_delay = 0, kick_amplitude=A0,
                induced_efield = Eprobe, damping_factor = 0.06, frequency_range = Omega_au)

        plt.subplot(2,1, 1)
        plt.plot(Wavelength, epsilon.imag, "-", label= "{:.2f}".format(i) )
        plt.ylabel("$Im(\epsilon)$", fontsize=fs)

        plt.subplot(2,1,2)
        plt.plot(Wavelength, epsilon.real, "-", label= "{:.2f}".format(i) )
        plt.ylabel("$Re(\epsilon)$", fontsize = fs)
        plt.xlabel("$Wavelength$ [nm]", fontsize=fs)

# 画出实验参考值
f_exp = np.loadtxt(exp_path[0], delimiter=",", skiprows=1, encoding="UTF-8")
eps1_exp, eps2_exp = op.n2eps(f_exp[:,1],f_exp[:,2])
Omega_exp = f_exp[:,0] # eV
Omega_exp_nm = cx.Energy2Wavelength(Omega_exp) # [nm]

plt.subplot(2,1, 1)
plt.plot(Omega_exp_nm, eps2_exp, "--b", label="exp JC")
plt.xlim([400,1000])
plt.ylim([0,10])

plt.subplot(2,1, 2)
plt.plot(Omega_exp_nm, eps1_exp, "--r", label="exp JC")
plt.legend(fontsize=fs-2)
plt.xlim([400,1000])
plt.ylim([-60,0])


plt.show()
