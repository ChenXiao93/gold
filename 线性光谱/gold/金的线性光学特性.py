import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op

##### 实验值读取对比
## JC1972
f_exp = np.loadtxt("./exp/JC-1972.csv", delimiter=",", skiprows=1, encoding="UTF-8")
eps1_exp, eps2_exp = op.n2eps(f_exp[:,1], f_exp[:,2])
Omega_exp = f_exp[:,0] # eV
Omega_exp_nm = cx.Energy2Wavelength(Omega_exp) # [nm]
#### Werner-2009DFT
f_wernerdft = np.loadtxt("./exp/Werner-2009DFT.csv", delimiter=",", encoding="UTF-8")
nm_wernerdft = f_wernerdft[:,0]*1e3 # micrometer to nanometer
eps1_wernerdft, eps2_wernerdft = op.n2eps(f_wernerdft[:,1], f_wernerdft[:,2])
#### Mag2019
f_mag = np.loadtxt("./exp/Magnozzi-25C-2019.csv", delimiter=",", encoding="UTF-8")
nm_mag = f_mag[:,0]*1e3 # micrometer to nanometer
eps1_mag, eps2_mag = op.n2eps(f_mag[:,1], f_mag[:,2])
#### Olmon-sc-2012
f_ol = np.loadtxt("./exp/Olmon-sc-2012.csv", delimiter=",", encoding="UTF-8")
nm_ol = f_ol[:,0]*1e3 # micrometer to nanometer
eps1_ol, eps2_ol = op.n2eps(f_ol[:,1], f_ol[:,2])
########设置画布
fig = plt.figure()
grid = plt.GridSpec(2, 1, wspace=0, hspace=0)
ax2 = fig.add_subplot(grid[0:1, 0:1])
ax1 = fig.add_subplot(grid[1:2, 0:1])
fs = 12
### 其它参数
V = 114.5822 # Bohr^3
A0  = 0.1
wavelength = np.arange(200,1000,1) # nm
energy_range_ev = cx.Wavelength2Energy(wavelength) # eV
energy_range_au = energy_range_ev / cx.hartree2ev # au
eta_ev_a = 0.15 # eV
eta_ev_j = 0.15 # eV
eta_au_a = eta_ev_a / cx.hartree2ev  # au
eta_au_j = eta_ev_j / cx.hartree2ev  # au
###################### 计算代码  ######################################
t,A,dAdt2 = cx.ReadData_xyz("./Akick/gauge_field",  1,2,8)
tj,j = cx.ReadData("./Astep/total_current",  1,2)
damp_a = cx.DampingMethod(t, eta_ev_a, method="gaussian")
damp_j = cx.DampingMethod(tj, eta_ev_j, method="gaussian")
plt.figure()
plt.plot(t,damp_a)
### 调用计算函数！
eps_j = cx.current_step(tj, A0, j,V, energy_range_ev, damp_j)  ## 电流法
eps_a = cx.gauge_kick(t, A0, A, energy_range_ev, damp_a) ## 规范场法
############## 画出计算值
## 电流法
ax1.plot(wavelength, eps_j.real, label = "$current$")
ax2.plot(wavelength, eps_j.imag, label = "$current$")
## 规范场法
ax1.plot(wavelength, eps_a.real,label = "$gauge$")
ax2.plot(wavelength, eps_a.imag,label = "$gauge$")
## 画出实验值
ax1.plot(nm_wernerdft, eps1_wernerdft, "k--", label="DFT(RPA):Werner-2009") # 实部
ax1.plot(Omega_exp_nm, eps1_exp, "g--", label="Exp:J&C-1972")
ax1.plot(nm_ol, eps1_ol, "r--", label="Exp:Olmon-2012")
ax1.plot(nm_mag, eps1_mag, "c--", label="Exp:Magnozzi-2019")

ax2.plot(nm_wernerdft, eps2_wernerdft, "k--", label="DFT(RPA):Werner-2009") # 虚部
ax2.plot(Omega_exp_nm, eps2_exp, "g--", label="Exp:J&C-1972")
ax2.plot(nm_ol, eps2_ol, "r--", label="Exp:Olmon-2012")
ax2.plot(nm_mag, eps2_mag, "c--", label="Exp:Magnozzi-2019")
########  其它做图参数
ax1.legend(fontsize = fs-2) # 标签
ax2.set_xticks([])
ax1.set_xlim([200, 800]) # 横坐标范围
ax2.set_xlim([200, 800])
ax2.set_ylim([0,10]) # 纵坐标范围
ax1.set_ylim([-27,2])
ax1.set_ylabel("$Re(\epsilon)$", fontsize = fs) # 轴标识
ax2.set_ylabel("$Im(\epsilon)$", fontsize = fs)
ax1.set_xlabel("$\lambda [nm]$", fontsize = fs)

plt.show()