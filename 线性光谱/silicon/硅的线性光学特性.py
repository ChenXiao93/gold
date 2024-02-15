import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op

##### 实验值读取对比
## JC1972
f_exp = np.loadtxt("./Aspnes.csv", delimiter=",", encoding="UTF-8")
eps1_exp, eps2_exp = op.n2eps(f_exp[:,1], f_exp[:,2])
wavelength_exp = f_exp[:,0] * 1000 # nm
########设置画布
fig = plt.figure()
grid = plt.GridSpec(2, 1, wspace=0, hspace=0)
ax2 = fig.add_subplot(grid[0:1, 0:1])
ax1 = fig.add_subplot(grid[1:2, 0:1])
fs = 12
### 其它参数
V = 263.7445 # Bohr^3
wavelength = np.arange(200,800,1) # nm
energy_range_ev = cx.Wavelength2Energy(wavelength) # eV
A0  = 0.1
######################  oct-conductivity 计算代码  ######################################
t,A,dAdt2 = cx.ReadData_xyz("./Akick/gauge_field",  1,2,8)
j = cx.ReadOne("./Astep/total_current", 2)
## 电流法
damp_j = cx.DampingMethod(t, eta = 0., method="gaussian")
eps_j = cx.current_step(t, A0, j,V, energy_range_ev, damp_j)
## 电场法
damp_a = cx.DampingMethod(t, eta = 0., method="gaussian")
eps_a = cx.gauge_kick(t, A0, A, energy_range_ev, damp_a)
############## 画出计算值
## 电流法
ax1.plot(wavelength, eps_j.real, label = "$current$")
ax2.plot(wavelength, eps_j.imag, label = "$current$")
# ## 规范场法
ax1.plot(wavelength, eps_a.real,label = "$gauge$")
ax2.plot(wavelength, eps_a.imag,label = "$gauge$")
## 画出实验值
ax1.plot(wavelength_exp, eps1_exp, "--", label="exp") # 实部
ax2.plot(wavelength_exp, eps2_exp, "--", label="exp") # 虚部
########  其它做图参数
ax1.legend(fontsize = fs-2) # 标签
ax2.set_xticks([])
ax1.set_xlim([200, 800]) # 横坐标范围
ax2.set_xlim([200, 800])
ax1.set_ylabel("$Re(\epsilon)$", fontsize = fs) # 轴标识
ax2.set_ylabel("$Im(\epsilon)$", fontsize = fs)
ax1.set_xlabel("$\lambda [nm]$", fontsize = fs)

plt.show()