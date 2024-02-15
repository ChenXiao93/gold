import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx



########  图片设定
fs = 14

fig0 = plt.figure(figsize=(6, 4))  ##  eps 比较
ax01 = fig0.add_subplot(211)
ax02 = fig0.add_subplot(212)

fig1 = plt.figure(figsize=(6, 4))  ##  dA/dt 和 damping 修正
ax1 = fig1.add_subplot(111)

start =0.1 # eV
end =2 # eV
w_ev = np.arange(start,end + 0.01,0.01) # eV
wavelen_nm = cx.Energy2Wavelength(w_ev) # nm
########
#### kick 方法计算的介电函数
t_cal, dAinddt_cal = cx.ReadData("./td.general/gauge_field", 1, 7)

etaev = 0. # eV
damping = cx.DampingMethod(t = t_cal, eta =etaev,  method = "exponential")
dAinddt_cal_damp = dAinddt_cal * damping


##### dA/dt 和 damping 修正
ax1.plot(t_cal * 0.02419, dAinddt_cal, label = "dA/dt")
ax1.plot(t_cal * 0.02419, damping * max(dAinddt_cal), label = "$exp(-\eta t)$, $\eta$=%.1feV" % etaev)
ax1.plot(t_cal * 0.02419, dAinddt_cal_damp, label = "dA/dt * $exp(-\eta t)$")

ax1.set_xlabel("t  [fs]", fontsize = fs)
ax1.set_ylabel(" dA/dt [a.u.]", fontsize = fs )
ax1.legend(fontsize =fs)

#### eps
eps_cal_linear =cx.gauge_kick_dAdt(t= t_cal, A0 = 0.1, dAdt = dAinddt_cal_damp, w_ev = w_ev)
### exp 实验值
wavelength_exp1, n_exp1, k_exp1 = cx.ReadData_xyz("./exp/Aguilar.txt", 0,1,2) #  微米
w_ev_exp1 = cx.Wavelength2Energy(wavelength_exp1 * 1E3) # eV
eps1_exp1, eps2_exp1 = cx.n2eps(n_exp1, k_exp1) # 无量纲

wavelength_exp2, n_exp2, k_exp2 = cx.ReadData_xyz("./exp/Stelling.txt", 0,1,2) #  微米
w_ev_exp2 = cx.Wavelength2Energy(wavelength_exp2 * 1E3) # eV
eps1_exp2, eps2_exp2 = cx.n2eps(n_exp2, k_exp2) # 无量纲
####
ax01.plot(w_ev, eps_cal_linear.real,"-")
ax01.plot(w_ev_exp1, eps1_exp1, "--")
ax01.plot(w_ev_exp2, eps1_exp2, "-.")

ax02.plot(w_ev, eps_cal_linear.imag,"-" ,label = "kick")
ax02.plot(w_ev_exp1, eps2_exp1,"--", label = "Exp:Aguilar")
ax02.plot(w_ev_exp2, eps2_exp2,"-.", label = "Exp:Stelling")

ax02.legend(loc = "upper left", fontsize = fs-2)
ax02.set_xlim([start,end])
ax01.set_xlim([start,end])
ax01.set_ylim([2.2, 4])
ax02.set_ylim([-0.05, 1])
ax02.set_xlabel("$\hbar \omega$ [eV]", fontsize  =fs)
ax01.set_ylabel("$\epsilon_1$", fontsize = fs)
ax02.set_ylabel("$\epsilon_2$", fontsize = fs)
plt.show()