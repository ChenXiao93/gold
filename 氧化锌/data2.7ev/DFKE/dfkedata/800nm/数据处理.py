import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx

########  图片设定
fs = 14

fig0 = plt.figure(figsize=(4, 6))  ##  eps 比较
ax01 = fig0.add_subplot(211)
ax02 = fig0.add_subplot(212)

fig1 = plt.figure(figsize=(6, 4))  ##  吸收系数比较
ax1 = fig1.add_subplot(111)

fig2 = plt.figure(figsize=(6, 4))  ##  吸收系数比较
ax2 = fig2.add_subplot(111)

start = 1  # eV
end = 10 # eV
w_ev = np.arange(start, end + 0.01, 0.1)  # eV
wavelen_nm = cx.Energy2Wavelength(w_ev)  # nm
t0 = 430  ## 800nm case, delay 430, 总传播时间600，  kick传播170
etaev_probe = 0.1 # ev
################# 基态的kick
#### kick 方法计算的介电函数
t_gs, dAinddt_gs = cx.ReadData("../linear_spectrum/td.general/gauge_field", 1, 7)

etaev_gs = 0.1 # eV
damping_gs = cx.DampingMethod(t = t_gs, eta =etaev_gs,  method = "exponential")
dAinddt_gs_damp = dAinddt_gs * damping_gs
eps_gs =cx.gauge_kick_dAdt(t= t_gs, A0 = 0.1, dAdt = dAinddt_gs_damp, w_ev = w_ev)
alpha_gs = cx.absorption_ev(w_ev, eps_gs)

ax01.plot(w_ev, eps_gs.real)
ax02.plot(w_ev, eps_gs.imag, label="gs")

ax1.plot(w_ev, alpha_gs, label="gs")
################# 激发态的kick
for amp in [5,10,15,20,25,30,55,90]:
    time, dAinddt_ref = cx.ReadData("./A%s/td.ref/gauge_field" % (amp), 1, 7)  # kick probe + pump
    dAinddt_general = cx.ReadOne("./A%s/td.general/gauge_field" % (amp), 7)  # pump
    dAinddt_kick = dAinddt_general - dAinddt_ref

    damping = cx.RatioDampingMethod(t=time, t0=t0, eta=etaev_probe, method="exponential")
    dAinddt_kick_damp = dAinddt_kick * damping

    # 修正dAinddt_kick_damp
    ax2.plot(time,dAinddt_kick_damp)

    eps_probe = cx.gauge_kick_delay(t=time, A0=0.1, dAdt=dAinddt_kick_damp, w_ev=w_ev, t0=t0)

    ax01.plot(w_ev, eps_probe.real)
    ax02.plot(w_ev, eps_probe.imag, label="A=%s" % amp)

    alpha_probe = cx.absorption_ev(w_ev, eps_probe)
    ax1.plot(w_ev, alpha_probe, label="A=%s" % amp)

ax02.legend(fontsize = fs-2)
ax02.set_xlim([start,end])
ax01.set_xlim([start,end])
ax02.set_xlabel("$\hbar \omega$ [eV]", fontsize  =fs)
ax01.set_ylabel("$\epsilon_1$", fontsize = fs)
ax02.set_ylabel("$\epsilon_2$", fontsize = fs)

ax1.legend(fontsize = fs-2)
ax1.set_xlim([start,end])
ax1.set_xlabel("$\hbar \omega$ [eV]", fontsize  =fs)
ax1.set_ylabel("$\\alpha$ [1/nm]", fontsize = fs)

plt.show()