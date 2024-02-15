import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx

########  图片设定
fs = 14
# fig0 = plt.figure(figsize=(6,4))  ##  eps_2 比较
# ax2 = fig0.add_subplot(111)

fig1 = plt.figure(figsize=(6, 4))  ##  吸收系数比较
plt.rcParams['xtick.direction'] = 'in'
ax1 = fig1.add_subplot(111)
#############  参数设定
start = 0.5 # eV
end = 4 # eV
w_ev = np.arange(start, end + 0.01, 0.01)  # eV
wavelen_nm = cx.Energy2Wavelength(w_ev)  # nm
t0 = 1200
wavelength = 3500 # nm
energy = cx.Wavelength2Energy(wavelength) # eV
omega_au = energy / cx.hartree2ev  # a.u.
eta = 0.3
etaev_probe = eta# ev
etaev_gs = eta # eV
################# 基态的kick方法计算的介电函数
t_gs, dAinddt_gs = cx.ReadData("/Users/chenxiao/Project/paper/氧化锌/data2.7ev/DFKE/dfkedata/linear_spectrum/td.general/gauge_field", 1, 7)
damping_gs = cx.DampingMethod(t = t_gs, eta =etaev_gs,  method = "exponential")
dAinddt_gs_damp = dAinddt_gs * damping_gs
eps_gs =cx.gauge_kick_dAdt(t= t_gs, A0 = 0.1, dAdt = dAinddt_gs_damp, w_ev = w_ev)
alpha_gs = cx.absorption_ev(w_ev, eps_gs)
################# 激发态的kick方法计算的介电函数
for amp in [143]:
    time, dAinddt_ref = cx.ReadData("./td.ref/gauge_field", 1, 7)  # kick probe + pump
    dAinddt_general = cx.ReadOne("./td.general/gauge_field", 7)  # pump
    dAinddt_kick = dAinddt_general - dAinddt_ref  # (pump+probe) - pump

    damping = cx.RatioDampingMethod(t=time, t0=t0, eta=etaev_probe, method="exponential")
    dAinddt_kick_damp = dAinddt_kick * damping

    eps_probe = cx.gauge_kick_delay(t=time, A0=0.1, dAdt=dAinddt_kick_damp, w_ev=w_ev, t0=t0)
    alpha_probe = cx.absorption_ev(w_ev, eps_probe)


    # ax2.plot(w_ev, eps_probe.imag,"--", label = "excited")
    # ax2.plot(w_ev, eps_gs.imag, "--",label = "gs")
    # ax21 = ax2.twinx()
    # ax21.plot(w_ev, (eps_probe.imag - eps_gs.imag), label = "diff")


    ax1.plot(w_ev, (alpha_probe)*1000, "--",color="orange",label = "excited")
    ax1.plot(w_ev, (alpha_gs) * 1000,"--",color="red",  label = "gs")
    ax11 = ax1.twinx()
    ax11.plot(w_ev, (alpha_probe - alpha_gs) * 1000,"-", color="#1f77b4", label = "diff")

    xx = np.arange(start,end+0.1, 0.1)
    yy = np.zeros_like(xx)
    ax1.plot(xx,yy,"--", color="black")
# ax2.legend(fontsize = fs-2)
# ax2.set_xlim([start,end])
# ax2.set_xlabel("$\hbar \omega$ [eV]", fontsize  =fs)
# ax2.set_ylabel("$\epsilon_2$", fontsize = fs)

# ax1.legend(fontsize = fs-2)
# ax11.legend(fontsize = fs-2)
ax11.set_ylim([-1.6, 1.6])
ax1.set_ylim([-10,10])


ax1.set_xlim([start,end])
ax1.set_xlabel("$\hbar \omega$ [eV]", fontsize  =fs)
ax1.set_ylabel("$\\alpha \ [1/\mu m]$ ", fontsize = fs)

plt.show()