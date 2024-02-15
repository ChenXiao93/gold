import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx

########  图片设定
fs = 14


fig0 = plt.figure(figsize= (6,9))  ## E_ext 和 E_ind
grid = plt.GridSpec(8, 1, wspace=0, hspace=1)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
ax01 = fig0.add_subplot(grid[0:2, 0:1])
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
ax02 = fig0.add_subplot(grid[2:4, 0:1])
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
ax1 = fig0.add_subplot(grid[4:8, 0:1])
#############  参数设定
start = 0 # eV
end = 6 # eV
w_ev = np.arange(start, end + 0.01, 0.01)  # eV
wavelen_nm = cx.Energy2Wavelength(w_ev)  # nm
t0 = 1200
wavelength = 3500 # nm
energy = cx.Wavelength2Energy(wavelength) # eV
omega_au = energy / cx.hartree2ev  # a.u.
eta = 0.3
etaev_probe = eta# ev
etaev_gs = eta # eV
################################
time, dAinddt_ref = cx.ReadData("./td.ref/gauge_field", 1, 7)  #  pump
dAinddt_general = cx.ReadOne("./td.general/gauge_field", 7)  # pump + probe
E_general =  (-1/cx.c_au) * dAinddt_general *  51.44  # 51.44= 27.2114/0.529  # V / Angstrom
E_ref = (-1/cx.c_au) * dAinddt_ref *  51.44
E_diff =  E_general - E_ref

ax01.plot(time * 0.02419, E_general,color="blue", label = "Pump+Probe")
ax01.set_ylabel("$E_{ind} \ [V/\mathring{A}]$", fontsize=fs)
ax01.legend(loc="lower left", fontsize=fs)
ax01.set_xlim([0,36.3])
ax01.set_xlabel("Time [fs]", fontsize=fs)
ax01.set_yticks([-0.5,0,0.5],[-0.5,0,0.5],fontsize = fs-2)

ax02.plot(time * 0.02419, E_diff, color="blue", label = "(Pump+probe) - Pump")
ax02.set_ylabel("$E_{ind} \ [V/\mathring{A}]$", fontsize=fs)
ax02.set_yticks([0,0.01],[0,0.01],fontsize = fs-2)
ax02.set_xlabel("Time [fs]", fontsize=fs)
ax02.legend(loc="upper left", fontsize=fs)
ax02.set_xlim([0,36.3])
ax02.set_xticks([0,10,20,30],[0,10,20,30],fontsize = fs-2)
ax01.set_xticks([0,10,20,30],[0,10,20,30],fontsize = fs-2)
################# 基态的kick方法计算的介电函数
t_gs, dAinddt_gs = cx.ReadData("/Users/chenxiao/Project/paper/氧化锌/data2.7ev/DFKE/dfkedata/linear_spectrum/td.general/gauge_field", 1, 7)
damping_gs = cx.DampingMethod(t = t_gs, eta =etaev_gs,  method = "gaussian")
dAinddt_gs_damp = dAinddt_gs * damping_gs
eps_gs =cx.gauge_kick_dAdt(t= t_gs, A0 = 0.1, dAdt = dAinddt_gs_damp, w_ev = w_ev)
alpha_gs = cx.absorption_ev(w_ev, eps_gs)
################# 激发态的kick方法计算的介电函数
for amp in [143]:
    time, dAinddt_ref = cx.ReadData("./td.ref/gauge_field", 1, 7)  # kick probe + pump
    dAinddt_general = cx.ReadOne("./td.general/gauge_field", 7)  # pump
    dAinddt_kick = dAinddt_general - dAinddt_ref  # (pump+probe) - pump

    damping = cx.RatioDampingMethod(t=time, t0=t0, eta=etaev_probe, method="gaussian")
    dAinddt_kick_damp = dAinddt_kick * damping

    eps_probe = cx.gauge_kick_delay(t=time, A0=0.1, dAdt=dAinddt_kick_damp, w_ev=w_ev, t0=t0)
    alpha_probe = cx.absorption_ev(w_ev, eps_probe)

    ax1.plot(w_ev, (alpha_probe)*1000, "--",color="green",label = "excited")
    ax1.plot(w_ev, (alpha_gs) * 1000,"--",color="red",  label = "gs")
    ax11 = ax1.twinx()
    ax11.plot(w_ev, (alpha_probe - alpha_gs) * 1000,"-", color="blue", label = "diff")
    ax11.set_yticks([-1, 0, 1], [-1, 0, 1], fontsize=fs - 2)
    xx = np.arange(start,end+0.1, 0.1)
    yy = np.zeros_like(xx)
    ax1.plot(xx,yy,"--", color="blue")


#################




ax1.legend(loc="upper left", fontsize=fs)
ax11.legend(loc="lower right", fontsize=fs)
ax11.set_ylim([-1.6, 1.6])
ax1.set_ylim([-20,20])
ax1.set_xlim([start,end])
ax1.set_yticks([-20,-10,0,10,20],[-20,-10,0,10,20],fontsize = fs-2)
ax1.set_xticks([1,2,3,4,5,6],[1,2,3,4,5,6],fontsize = fs-2)
ax1.set_xlabel("$\hbar \omega$ [eV]", fontsize  =fs)
ax1.set_ylabel("$\\alpha \ [1 / \mu m]$ ", fontsize = fs)
plt.savefig("dfke.eps", dpi = 900, bbox_inches='tight')
plt.show()