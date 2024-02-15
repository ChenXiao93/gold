import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
########  图片设定
fs = 14

fig0 = plt.figure(figsize= (6, 9))
grid = plt.GridSpec(100, 100, wspace=0, hspace=10)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
ax01 = fig0.add_subplot(grid[0:25, 0:100]) ### Eind:  pump
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
ax02 = fig0.add_subplot(grid[60:100, 0:100]) ### Eind:  (pump+probe) - pump
ax03 = fig0.add_subplot(grid[25:50, 0:100]) ### Eind:  (pump+probe) - pump
#######################################
delay_au = [0, 200,300,400,500, 1100, 1200]

dt = 0.05
start,end,dw = 3, 18, 0.1   # eV
w_ev = np.arange(start, end + dw, dw)  # eV
omega_au = w_ev / cx.hartree2ev  # a.u.
length = 500 # a.u.  最大500
eta = 0.3 # eV
method = "gaussian"     # "gaussian" # "exponential"
######## 纯 pump
text_pump, Aext_pump = cx.ReadData("./pump/IRpump",0,1)
Eext_pump = cx.A2E(text_pump, Aext_pump)
ax01.plot(text_pump, Eext_pump, "b", label = "pump") #  * cx.timeau2fs
ax01.legend(fontsize = fs-2)
ax01.set_xlim([0,1300])
ax01.set_ylim([-3,3])
ax01.set_ylabel("$E_{ind}$ [$V/\mathring{A}$]", fontsize = fs)
ax01.set_xticks([0,300,600,900,1200],[],fontsize = fs-2)
ax01.set_yticks([-2,0,2],[-2,0,2],fontsize = fs-2)

t_pump,dA_pump = cx.ReadData("./pump/gauge_field",1,5)
#######  gs 吸收谱 ########
t_probe, dA_probe = cx.ReadData("./gs/gauge_field", 1, 5)
t_probe = t_probe[0:int(length/0.05)] ###  截断
dA_probe = dA_probe[0:int(length/0.05)] ###  截断
damping_probe = cx.RatioDampingMethod(t=t_probe, t0=0, eta=eta, method=method) ###  衰减
dA_probe_damp = dA_probe * damping_probe ###  衰减
eps_gs = cx.gauge_kick_delay(t=t_probe, A0=0.22, dAdt=dA_probe_damp, w_ev=w_ev, t0=0) # 复数
alpha_gs = cx.absorption_ev(w_ev, eps_gs)  * 1000 # [1/mu m]
#######  不同延迟的 probe   ##################################
for i in range(len(delay_au)):
    t_delay, dA_dely = cx.ReadData("./%s/gauge_field" % delay_au[i],1, 5)
    t_delay = t_delay[0:int( (length +  delay_au[i]) /dt)]  ###  截断
    dA_dely = dA_dely[0:int( (length +  delay_au[i]) / dt)]  ###  截断  pump+probe
    dA_pump_cutoff  =dA_pump[0:int( (length +  delay_au[i]) / dt)] ###  截断   pump
    dAkick = dA_dely - dA_pump_cutoff  ##  pump+probe - pump
    damping = cx.RatioDampingMethod(t=t_delay, t0=delay_au[i], eta=eta, method=method)
    dAkick_damp = dAkick * damping   ## 衰减后的 pump+probe - pump
    ####  画吸收谱
    eps_delay = cx.gauge_kick_delay(t=t_delay, A0=0.22, dAdt=dAkick_damp, w_ev=w_ev, t0=delay_au[i])
    alpha_delay = cx.absorption_ev(w_ev, eps_delay) * 1000  # [1/nm]
    diff = (alpha_delay - alpha_gs)
    ax03.plot(t_delay, dAkick_damp)
    ax02.plot(w_ev, alpha_delay, label = delay_au[i])


ax02.plot(w_ev, alpha_gs, "k--", label = "gs")
ax02.legend()
ax02.set_ylabel("$\\alpha [1/ \mu m]$]", fontsize = fs)
ax02.set_xlabel("$\hbar\omega [eV]$", fontsize = fs)
ax02.set_xlim([start,end])

ax03.set_xlim([0,1300])
ax03.set_xticks([0,300,600,900,1200],[0,300,600,900,1200],fontsize = fs-2)
ax03.set_xlabel("t [a.u.]", fontsize = fs)
ax03.set_ylabel("dA/dt [a.u.]", fontsize = fs)


plt.show()