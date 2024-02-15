import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import matplotlib.colors as colors
########  图片设定
fs = 14

fig0 = plt.figure(figsize= (9, 6))
grid = plt.GridSpec(6, 100, wspace=0, hspace=1.5)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
ax01 = fig0.add_subplot(grid[0:2, 0:80]) ### Eind:  pump
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
ax02 = fig0.add_subplot(grid[2:6, 0:100]) ### Eind:  (pump+probe) - pump
#######################################
dt = 0.05
start,end,dw = 1, 9, 0.1   # eV
w_ev = np.arange(start, end + dw, dw)  # eV
omega_au = w_ev / cx.hartree2ev  # a.u.
length = 500 # a.u.  最大500
eta = 0.3 # eV
method = "gaussian"      # "gaussian" # "exponential"
######## 纯 pump
text_pump, Aext_pump = cx.ReadData("./pump/IRpump",0,1)
Eext_pump = cx.A2E(text_pump, Aext_pump)
ax01.plot(text_pump * cx.timeau2fs, Eext_pump, "b", label = "pump")
ax01.legend(fontsize = fs-2)
ax01.set_xlim([0,30])
ax01.set_ylim([-3,3])
ax01.set_ylabel("$E_{ind}$ [$V/\mathring{A}$]", fontsize = fs)
ax01.set_xlabel("t [fs]", fontsize = fs)
ax01.set_yticks([-2,0,2],[-2,0,2],fontsize = fs-2)

t_pump,dA_pump = cx.ReadData("./pump/gauge_field",1,5)
#######  gs 吸收谱 ########
t_probe, dA_probe = cx.ReadData("./gs/gauge_field", 1, 5)
t_probe = t_probe[0:int(length/0.05)] ###  截断
dA_probe = dA_probe[0:int(length/0.05)] ###  截断
damping_probe = cx.RatioDampingMethod(t=t_probe, t0=0, eta=eta, method=method) ###  衰减
dA_probe_damp = dA_probe * damping_probe ###  衰减
eps_probe = cx.gauge_kick_delay(t=t_probe, A0=0.22, dAdt=dA_probe_damp, w_ev=w_ev, t0=0) # 复数
alpha_probe = cx.absorption_ev(w_ev, eps_probe)  * 1000 # [1/mu m]
#######  不同延迟的 probe   ##################################
###  还有 100， 700，800，900，1000 这几个延迟失败了，正在重新计算
## 600 也不可以！！！！！

delay_au = [1000]
for i in range(len(delay_au)):
    t_delay, dA_dely = cx.ReadData("./%s/gauge_field" % delay_au[i],1, 5)


    t_delay = t_delay[0:int( (length +  delay_au[i]) / 0.05)]  ###  截断
    dA_dely = dA_dely[0:int( (length +  delay_au[i]) / 0.05)]  ###  截断
    dA_pump_cutoff  =dA_pump[0:int( (length +  delay_au[i]) / 0.05)]

    dAkick = dA_dely - dA_pump_cutoff

    ax02.plot(t_delay, dAkick, label = delay_au[i])



plt.show()