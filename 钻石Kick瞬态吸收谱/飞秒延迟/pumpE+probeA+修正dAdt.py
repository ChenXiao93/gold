import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx

########  图片设定
fs = 14
fig0 = plt.figure(figsize= (9,6))
grid = plt.GridSpec(6, 1, wspace=0, hspace=0.1)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
ax01 = fig0.add_subplot(grid[0:2, 0:1]) ### Eind:  pump
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
ax02 = fig0.add_subplot(grid[2:4, 0:1]) ### Aind:  (pump+probe) - pump
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
ax03 = fig0.add_subplot(grid[4:6, 0:1]) ### dAind/dt:  (pump+probe) - pump
#######################################
method = "gaussian"  #  "exponential"
eta = 0.3 # eV
length = 500 # a.u.

dt = 0.05
delay = [0,40,80,108,120, 160, 200,240,280, 320, 360,400,440,480, 518,559,600,641,682]
largekick = [160, 200,240,280, 320, 360,400,440,480]
longkick = [760,800,840,880,920,960,1000]
#### [0,40,80,120, 518,559,600,641,682,723]   kick = 0.22
### 160, 200,240,280, 320, 360,400,440,480 kick = 0.43925
################# 纯pump #####################
time, Aref = cx.ReadData("./pump/gauge_field",1,2)
dAref = cx.ReadOne("./pump/gauge_field",5)
##########   外电场
laserA = cx.ReadOne("./pump/laser",2)
laserE = cx.A2E(time,laserA)
ax01.plot(time,laserE, "b", label = "pump")
##############  纯probe ################
t_probe, A_probe = cx.ReadData("./probe/gauge_field", 1,2)
dA_probe = cx.ReadOne("./probe/gauge_field", 5)
t_probe = t_probe[0:int(length/0.05)] ###  截断
A_probe = A_probe[0:int(length/0.05)] ###  截断
dA_probe = dA_probe[0:int(length/0.05)] ###  截断
damping_probe = cx.RatioDampingMethod(t=t_probe, t0=0, eta=eta, method=method) ###  衰减
A_probe_damp = A_probe * damping_probe ###  衰减
dA_probe_damp = dA_probe * damping_probe ###  衰减

ax02.plot(t_probe, A_probe, "--r")  # 截断
ax03.plot(t_probe, dA_probe_damp, "--r")  #  衰减
##########  提取probe = (pump+probe) - pump  ####################
for i in range(len(delay)):
    Atotal = cx.ReadOne("./%s/gauge_field" % delay[i], 2)
    dAtotal = cx.ReadOne("./%s/gauge_field" % delay[i], 5)
    Akick = Atotal - Aref
    dAkick = dAtotal -  dAref

    ax02.plot(time, Akick)    # 原始数据
    # 截断
    cutoff = int((delay[i] + length) / dt)
    t = time[0:cutoff]
    Acutoff = Akick[0:cutoff]
    dAcutoff = dAkick[0:cutoff]
    # 衰减
    damping = cx.RatioDampingMethod(t=t, t0=delay[i], eta=eta, method=method)
    Acutoff_damp = Acutoff * damping
    dAcutoff_damp = dAcutoff * damping
    # 画图: 截断+衰减
    ax03.plot(t, dAcutoff_damp)  # 截断+衰减
    ax03.plot(t, max(dAcutoff) * damping)  # 截断+衰减包络

########################## 图形设置
ax01.legend(fontsize = fs-2)
ax01.set_xlim([0,1200])
ax01.set_ylim([-3,3])
ax01.set_ylabel("$E_{ind}$ [$V/\mathring{A}$]", fontsize = fs)
ax01.set_xticks([0,200,400,600,800,1000,1200],["","","","","","",""],fontsize = fs)
ax01.set_yticks([-2,0,2],[-2,0,2],fontsize = fs-2)

ax02.set_ylabel("$A_{ind}$ [$Ha$]", fontsize = fs)
ax02.set_xlim([0,1200])
ax02.set_xticks([0,200,400,600,800,1000,1200],["","","","","","",""],fontsize = fs)

ax03.set_xticks([0,200,400,600,800,1000,1200],[0,200,400,600,800,1000,1200],fontsize = fs)
ax03.set_ylabel("$dA_{ind}/dt$ [$a.u.$]", fontsize = fs)
ax03.set_xlabel("Time [$a.u.$]", fontsize = fs)
ax03.set_xlim([0, 1200])

plt.show()