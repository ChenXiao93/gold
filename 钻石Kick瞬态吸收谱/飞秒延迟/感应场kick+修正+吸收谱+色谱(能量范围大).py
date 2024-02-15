import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import matplotlib.colors as colors
########  图片设定
fs = 14
fig0 = plt.figure(figsize= (6,9))
grid = plt.GridSpec(9, 1, wspace=0, hspace=1)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
ax01 = fig0.add_subplot(grid[0:2, 0:1]) ### Eind:  pump
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
ax02 = fig0.add_subplot(grid[2:4, 0:1]) ### Eind:  (pump+probe) - pump
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
ax03 = fig0.add_subplot(grid[4:6, 0:1]) ### dA/dt:  (pump+probe) - pump
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
ax04 = fig0.add_subplot(grid[6:9, 0:1])  ### 画吸收谱！
######################################
# delay = [0,40,80,518,559,600,641,682,723]
delay = [0,40,80,120,160, 200,240,280, 320, 360,400,440,480, 518,559,600,641,682]
### delay44 =
### delay44 =
delaylist = delay
dt = 0.05
### 能量范围设置
start,end,dw = 1,10, 0.1   # eV
w_ev = np.arange(start, end + dw, dw)  # eV
omega_au = w_ev / cx.hartree2ev  # a.u.
###  衰减修正参数设置
method = "gaussian" # "exponential"
length = 500 # a.u.
eta = 0.3 # eV
################# 纯pump #####################
time, dAref = cx.ReadData("./pump/gauge_field",1,5)
Aref = cx.ReadOne("./pump/gauge_field",2)
laserA = cx.ReadOne("./pump/laser",2)
laserE = cx.A2E(time,laserA)
ax01.plot(time ,laserE, "b", label = "pump")
ax01.plot(np.ones_like(np.arange(-3,4,1)) * 600, np.arange(-3,4,1), "--r") #红色竖虚线
ax01.legend(fontsize = fs-2)
ax01.set_xlim([0,1200])
ax01.set_ylim([-3,3])
ax01.set_ylabel("$E_{ind}$ [$V/\mathring{A}$]", fontsize = fs)
ax01.set_xlabel("t [fs]", fontsize = fs)
ax01.set_yticks([-2,0,2],[-2,0,2],fontsize = fs-2)
ax01.set_xticks([0,200,400,600,800,1000,1200],["","","","","","",""],fontsize = fs)
######  纯probe   kick = 0.22  ########
t_probe, dAdt_probe = cx.ReadData("./probe/gauge_field", 1,5)
t_probe = t_probe[0:int(length/0.05)] ###  截断
dAdt_probe = dAdt_probe[0:int(length/0.05)] ###  截断
damping_probe = cx.RatioDampingMethod(t=t_probe, t0=0, eta=eta, method=method) ###  衰减
dAdt_probe_damp = dAdt_probe * damping_probe ###  衰减
eps_pureprobe = cx.gauge_kick_delay(t=t_probe, A0=0.22, dAdt=dAdt_probe_damp, w_ev=w_ev, t0=0) # 复数
alpha_pureprobe = cx.absorption_ev(w_ev, eps_pureprobe)  * 1000 # [1/mu m]

ax02.plot(t_probe, dAdt_probe_damp, "r", label = "gs probe")
ax03.plot(t_probe, max(dAdt_probe) * damping_probe, "r-")
ax03.plot(t_probe, dAdt_probe_damp, "r-")
ax04.plot(w_ev, alpha_pureprobe, "r-", label = "gs")
######  纯probe kick = 0.44 ########
t_probe1, dAdt_probe1 = cx.ReadData("./probe0.44/gauge_field", 1,5)
t_probe1 = t_probe1[0:int(length/0.05)] ###  截断
dAdt_probe1 = dAdt_probe1[0:int(length/0.05)] ###  截断
damping_probe1 = cx.RatioDampingMethod(t=t_probe1, t0=0, eta=eta, method=method) ###  衰减
dAdt_probe_damp1 = dAdt_probe1 * damping_probe1 ###  衰减
eps_pureprobe1 = cx.gauge_kick_delay(t=t_probe1, A0=0.44, dAdt=dAdt_probe_damp1, w_ev=w_ev, t0=0) # 复数
alpha_pureprobe1 = cx.absorption_ev(w_ev, eps_pureprobe1)  * 1000 # [1/mu m]

ax02.plot(t_probe1, dAdt_probe_damp1,  label = "gs probe1")
ax03.plot(t_probe1, max(dAdt_probe1) * damping_probe1)
ax03.plot(t_probe1, dAdt_probe_damp1)
ax04.plot(w_ev, alpha_pureprobe1,  label = "gs1")
############ 提取probe = (pump+probe) - pump ##################
Matrix1D = []  # 准备热力图数据
for i in range(len(delay)):
    Atotal = cx.ReadOne("./%s/gauge_field" % delay[i], 2)
    Akick = Atotal - Aref

    dAtotal = cx.ReadOne("./%s/gauge_field" % delay[i], 5)
    dAkick = dAtotal - dAref

    ax02.plot(time, dAkick)
    # 截断
    cutoff = int((delay[i] + length) / dt)
    t = time[0:cutoff]
    dAcutoff = dAkick[0:cutoff]
    Acutoff = Akick[0:cutoff]
    # 衰减
    damping = cx.RatioDampingMethod(t=t, t0=delay[i], eta=eta, method=method)
    dAcutoff_damp = dAcutoff * damping
    Acutoff_damp = Akick[0:cutoff]* damping
    kick= np.abs(max(Acutoff_damp)) ##  注意kick大小不一样

    ax03.plot(t, dAcutoff_damp)

    eps_probe = cx.gauge_kick_delay(t=t, A0=kick, dAdt=dAcutoff_damp, w_ev=w_ev, t0=delay[i])
    alpha_probe = cx.absorption_ev(w_ev, eps_probe) * 1000  # [1/nm]

    ax04.plot(w_ev, alpha_probe , "--", label = "%.f" % delay[i])

    diff = alpha_probe - alpha_pureprobe
    Matrix1D = np.concatenate((Matrix1D, diff), axis=0)  # 将各行数据拼接起来

ax02.plot(np.ones_like(np.arange(-0.15,0.16,0.01)) * 600, np.arange(-0.15,0.16,0.01), "--r")  #红色竖虚线
########################## 图形设置
ax02.legend(loc = "upper right", fontsize = fs-2)
ax02.set_ylabel("$dA_{ind}/dt$ [a.u.]", fontsize = fs)
ax02.set_xlim([0, 1200])
ax02.set_ylim([-0.3,0.3])
ax02.set_yticks([-0.1,0,0.1],[-0.1,0,0.1],fontsize = fs-2)
ax02.set_xticks([0,200,400,600,800,1000,1200],["","","","","","",""],fontsize = fs)

ax03.set_xlim([0, 1200])
ax03.set_ylabel("$dA_{ind}/dt$ [a.u.]", fontsize=fs)
ax03.set_xlabel("$t$ [a.u.]", fontsize=fs-2)
ax03.set_xticks([0,200,400,600,800,1000,1200],["","","","","","",""],fontsize = fs)

ax04.legend(loc = "lower left", fontsize = fs-2)
ax04.set_xlim([start, end])
ax04.set_xlabel("$\hbar\omega \ [eV]$",fontsize=fs)
ax04.set_ylabel("$\\alpha [1/\mu m]$",fontsize=fs)

####
def colormap(x,y,z):
    X, Y = np.meshgrid(x, y)
    wstart, wend = y[0], y[-1]
    mid =  "{:.1f}".format((wstart+wend)/2)
    Matrix2D_amp = z.reshape(len(list(x)), len(y)).T # 将拼接起来的数据重新排布成和画布相适应的矩阵
    fig = plt.figure(figsize = (4,9))
    ax = fig.add_subplot(111)
    heatmap = ax.imshow(Matrix2D_amp,
                        interpolation='bilinear',  # 'bilinear' , "none"
                        cmap='RdBu',
                        norm=colors.CenteredNorm(),
                        aspect='auto',
                        extent=[x[0], x[-1], wend, wstart])
    ax.invert_yaxis()
    ax.set_xlabel('$\\tau $[fs]', fontsize=fs)
    ax.set_ylabel('$\\omega $[eV]', fontsize=fs)
    ax.set_xticks(list(x), labels= list(x),fontsize=fs)  # 设置x轴标签
    ax.set_yticks([wstart, (wstart+wend)/2,  wend], labels= [ "{:.1f}".format(wstart), mid,  "{:.1f}".format(wend)],fontsize=fs)  # 设置y轴标签
    cbar = plt.colorbar(heatmap, ax=ax)
    cbar.ax.tick_params(labelsize=fs)
    return()

colormap(x =np.array(delaylist), y = w_ev, z= Matrix1D)

plt.show()