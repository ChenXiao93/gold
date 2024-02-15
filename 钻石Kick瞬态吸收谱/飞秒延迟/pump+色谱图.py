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
ax01 = fig0.add_subplot(grid[0:2, 0:80]) ### pump
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
ax02 = fig0.add_subplot(grid[2:6, 0:100]) ### 色谱图
#######################################
delay = [0,40,80,108,120, 518,559,600,641,682,723]
largekick = [160, 200,240,280, 320, 360,400,440,480]
delaylist = delay
dt = 0.05
kick1 = 0.22
kick2 = 0.44 # 0.43925
### 能量范围设置
start,end,dw = 1, 10, 0.1   # eV
w_ev = np.arange(start, end + dw, dw)  # eV
omega_au = w_ev / cx.hartree2ev  # a.u.
###  衰减修正参数设置
method = "gaussian"      # "gaussian" # "exponential"
length = 500 # a.u.
eta = 0.3 # eV
################# 纯pump #####################
time, dAref = cx.ReadData("./pump/gauge_field",1,5)
laserA = cx.ReadOne("./pump/laser",2)
laserE = cx.A2E(time,laserA)
ax01.plot(time ,laserE, "b", label = "pump")
ax01.legend(fontsize = fs-2)
ax01.set_xlim([0,1200])
ax01.set_ylim([-3,3])
ax01.set_ylabel("$E_{ind}$ [$V/\mathring{A}$]", fontsize = fs)
ax01.set_xlabel("t [fs]", fontsize = fs)
ax01.set_yticks([-2,0,2],[-2,0,2],fontsize = fs-2)
#################  纯probe ######################
t_probe, dAdt_probe = cx.ReadData("./probe/gauge_field", 1,5)
E_probe = dAdt_probe * (-1/cx.c_au)  * 51.44
t_probe = t_probe[0:int(length/0.05)] ###  截断
dAdt_probe = dAdt_probe[0:int(length/0.05)] ###  截断
damping_probe = cx.RatioDampingMethod(t=t_probe, t0=0, eta=eta, method=method) ###  衰减包络
dAdt_probe_damp = dAdt_probe * damping_probe ###  衰减
eps_pureprobe = cx.gauge_kick_delay(t=t_probe, A0=kick1, dAdt=dAdt_probe_damp, w_ev=w_ev, t0=0) # 复数
alpha_pureprobe = cx.absorption_ev(w_ev, eps_pureprobe)  * 1000 # [1/mu m]  吸收系数
############  提取probe = (pump+probe) - pump ##################
Matrix1D = []  # 准备热力图数据
for i in range(len(delay)):
    dAtotal = cx.ReadOne("./%s/gauge_field" % delay[i], 5)
    dAkick = dAtotal - dAref
    # 截断
    cutoff = int((delay[i] + length) / dt)
    t = time[0:cutoff]
    dAcutoff = dAkick[0:cutoff]
    # 衰减
    damping = cx.RatioDampingMethod(t=t, t0=delay[i], eta=eta, method=method)
    dAcutoff_damp = dAcutoff * damping
    ####  注意kick大小不一样
    if delay[i] in largekick:
        kick = kick2
    else:
        kick = kick1
    ####  计算介电函数
    eps_probe = cx.gauge_kick_delay(t=t, A0=kick, dAdt=dAcutoff_damp, w_ev=w_ev, t0=delay[i])
    alpha_probe = cx.absorption_ev(w_ev, eps_probe) * 1000  # [1/nm]
    diff = (alpha_probe - alpha_pureprobe)
    Matrix1D = np.concatenate((Matrix1D, diff), axis=0)  # 将各行数据拼接起来

######## 色谱图函数
def colormap(x,y,z,ax):
    X, Y = np.meshgrid(x, y)
    wstart, wend = y[0], y[-1]
    mid =  "{:.1f}".format((wstart+wend)/2)
    Matrix2D_amp = z.reshape(len(list(x)), len(y)).T # 将拼接起来的数据重新排布成和画布相适应的矩阵
    heatmap = ax.imshow(Matrix2D_amp,
                        interpolation='bilinear',  # 'bilinear' , "none"
                        cmap='RdBu',
                        norm=colors.CenteredNorm(),
                        aspect='auto',
                        extent=[x[0], x[-1],wend, wstart])
    ax.invert_yaxis()
    ax.set_xlabel('$\\tau $[fs]', fontsize=fs)
    ax.set_ylabel('$\\omega $[eV]', fontsize=fs)
    # ax.set_xticks(list(x), labels= list(x))  # 设置x轴标签
    # ax.set_yticks([wstart, (wstart+wend)/2,  wend], labels= [ "{:.1f}".format(wstart), mid,  "{:.1f}".format(wend)],fontsize=fs)  # 设置y轴标签
    cbar = plt.colorbar(heatmap, ax=ax)
    cbar.ax.tick_params(labelsize=fs-2)
    cbar.set_label("$\Delta \\alpha(\omega) [1/\mu m]$", fontsize = fs)
    return()

##  画色谱图
colormap(x =np.array(delaylist), y = w_ev, z= Matrix1D, ax=ax02)
# 画直接/间接带隙
exp_idgap = 4.96 * np.ones_like(delaylist)
exp_dgap = 6.42 * np.ones_like(delaylist)
ax02.plot(delaylist, exp_idgap,"k--")
ax02.plot(delaylist, exp_dgap,"k-")

plt.show()