import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import matplotlib.colors as colors
### 色谱图函数设置
def colormap(x,y,z,ax):
    X, Y = np.meshgrid(x, y)
    wstart, wend = y[0], y[-1]
    mid =  "{:.1f}".format((wstart+wend)/2)
    Matrix2D_amp = z.reshape(len(list(x)), len(y)).T # 将拼接起来的数据重新排布成和画布相适应的矩阵
    heatmap = ax.imshow(Matrix2D_amp,
                        interpolation='bilinear',                   # 'none',  'bilinear'
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
########  图片设定
fs = 14

fig0 = plt.figure(figsize= (6,4))
grid = plt.GridSpec(100, 100, wspace=0, hspace=0)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
ax01 = fig0.add_subplot(grid[0:100, 0:100])

fig2 = plt.figure(figsize= (8,6))
grid = plt.GridSpec(100, 100, wspace=0, hspace=10)
ax1 = fig2.add_subplot(grid[0:30, 0:80])
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
ax2 = fig2.add_subplot(grid[40:100, 0:100])
#######################################
delay_au = [0, 200,300,400,500, 1100, 1200]

dt = 0.05
start,end,dw = 2,10, 0.1   # eV
# start,end,dw = 4, 8, 0.01   # eV
w_ev = np.arange(start, end + dw, dw)  # eV
omega_au = w_ev / cx.hartree2ev  # a.u.
length = 500 # a.u.  最大500
eta = 0.3 # eV
method = "gaussian"    # "gaussian" # "exponential"
######## 纯 pump
text_pump, Aext_pump = cx.ReadData("./pump/IRpump",0,1)
t_pump,dA_pump = cx.ReadData("./pump/gauge_field",1,5)
Eext_pump = cx.A2E(text_pump, Aext_pump)
ax1.plot(text_pump, Eext_pump, "b", label = "pump")
ax1.legend(fontsize = fs-2)
ax1.set_xlim([0,1200])
ax1.set_ylim([-3,3])
ax1.set_ylabel("$E_{ind}$ [$V/\mathring{A}$]", fontsize = fs-2)
ax1.set_xlabel("t [a.u.]", fontsize = fs-2)
ax1.set_yticks([-2,0,2],[-2,0,2],fontsize = fs-2)
#######  gs 吸收谱 ########
t_probe, dA_probe = cx.ReadData("./gs/gauge_field", 1, 5)
t_probe = t_probe[0:int(length/0.05)] ###  截断
dA_probe = dA_probe[0:int(length/0.05)] ###  截断
damping_probe = cx.RatioDampingMethod(t=t_probe, t0=0, eta=eta, method=method) ###  衰减
dA_probe_damp = dA_probe * damping_probe ###  衰减
eps_gs = cx.gauge_kick_delay(t=t_probe, A0=0.22, dAdt=dA_probe_damp, w_ev=w_ev, t0=0) # 复数
alpha_gs = cx.absorption_ev(w_ev, eps_gs)  * 1000 # [1/mu m]
#######
delayref = 500
t_delay, dA_dely = cx.ReadData("./%s/gauge_field" % delayref, 1, 5)
t_delay = t_delay[0:int((length + delayref) / dt)]  ###  截断
dA_dely = dA_dely[0:int((length + delayref) / dt)]  ###  截断  pump+probe
dA_pump_cutoff = dA_pump[0:int((length + delayref) / dt)]  ###  截断   pump
dAkick = dA_dely - dA_pump_cutoff  ##  pump+probe - pump
damping = cx.RatioDampingMethod(t=t_delay, t0=delayref, eta=eta, method=method)
dAkick_damp = dAkick * damping  ## 衰减后的 pump+probe - pump
eps_delay = cx.gauge_kick_delay(t=t_delay, A0=0.22, dAdt=dAkick_damp, w_ev=w_ev, t0=delayref)
alpha_delayref = cx.absorption_ev(w_ev, eps_delay) * 1000  # [1/nm]
#######  不同延迟的 probe   ##################################
Matrix1D = []  # 准备热力图数据
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
    # diff = alpha_delay - alpha_delayref
    ax01.plot(w_ev, diff, label = delay_au[i])

    Matrix1D = np.concatenate((Matrix1D, diff), axis=0)  # 将各行数据拼接起来

ax01.plot(w_ev, alpha_gs - alpha_gs , "k--", label = "gs")
ax01.legend()
ax01.set_ylabel("$\\alpha_{kick} - \\alpha_{%s} [1/ \mu m]$]" % delayref, fontsize = fs)
ax01.set_xlabel("$\hbar\omega [eV]$", fontsize = fs)
ax01.set_xlim([start,end])



##########################画色谱图
colormap(x =np.array(delay_au), y = w_ev, z= Matrix1D, ax=ax2)
exp_idgap = 4.94 * np.ones_like(delay_au)
exp_dgap = 6.42 * np.ones_like(delay_au)
ax2.plot(delay_au , exp_idgap ,"k--")
ax2.plot(delay_au , exp_dgap ,"k-")

plt.show()