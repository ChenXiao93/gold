import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx

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
# delay_au = list(range(0,1300,100))
delay_au = [0,200,300,400,500, 1100, 1200]
delay_fs = np.around( np.array(delay_au) * cx.timeau2fs, 1)

dt = 0.05
start,end,dw = 1, 10, 0.1   # eV
w_ev = np.arange(start, end + dw, dw)  # eV
omega_au = w_ev / cx.hartree2ev  # a.u.
length = 500 # a.u.  最大500
eta = 0.3 # eV
method = "gaussian"     # "gaussian" # "exponential"
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
eps_gs = cx.gauge_kick_delay(t=t_probe, A0=0.22, dAdt=dA_probe_damp, w_ev=w_ev, t0=0) # 复数
alpha_gs = cx.absorption_ev(w_ev, eps_gs)  * 1000 # [1/mu m]
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
    Matrix1D = np.concatenate((Matrix1D, diff), axis=0)  # 将各行数据拼接起来

########################## 图形设置
def colormap(x,y,z,ax):
    X, Y = np.meshgrid(x, y)
    wstart, wend = y[0], y[-1]
    mid =  "{:.1f}".format((wstart+wend)/2)
    Matrix2D_amp = z.reshape(len(list(x)), len(y)).T # 将拼接起来的数据重新排布成和画布相适应的矩阵
    heatmap = ax.imshow(Matrix2D_amp,
                        interpolation='none',  # 'bilinear'
                        cmap='RdBu',
                        # norm=colors.CenteredNorm(),
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

colormap(x =np.array(delay_fs ), y = w_ev, z= Matrix1D, ax=ax02)
exp_gap = 4.94 * np.ones_like(delay_fs )
ax02.plot(delay_fs , exp_gap,"k-")


plt.show()