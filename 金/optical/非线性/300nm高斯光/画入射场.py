import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.Laser as la


RATIO = [0.2, 0.5, 0.5, 0.5, 0.58]
AST = [1E-8,1E-8,1E-8,1E-8,1E-8]# Ha
amp = [1,4,7,10]
T= [500,1500,3000,5000,10000] # a.u.
DT = [0.2, 0.2, 0.2, 0.2, 0.05]
extratime = 10000 - np.array(T)

fig = plt.figure(figsize=(6,7))
fs  =12
grid = plt.GridSpec(5, 1, wspace=0, hspace=0.4)
plt.rcParams['xtick.direction'] = 'in'



for i in range(5):
    pulse_duration = T[i]
    dt = DT[i]
    A_standard = AST[i]
    ratio = RATIO[i]
    Amplitude = 10
    time_ext, A_ext = cx.ReadData("./laser/t%s/A%s" % (pulse_duration, Amplitude), 0, 1)
    time_ind, A_ind = cx.ReadData("./data/t%s/A%s/gauge_field" % (pulse_duration,Amplitude), 1, 2)
    start = 0
    end = int((pulse_duration+100)/dt)  # 把ext场和ind场调整到长度一致. 100是额外添加的尾部，用于平复感应场的滞后波动
    time_ext,  A_ext = cx.shortten(time_ext, A_ext, start, end)
    time_ind, A_ind = cx.shortten(time_ind, A_ind, start, end)

    #### 加上自动化调控衰减！ ####
    eta_au = np.log( np.abs(A_ind[-1] / A_standard))  / time_ind[-1]  # atom unit
    time, envy = la.env_exp(mid=(pulse_duration*(1-ratio)), tail=(pulse_duration*ratio) + 100, dt=dt, gamma=eta_au)
    A_ind = A_ind * envy
    A_ext = A_ext * envy
    A_tot = A_ext + A_ind


    A_ext = A_tot

    ax = fig.add_subplot(grid[i:i+1, 0:1])
    ax.plot(time_ext * 0.02419, A_ext)
    ax.plot(time_ext * 0.02419, envy * max(A_ext))
    end = (max(time_ext)*0.02419)
    ax.set_xlim([0,end])


ax.set_ylabel("A [Ha]", fontsize = fs)
ax.set_xlabel("Propagation time t [fs]", fontsize = fs)

plt.show()