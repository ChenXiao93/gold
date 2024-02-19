import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.Laser as la


#########################################
# amp = [1,4,7,10] # Ha
# T= [500, 1500, 3000, 5000, 10000] # a.u.
amp = [10] # Ha
T= [1500]
Ts= np.array(T) * cx.timeau2s * 1e12 * 2/5 # ps  # FWHM, 因为高斯光产生的时候设置PulseDuration=int(t*2/5)
RATIO = [0.1,0.1,0.1,0.1,0.1]
STD = [1E-8,1E-8,1E-8,1E-8,1E-8]# Ha
CHI3R,CHI3I = [],[]
dt = 0.2 # a.u.
wavelegth = 300 # nm

for t in range(len(T)):
    print("eps1 = -11.980  eps2 = 1.1524 (ex)")  # 实验值
    eps1_current, eps2_current, Intensity, N,K = [], [], [],[],[] # 初始化
    for i in range(len(amp)):
        # 通过电流计算介电函数
        time_ind, A_ind = cx.ReadData("./data/t%s/A%s/gauge_field" % (T[t], amp[i]), 1, 2)

        time_ext, A_ext, envy = la.GaussianPulse(PropTime=T[t] + dt,
                                                PulseDuration=int(T[t] * 2 / 5),
                                                Wavelength=wavelegth ,
                                                Amplitude=amp[i], TimeStep=dt, Phase=0)

        J_cal = cx.current_cal(t_ind=time_ind, A_ind=A_ind, V=114.5822) # 感应电流
        J_oct = cx.ReadOne("./data/t%s/A%s/total_current" % (T[t], amp[i]), 2)

        plt.plot(time_ind, J_cal, label = "cal")
        plt.plot(time_ind, J_oct, label = "oct")


plt.legend()
plt.show()
