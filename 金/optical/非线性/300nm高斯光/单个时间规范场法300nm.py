import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
####################
wavelength = 300 # nm
omega = cx.Wavelength2Energy(wavelength) # eV
amp = [1,4,7,10] # Ha
T= [1500] # a.u.
print("eps1 = -1.237  eps2 = 6.031 (exp)")
print("T = %s a.u." % T[0])
for t in range(len(T)):
    for i in range(len(amp)):
        ######  A_ext 外场 and A_ind 感应场 ###############
        time_ext, A_ext = cx.ReadData("./laser/t%s/A%s" % (T[t], amp[i]), 0, 1)
        A_ind = cx.ReadOne("./data/t%s/A%s/gauge_field" % (T[t], amp[i]), 2)
        # 把ext场和ind场调整到长度一致
        dt = time_ext[2] - time_ext[1]
        start = 0
        end = int(T[t] / dt)
        time_ext = time_ext[start:end]
        A_ext = A_ext[start:end]
        A_ind = A_ind[start:end]
        ## 衰减函数配置
        damp = np.ones_like(time_ext)
        ###### 电场法
        epsilon_a =  cx.gauge_mono(time_ext, A_ext, A_ind, omega,damp)
        print("eps1_a=%.4f,eps2_a=%.4f, A =%s Ha (电场法)" %(epsilon_a.real, epsilon_a.imag, amp[i]))

plt.show()
