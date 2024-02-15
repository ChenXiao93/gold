import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op
import OctopusTools.Laser as la
import OctopusTools.chi3 as ch

def env(x,y, pulse_duration):
    #  sin2 attenuation
    tail = 500
    envy =   max(y) * (np.heaviside(x - pulse_duration, 1)   * \
             np.square( np.sin( (x - pulse_duration)  *  np.pi / tail /2    + np.pi/2)) + \
             np.heaviside(pulse_duration - x - 0.01, 1) ) *  np.heaviside( pulse_duration + tail - x, 1)
    newy = envy * y / max(envy)
    return(newy)



def N_Cal(Amplitude, pulse_duration):
    ######   A_ext 外场 and A_ind 感应场 ###############
    time_ext, A_ext = cx.ReadData("./laser/A%s" % Amplitude, 0, 1)
    time_ind, A_ind = cx.ReadData("./gauge_field/A%s" % Amplitude, 1, 2)
    start = 0
    end = int((500+500+pulse_duration)/0.2)  # 把ext场和ind场调整到长度一致
    time_ext,  A_ext= cx.shortten(time_ext, A_ext, start, end)
    time_ind, A_ind = cx.shortten(time_ind, A_ind, start, end)
    # plt.figure()
    # plt.plot(time_ext * 0.02419,  A_ext, "blue", label = "ext")
    # plt.plot(time_ind*0.02419 , A_ind, "orange", label = "ind")

    ######## Add the envelope to cut it 加包络 切断时间 ######################
    head = 500
    time, envy = la.env_exp(mid=pulse_duration+head, tail=500, dt=0.2, gamma=0.0)
    A_ext = envy * A_ext
    A_ind = envy * A_ind
    # A_ext = env(time_ext, A_ext, pulse_duration+head)
    # A_ind = env(time_ind, A_ind, pulse_duration+head)
    # plt.figure()
    # plt.plot(time_ext * 0.02419,  A_ext, "blue", label = "ext")
    # plt.plot(time_ind*0.02419 , A_ind, "orange", label = "ind")

    ######## 求介电函数 ####################
    wavelength = 300  # nm
    Omega = cx.Wavelength2Energy(wavelength) # eV
    epsilon_real, epsilon_imag = op.Get_Epsilon_Gauge(time_ext, A_ext, time_ind, A_ind, Omega, eta = 0.00)
    print("eps1=%.4f  eps2=%.4f (A=%s)" % (epsilon_real,epsilon_imag,Amplitude ))
    return()


def chi3_cal(pulse_duration):
    for a in [10]:
        N_Cal(a,pulse_duration)
    return()

print("eps1=-1.2360  eps2=5.7608 (exp)")
# Time= np.arange(0,9000,8000)
Time= [3000]
for pulse_duration in Time:
    chi3_cal(pulse_duration)

plt.ylabel("A [Ha]")
plt.xlabel("t [fs]")
plt.legend()
plt.show()