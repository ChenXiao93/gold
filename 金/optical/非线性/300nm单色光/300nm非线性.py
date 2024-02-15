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
    end = min(len(A_ind),len(A_ext))  # 把ext场和ind场调整到长度一致
    time_ext,  A_ext= cx.shortten(time_ext, A_ext, start, end)
    time_ind, A_ind = cx.shortten(time_ind, A_ind, start, end)
    ######## Add the envelope to cut it 加包络 切断时间 ######################
    head = 500
    A_ext = env(time_ext, A_ext, pulse_duration+head)
    A_ind = env(time_ind, A_ind, pulse_duration+head)
    ######## 求介电函数 ####################
    wavelength = 300  # nm
    Omega = cx.Wavelength2Energy(wavelength) # eV
    epsilon_real, epsilon_imag = op.Get_Epsilon_Gauge(time_ext, A_ext, time_ind, A_ind, Omega, eta = 0.0)
    print("eps1=%.4f  eps2=%.4f" % (epsilon_real,epsilon_imag))
    # Johnson and Christy 1972   eps1 = -1.2360 eps2 = 5.7608
    ###### Refractive Index  求折射率和光强 ##################################
    # modification = np.sqrt(epsilon_real**2 + epsilon_imag**2)   # 修正1  ｜\epsilon｜
    modification = np.abs(epsilon_real)  # 修正2  ｜\epsilon_{real}｜
    n_real, n_imag = op.refractive_index_positive(epsilon_real, epsilon_imag)
    Emax = max(cx.A2E(time_ext, A_ext))  / modification  # epsilon modification
    # 真空中的光强！才可以用epsilon修正 来近似 介质内部的光强！  cx.E2I 里面是I= 1/2 xxx， 文献里是 I= 2 xxx， 所以有个4的倍数关系
    I = cx.E2I(E = Emax, Refractive_index = 1) * 4  # [W/cm^2]  # cited from the 'n0n2k0k2'
    # I = cx.E2I(E = Emax, Refractive_index = n_real) * 4  # 介质内部的光强！ 不要epsilon修正
    return(I, n_real, n_imag)


def chi3_cal(pulse_duration):
    I, n_real, n_imag = [], [], []
    for a in [10,20,30,40]:
        I_a, n_real_a, n_imag_a = N_Cal(a,pulse_duration)
        I.append(I_a)
        n_real.append(n_real_a)
        n_imag.append(n_imag_a)
    #######   线性拟合 #######
    fitted_yr, pr = ch.nfit(I, n_real, 1)
    fitted_yi, pi = ch.nfit(I, n_imag, 1)
    n0 = pr[0]
    n2 = pr[1]*1E-4  # [m^2/W]
    k0 =  pi[0]
    k2 = pi[1]*1E-4 # [m^2/W]

    chi3r = ( 4 * n0 * cx.Vaccum_permitivity  * cx.Light_speed / 3) * (n0*n2 - k0*k2)  # [m^2/V^2]
    chi3i = ( 4 * n0 * cx.Vaccum_permitivity  * cx.Light_speed / 3) * (n0*k2 + k0*n2)  # [m^2/V^2]
    return(chi3r,chi3i)

print("eps1=-1.2360  eps2=5.7608  Johnson and Christy 1972")
Time= np.arange(0,8000,1000)
CHI3R, CHI3I = [],[]
for pulse_duration in Time:
    chi3r,chi3i = chi3_cal(pulse_duration)
    CHI3R.append(chi3r)
    CHI3I.append(chi3i)


# print(CHI3R,CHI3I)
plt.subplot(1,2,1)
plt.loglog(Time * 0.02419 * 1E-15, np.abs(CHI3R), "--o")

plt.subplot(1,2,2)
plt.loglog(Time * 0.02419 * 1E-15, np.abs(CHI3I), "--o")

plt.show()