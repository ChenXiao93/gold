import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as b

def pump_pulse(wavelength, intensity, pulse_duration):
    #### 输入参数 #####
    # wavelength      [nm]  波长
    # intensity       [W/cm^2] 真空/空气中激光强度
    # pulse duration  [fs] 脉冲持续时间
    #### 输出参数 ####
    #  out_t    [fs]  实数
    #  out_E    [V/Angstrom]  复数
    #  out_A    [Ha]    复数
    #################
    dt=0.0015* b.timeau2fs    # time step magnitude  [au]  to [fs]
    t = np.arange(-60,10,dt)   #  fs
    omega_s = b.Wavelength2Omega(wavelength)   # [1/s]
    omega_fs = omega_s * 1e-15 # [1/s] -- [1/fs]

    tau=5
    phi=0   # carrier-envelope phase (CEP) is not stabilised in experiment!
    E0 = b.I2E(I= intensity, Refractive_index=1)  # 外界空气折射率为1   E[V/Angstrom]

    epsilont=E0*np.exp(-2*np.log(2)*np.square(t/pulse_duration))        # log(2) 以e为底
    out_t = np.arange(0, 70, dt)  # [fs]
    Et_complex = epsilont*np.exp((t*omega_fs+phi)*1j)  # 复数
    out_E_complex = np.flipud(Et_complex)   ## np.flipud()方法用于翻转列表，将矩阵进行上下翻转 out_E[V/Angstrom]
    out_A_complex = b.E2A(out_t / b.timeau2fs, out_E_complex)

    return(out_t, out_E_complex, out_A_complex)

def probe_pulse(wavelength, intensity, pulse_duration):
    #### 输入参数 #####
    # wavelength      [nm]  波长
    # intensity       [W/cm^2] 真空/空气中激光强度
    # pulse duration  [fs] 脉冲持续时间
    #### 输出参数 ####

    #################
    dt=0.1* b.timeau2fs # time step magnitude  from [hbar/Hartree] to [fs]
    t = np.arange(-60, 10, dt)  # fs
    omega_s = b.Wavelength2Omega(wavelength)   # [1/s]
    omega_fs = omega_s * 1e-15 # [1/s] -- [1/fs]
    tau = pulse_duration       #  [fs]
    I=intensity
    chrip=10     #[fs^2]
    impedance= b.Light_speed*b.Vacuum_permeability     # vacuum impedance Z [H/s]  真空阻抗
    sigma0= tau/(2*np.sqrt(np.log(2)))    #  1.8[fs]   高斯分布的标准差
    sigma=np.sqrt(np.square(sigma0)+np.square(chrip)/np.square(sigma0))    #  5.84[fs]
    G2=0.5*chrip/(np.power(sigma0,4)+np.square(chrip))             # 是啥？  0.045 [1/fs^2]
    E0=np.sqrt(2*impedance*I)*(1E-10)  # 0.00087 [V/A]
    E1=E0*np.sqrt(sigma0/sigma)        # 0.00048 [V/A]
    k=omega*(1e15)/b.Light_speed                  # 2.6E7[1/m]

    v_phase = omega/k    #  [m/fs]
    ceo=0        # 待定=0 carrier envelope offset 载波包络偏移
    TOD=9        # [fs^3] P88页，不知是石英玻璃的TOD还是钻石的TOD，待定 third order dispersion
    w0=1        # 1微米= 10-6 m  估计，待定
    epsilon=7.1  # experiment data 7.1, octopus 7.5

    E_U=(1/np.sqrt(np.square(epsilon)+1)) *  \
        E0*np.exp(-0.5*np.square(t/sigma0)) * \
        np.cos(omega*t+ceo+G2*np.square(t))

    E_V=(1/np.sqrt(np.square(epsilon)+1)) *  \
        E1*np.exp(-0.5*np.square(t/sigma)) * \
        np.cos(omega * t+ceo+G2*np.square(t))

    new_EU=np.flipud(E_U*1.77)  # 1.77 is np.sqrt(pi), is the integral for y direction.
    new_EV=np.flipud(E_V*1.77)
    new_t = np.arange(0, 70, dt)  # fs

    plt.plot(new_t, new_EU)
    plt.plot(new_t,new_EV)

    return()


def Gaussian_pulse(Intensity, wavelength, pulse_width, pulse_mid_time):
    # Intensity [W/cm^2], Wavelength [nm] pulse_width[ps]
    # E(t) = E_0 exp{-(t/t_{FWHM})^2} exp{i w_0 t}
    dt = 0.01   # 1ps
    t = np.arange(0,30,dt) # ps
    E_0 = my.I2E(Intensity, 1)  #I [W/cm^2] -->  E [V/Angstrom]
    omega = my.Wavelength2Omega(wavelength)  #  Wavelength [nm]  to  Omega[1/s]
    omega = omega * 1E12  # [1/s] --> [1/ps]
    E_real = E_0  * np.cos(omega * t ) * np.exp(-np.square((t-pulse_mid_time)/pulse_width))  # 1fs =   1E-15[s]
    E_imag = E_0 * np.cos(omega * t )  * np.exp(-np.square((t-pulse_mid_time)/pulse_width))# 1fs =   1E-15[s]
    return(t,E_real)

#####################
def test():
    #### 测试pump pulse ##
    # t, E, A = pump_pulse(wavelength = 720,  intensity=1.7E12, pulse_duration = 5)
    # plt.subplot(1,2,1)
    # plt.plot(t, E.real)
    # plt.plot(t, E.imag)
    # plt.subplot(1,2,2)
    # plt.plot(t, A.real)
    # plt.plot(t, A.imag)
    # print(len(t))
    # 测试 probe pulse U
    probe_pulse(wavelength = 240,  intensity = 1E7, pulse_duration = 3.4 )
    # 测试 probe pulse V
    probe_pulse(wavelength = 240,  intensity = 1E7, pulse_duration = 10 )



    plt.show()

    return()

if __name__ == "__main__" :
    test()



