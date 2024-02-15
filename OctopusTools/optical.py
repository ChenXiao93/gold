import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx

def damping(prop_time,  kick_delay, damping_factor):
    # time,t0 [hbar/Ha]; eta [eV]
    eta_au = damping_factor / 27.2114  #  [Ha/hbar]
    f_damp = []
    for i in range(len(prop_time)):
        if prop_time[i] < kick_delay:
            f_damp.append(1)
        else:
            value = np.exp(- eta_au * (prop_time[i] -  kick_delay))
            f_damp.append(value)
    f_damp = np.array(f_damp)
    return(f_damp)

def gauge_field_kick(prop_time, kick_delay, kick_amplitude, dAprobe_dt, damping_factor, frequency_range):
    ############################################
    # 使用kick法，利用规范场计算介电函数
    # 输入参数 prop_time # [hbar/Ha]   总传播时间Get_Epsilon_GaugeGet_Epsilon_GaugeGet_Epsilon_GaugeGet_Epsilon_GaugeGet_Epsilon_GaugeGet_Epsilon_Gauge
    # kick_delay # [hbar/Ha]  kick延迟时间
    # kick_amplitude # [Ha]   kick大小
    # induced_Afield # [Ha]   感应规范场
    # induced_efield = -1 * cx.ReadOne("./gauge_field", 5) / cx.c_au  # E= -1/c  dA/dt  原子单位c=137
    # induced_efield   # Eprobe 感应电场，用原子单位  E =（-1/c）*（dA/dt）  原子单位c=137
    # Eind_si = Eprobe * 51.44  # 1 [Ha] = 27.2114 [eV], 1 [Bohr] = 0.529 [Angstrom]   27.2114/0.529 =51.44
    # damping_factor   #  衰减因子  eta [eV]
    # frequency_range  # [eV]
    # 输出： epsilon   # 无量纲复数
    #############################################
    frequency_range = frequency_range / 27.2114  # 从 eV 改到 原子单位
    dt = prop_time[2] - prop_time[1]
    fdamp = damping(prop_time, kick_delay, damping_factor)
    dAdt = dAprobe_dt * fdamp
    epsilon = np.zeros_like(frequency_range, dtype=complex)
    for w in range(len(frequency_range)):
        real, imag = 0, 0  # initialize the real and imag of FT
        for i in range(len(prop_time)):
            real += dAdt[i] * dt * np.cos(frequency_range[w] * prop_time[i])
            imag += dAdt[i] * dt * np.sin(frequency_range[w] * prop_time[i])
        eps_inv = 1 + (real + 1j * imag) / (kick_amplitude * np.exp(1j * frequency_range[w] * kick_delay))
        epsilon[w] = 1 / eps_inv
    return(epsilon)

def gauge_field_kick_mono(prop_time, kick_delay, kick_amplitude, dAprobe_dt, damping_factor, frequency):
    ############################################
    # 使用kick法，利用规范场计算介电函数
    # 输入参数 prop_time # [hbar/Ha]   总传播时间
    # kick_delay # [hbar/Ha]  kick延迟时间
    # kick_amplitude # [Ha]   kick大小
    # induced_Afield # [Ha]   感应规范场
    # induced_efield = -1 * cx.ReadOne("./gauge_field", 5) / cx.c_au  # E= -1/c  dA/dt  原子单位c=137
    # induced_efield   # Eprobe 感应电场，用原子单位  E =（-1/c）*（dA/dt）  原子单位c=137
    # Eind_si = Eprobe * 51.44  # 1 [Ha] = 27.2114 [eV], 1 [Bohr] = 0.529 [Angstrom]   27.2114/0.529 =51.44
    # damping_factor   #  衰减因子  eta [eV]
    # frequency_range  # [eV]  是一个值，不是列表
    # 输出： epsilon   # 无量纲复数
    #############################################
    frequency = frequency / 27.2114  # 从 eV 改到 原子单位
    dt = prop_time[2] - prop_time[1]
    fdamp = damping(prop_time, kick_delay, damping_factor)
    dAdt = dAprobe_dt * fdamp
    real, imag = 0, 0  # initialize the real and imag of FT
    for i in range(len(prop_time)):
        real += dAdt[i] * dt * np.cos(frequency * prop_time[i])
        imag += dAdt[i] * dt * np.sin(frequency * prop_time[i])
    eps_inv = 1 + (real + 1j * imag) / (kick_amplitude * np.exp(1j * frequency * kick_delay))
    epsilon = 1 / eps_inv

    return(epsilon)

def epsilon_kick(prop_time, kick_delay, kick_amplitude, induced_efield, damping_factor, frequency_range):
    ############################################
    # 使用kick法，利用规范场计算介电函数
    # 输入参数 prop_time # [hbar/Ha]   总传播时间
    # kick_delay # [hbar/Ha]  kick延迟时间
    # kick_amplitude # [Ha]   kick大小
    # induced_efield   # Eprobe 感应电场，用原子单位  E =（-1/c）*（dA/dt）  原子单位c=137
    # Eind_si = Eprobe * 51.44  # 1 [Ha] = 27.2114 [eV], 1 [Bohr] = 0.529 [Angstrom]   27.2114/0.529 =51.44
    # damping_factor   #  衰减因子  eta [eV]
    # frequency_range  # [eV]
    # 输出： epsilon   # 无量纲复数
    #############################################
    frequency_range = frequency_range / 27.2114  # 从 eV 改到 原子单位
    dt = prop_time[2] - prop_time[1]
    fdamp = damping(prop_time, kick_delay, damping_factor)
    dAdt = induced_efield * (-cx.c_au) * fdamp
    epsilon = np.zeros_like(frequency_range, dtype=complex)
    for w in range(len(frequency_range)):
        real, imag = 0, 0  # initialize the real and imag of FT
        for i in range(len(prop_time)):
            real += dAdt[i] * dt * np.cos(frequency_range[w] * prop_time[i])
            imag += dAdt[i] * dt * np.sin(frequency_range[w] * prop_time[i])
        eps_inv = 1 + (real + 1j * imag) / (kick_amplitude * np.exp(1j * frequency_range[w] * kick_delay))
        epsilon[w] = 1 / eps_inv
    return(epsilon)

def Get_Epsilon_Gauge(t_ext, A_ext, t_ind, A_ind, Omega,eta):
    ###########################
    # calculate complex dielecreic function epsilon by using gauge field method!
    ## Input parameters:
    #  t_ext  # [hbar/Ha]  time of external laser field
    #  A_ext  # [Ha]     gauge field amplitude of external laser field 
    #  t_ind  # [habr/Ha]  time of induced field 
    #  A_ind  # [Ha]     gauge field amplitude of induced field 
    #  Omega  # eV, can be a value 数 or  a list 列表
    #  eta  # eV,  a value  数
    ## Output parameters
    #   epsilon_real  # dimensionless
    #   epsilon_imag  # dimensionless
    ######################
    # 计算傅里叶变换
    dt = t_ext[2] - t_ext[1]
    dAdt_ext = cx.dA_over_dt(t_ext, A_ext)
    dAdt_ind = cx.dA_over_dt(t_ind, A_ind)
    omega_au = Omega * cx.timeau2s / cx.hbar  # omega_au  [Ha/hbar] * cx.timeau2s / cx.hbar = 1 / 27.2114
    eta_au = eta * cx.timeau2s / cx.hbar  # [Ha/hbar]
    ind_real, ind_imag = cx.FT_eta(t_ind, dAdt_ind, eta_au, omega_au)
    ext_real, ext_imag = cx.FT_eta(t_ind, dAdt_ext, eta_au, omega_au)
    # 计算介电函数 epsilon
    Denominator = np.square(ind_real + ext_real) + np.square(ind_imag + ext_imag)
    epsilon_real = (ext_real * (ind_real + ext_real) + ext_imag * (ind_imag + ext_imag)) / Denominator
    epsilon_imag = (ind_real * ext_imag - ext_real * ind_imag) / Denominator

    return(epsilon_real, epsilon_imag)

def Get_Epsilon_Current(time, A_ext, current, omega):
    ###########################
    # 第二种： 从电流出发
    # calculate complex dielecreic function epsilon by using current method!
    ## Input parameters:
    #  time  # [hbar/Ha]  propagation time
    #  A_ext # [Ha]  external field
    #  current # [au]    induced current  j(t)
    #  Omega  # eV,   can be a value or  a list
    ## Output parameters
    #epsilon_real  # dimensionless
    #epsilon_imag  # dimensionless
    ###########################
    dt = time[2] - time[1]  # dt [hbar / Ha]
    # octopus输出的电流J 需要乘 - 1 !
    current = current * (-1)
    omega_au = omega * cx.timeau2s / cx.hbar
    # J(t) FT J(w)
    Jw_real,Jw_imag = 0, 0
    for i in range(len(current)):
        Jw_real += current[i] * dt * np.cos(omega_au * time[i])  # J_w [au]
        Jw_imag += current[i] * dt * np.sin(omega_au * time[i])
    # A(t) FT A(w)
    Aw_real, Aw_imag = 0, 0  # initialize for FT
    for i in range(len(A_ext)):
        Aw_real += A_ext[i] * dt * np.cos(omega_au * time[i] )
        Aw_imag += A_ext[i] * dt * np.sin(omega_au * time[i] )
    # A(w) ->  E(w)
    E_imag = omega_au * Aw_real / 137
    E_real = -omega_au * Aw_imag / 137
    #conductivity 电导  sigma
    denominator = np.square(E_real)+np.square(E_imag)
    sigma_real = ( Jw_real* E_real + Jw_imag * E_imag ) / denominator
    sigma_imag = ( Jw_imag* E_real - Jw_real * E_imag ) / denominator
    #从电导求介电函数eps, 这里使用了原子单位没有4pi, 但cgs单位制多个4pi
    epsilon_real = 1 -  sigma_imag / omega_au
    epsilon_imag =   sigma_real / omega_au

    return(epsilon_real,epsilon_imag)


def Get_Epsilon_Current_damping(time, A_ext, current, omega, eta):
    ###########################
    # 第二种： 从电流出发
    # calculate complex dielecreic function epsilon by using current method!
    ## Input parameters:
    #  time  # [hbar/Ha]  propagation time
    #  A_ext # [Ha]  external field
    #  current # [au]    induced current  j(t)
    #  Omega  # eV,   can be a value or  a list
    #  eta # [Ha/hbar] damping parameter    eta = eta_ev / cx.hbar *  cx.timeau2s
    ## Output parameters
    #epsilon_real  # dimensionless
    #epsilon_imag  # dimensionless
    ###########################
    dt = time[2] - time[1]  # dt [hbar / Ha]
    # octopus输出的电流J 需要乘 - 1 !
    current = current * (-1)
    omega_au = omega * cx.timeau2s / cx.hbar
    # J(t) FT J(w)
    Jw_real,Jw_imag = 0, 0
    for i in range(len(current)):
        Jw_real += current[i] * dt * np.cos(omega_au * time[i]) * np.exp( - eta * time[i])  # J_w [au]
        Jw_imag += current[i] * dt * np.sin(omega_au * time[i]) * np.exp( - eta * time[i])
    # A(t) FT A(w)
    Aw_real, Aw_imag = 0, 0  # initialize for FT
    for i in range(len(A_ext)):
        Aw_real += A_ext[i] * dt * np.cos(omega_au * time[i] ) * np.exp( - eta * time[i])
        Aw_imag += A_ext[i] * dt * np.sin(omega_au * time[i] ) * np.exp( - eta * time[i])
    # A(w) ->  E(w)
    E_imag = omega_au * Aw_real / 137
    E_real = -omega_au * Aw_imag / 137
    #conductivity 电导  sigma
    denominator = np.square(E_real)+np.square(E_imag)
    sigma_real = ( Jw_real* E_real + Jw_imag * E_imag ) / denominator
    sigma_imag = ( Jw_imag* E_real - Jw_real * E_imag ) / denominator
    #从电导求介电函数eps, 这里使用了原子单位没有4pi, 但cgs单位制多个4pi
    epsilon_real = 1 -  sigma_imag / omega_au
    epsilon_imag =   sigma_real / omega_au

    return(epsilon_real,epsilon_imag)
####################  Optical relations ###################################
def absorption(energy,epsilon_real,epsilon_imag):
    # energy [eV]  hbar
    omega = energy / cx.hbar   # [1/s]
    absorption_coefficient = 1e-9 * np.sqrt(2) * omega / cx.Light_speed *  \
        np.sqrt(np.sqrt(np.square(epsilon_real) + np.square(epsilon_imag))-epsilon_real)
    return(absorption_coefficient) # [1/nm]


def BeerLambert(omega_si, n_imag, distance):
    #Input:  omega_si [1/s],  m_imag dimensionless, distance [m]
    #Output: 吸收系数 alpha [1/m]
    alpha = 2 * omega_si * n_imag / cx.Light_speed
    # 透射率的大小 T = I_out / I_in
    Transmission = np.exp(-1 * alpha * distance)
    return(alpha)

def Transmission(Omega, epsilon_real,epsilon_imag, distance):
    #Input:  Omega [eV],  m_imag dimensionless, distance [m]
    #Output: 吸收系数 alpha [1/m]
    omega_si = Omega / cx.hbar  # omega_si [1/s]
    n_real, n_imag = refractive_index(epsilon_real, epsilon_imag)
    alpha = 2 * omega_si * n_imag / cx.Light_speed
    # 透射率的大小 T = I_out / I_in
    T = np.exp(-1 * alpha * distance) # T = Iout/Iin = exp(-alpha d)
    return(T)

def n2eps(n,k):
    eps1 = np.square(n) - np.square(k)
    eps2 = 2 * n * k
    return(eps1,eps2)

def refractive_index(epsilon_real, epsilon_imag):
    # epsilon is np.array()
    #n_imag的正负号和epsilon_imag一致
    sign_epsilon_imag = []
    for i in range(len(epsilon_imag)):
        if (epsilon_imag[i] > 0) or (epsilon_imag[i] == 0):
            sign_epsilon_imag.append(1)
        else:
            sign_epsilon_imag.append(-1)
    sign_epsilon_imag = np.array(sign_epsilon_imag)

    modulus = np.sqrt(np.square(epsilon_real) + np.square(epsilon_imag))
    n_real = np.sqrt((modulus + epsilon_real) / 2)
    # n_imag的正负号和epsilon_imag一致
    n_imag = np.sqrt((modulus - epsilon_real) / 2) * sign_epsilon_imag

    return(n_real, n_imag)

def refractive_index_positive(epsilon_real, epsilon_imag):
    # epsilon is np.array().  n,k 全正
    modulus = np.sqrt(np.square(epsilon_real) + np.square(epsilon_imag))
    n_real = np.sqrt((modulus + epsilon_real) / 2)
    n_imag = np.sqrt((modulus - epsilon_real) / 2)
    return(n_real, n_imag)

def eps2n_fomo_630nmgold(eps1, eps2):
    # epsilon is np.array().  eps1<0, eps2>0, n>0, k>0
    modulus = np.sqrt(1 + np.square(eps1) + np.square(eps2))
    n = np.sqrt( (1 - modulus) * eps1 / 2)
    k = np.sqrt( (-1 - modulus) * eps1 / 2)
    return(n, k)



def simple_refractive_index(eps1, eps2):
    # epsilon is a value
    modulus = np.sqrt(np.square(eps1) + np.square(eps2))
    n_real = np.sqrt((modulus + eps1) / 2)
    # n_imag的正负号和epsilon_imag一致
    if eps2 > 0 :
        n_imag = np.sqrt((modulus - eps1) / 2)
    else:
        n_imag = np.sqrt((modulus - eps1) / 2) * (-1)
    return(n_real, n_imag)

def reflection(n,k):
    R = (np.square(n-1) + np.square(k))/ (np.square(n+1) + np.square(k)  )
    return(R)

#########################  画成热力图 ######################
def colormap(x,y,z):
    X, Y = np.meshgrid(x, y)
    wstart, wend = y[0], y[-1]
    mid =  "{:.1f}".format((wstart+wend)/2)
    Matrix2D_amp = z.reshape(len(list(x)), len(y)).T # 将拼接起来的数据重新排布成和画布相适应的矩阵
    fig = plt.figure(figsize = (4,9))
    ax = fig.add_subplot()
    heatmap = ax.imshow(Matrix2D_amp, interpolation='bilinear', cmap='jet', aspect='auto', extent=[x[0], x[-1],wend, wstart])
    ax.invert_yaxis()
    fs = 20
    ax.set_xlabel('$\\tau $[fs]', fontsize=fs)
    ax.set_xticks(list(x), labels= list(x),fontsize=fs)  # 设置x轴标签
    ax.set_yticks([wstart, (wstart+wend)/2,  wend], labels= [ "{:.1f}".format(wstart), mid,  "{:.1f}".format(wend)],fontsize=fs)  # 设置y轴标签
    cbar = plt.colorbar(heatmap, ax=ax)
    cbar.ax.tick_params(labelsize=fs)

    return()

####################################################################
def test():

    pass
    return()

if __name__ == "__main__" :
    test()
