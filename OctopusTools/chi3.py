import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.Laser as la
#### N-th order polynomial fitting  ##################################
def nfit(x, y, n):
    x = np.array(x)
    fitting = np.polyfit(x, y, n)
    parameters = np.poly1d(fitting)
    fitted_y = 0
    for i in range(n + 1):
        fitted_y += np.power(x, i) * parameters[i]

    return (fitted_y, parameters)
##########################################################################
def Get_Epsilon_Gauge(time, A_ext, A_ind, Omega, eta):
    ###########################
    # calculate complex dielecreic function epsilon by using gauge field method!
    ## Input parameters:
    #  time # [hbar/Ha]  time
    #  A_ext  # [Ha]     gauge field amplitude of external laser field
    #  A_ind  # [Ha]     gauge field amplitude of induced field
    #  Omega  # eV, can be a value or  a list
    ## Output parameters
    #   epsilon_real  # dimensionless
    #   epsilon_imag  # dimensionless
    ######################
    # 计算傅里叶变换
    dt = time[2] - time[1]
    dAdt_ext = cx.dA_over_dt(time, A_ext)
    dAdt_ind = cx.dA_over_dt(time, A_ind)
    omega_au = Omega * cx.timeau2s / cx.hbar  # omega_au  [Ha/hbar]
    ind_real, ind_imag = cx.FT_eta(time, dAdt_ind, eta, omega_au)
    ext_real, ext_imag = cx.FT_eta(time, dAdt_ext, eta, omega_au)
    # 计算介电函数 epsilon
    Denominator = np.square(ind_real + ext_real) + np.square(ind_imag + ext_imag)
    epsilon_real = (ext_real * (ind_real + ext_real) + ext_imag * (ind_imag + ext_imag)) / Denominator
    epsilon_imag = (ind_real * ext_imag - ext_real * ind_imag) / Denominator

    return(epsilon_real, epsilon_imag)
###########################################################################################
def Get_Epsilon_Current(time, A_ext, current, omega, eta):
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
        Jw_real += current[i] * dt * np.cos(omega_au * time[i]) * np.exp( -1 *  eta * time[i]) # J_w [au]
        Jw_imag += current[i] * dt * np.sin(omega_au * time[i]) * np.exp( -1 *  eta * time[i])
    # A(t) FT A(w)
    Aw_real, Aw_imag = 0, 0  # initialize for FT
    for i in range(len(A_ext)):
        Aw_real += A_ext[i] * dt * np.cos(omega_au * time[i] ) * np.exp( -1 *  eta * time[i])
        Aw_imag += A_ext[i] * dt * np.sin(omega_au * time[i] ) * np.exp( -1 *  eta * time[i])
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
########################################################################
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
#########################################################################
def nIcurve(path, lasermode, Timestep, Propagationtime, Tail, wavelength, n00, Env, IntensityMode, plot):
    omega = cx.Wavelength2Omega(wavelength) * cx.timeau2s  # [1/s] --> 1 / [hbar/Ha]
    I = []
    n = []
    n_im = []
    ep1 = []
    ep2 = []
    e0 = cx.Vaccum_permitivity  # [F/m]

    if plot  == True :
        plt.figure()
    else:
        pass
    ##########################################################################################
    for i in range(len(path)):
        ## 读取数据
        if lasermode == 1:
            t_ext, A_ext = cx.ReadData(path[i] + "/laser", 0, 1)
        else:
            t_ext, A_ext = cx.ReadData(path[i] + "/laser", 1, 2)

        t_ind, A_ind = cx.ReadData(path[i] + "/gauge_field", 1, 2)
        # 特定时间的 截断
        start = 0
        end = int(Propagationtime/Timestep)
        t_ext, A_ext = cx.shortten(t_ext, A_ext, start, end)
        t_ind, A_ind = cx.shortten(t_ind, A_ind, start, end)

        ####### 画图: 加载包络之前的脉冲波形 #######################################################
        if plot == True:
            plt.subplot( len(path),2,  i*2+1)
            plt.plot(t_ext, A_ext, label = "A%.2f:Origin" % max(A_ext))
            plt.legend(fontsize=12, loc = 2 )
            if i == (len(path) - 1):
                plt.xlabel("Propagation time [$\hbar/Ha$]", fontsize=14)
                plt.ylabel(" $A_{ext}$ ", fontsize=16)
            else: pass

            plt.subplot(len(path), 2, i*2+2)
            plt.plot(t_ind * 0.02419, A_ind)
            if i == (len(path) - 1):
                plt.xlabel("Propagation time [fs]", fontsize=14)
                plt.ylabel(" $A_{ind}$", fontsize=16)
            else: pass
        ##################################################################################

        # 包络, 单位幅度大小为 1
        time, envy = la.env_sin2( mid=(Propagationtime - Tail), tail=Tail, dt=Timestep)
        # 判定是否加载包络
        if Env == True :
            A_ext = A_ext * envy
            A_ind = A_ind * envy
        else:
            A_ext = A_ext
            A_ind = A_ind
            time = t_ext
        ####### 画图: 加载包络之后的脉冲波形 #######################################################
        if plot == True:
            plt.subplot( len(path),2,  i*2+1)
            plt.plot(time, A_ext, label = "Enveloped")
            plt.legend(fontsize=12, loc = 2 )
            if i == (len(path) - 1):
                plt.xlabel("Propagation time [$\hbar/Ha$]", fontsize=14)
                plt.ylabel(" $A_{ext}$ ", fontsize=14)
            else: pass

            plt.subplot(len(path), 2, i*2+2)
            plt.plot(time * 0.02419, A_ind)
            if i == (len(path) - 1):
                plt.xlabel("Propagation time [fs]", fontsize=14)
                plt.ylabel(" $A_{ind}$", fontsize=16)
            else: pass
        ##############################################################################

        #  计算 dA/dt
        dAdt_ext = cx.dA_over_dt(time, A_ext)
        dAdt_ind = cx.dA_over_dt(time, A_ind)

        #  计算傅里叶变换
        ext1, ext2 = cx.FT(time, dAdt_ext, omega)
        ind1, ind2 = cx.FT(time, dAdt_ind, omega)

        #   计算介电函数
        Denominator = np.square(ind1 + ext1) + np.square(ind2 + ext2)
        epsilon_real = (ext1 * (ind1 + ext1) + ext2 * (ind2 + ext2)) / Denominator
        epsilon_imag = (ind1 * ext2 - ext1 * ind2) / Denominator

        # #   计算折射率
        n_real, n_imag = cx.refractive_index(epsilon_real, epsilon_imag)

        ############### 计算光强！！！  #####################
        #   计算光强  方法 1
        Emax = max(cx.A2E(time, A_ext))
        I1 = cx.E2I(Emax, 1)

        #   计算光强  方法2
        c = 137.03604
        Iau = max(A_ext) * omega / c
        I2 = np.square(Iau * np.sqrt(3.509470 * 1e16))

        #  计算光强  方法3
        I3 = cx.E2I(Emax, n00)

        #  计算光强  方法4
        modification =  epsilon_real
        Emax = max(cx.A2E(time, A_ext)) / modification  # epsilon modification
        I4= cx.E2I(Emax, n_real)   # [W/cm^2]  # cited from the 'n0n2k0k2' paper

        #  计算光强  方法5
        modification =  - 12
        Emax = max(cx.A2E(time, A_ext)) / modification  # epsilon modification
        I5= cx.E2I(Emax, n_real)   # [W/cm^2]  # cited from the 'n0n2k0k2' paper

        if IntensityMode == 1:
            I.append(I1)
        elif IntensityMode == 2:
            I.append(I2)
        elif IntensityMode == 3:
            I.append(I3)
        elif IntensityMode == 4:
            I.append(I4)
        elif IntensityMode == 5:
            I.append(I5)
        else:
            print("请选择正确的光强计算方式")

        n.append(n_real)
        n_im.append(n_imag)
        ep1.append(epsilon_real)
        ep2.append(epsilon_imag)

        print("I is %.2e W/cm^2,  epsilon1 is %.4f, epsilon2 is %.4f,  n = %.4f, k = %.4f"
              % (I[i], epsilon_real, epsilon_imag, n_real, n_imag))

    return (I, n, n_im, ep1, ep2 )


def GetChi3(path, lasermode, Timestep, Propagationtime, Tail, wavelength, n00, Env, IntensityMode, plot):
    I, n, n_im,  ep1, ep2 = nIcurve(path, lasermode, Timestep, Propagationtime, Tail, wavelength, n00, Env, IntensityMode, plot)
    # 线性拟合
    fitted_yr, pr = nfit(I, n, 1)
    fitted_yi, pi = nfit(I, n_im, 1)
    fitted_ep1, pep1 = nfit(I, ep1, 1)
    fitted_ep2, pep2 = nfit(I, ep2, 1)
    # 提取拟合参数， 斜率 n2， k2 和 截距 n0， k0
    n0 = pr[0]
    n2 = pr[1] * 1E-4  # [m^2/W]

    k0 = pi[0]
    k2 = pi[1] * 1E-4  # [m^2/W]
    #  计算三阶电极化率
    chi3r = ( 1 * n0 * cx.Vaccum_permitivity  * cx.Light_speed / 3) * (n0*n2 - k0*k2)  # [m^2/V^2]
    chi3i = ( 1 * n0 * cx.Vaccum_permitivity  * cx.Light_speed / 3) * (n0*k2 + k0*n2)  # [m^2/V^2]
    # 输出结果
    print("n0 = %.2f , n2 = %.2e m^2/W" % (n0, n2))
    print("k0 = %.4f, k2 = %.2e m^2/W" % (k0, k2))
    print("epreal0 = %.4f, epimag0 = %.4f " % (pep1[0], pep2[0]))
    print("chi3r = %.4e [m^2/V^2], chi3i = %.4e [m^2/V^2]" % (chi3r, chi3i))


    if plot == True:
        # 画图: epsilon-I 原始点 和 线性拟合
        plt.figure()
        plt.subplot(1, 2, 1)
        plt.plot(I, ep1, "o")
        plt.plot(I, fitted_ep1, "-")
        plt.xlabel("Intensity [$W/cm^2$]", fontsize = 14)
        plt.ylabel(" $\epsilon_{real}$ ", fontsize=16)

        plt.subplot(1, 2, 2)
        plt.plot(I, ep2, "o")
        plt.plot(I, fitted_ep2, "-")
        plt.xlabel("Intensity [$W/cm^2$]", fontsize = 14)
        plt.ylabel(" $\epsilon_{imag}$", fontsize=16)

        # 画图: n-I 原始点 和 线性拟合
        plt.figure()
        plt.subplot(1, 2, 1)
        plt.plot(I, n, "o")
        plt.plot(I, fitted_yr, "-")
        plt.xlabel("Intensity [$W/cm^2$]", fontsize = 14)
        plt.ylabel(" n ", fontsize=16)

        plt.subplot(1, 2, 2)
        plt.plot(I, n_im, "o")
        plt.plot(I, fitted_yi, "-")
        plt.xlabel("Intensity [$W/cm^2$]", fontsize = 14)
        plt.ylabel(" k ", fontsize=16)

        plt.show()
    else: pass

    return(n0,n2,k0,k2,chi3i, chi3r)

###########
def TimeEnvolution(path, lasermode, Timestep, Time, Tail, wavelength, n00, Env, IntensityMode, plot):

    N0,N2,K0,K2,CHI3R,CHI3I = [],[],[],[],[],[]

    for i in range(len(Time)):

        I, n, n_im, ep1, ep2 = nIcurve(path, lasermode, Timestep, Time[i], Tail, wavelength, n00, Env,
                                       IntensityMode, plot)
        # 线性拟合
        fitted_yr, pr = nfit(I, n, 1)
        fitted_yi, pi = nfit(I, n_im, 1)
        fitted_ep1, pep1 = nfit(I, ep1, 1)
        fitted_ep2, pep2 = nfit(I, ep2, 1)
        # 提取拟合参数， 斜率 n2， k2 和 截距 n0， k0
        n0 = pr[0]
        n2 = pr[1] * 1E-4  # [m^2/W]
        # n2 = np.abs(n2)  # !!!!! 需要修改并说明 ！！！！
        k0 = pi[0]
        k2 = pi[1] * 1E-4  # [m^2/W]
        #  计算三阶电极化率
        chi3r = ( 1 * n0 * cx.Vaccum_permitivity  * cx.Light_speed / 3) * (n0*n2 - k0*k2)  # [m^2/V^2]
        chi3i = ( 1 * n0 * cx.Vaccum_permitivity  * cx.Light_speed / 3) * (n0*k2 + k0*n2)  # [m^2/V^2]
        N0.append(n0)
        N2.append(n2)
        K0.append(k0)
        K2.append(k2)
        CHI3R.append(chi3r)
        CHI3I.append(chi3i)

    return(N0,N2,K0,K2,CHI3R,CHI3I)




#################################################################################################
def test():

    pass
    return()

if __name__ == "__main__" :
    test()
