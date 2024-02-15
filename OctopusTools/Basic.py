import numpy as np
import matplotlib.pyplot as plt
######### Basic Physical constants  #################################
Vaccum_permitivity = 8.8542E-12   # [Farad/meter]
Vacuum_permeability = 4*np.pi*1E-7     # 真空磁导率 [H/m], [V·s/(A·m)]
Light_speed = 299792458 # [Meter/Second]
Impedance= Light_speed * Vacuum_permeability     # vacuum impedance Z [H/s]  真空阻抗
h = 4.1356676969 * 1e-15  # planck constant [eV·s]
hbar = h / (2 * np.pi)  # reduced planck constant  [eV·s]
e = 1.60217662E-19  # [Coulombs] elements charges
Na = 6.02214076E23   # Avogadro constant
electron_mass = 9.10938356E-31  # kg
################ Unit Conversion constant###############################
bohr2m = 5.2917721067e-11   # 1 [Bohr]    = 5.29E-11 [m]
timeau2s = 2.41888435E-17   # 1 [hbar/Hartree] =  2.41888435E-17 [s]
timeau2fs = 2.41888435E-2    # 1 [hbar/Ha] = 0.0241888435 [fs]
c_au = 137.036 # light speed in atom unit
hartree2ev = hbar/timeau2s      # 1 [Hartree] = 27.211386 [eV]   =   cx.hbar /cx.timeau2s = 27.211385982430674
################ Read&Save data from file  ################################
def ReadData(filename,column1,column2):
    f = np.loadtxt(filename, encoding="UTF-8")
    x = f[:,column1]
    y = f[:,column2]
    return(x,y)

def ReadOne(filename,column):
    f = np.loadtxt(filename, encoding="UTF-8")
    x = f[:,column]
    return(x)

def ReadData_xyz(filename,column1,column2,column3, ):
    f = np.loadtxt(filename, encoding="UTF-8")
    x = f[:,column1]
    y = f[:,column2]
    z = f[:,column3]
    return(x,y,z)

def save(filename,x,y, FirstLineComment):
    f = open(filename,"w", encoding="UTF-8")
    f.write("#####" + FirstLineComment + "\n")
    for i in range(len(x)):
        f.write("%.4f            %.8f    \n" % (x[i],y[i]))
    f.close()
    return()

def save_xyz(filename,x,y,z, FirstLineComment):
    f = open(filename,"w", encoding="UTF-8")
    f.write("#####" + FirstLineComment + "\n")
    for i in range(len(x)):
        f.write("%.4f                %.4f        %.4f\n" % (x[i],y[i],z[i]))
    f.close()
    return()

def Align(x1,y1,x2,y2):
    #Align two sets of data {x1,y1} and {x2,y2} from start
    count = min(len(x1),len(x2))
    newx1 = x1[0:count]
    newy1 = y1[0:count]
    newx2 = x2[0:count]
    newy2 = y2[0:count]
    return(newx1, newy1, newx2, newy2)

def shortten(x1,y1,start,end):
    x2 = x1[start:end]
    y2 = y1[start:end]
    return(x2,y2)

def GetIndex(value, list):
# In a monotonically increasing "list", to find the index of "value".
    count = 0
    for i in range(len(list)):
        if value > list[i] :
            count += 1
        else:
            break
    index = count
    return index

def sgn_list(list_x):
# sign function. x>0, y=1; x<=0, y=-1
    list_y = []
    for i in range(len(list_x)):
        if list_x[i] > 0 :
            list_y.append(1)
        else:
            list_y.append(-1)

    return(np.array(list_y))


def sgn(x):
# sign function. x>0, y=1; x<=0, y=-1
    if x > 0 :
        return(1)
    else:
        return(-1)

def sgn_list_2(t0,tau0,t):
# sign function: If abs(t-t0)>tau0 then f(t) = 0
    sgn_function = []
    for i in range(len(t)):
        if abs(t[i]-t0) > tau0 :
            sgn_function.append(0)
        else:
            sgn_function.append(1)
    sgn_function = np.array(sgn_function)
    return(sgn_function)
######### Some conversions: wavelength&omega, wavelength&photonenergy, E&I, E&A  ###############
def Wavelength2Omega(wavelength):
    # Wavelength [nm]  to  Omega[1/s]
    return 2 * np.pi * Light_speed * 1E9   / wavelength

def Omega2Wavelength(omega):
    # Omega[1/s]  to  Wavelength [nm]
    return  2 * np.pi * Light_speed * 1E9   / omega

def Energy2Wavelength(energy):
    # Energy [eV]  to Wavelength [nm]
    return 2 * np.pi * Light_speed * 1E9 * hbar  / energy

def Wavelength2Energy(wavelength):
    # wavelength [nm] to  Energy [eV]
    energy = 2 * np.pi * Light_speed * 1E9 * hbar  / wavelength
    return energy

def E2I(E, Refractive_index):
    # E [V/Angstrom]  to  I [W/cm^2]
    I = 0.5 * Refractive_index * Vaccum_permitivity * Light_speed * np.square(E)  * 1E16
    return I

def I2E(I, Refractive_index):
    # I [W/cm^2]  to  E [V/Angstrom]      ,   1[V/cm] = 1E-8 [V/Angstrom]
    E = np.sqrt(2 * I / (Refractive_index * Vaccum_permitivity * Light_speed)) *1E-8
    return(E)

def A2E(t,A):
    # t [hbar/Hartree]  A [Hartree] -->  t [hbar/Hartree]  E [V/Angstrom]  $\mathring{A}$
    # E = -1/c * dA/dt
    dt = t[2] - t[1]   # 1 [Ha] = 27.2114 [eV], 1 [Bohr] = 0.529 [Angstrom]
    E=[0]
    for i in range(len(A)-1):
        dA = A[i+1] - A[i]
        e =  (-1/c_au) * dA / dt  * 51.44
        E.append(e)
    E = np.array(E)   # V/Angstrom
    return(E)

def A2Eau(t,A):
    # t [hbar/Hartree]  A [Hartree] -->  t [hbar/Hartree]  E [a.u.]
    # E = -1/c * dA/dt  # 1 [Ha] = 27.2114 [eV], 1 [Bohr] = 0.529 [Angstrom]
    dt = t[2] - t[1]   # a.u.
    E=[0]
    for i in range(len(A)-1):
        dA = A[i+1] - A[i]
        e =  (-1/c_au) * dA / dt
        E.append(e)
    Eau = np.array(E)   # a.u.
    return(Eau)

def E2A(time,E,A0):
    # t [hbar/Hartree]  E [V/Angstrom]  -->  t [hbar/Hartree]  A [Hartree]
    # A = -137.131/51.44* \int E dt
    A = []
    A_field = 0
    dt = time[1] - time[0]
    for i in range(len(time)):
        A_field += -137.131/51.44*dt*E[i]   # 51.44= 27.2114/0.529
        A.append(A_field)
    A = np.array(A)
    A = A + A0
    return(A)

def dA_over_dt(t,A):
    # t [hbar/Hartree]  A [Hartree]
    dt = t[2] - t[1]
    dAdt=[0]
    for i in range(len(A)-1):
        dA = A[i+1] - A[i]
        e =  dA / dt
        dAdt.append(e)
    dAdt= np.array(dAdt)   # Ha^2 / hbar
    return(dAdt)

def current_cal(t_ind, A_ind, V):
    # V [b^3]
    dAdt = dA_over_dt(t_ind, A_ind)
    dAdt2 = dA_over_dt(t_ind, dAdt)
    J_cal = dAdt2 * V / (4 * np.pi * c_au)
    return(J_cal)

def positivecomplex(complex1):
    return np.abs(complex1.real) + 1j * np.abs(complex1.imag)


def RatioDampingMethod(t, t0, eta, method):
    ## t [a.u.]  eta [eV]  method: exponential, gaussian, polynomial
    eta_au = eta / hartree2ev # a.u.
    f = []
    for i in range(len(t)):
        if t[i]<t0:
            f.append(1)
        else:
            if method == "exponential":
                f0 = np.exp( -eta_au * (t[i] - t0) )
            elif method == "gaussian":
                f0 = np.exp( - np.square(eta_au * (t[i] - t0)) )
            elif method == "polynomial":
                f0 = 1 - 3 * np.square(eta_au * (t[i] - t0)) + 2 * np.power(eta_au * (t[i] - t0),3)
            else:
                f0 = 1
                print("Please input the right daming method: exponential, gaussian or polynomial")
            f.append(f0)
    f = np.array(f)
    return(f)

def DampingMethod(t, eta, method):
    ## t [a.u.]  eta [eV]  method: exponential, gaussian, polynomial
    eta_au = eta / hartree2ev # a.u.

    if method == "exponential":
        f = np.exp( -eta_au * t )
    elif method == "gaussian":
        f = np.exp( - np.square(eta_au * t) )
    elif method == "polynomial":
        f = 1 - 3 * np.square(eta_au * t) + 2 * np.power(eta_au * t,3)
    else:
        f = np.ones_like(t)
        print("Please input the right daming method: exponential, gaussian or polynomial")
    return(f)

def gauge_kick(t, A0, A_ind, w_ev, damp):   # 电场法 gauge field kick场
    w_au = w_ev / hartree2ev  # 从 eV 改到 原子单位
    dt = t[2] - t[1]
    dAdt_ind = dA_over_dt(t, A_ind)
    Indw = 0 # 初始化
    for i in range(len(dAdt_ind)):
        Indw += dAdt_ind[i] * dt * np.exp(w_au * t[i] * 1j) * damp[i]
    ### 计算epsilon
    inv_epsilon = 1 +  Indw / A0
    epsilon = 1 / inv_epsilon
    return(epsilon)


def gauge_kick_dAdt(t, A0, dAdt, w_ev):   # 电场法 gauge field kick场
    w_au = w_ev / hartree2ev  # 从 eV 改到 原子单位
    dt = t[2] - t[1]
    Indw = 0 # 初始化
    for i in range(len(dAdt)):
        Indw += dAdt[i] * dt * np.exp(w_au * t[i] * 1j)
    ### 计算epsilon
    inv_epsilon = 1 +  Indw / A0
    epsilon = 1 / inv_epsilon
    return(epsilon)

def gauge_kick_delay(t, A0, dAdt, w_ev, t0):   # 电场法 gauge field kick场, 延迟t0
    w_au = w_ev / hartree2ev  # 从 eV 改到 原子单位
    dt = t[2] - t[1]
    Indw = 0 # 初始化
    for i in range(len(dAdt)):
        Indw += dAdt[i] * dt * np.exp(w_au * t[i] * 1j)
    ### 计算epsilon
    inv_epsilon = 1 +  Indw / (A0 * np.exp(w_au * t0 * 1j ))
    epsilon = 1 / inv_epsilon
    return(epsilon)


def gauge_mono(t, A_ext, A_ind, w_ev,damp):  # 电场法 激光单色场
    dt = t[2] - t[1]
    w_au = w_ev / hartree2ev
    dAdt_ext = dA_over_dt(t, A_ext)
    dAdt_ind = dA_over_dt(t, A_ind)
    Extw, Indw = 0, 0  # initialize for Fourier Transform
    for i in range(len(A_ind)):
        Extw += dAdt_ext[i] * dt * np.exp(w_au * t[i] * 1j)
        Indw += dAdt_ind[i] * dt * np.exp(w_au * t[i] * 1j) * damp[i]
    inv_epsilon = 1 +  Indw / Extw
    epsilon = 1 / inv_epsilon
    return(epsilon)

def current_step(t, A0, j,V, w_ev, damp):  # 电流法  Ekick = Astep场
    # time[a.u.] A0[Ha] j [a.u.] V [bohr^3] omega_ev [eV] damp f(t)
    dt = t[2] - t[1]
    w_au = w_ev /hartree2ev # a.u.
    Jw = 0  # initialize for Fourier Transform
    for i in range(len(j)):
        Jw += j[i] * dt * np.exp(w_au * t[i] * 1j) * damp[i]
    sigma = Jw * c_au / (A0) # 负号未引入， 因为 E0 = -A0/c
    epsilon = (1 + 4 * np.pi * 1j * sigma / (w_au * V))
    return(epsilon)

def current_mono(t, A, J, w_ev, V, damp): # 电流法 激光单色场
    dt = t[2] - t[1]
    w_au = w_ev / hartree2ev  # a.u.
    Et = A2Eau(t,A) # a.u.
    Ew, Jw = 0,0 # initialize for Fourier Transform
    for i in range(len(A)):
        Ew += Et[i] * dt * np.exp(w_au * t[i] * 1j)  # 外场无需加修正
        Jw += J[i] * dt * np.exp(w_au * t[i] * 1j) * damp[i] # 感应场杂峰多，加修正
    sigma = Jw / (Ew * V)
    return(sigma)


def n2eps(n,k):
    eps1 = np.square(n) - np.square(k)
    eps2 = 2 * n * k
    return(eps1,eps2)

def eps2n(eps):
    # epsilon is np.array().  n,k 全正
    n_real = np.sqrt(( np.abs(eps) + eps.real) / 2)
    n_imag = np.sqrt(( np.abs(eps) - eps.real) / 2)
    return(n_real, n_imag)

def absorption(wavelength, eps):
    # wavelength [nm]
    absorption_coefficient = 4 * np.pi * np.sqrt((np.abs(eps)-eps.real) / 2) / (wavelength)
    return(absorption_coefficient) # [1/nm]

def absorption_ev(w_ev, eps):
    # wavelength [nm]
    wavelength = Energy2Wavelength(w_ev)  # nm
    absorption_coefficient = 4 * np.pi * np.sqrt((np.abs(eps)-eps.real) / 2) / (wavelength)
    return(absorption_coefficient) # [1/nm]

##################### Fourier Transform ###################################################
def FT(t, dAdt,omega):
# omega can be one value, or even list.   [Ha/hbar]
# Eta is for the imaginary part of response
# t [hbar/Ha]; dA [Hartree]
# return the real and imagnary part of FT result.
    dt = t[2] - t[1]  # dt [hbar / Ha]
    a = 0  # initialize the real of FT
    b = 0  # initialize the imag of FT
    for i in range(len(t)):
        a += dAdt[i] * dt  * np.cos(omega * t[i])
        b += dAdt[i] * dt  * np.sin(omega * t[i])
    return(a,b)

def FT_eta(t, dAdt,Eta,omega):
# omega can be one value, or even list.   [Ha/hbar]
# Eta is for the imaginary part of response [Ha/hbar]
# t [hbar/Ha]; dA [Hartree]
# return the real and imagnary part of FT result.
    dt = t[2] - t[1]  # dt [hbar / Ha]
    a = 0  # initialize the real of FT
    b = 0  # initialize the imag of FT
    for i in range(len(t)):
        a += dAdt[i] * dt * np.exp(-Eta * t[i]) * np.cos(omega * t[i])
        b += dAdt[i] * dt * np.exp(-Eta * t[i]) * np.sin(omega * t[i])
    return(a,b)

####################### Vectors Calculation #######################      
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in degree between vectors 'v1' and 'v2'
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)) * 180 / np.pi
########################## 字典排序 #####################################
def dict_sort_by_key(dict):
    #按键(key)排序
    key = dict.keys()
    value = dict.values()

    new_dict = {}
    new_key = sorted(key)

    for i in range(len(new_key)):
        new_dict[new_key[i]] = dict[new_key[i]]

    new_dict_key = new_dict.keys()
    new_dict_value = new_dict.values()

    return(new_dict,new_dict_key,new_dict_value)

def TwoList_sort_by_list1(list1,list2):
    #按键(key)排序
    dictionary = dict(zip(list1,list2)) # connect two lists to one dictionary, to fix the relation

    key = dictionary.keys()
    value = dictionary.values()

    new_dict = {}
    new_key = sorted(key)

    for i in range(len(new_key)):
        new_dict[new_key[i]] = dictionary[new_key[i]]

    new_list1 = list(new_dict.keys())
    new_list2 = list(new_dict.values())

    return(new_list1, new_list2)

def dict_sort_by_value(dict):
    #按值(value)排序, 值不重复
    key = list(dict.keys())
    value = list(dict.values())

    new_dict = {}
    new_value = sorted(value)

    for i in range(len(new_value)):
        index = value.index(new_value[i])
        new_dict[key[index]] = new_value[i]

    new_dict_key = new_dict.keys()
    new_dict_value = new_dict.values()

    return(new_dict,new_dict_key,new_dict_value)

#### N-th order polynomial fitting  ##################################
def nfit(x, y, n):

    if n == 0:
        fitted_y = y
    else:
        x = np.array(x)
        fitting = np.polyfit(x, y, n)
        parameters = np.poly1d(fitting)
        fitted_y = 0
        for i in range(n + 1):
            fitted_y += np.power(x, i) * parameters[i]

    return (parameters, fitted_y)


def spacing2cutoff(spacing):
    # spacing [bohr]   electron_mass [kg]
    energy_cutoff = np.square(hbar) * e  * np.pi / (2 * electron_mass * spacing * bohr2m)
    return(energy_cutoff)





#  打印变量名, but 当多个变量名指向同一个值的时候就会发错名字混乱
def var_name(var, all_var=locals()):
    return [var_name for var_name in all_var if all_var[var_name] is var][0]


#########################################################################
def test():
    print(spacing2cutoff(spacing = 0.2))

    pass
    return()

if __name__ == "__main__" :
    test()

