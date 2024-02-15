import numpy as np
import matplotlib.pyplot as plt
import MyPackage.Basic_Function  as my


h = 4.1356676969 * 1e-15  # planck constant eV·s
hbar = h / (2 * np.pi)  # reduced planck constant  eV·s
n0 = 2.4  # linear refractive index of diamond
# Eig = 5.54 # [eV] indirect band gap energy of diamond
Light_speed = 3E8 # m/s   meter per second

def F_function(i,x):
    # i = 0,1,2 ,  x is value, not list
    F = np.power((2*x-1),(i+2)) / np.power(2*x,5)

    return(F)


def TPA_coefficient_beta(wavelength,Eig = 5.54):
    # wavelength  [nm] is a value, not a list
    omega = my.Wavelength2Omega(wavelength)  # [1/s]
    if omega * hbar / Eig > 0.5 :

        k0 = 62.86E-9  # m* eV^3/W
        k1 = 5.71E-9  # m* eV^3/W
        k2 = 0    # m* eV^3/W
        n = hbar * omega / Eig   # list[]
        beta = k0 / ( np.power(n0,2) * np.power(Eig,3) ) * F_function(0, n) + \
               k1 / (np.power(n0, 2) * np.power(Eig, 3)) * F_function(1, n) + \
               k2 / (np.power(n0, 2) * np.power(Eig, 3)) * F_function(2, n)
        # beta  [m/W]
        beta = beta * 1E11 # [cm/GW]
    else:
        beta = 0
    return(beta)


def TPA_coefficient_beta_change_K(wavelength,Eig = 5.54):
    # increasing factor 3.5, due to the use of the degenerated 2PA at the mean frequency, instead of the nondegenerated one!!!
    # wavelength  [nm] is a value, not a list
    omega = my.Wavelength2Omega(wavelength)  # [1/s]
    if omega * hbar / Eig > 0.5 :

        k0 = 62.86E-9  # m* eV^3/W
        k0 = k0 * 3.5  # modification
        k1 = 5.71E-9  # m* eV^3/W
        k1 = k1 * 3.5  # modification
        k2 = 0    # m* eV^3/W
        n = hbar * omega / Eig   # list[]
        beta = k0 / ( np.power(n0,2) * np.power(Eig,3) ) * F_function(0, n) + \
               k1 / (np.power(n0, 2) * np.power(Eig, 3)) * F_function(1, n) + \
               k2 / (np.power(n0, 2) * np.power(Eig, 3)) * F_function(2, n)
    else:
        beta = 0

    return(beta)  # beta  [m/W]

def Kerr_n2(wavelength):
    # wavelength   [nm]
    n2 = 0  # initialization
    d = 1 # [nm]
    wave = np.arange(1,450,d)     # [nm]

    for w in wave:
        n2 += TPA_coefficient_beta_change_K(w) * d / (w-wavelength)
        # n2 += TPA_coefficient_beta_change_K(w)
    n2 = n2 * (1/np.pi)

    return(n2)


#########################################################################
def test():
    # wavelength = np.arange(450,501,1)  # [nm]
    # y = []
    # for wave in wavelength:
    #     y.append(TPA_coefficient_beta_change_K(wave,Eig = 5.54))
    # print(y)
    # plt.plot(wavelength, y, "--o")
    # plt.show()
    wavelength = 400
    y = Kerr_n2(wavelength)
    print("%.2e" % y)

    return()

if __name__ == "__main__" :
    test()

