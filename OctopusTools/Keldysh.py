import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx

###############  Keldysh Model  #######################################################
def Ionization_potential(e_kj_mol):
    # e [kJ/mol] to  E [eV]
    E = e_kj_mol / cx.Na / cx.e * 1E3
    return(E)

def Ponderomotive_potential(I,lamda):
    # Up[eV]  I [W/cm^2] lamda[m]
    Up = 0.09337 * I * np.square(lamda)
    return(Up)

def Precision(Hartree,atom_number):
    # precision [meV/atom]
    precision = Hartree * cx.hartree2ev * 1000 / atom_number
    return(precision)


def Keldysh1(Intensity,photon_energy):
    # Kp = np.sqrt(E_ip / (2*Up))  ZnO
    wavelength = cx.Energy2Wavelength(photon_energy)  # nm
    Up = Ponderomotive_potential(Intensity,wavelength*1E-9)  # eV
    Ip = 7.5  # eV  ionization potential
    Kp1 = np.sqrt(Ip/(2*Up))  #  Keldysh parameter
    return(Kp1)

def Keldysh2(Intensity,photon_energy,Refractive_index, reduced_electron_mass, bandgap):
    # Kp = omega/(e*E_0) * np.sqrt(m* * bandgap)
    Electric_field = cx.I2E(Intensity, Refractive_index) * 1E10  # V/m
    omega = photon_energy / cx.hbar  #  1/s
    Kp2 =  omega / (cx.e * Electric_field) * np.sqrt(reduced_electron_mass * bandgap * cx.e )   #  Keldysh parameter
    return(Kp2)


def gamma1(gamma):

    gamma1 = gamma / np.sqrt(1 + np.square(gamma))

    return(gamma1)


def gamma2(gamma):

    gamma2 = 1 / np.sqrt(1 + np.square(gamma))

    return (gamma2)


def  Elliptic_Integral_1_list(x):
# 第一类完全椭圆积分  Complete elliptic integral of the first kind
    k1 = []
    for i in range(len(x)):
        x = list(x)
        dtheta = 0.01
        theta = np.arange(0,np.pi/2,dtheta)
        kk = 0
        for j in range(len(theta)):
            kk += dtheta / np.sqrt( 1 - np.square(x[i]) * np.square( np.sin(theta[j])))
        k1.append(kk)

    k1 = np.array(k1)
    return(k1)

def  Elliptic_Integral_2_list(x):
# 第二类完全椭圆积分 	Complete elliptic integral of the second kind
    k2 = []
    for i in range(len(x)):
        dtheta = 0.01
        theta = np.arange(0,np.pi/2,dtheta)
        kk = 0
        for j in range(len(theta)):
            kk += dtheta * np.sqrt( 1 - np.square(x[i]) * np.square( np.sin(theta[j])))
        k2.append(kk)

    k2 = np.array(k2)
    return(k2)


def  Elliptic_Integral_1_value(x):
# 第一类完全椭圆积分  Complete elliptic integral of the first kind
    dtheta = 0.001
    theta = np.arange(0,np.pi/2,dtheta)
    k = 0
    for j in range(len(theta)):
        k += dtheta / np.sqrt( 1 - np.square(x) * np.square( np.sin(theta[j])))
    return(k)

def  Elliptic_Integral_2_value(x):
# 第二类完全椭圆积分 	Complete elliptic integral of the second kind
    dtheta = 0.001
    theta = np.arange(0,np.pi/2,dtheta)
    k = 0
    for j in range(len(theta)):
        k += dtheta * np.sqrt( 1 - np.square(x) * np.square( np.sin(theta[j])))
    return(k)

def PHI_X(x):
    # \Phi(x) = \int^x_0 exp(y^2 - x^2) dy.   x is a value not list
    dy = x/1000
    y = np.arange(0,x,dy)
    PHI = 0
    for i in range(len(y)):
        PHI += dy * np.exp(np.square(y[i])-np.square(x))
    return(PHI)


def Q(gamma,x):

    k1g1 = Elliptic_Integral_1_value(gamma1(gamma))
    k1g2 = Elliptic_Integral_1_value(gamma2(gamma))
    k2g1 = Elliptic_Integral_2_value(gamma1(gamma))
    k2g2 = Elliptic_Integral_2_value(gamma2(gamma))
    sum = 0
    cutoff = 200 # integer, cutoff of the sum fo n.
    for n in np.arange(0,cutoff): # n should be 0 to infinity, here there is a cutoff
        sum += np.exp(-np.pi * n * (k1g1 - k2g1) / k2g2) * \
            PHI_X(np.sqrt( (np.square(np.pi) *(2 * int(x+1) - 2*x + n)) \
            / (2*k1g2*k2g2) ) )
    Q = sum * np.sqrt(np.pi /(2 *k1g2))

    return(Q)

def Kelydsh_W(Intensity,photon_energy,band_gap,RefractiveIndex):
    ###################### Initial parameter setting ############
    #Intensity = 2.84E11  # W/cm^2
    #photon_energy = 1.55  # 0.41328  eV  ~~ 3000nm;  1.55eV ~ 800nm
    #band_gap = 2.8  # V
    ############################################################
    wavelength = cx.Energy2Wavelength(photon_energy)  # nm
    omega = cx.Wavelength2Omega(wavelength) # 1/s
    # reduced_electron_mass = 1.88 * cx.electron_mass  # kg    effective mass calculation, can be get from vaspkit : Gamma to L point
    reduced_electron_mass = 1.88 * cx.electron_mass  # kg
    gamma = Keldysh2(Intensity = Intensity, photon_energy = photon_energy, Refractive_index = RefractiveIndex, reduced_electron_mass = reduced_electron_mass ,bandgap = band_gap)
    g1 = gamma1(gamma)
    g2 = gamma2(gamma)
    k1g1 = Elliptic_Integral_1_value(g1)
    k2g1 = Elliptic_Integral_2_value(g1)
    k2g2 = Elliptic_Integral_2_value(g2)
    eff_ion_pot = (2 / np.pi) * (band_gap / g1) * k2g2

    # k-photon absorption   2 for octopus bandgap 2.8eV or 3 for experiment data 3.3eV
    k = int(1 + eff_ion_pot / photon_energy)
    finit_k = 1 + eff_ion_pot / photon_energy
    print(k)
    print("%.2f-photon absorption" % finit_k)
    first_term = 2 * omega / (9 * np.pi) # [1/s]
    coefficient =  reduced_electron_mass * omega / (g1 * cx.hbar * cx.e)
    second_term = np.power(coefficient,3/2)
    third_term = Q(gamma, eff_ion_pot / band_gap)  # dimensionless
    forth_term = np.exp(-np.pi * k * (k1g1 - k2g1) / k2g2)  # dimensionless

    W = first_term * second_term * third_term * forth_term   # [m^{-3}*s^{-1}]
    return(gamma,W)
#########################################################################
def test():

    pass
    return()

if __name__ == "__main__" :
    test()