import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as my
from mpl_toolkits.mplot3d import Axes3D

########################################################
def n_Level_Legendre_Polynomial(x,n):
# 勒让德多项式
    P = np.zeros_like(x)
    for k in range(0,int(n/2)+1):
        coefficient = np.power(-1,k) * np.math.factorial(2*n-2*k) / (np.power(2,n) \
                    * np.math.factorial(k) * np.math.factorial(n-k) \
                    * np.math.factorial(n- 2*k))
        P +=  coefficient * np.power(x, (n-2*k))

    return(P)

def test_n_Level_Legendre_Polynomial():
    x = np.arange(-1,1,0.01)
    for n in [0,1,2,3,4,5]:
        P = n_Level_Legendre_Polynomial(x, n)
        plt.plot(x,P, label = n)

    plt.legend()
    plt.show()
    return()



def Associated_Legendre_Polynomials(x,l,m):
#  伴随(缔合/连带/关联)勒让德多项式
# 仅当ℓ和m均为整数且满足0≤m≤ℓ时,才在区间[−1, 1]上有非奇异解，
    P0 = n_Level_Legendre_Polynomial(x,l) # m=0
    for i in range(abs(m)):
        x,y,P0 = derivative_1(x, P0)
    P0 = P0 * np.power(1-np.square(x),np.abs(m/2))
    return(x,P0)


def test_Associated_Legendre_Polynomials():
    x = np.arange(-1,1.01,0.01)
    l = 5
    for m in [0,2,4]:
        Normalize_factor = np.sqrt(np.math.factorial(l - abs(m)) * \
            (2 * l + 1) / (4 * np.pi * np.math.factorial(l + abs(m))))

        xx,Pml = Associated_Legendre_Polynomials(x, l, m)
        if m%2 !=0 :
            Pml = Pml * (-1)

        plt.plot(xx,Pml * Normalize_factor)

    print(len(x))
    print(len(xx))

    plt.xlim([-1,1])
    plt.ylim([-1, 1])
    plt.grid(linestyle ='--')
    plt.show()
    return()

##############################################################
def derivative_1(x,y):
# Get the first derivative, to avoid the singular point at the end of list, I delete the last data of x and y.
# So the length of x and y will decrease for 1!
    dx = x[2] - x[1]
    dydx = []

    for i in range(len(x)-1):
        dy = y[i+1] - y[i]
        dydx.append(dy / dx)

    dydx= np.array(dydx)
    y = np.delete(y, -1)
    x = np.delete(x, -1)

    return(x,y,dydx)

def derivative_2(x,y):
    dx = x[2] - x[1]
    dydx = []
    dydx.append((y[2] - y[1]) / dx)

    for i in range(1, len(x)-1):
        dy = y[i+1] - y[i-1]
        dydx.append(dy / (2 * dx))

    dydx.append((y[-1] - y[-2]) / dx)
    dydx= np.array(dydx)

    return(dydx)

def Hermit(x,y,n):
# Hermit  厄米多项式
    for i in range(n):
        x,y,dydx = derivative(x,y)
        y = dydx

    Hn = np.power(-1,n) * np.exp(np.square(x)) * y

    return(x,Hn)

def test_1D_harmonic_oscillator(n):
#The wave function solution of one dimensional harmonic oscillator
    x = np.arange(-10,10,0.01)
    y = np.exp(-1 * np.square(x) )

    x,Hn = Hermit(x,y,n)
    solution = Hn * np.exp(-0.5 * np.square(x))

    plt.plot(x,solution)
    plt.title("1d Harmonic Oscillator n = %i" % n, fontsize = 16)
    plt.show()

    return()

#############################################################
def test_1D_inf_potential_well(n):

    x = np.arange(-0.5,0.5,0.01)
    if n%2 == 0 :  # even
        y = np.sqrt(2) * np.sin( n * np.pi * x)
    else:  # odd
        y = np.sqrt(2) * np.cos(n * np.pi * x)

    plt.plot(x,y)
    plt.title("1d infinite potential well n = %i" % n, fontsize = 16)
    plt.show()
    return()
######################################################
def Spherical_harmonics(m,l):
    # 球谐函数
    # for each l, m = -l, ..., 0, ..., l

    theta = np.arange(0, np.pi, 0.01)
    phi = np.arange(0, 2 * np.pi, 0.01)

    Normalize_factor = np.sqrt(np.math.factorial(l-abs(m)) * \
       (2 * l +1) / (4*np.pi* np.math.factorial(l+abs(m))) )

    cos_theta,Pml = Associated_Legendre_Polynomials(np.cos(theta), l, m)

    Ylm = Normalize_factor * Pml * np.exp(1j * m * phi)


    return(Ylm)



def plot_Ylm():

    l = 1
    m = 0


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')



#########  TEST  #############################################
if __name__ == "__main__" :

    # test_n_Level_Legendre_Polynomial() #  n阶勒让德多项式
    # test_1D_harmonic_oscillator(4)  #  一维简谐振子
    # test_1D_inf_potential_well(4)  #  一维无限深势阱
    # test_Associated_Legendre_Polynomials()  # 连带勒让德多项式
    # Spherical_harmonics(1, 2)  # 球谐函数
    # plot_Ylm()
########  Parameters ############################################
    mass_electron =  9.1094E-31  # [kg]
    Intensity = 5 * 1E13   # [W/cm^2]
    reduced_mass = 0.95 * mass_electron  #  should be calculated !!!
    Band_gap =  9.3 # [eV]
    refractive_index = 1.453  # dimensionless
    wavelength = 800
    frequency = my.Wavelength2Omega(wavelength) # [800 nm] -> [1/s]
#################################################################


# Quantitative identification of differentstrong-field ionization channels in the transitionregime


# Theories of photoelectron correlation in laser-driven multiple atomic ionization