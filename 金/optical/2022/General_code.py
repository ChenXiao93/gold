import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx

def FT(t, dAdt, omega):
# omega can be one value, or even list.   [Ha/hbar]
# t [hbar/Ha]; dA [Hartree]
# return the real and imagnary part of FT result.
    dt = t[2] - t[1]  # dt [hbar / Ha]
    a = 0  # initialize the real of FT
    b = 0  # initialize the imag of FT
    for i in range(len(t)):
        a += dAdt[i] * dt  * np.cos(omega * t[i])
        b += dAdt[i] * dt  * np.sin(omega * t[i])
    return(a,b)

def nfit(x, y, n):
    x = np.array(x)
    fitting = np.polyfit(x, y, n)
    parameters = np.poly1d(fitting)
    fitted_y = 0
    for i in range(n + 1):
        fitted_y += np.power(x, i) * parameters[i]

    return (fitted_y, parameters)

def env(x,y, pulse_duration):
    #  sin2 attenuation
    envy =   max(y) * (np.heaviside(x - pulse_duration, 1)   * \
             np.square( np.sin( (x - pulse_duration)  *  np.pi / tail /2    + np.pi/2)) + \
             np.heaviside(pulse_duration - x - 0.01, 1) ) *  np.heaviside( pulse_duration + tail - x, 1)
    newy = envy * y / max(envy)
    return(newy)



rotenberg = np.loadtxt("./exp/Rotenberg630real.txt", delimiter=',')
rotenbergimag = np.loadtxt("./exp/Rotenberg630imag.txt", delimiter=',')
chi3 = np.loadtxt("./exp/test")
time  = (np.arange(0,3300,300) + 2000) * 0.02419 * 1E-15
chi3r = chi3[:,0]
chi3i = chi3[:,1] 
#############################
plt.subplot(1,2,1)
plt.loglog(time, np.abs(chi3r), "--o")
plt.loglog(rotenberg[:,0],  rotenberg[:,1],"--o", label = "exp", markersize = 5)
plt.subplot(1,2,2)
plt.loglog(time, np.abs(chi3i), "--o")
plt.loglog(rotenbergimag[:,0],  rotenbergimag[:,1],"--o", label = "exp", markersize = 5)

plt.show()




