import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx

def pump_pulse_3500(a):

    dt= 0.015   # time step magnitude  [hbar/Hartree]
    T = 1300
    T1 = 200
    t = np.arange(0,T,dt)   #  [hbar/Hartree] normal propagation time
    t1 = np.arange(T, T+T1, dt)  # add a flat tail
    tt = np.append(t,t1)  # total propagation time

    A = a * np.sin( omega_au * t )  * np.square(np.sin(np.pi * t/T))
    A1 = np.zeros_like(t1)
    AA = np.append(A,A1)


    return(tt,AA)

##############################
wavelength = 3500  # nm
energy = cx.Wavelength2Energy(wavelength) # eV
omega_au = energy / cx.hartree2ev  # a.u.
Intenisty = [1e11,3e11,6e11,1e12,3e12,6e12,1e13]
eps1_3500 = 3.3
for i in Intenisty:
    Iau = np.sqrt(i) / np.sqrt(3.509470 * 1E16)
    a = Iau/omega_au*cx.c_au
    a_modified = a * eps1_3500
    print(a)
    t,A = pump_pulse_3500(a_modified)

    plt.plot(t,A)
    plt.xlabel("Propagation time [$\hbar/Ha$]", fontsize=14)
    plt.ylabel("Gauge field A [Ha]", fontsize=14)
    cx.save("./A%i" % a_modified, t, A, "t [hbar/Ha]  A [Ha]")
plt.show()

# Tpulse = 1300
# Iau  = (sqrt(1*10^11)/sqrt(3.509470*10^16)) #Conservion in a.u. from W/cm^2
# omega = 0.35424*eV
# phi = 0.5*pi
#
# %TDExternalFields
#   vector_potential | 1.0 | 0.0 | 0.0 | omega | "envelope_sin2" | "phase"
# %
#
# %TDFunctions
#   "envelope_sin2" | tdf_from_expr | "-Iau/omega*c*sin(pi*t/Tpulse)^2*(1-step(t-Tpulse))"
#    "phase" | tdf_from_expr  | "phi"