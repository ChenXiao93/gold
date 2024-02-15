import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx

##########  Envelope ###################################################
def env_sin2(mid, tail, dt):
    # flat line + sin2 attenuation
    time = np.arange(0, mid+tail, dt)
    envy =   (np.heaviside(time - mid, 1)   * \
             np.square( np.sin( (time - mid)  *  np.pi / tail /2    + np.pi/2)) + \
             np.heaviside(mid - time - 0.01, 1) ) *  np.heaviside( mid + tail - time, 1)
    return(time, envy)

def env_sin2_full(head, mid, tail, dt):
    # head sin2 increase + flat line + sin2 attenuation
    time = np.arange(0, head+mid+tail, dt)
    envy =   (np.heaviside(time - mid - head, 1)   * \
             np.square( np.sin( (time - mid - head)  *  np.pi / tail /2    + np.pi/2))  + \
             np.heaviside(head + mid - time - 0.01, 1)    \
             +  np.heaviside(head  - time, 1)  * \
              ((np.square(np.sin( (time - head)  *  np.pi / head /2 + np.pi/2)) -1) ) )

    return(time, envy)

def env_exp(mid, tail, dt, gamma):
    #  flat line +  exp(-gamma*t) attenuation
    time = np.arange(0, mid+tail, dt)
    envy = np.exp(-gamma * (time - mid) ) * np.heaviside(time - mid-0.01,1)   \
           + np.heaviside( mid - time ,1)
    return(time, envy)

def env_exp_full(head, mid, tail, dt, gamma_head, gamma_tail):
    # head exp(gamma*t) increase + flat line +  exp(-gamma*t) attenuation
    time = np.arange(0, head+mid+tail, dt)
    envy = np.exp(-gamma_tail * (time - mid - head) ) * np.heaviside(time - mid - head -0.01,1)   \
           + np.heaviside( head + mid - time ,1)  + \
           np.heaviside(head - time, 1) * ( np.exp(gamma_head * (time - head))  -  1)
    return(time, envy)


##############  Generate monochromatic pulse #################################################################
def Env_Gaussian(time, ratio, eta_ev ):
    ################# Input parameters###########################
    envelope = np.ones_like(time)
    t0 = time[-1] * (1 - ratio)
    eta_au = eta_ev / 27.2114
    T = 1 / eta_au
    for i in range(len(time)):
        if time[i] > t0 :
            envelope[i] = np.exp(-np.square((time[i] - t0) / (T/2) ))
        else:
            pass
    return ( envelope)   # t[hbar/Ha] A [Hartree/e]


def Env_exponential(time, ratio, eta_ev):
    ################# Input parameters###########################
    envelope = np.ones_like(time)
    t0 = time[-1] * (1 - ratio)
    eta_au = eta_ev / 27.2114
    for i in range(len(time)):
        if time[i] > t0 :
            envelope[i] = np.exp(-(eta_au * (time[i] - t0)))
        else:
            pass
    return ( envelope)   # t[hbar/Ha] A [Hartree/e]

def Env_exponential_ratio_std(time, ratio, eta_ev):
    ################# Input parameters###########################
    envelope = np.ones_like(time)
    t0 = time[-1] * (1 - ratio)
    eta_au = eta_ev / 27.2114
    for i in range(len(time)):
        if time[i] > t0 :
            envelope[i] = np.exp(-(eta_au * (time[i] - t0)))
        else:
            pass

    return ( envelope)   # t[hbar/Ha] A [Hartree/e]


def GaussianPulse(PropTime, PulseDuration, Wavelength, Amplitude, TimeStep, Phase):
    ################# Input parameters###########################
    # PropTime is the total propagation time [hbar/Ha]
    # PulseDuration is the time range when the intensity decrease to 1/e [habr/Ha]
    PeakPosition = PropTime / 2  # Pulse peak position tau [fs]
    omega_si =  cx.Wavelength2Omega(Wavelength)  # Wavelength [nm]  to  Omega[1/s]
    omega_au = omega_si  * cx.timeau2s  #
    energy = cx.Wavelength2Energy(Wavelength) # wavelength [nm] to  Energy [eV]

    time = np.arange(0, PropTime, TimeStep)
    envelope = Amplitude * np.exp(-np.square((time - PeakPosition) / (PulseDuration/2) ))
    gauge_field = envelope * np.cos(omega_au * (time-PeakPosition) + Phase)
    return (time, gauge_field, envelope)   # t[hbar/Ha] A [Hartree/e]

#################################################################################################
def test():
    # time, envy = env_exp(mid = 1000, tail = 500, dt = 0.1, gamma = 0.02)
    # time, envy = env_sin2(mid=10, tail=20, dt=0.1)
    # time, envy = env_sin2_full(head = 500, mid=9000, tail=500, dt=0.2)
    # time, envy = env_exp_full(head = 0, mid = 100, tail = 100, dt = 0.1, gamma_head = 0.02, gamma_tail = 0.1)
    # plt.plot(time, envy)
    x = np.arange(0,1000,0.05)
    # y = Env_Gaussian(time = x, ratio = 0.5, eta_ev = 0.05)
    y = Env_exponential(time = x, ratio = 0.2, eta_ev = 0.2)
    # time, gauge_field, envy = GaussianPulse(PropTime = 500.05, PulseDuration= 200,  Wavelength = 800, Amplitude = 1, TimeStep =0.1, Phase = 0 )
    # cx.save("A1",time, gauge_field, "### t [hbar/Ha]    A [Ha/e]")

    # plt.plot(time * 0.02419, gauge_field)
    plt.plot(x,y)
    plt.xlabel("t [fs]")
    plt.ylabel("A [Ha]")

    plt.show()

    return()

if __name__ == "__main__" :
    test()
