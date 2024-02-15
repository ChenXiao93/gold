# 对于不同的T，E（w）有正有负， 而且是对于外场而言的，不应该啊，检查一下，应该都是正的！！
import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op
import OctopusTools.Laser as la
import OctopusTools.chi3 as ch
from scipy import fftpack

def compare_current(pulse_duration, Amplitude, dt):
    time_ext, A_ext, envy = la.GaussianPulse(PropTime=pulse_duration + dt,
                                    PulseDuration=int(pulse_duration * 2 / 5),
                                    Wavelength=630,
                                    Amplitude=Amplitude, TimeStep=dt, Phase=0)
    return(time_ext, A_ext)

def positivecomplex(complex1):
    return np.abs(complex1.real) + 1j * np.abs(complex1.imag)



def At2Ew(omega,t,A,J):
    # omega [eV]
    dt = t[2] - t[1]
    w0 = omega * cx.timeau2s / cx.hbar # a.u.
    Aw, Jw = 0,0 # initialize for Fourier Transform

    for i in range(len(A)):
        Aw += A[i] * dt * np.exp(w0 * t[i] * 1j)
        Jw += J[i] * dt * np.exp(w0 * t[i] * 1j)
    Ew =  1j * w0 * Aw / 137

    # Ew = positivecomplex(Ew)
    # Jw = positivecomplex(Jw)
    return(Ew, Jw)

amp = [10] # Ha
T= [1000] # a.u.
dt = 0.05  # a.u.
omega = np.arange(1.94,2.01,0.002)  # eV
omega_au = omega / 27.2114
omega630 = cx.Wavelength2Energy(630) # eV
omega630_au = omega630 / 27.2114
print("eps1 = -11.980  eps2 = 1.1524 (exp)")
n0exp, k0exp = 0.16629, 3.4652
print("n = %.4f, k = %.4f (exp)" % (n0exp, k0exp))

for t in range(len(T)):
    for i in range(len(amp)):
        time_ind, A_ind = cx.ReadData("./induced/t%s/A%s/gauge_field" % (T[t], amp[i]), 1, 2)
        time_ext, A_ext = compare_current(pulse_duration = T[t], Amplitude = amp[i], dt = dt)
        J_cal = cx.current_cal(t_ind=time_ind, A_ind=A_ind, V=114.5822)
        # J_cal= -J_cal
        ## 单色630nm的结果
        Ew630, Jw630 = At2Ew(omega630, time_ext, A_ext,J_cal)

        sigma630 = Jw630/Ew630
        sigma630 = positivecomplex(sigma630)
        epsilon630 = 1 + 1j * sigma630 / omega630_au
        print("630nm = %.4feV: Ereal = %.2f, Eimag = %.2f, |Ew| = %.4f" % (omega630, (Ew630.real), (Ew630.imag), np.abs(Ew630)))
        print("630nm = %.4feV: Jreal = %.2f, Jimag = %.2f, |Jw| = %.4f" % (omega630, (Jw630.real), (Jw630.imag), np.abs(Jw630)))
        print("sigma1 = %.4f, sigma2 = %.4f" % (sigma630.real, sigma630.imag))
        print("epsilon1 = %.4f, epsilon2 = %.4f" % (epsilon630.real, epsilon630.imag))
        # 多色630nm附近的结果
        Ew, Jw = [], []
        for w0 in omega:
            E1, J1 = At2Ew(w0, time_ext, A_ext, J_cal)
            Ew.append(E1)
            Jw.append(J1)
        Ew = np.array(Ew)
        Jw = np.array(Jw)
        sigmaw = Jw/Ew
        sigmaw = positivecomplex(sigmaw)
        epsilonw = 1 + 1j * sigmaw / omega_au

        ## plotting
        plt.subplot(4, 1, 1)
        plt.plot(omega,Ew.real, label = "$E(w)_{real}$")
        plt.plot(omega,Ew.imag, label = "$E(w)_{imag}$")
        plt.plot(omega, np.abs(Ew), label="$|E(w)|$")

        Ewlist = list(np.abs(Ew))
        Ewlistmax = max(Ewlist)
        Ewlistmaxindex = Ewlist.index(Ewlistmax)
        plt.scatter(omega[Ewlistmaxindex], np.abs(Ew)[Ewlistmaxindex])
        print("EwMax: omega = %.4f, |Ew| = %.4f" % (omega[Ewlistmaxindex], np.abs(Ew)[Ewlistmaxindex]))
        Eh = np.arange(min(np.abs(Ew)), max(np.abs(Ew)), 0.05)
        Wh = np.ones_like(Eh) * omega[Ewlistmaxindex]
        plt.plot(Wh, Eh)
        plt.legend(loc = "upper right")

        plt.subplot(4,1,2)
        plt.plot(omega,Jw.real, label = "$J(w)_{real}$")
        plt.plot(omega,Jw.imag, label = "$J(w)_{imag}$")
        plt.plot(omega, np.abs(Jw), label="$|J(w)|$")
        plt.legend(loc = "upper right")

        Jwlist = list(np.abs(Jw))
        Jwlistmax = max(Jwlist)
        Jwlistmaxindex = Jwlist.index(Jwlistmax)
        plt.scatter(omega[Jwlistmaxindex], np.abs(Jw)[Jwlistmaxindex])
        print("JwMax: omega = %.4f, |Jw| = %.4f" % (omega[Jwlistmaxindex], np.abs(Jw)[Jwlistmaxindex]))
        Jh = np.arange(min(np.abs(Jw)), max(np.abs(Jw)), 0.05)
        Wh = np.ones_like(Jh) * omega[Jwlistmaxindex]
        plt.plot(Wh, Jh)

        # if omega[Ewlistmaxindex] < omega[Jwlistmaxindex]:
        #     print("w(Emax) < w(Jmax) : blue shift")
        # else:
        #     print("w(Emax) > w(Jmax) : red shift")
        #
        # ## 验证一个猜想！sigma = Jwmax/Ewmax
        # sigma630mod = Jw[Jwlistmaxindex]/Ew[Ewlistmaxindex]
        # sigma630mod = positivecomplex(sigma630mod)
        # epsilon630mod = 1 + 1j * sigma630mod / (omega[Ewlistmaxindex] * cx.timeau2s / cx.hbar)
        # print("Jwmax/Ewmax:  sigma1 = %.4f, sigma2 = %.4f" % (sigma630mod.real, sigma630mod.imag))
        # print("Jwmax/Ewmax:  epsilon1 = %.4f, epsilon2 = %.4f" % (epsilon630mod.real, epsilon630mod.imag))
        #
        # plt.subplot(4, 1, 3)
        # plt.plot(omega,sigmaw.real, label = "$\sigma(w)_{real}$")
        # plt.plot(omega,sigmaw.imag, label = "$\sigma(w)_{imag}$")
        # plt.legend(loc = "upper right")
        # plt.xlim([1.95,2.0])
        # plt.ylim([0,1])
        #
        #
        # plt.subplot(4, 1, 4)
        # plt.plot(omega,epsilonw.real, label = "$\epsilon(w)_{real}$")
        # plt.plot(omega,epsilonw.imag, label = "$\epsilon(w)_{imag}$")
        # plt.legend(loc = "upper right")
        # plt.xlabel("$\hbar\omega [eV]$")
        # plt.xlim([1.95,2.0])
        # plt.ylim([-10,12])

plt.show()