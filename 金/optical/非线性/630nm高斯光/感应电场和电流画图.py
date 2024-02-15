import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op
import OctopusTools.Laser as la
import OctopusTools.chi3 as ch

def compare_current(pulse_duration, Amplitude, dt, index):
    time_ind, A_ind = cx.ReadData("./induced/t%s/A%s/gauge_field" % (pulse_duration, Amplitude), 1, 2)
    time_ext, A_ext, envy = la.GaussianPulse(PropTime=pulse_duration + dt,
                                    PulseDuration=int(pulse_duration * 2 / 5),
                                    Wavelength=630,
                                    Amplitude=Amplitude, TimeStep=dt, Phase=0)
    A_tot = A_ind + A_ext
    J_oct = cx.ReadOne("./induced/t%s/A%s/total_current" % (pulse_duration, Amplitude), 2)
    J_cal = cx.current_cal(t_ind = time_ind, A_ind= A_ind, V = 114.5822 )
    # 画 规范场A
    ax1[index].plot(time_ind * 0.02419, A_ind, label="A-ind")
    # ax1[index].plot(time_ext * 0.02419, A_ext, label="A-ext")
    # ax1[index].plot(time_ext * 0.02419, A_tot, label="A-tot")
    # 画电流
    # ax2[index].plot(time_ind * 0.02419, J_oct * (-1), label="-J-Oct")
    # ax2[index].plot(time_ind * 0.02419, J_cal, label="J-Cal")
    return()
#########################################
fig1 = plt.figure()
ax11 = fig1.add_subplot(221)
ax12 = fig1.add_subplot(222)
ax13 = fig1.add_subplot(223)
ax14 = fig1.add_subplot(224)
ax1 = [ax11,ax12,ax13,ax14]

fig2 = plt.figure()
ax21 = fig2.add_subplot(221)
ax22 = fig2.add_subplot(222)
ax23 = fig2.add_subplot(223)
ax24 = fig2.add_subplot(224)
ax2 = [ax21,ax22,ax23,ax24]
####################
wavelength = 630
omega = cx.Wavelength2Energy(wavelength) # eV
omega_au = omega/27.2114 # a.u.
amp = [1,4,7,10] # Ha
T= [3000] # a.u.

for t in range(len(T)):
    for i in range(len(amp)):
        pulse_duration = T[t]
        Amplitude = amp[i]
        dt = 0.05
        index= i
        time_ind, A_ind = cx.ReadData("./induced/t%s/A%s/gauge_field" % (pulse_duration, Amplitude), 1, 2)
        time_ext, A_ext, envy = la.GaussianPulse(PropTime=pulse_duration + dt,
                                                PulseDuration=int(pulse_duration * 2 / 5),
                                                Wavelength=630,
                                                Amplitude=Amplitude, TimeStep=dt, Phase=0)
        A_tot = A_ind + A_ext
        J_oct = cx.ReadOne("./induced/t%s/A%s/total_current" % (pulse_duration, Amplitude), 2)
        J_cal = cx.current_cal(t_ind=time_ind, A_ind=A_ind, V=114.5822)
        # 画 规范场Aind
        # ax1[index].plot(time_ind * 0.02419, A_ind, label="Aind")
        ax1[index].plot(time_ext * 0.02419, A_ext, label="A-ext")
        # ax1[index].plot(time_ext * 0.02419, A_tot, label="A-tot")
        # 画电流 Jind
        ax2[index].plot(time_ind * 0.02419, J_oct, label="-Jind")
        # ax2[index].plot(time_ind * 0.02419, J_cal, label="J-Cal")


    ax14.legend()
    ax24.legend()
plt.show()
