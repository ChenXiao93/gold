import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx

########  图片设定
fs = 14

fig0 = plt.figure(figsize=(6, 4))  ##  eps 比较
ax01 = fig0.add_subplot(211)
ax02 = fig0.add_subplot(212)

fig1 = plt.figure(figsize=(8, 4))  ##  dA/dt 和 damping 修正
ax1 = fig1.add_subplot(111)

fig2 = plt.figure(figsize=(6, 4))  ##  dA/dt 和 damping 修正
ax2 = fig2.add_subplot(111)

start =0.5 # eV
end =4 # eV
w_ev = np.arange(start,end + 0.01,0.01) # eV
wavelen_nm = cx.Energy2Wavelength(w_ev) # nm
#################
for amp in [143]:
    time, dAinddt_ref = cx.ReadData("./A%s/td.ref/gauge_field" % amp, 1, 7)  # kick probe + pump
    dAinddt_general = cx.ReadOne("./A%s/td.general/gauge_field" % amp, 7)    # pump
    dAinddt_kick = dAinddt_general - dAinddt_ref

    ax2.plot(time, dAinddt_general, label = "pump+probe")
    ax2.legend(loc="lower left", fontsize=fs)
    ax2.set_xlabel("t  [a.u.]", fontsize=fs)
    ax2.set_ylabel(" dA/dt [a.u.]", fontsize=fs)

    etaev = 0.2 # eV
    ax1.plot(time * 0.02419, dAinddt_kick, label = "dA/dt")
    damping = cx.RatioDampingMethod(t = time, t0 = 1200, eta = etaev, method = "exponential")
    ax1.plot(time * 0.02419,  max(dAinddt_kick) * damping, label = "$exp(-\eta t)$, $\eta$=%.1feV" % etaev)
    dAinddt_kick_damp = dAinddt_kick * damping
    ax1.plot(time * 0.02419, dAinddt_kick_damp, label = "dA/dt * $exp(-\eta t)$")

    eps_probe = cx.gauge_kick_delay(t=time, A0=0.1, dAdt=dAinddt_kick_damp, w_ev=w_ev, t0 = 1200)
    
    ax01.plot(w_ev, eps_probe.real, "-")
    ax02.plot(w_ev, eps_probe.imag, "-", label="probe")

ax1.legend(loc = "lower left", fontsize =fs)
ax1.set_xlabel("t  [fs]", fontsize = fs)
ax1.set_ylabel(" dA/dt [a.u.]", fontsize = fs )
ax1.set_xlim([1000 * 0.02419,1500 * 0.02419])


plt.show()