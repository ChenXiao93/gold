import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx

########  图片设定
fs = 14

fig0 = plt.figure(figsize= (6,5))
grid = plt.GridSpec(4, 1, wspace=0, hspace=0)
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
ax01 = fig0.add_subplot(grid[0:2, 0:1])
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
ax02 = fig0.add_subplot(grid[2:4, 0:1])
#######################################
start,end,dw = 0.1, 20, 0.01   # eV
w_ev = np.arange(start, end + dw, dw)  # eV
omega_au = w_ev / cx.hartree2ev  # a.u.
eta = 0.3 # eV
delay = [600]
time = cx.ReadOne("./pump/gauge_field",1)
dAref = cx.ReadOne("./pump/gauge_field",5)
dt = 0.05

d = 250 # nm
length = 600 # a.u.
######  纯probe ########
t_probe, dAdt_probe = cx.ReadData("./probe/gauge_field", 1,5)
E_probe = dAdt_probe * (-1/cx.c_au)  * 51.44
t_probe = t_probe[0:int(length/0.05)] ###  截断
dAdt_probe = dAdt_probe[0:int(length/0.05)] ###  截断
damping_probe = cx.RatioDampingMethod(t=t_probe, t0=0, eta=eta, method="gaussian") ###  衰减
dAdt_probe_damp = dAdt_probe * damping_probe ###  衰减
eps_gs = cx.gauge_kick_delay(t=t_probe, A0=0.22, dAdt=dAdt_probe_damp, w_ev=w_ev, t0=0) # 复数

ax01.plot(w_ev, eps_gs.real, "g--", label="gs")
ax02.plot(w_ev, eps_gs.imag, "g--", label="gs")

##############################
for i in range(len(delay)):
    dAtotal = cx.ReadOne("./%s/gauge_field" % delay[i], 5)
    dAkick = dAtotal - dAref
    Ekick = dAkick * (-1/cx.c_au)  * 51.44

    cutoff = int((delay[i] + length) / dt)
    t = time[0:cutoff]
    dAcutoff = dAkick[0:cutoff]
    ####################
    damping = cx.RatioDampingMethod(t=t, t0=delay[i], eta=eta, method="gaussian")
    dAcutoff_damp = dAcutoff * damping

    ####  画吸收谱
    eps_probe = cx.gauge_kick_delay(t=t, A0=0.22, dAdt=dAcutoff_damp, w_ev=w_ev, t0=delay[i])
    ax01.plot(w_ev, eps_probe.real, "r-",label = "excited: code")
    ax02.plot(w_ev, eps_probe.imag, "r-",label = "excited: code")


######  octopus utility:  oct-dielectric-function
# PropagationSpectrumMaxEnergy = 20*eV
# PropagationSpectrumDampMode = gaussian
# PropagationSpectrumDampFactor = 0.3*eV
# TransientAbsorptionReference = "../ref/td.general"

w_u, eps1_u, eps2_u = cx.ReadData_xyz("./600/dielectric_function", 0,1,2)
w_u = w_u * 27.2114
ax01.plot(w_u, eps1_u, "b--",label = "excited: utility")
ax02.plot(w_u, eps2_u,  "b--",label = "excited: utility")

############################################
ax02.legend(loc="upper left", fontsize=fs)

ax01.set_xlim([start,end])
ax02.set_xlim([start,end])

ax01.set_ylabel("$\epsilon_1$", fontsize = fs+2 )
ax02.set_ylabel("$\epsilon_2$", fontsize = fs+2 )
ax02.set_xlabel("$\hbar \omega$ [eV]", fontsize = fs+2 )


ax01.set_xticks([0,5,10,15,20],["","","","",""],fontsize = fs)
ax02.set_xticks([0,5,10,15,20],[0,5,10,15,20],fontsize = fs)


plt.show()