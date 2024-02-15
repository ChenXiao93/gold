import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx

########  图片设定
fs = 12


fig0 = plt.figure(figsize= (6,4))  ## E_ext 和 E_ind
grid = plt.GridSpec(2, 1, wspace=0, hspace=0)
plt.rcParams['xtick.direction'] = 'in'
ax01 = fig0.add_subplot(grid[0:1, 0:1])
plt.rcParams['xtick.direction'] = 'in'
ax02 = fig0.add_subplot(grid[1:2, 0:1])


################################
time, dAinddt_ref = cx.ReadData("./td.ref/gauge_field", 1, 7)  #  pump
dAinddt_general = cx.ReadOne("./td.general/gauge_field", 7)  # pump + probe
E_general =  (-1/cx.c_au) * dAinddt_general *  51.44  # 51.44= 27.2114/0.529  # V / Angstrom
E_ref = (-1/cx.c_au) * dAinddt_ref *  51.44
E_diff =  E_general - E_ref

ax01.plot(time * 0.02419, E_general, label = "Pump+Probe")
ax01.set_ylabel("$E_{ind} \ [V/\mathring{A}]$", fontsize=fs)
ax01.legend(loc="lower left", fontsize=fs)
ax01.set_xlim([0,37])
ax01.set_xticks([0,5,10,15,20,25,30,35],[],fontsize = fs)

ax02.plot(time * 0.02419, E_diff, label = "(Pump+probe) - Pump")
ax02.set_ylabel("$E_{ind} \ [V/\mathring{A}]$", fontsize=fs)
ax02.set_xlabel("Time [fs]", fontsize=fs)
ax02.legend(loc="upper left", fontsize=fs)
ax02.set_xlim([0,37])

plt.show()