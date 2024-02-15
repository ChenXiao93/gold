import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
#######################################
t1,n1 = cx.ReadData("./n_ex/n_A25",1,2)
t2,n2 = cx.ReadData("./n_ex/n_A35_no",1,2)
t3,Aext = cx.ReadData("./laser/A25",0,1)
t4,Aind = cx.ReadData("./gauge_field/gaugefield_A25",1,4)
t5,Ano = cx.ReadData("./laser/A35_no_modify",0,1)
time,energy = cx.ReadData("./energy/e_A25", 1,2)

t3 = t3[0:len(t4)]
Aext = Aext[0:len(Aind)]

Atot = Aext + Aind
Eext = cx.A2E(t3,Aext)
Etot = cx.A2E(t4,Atot)

##########################################
fs = 14
fig = plt.figure(figsize= (9,8))
grid = plt.GridSpec(3, 1, wspace=0, hspace=0)
plt.rcParams['xtick.direction'] = 'in'
ax1 = fig.add_subplot(grid[0:1, 0:1])
plt.rcParams['xtick.direction'] = 'in'
ax2 = fig.add_subplot(grid[1:2, 0:1])
ax3 = fig.add_subplot(grid[2:3, 0:1])

ax1.plot(t3*0.02419, Eext, "--b", label = "ext")
ax1.plot(t4*0.02419, Etot, "-r", label = "tot")


ax1.legend(loc = "upper right", fontsize = fs)
ax1.set_xlim([0,14])
ax1.set_xticks([0,2,4,6,8,10,12,14],[],fontsize = fs)
ax1.set_ylabel("E [$V/\mathring{A}$]", fontsize = fs+2)
ax1.set_yticks([-1,0,1],[-1,0,1],fontsize = fs)
ax1.legend(loc = "lower right", fontsize = fs-1)

ax2.plot(t1 * 0.02419,n1, "--b", label = "$\epsilon$ modification")
ax2.plot(t2 * 0.02419,n2,"-r", label = "no modification")
ax2.set_xlim([0,14])
ax2.set_ylim([0,1.5])
ax2.set_xticks([0,2,4,6,8,10,12,14],[],fontsize = fs)
ax2.set_yticks([0,0.5,1],[0,0.5,1],fontsize = fs)
ax2.set_ylabel("$n_{ex}$ [1/cell]", fontsize = fs+2)

ax2.legend(loc = "upper right", fontsize = fs)


ax3.plot(time * 0.02419, energy * 27.2114 + 12276, "--b", label = "$E_{tot} + 12276$")
ax3.set_xlim([0,14])
ax3.set_xticks([0,2,4,6,8,10,12,14],[0,2,4,6,8,10,12,14],fontsize = fs)
ax3.legend(loc = "upper right", fontsize = fs)
ax3.set_ylabel("$E_{tot}$ [eV]", fontsize = fs)
ax3.set_xlabel("Time [fs]", fontsize = fs+2)
ax3.set_yticks([-0.6, -0.4, -0.2],[-0.6, -0.4, -0.2],fontsize = fs)
fig.savefig("fig1.eps", dpi = 900, bbox_inches='tight')

plt.show()