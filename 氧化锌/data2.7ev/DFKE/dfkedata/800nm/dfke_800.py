import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx

########  图片设定
fs = 14



fig1 = plt.figure(figsize=(8, 4))  ##  dA/dt 和 damping 修正
ax1 = fig1.add_subplot(111)

start =0.5 # eV
end =4 # eV
w_ev = np.arange(start,end + 0.01,0.1) # eV
wavelen_nm = cx.Energy2Wavelength(w_ev) # nm

for amp in [5]:
    time, dAinddt_ref = cx.ReadData("./A%s/td.ref/gauge_field" % amp, 1, 7)  # kick probe + pump
    dAinddt_general = cx.ReadOne("./A%s/td.general/gauge_field" % amp, 7)    # pump
    dAinddt_kick = dAinddt_general - dAinddt_ref



    # ax1.plot(time, dAinddt_ref)
    # ax1.plot(time, dAinddt_general)
    ax1.plot(time,dAinddt_kick)

plt.show()