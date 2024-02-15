import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx


# for amp in [5,10,15,20,25,30,55,90]:
for amp in [5,10,15,20,25,30,55,90]:
    # f = open("./A%s" % amp, "a")
    # for t in np.arange(600.005,800.01,0.015):
    #     f.write("%.4f                %.4f\n" % (t,0))
    # f.close()
    t,A = cx.ReadData("./A%s" % amp, 0,1)
    plt.plot(t,A)
    print(max(A))


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
# %