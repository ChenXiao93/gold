import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op

t,dAdt = cx.ReadData("./gauge_field",1,5)  # 纯probe场

plt.plot(t,dAdt)

plt.show()