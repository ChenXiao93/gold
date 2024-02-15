import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.band as band


x1 = [0, 0.16489717, 0.31032280, 0.60117407,0.68967720,0.85457437, 1]
x2 = [r'$\Gamma$', 'M', 'K', r'$\Gamma$', 'A', 'L', 'H']
band.band(file = "bandstructure",valence_band = 26, conduction_band = 8,\
    efermi = -0.119379, xticks1 = x1, xticks2 = x2)



