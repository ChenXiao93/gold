import numpy as np
band=np.loadtxt("./static/bandstructure")
k_octopus=band[:,0]
valence_band = 26
conduction_band = 8
CDMI = min(band[:,valence_band+4]) * 27.2114 # conduction band minimum point
VBMX = max(band[:,valence_band+4-1]) * 27.2114 # valence band maximum point
band_gap = CDMI - VBMX

print(band_gap)

