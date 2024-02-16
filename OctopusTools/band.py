import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx

#！！！添加测试


def band(file,valence_band,conduction_band, efermi, xticks1, xticks2):

    band=np.loadtxt(file)
    k = band[:,0]
    CDMI = min(band[:,valence_band+4]) * 27.2114 # conduction band minimum point
    VBMX = max(band[:,valence_band+4-1]) * 27.2114 # valence band maximum point
    bandgap = CDMI - VBMX

    band_number = 4
    for i in range(valence_band-band_number,valence_band+band_number):
        each_band = (band[:, i + 4]) * 27.2114 - VBMX
        plt.plot(k, each_band, "blue")

    xx = np.arange(0,1.1,0.1)
    yy1 = np.ones_like(xx) * (0)
    yy2 =  np.ones_like(xx) * (CDMI - VBMX)
    plt.plot(xx,yy1,"--", color = "red")
    plt.plot(xx, yy2, "--", color="red")

    fs = 16
    plt.ylabel("$E-E_f$  [eV]", fontsize=fs)
    plt.xlabel("Wave vector $ \\vec{k}$", fontsize=fs)
    plt.xticks([0, 0.16489717, 0.31032280, 0.60117407,0.68967720,0.85457437, 1], [r'$\Gamma$','M','K', r'$\Gamma$','A','L','H'], fontsize=fs)

    plt.text(0.2,bandgap/3,"Direct Bandgap = %.2f eV" % bandgap, fontsize = fs)
    plt.title("ZnO Band Structure", fontsize = fs+2)
    plt.xlim(0, 1)
    plt.show()

    return()


def testband():
    band()
    return()

if __name__ == "__main__" :
    testband()
