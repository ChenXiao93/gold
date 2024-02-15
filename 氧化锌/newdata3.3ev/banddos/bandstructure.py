import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx

def band(file,valence_band,conduction_band, efermi, xticks1, xticks2):

    band=np.loadtxt(file)
    k = band[:,0]
    CDMI = min(band[:,valence_band+4]) * 27.2114 # conduction band minimum point
    VBMX = max(band[:,valence_band+4-1]) * 27.2114 # valence band maximum point
    bandgap = CDMI - VBMX

    band_number = 6
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
    plt.xticks(xticks1, xticks2, fontsize=fs)
    plt.text(0.2,bandgap/3,"Direct Bandgap = %.3f eV" % bandgap, fontsize = fs)
    plt.title("ZnO Band Structure", fontsize = fs+2)
    plt.xlim(0, 1)
    plt.show()

    return()

def BandDos(file1,valence_band,conduction_band, efermi, xticks1, xticks2, file2):

    band=np.loadtxt(file1)
    k = band[:,0]
    CDMI = min(band[:,valence_band+4]) * 27.2114 # conduction band minimum point
    VBMX = max(band[:,valence_band+4-1]) * 27.2114 # valence band maximum point
    bandgap = CDMI - VBMX

    grid = plt.GridSpec(1, 3, wspace=0.3, hspace=0.5)
    # 第一个子图
    plt.subplot(grid[0:1, 0:2])

    band_number = 6
    for i in range(valence_band - band_number, valence_band + band_number):
        each_band = (band[:, i + 4]) * 27.2114 - VBMX
        plt.plot(k, each_band, "blue")
    # 标记出带隙边界的上下两条线
    xx = np.arange(0,1.1,0.1)
    yy1 = np.ones_like(xx) * (0)
    yy2 =  np.ones_like(xx) * (CDMI - VBMX)
    plt.plot(xx,yy1,"--", color = "red")
    plt.plot(xx, yy2, "--", color="red")

    fs = 14
    plt.ylabel("$E-E_f$  [eV]", fontsize=fs)
    plt.xlabel("Wave vector $ \\vec{k}$", fontsize=fs)
    plt.xticks(xticks1, xticks2, fontsize=fs)
    plt.text(0.2,bandgap/3,"Direct Bandgap = %.3f eV" % bandgap, fontsize = fs)
    plt.title("ZnO Band Structure", fontsize = fs+2)
    plt.xlim(0, 1)
    plt.ylim(-5, 15)
    # 第二个子图
    plt.subplot(grid[0:1, 2:3])

    f = np.loadtxt(file2)

    energy = (f[:,0] - 0.224918)*27.2114
    dos = f[:,1]

    plt.plot(dos, energy, "b")
    plt.ylim(-5, 15)
    plt.xlim(0,45)
    ax = plt.gca()
    # ax.set_yticks([])  # 关闭x轴坐标刻度
    plt.xlabel("$DOS (arb.u.)$", fontsize=fs)
    # plt.xticks([5, 10], [5,10], fontsize=fs)

    #画出两条分割能隙的边界线
    x_bound = np.arange(0,45.1, 0.1)
    y1_bound = np.ones_like(x_bound) * 0
    y2_bound = np.ones_like(x_bound) * ( CDMI - efermi * 27.2114 )
    plt.plot(x_bound, y1_bound, "--", color = "k")
    plt.plot(x_bound, y2_bound, "--", color="k")

    plt.show()

    return()


def testband():
    band(file = "./bandstructure",valence_band = 34-8,conduction_band = 8, efermi = -0.135515,\
        xticks1 = [0,0.16664393,0.26285584,0.45527968,0.54472032,0.71136425,0.80757617,1], \
        xticks2 = ["$\Gamma$", "M", "K","$\Gamma$", "A", "L", "H", "A"])

    BandDos(file1 = "./bandstructure",valence_band = 34-8,conduction_band = 8, efermi = -0.135515,\
        xticks1 = [0,0.16664393,0.26285584,0.45527968,0.54472032,0.71136425,0.80757617,1], \
        xticks2 = ["$\Gamma$", "M", "K","$\Gamma$", "A", "L", "H", "A"], \
        file2 = "./total-dos.dat")

    return()

if __name__ == "__main__" :
    testband()
