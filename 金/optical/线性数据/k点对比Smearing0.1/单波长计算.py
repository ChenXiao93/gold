# Smearing = 0.3*eV
# SmearingFunction = fermi_dirac
# RestartFixedOccupations = yes
# ExtraStates = 5
# a = 4.08*angstrom
# nk=10，20，40  对比不同k值
# prop = 500
import numpy as np
import matplotlib.pyplot as plt
import OctopusTools.Basic as cx
import OctopusTools.optical as op

exp_path = ["../exp/JC.csv"]  # 实验路径
# path = ["k10", "k20","k40"]   # 模拟路径


# 画出tddft结果
start, end = 1, 6  # eV

Wavelength = 300 #  Energy [eV]  to Wavelength [nm]
Omega = cx.Wavelength2Energy(Wavelength)  # eV
A0 = 1   # Ha
damp = 0.2# eV

time = cx.ReadOne("k20", 1)
epsilon = op.gauge_field_kick_mono(prop_time = time, kick_delay = 0, kick_amplitude=A0,
        dAprobe_dt = cx.ReadOne("k20", 5), damping_factor = damp, frequency = Omega)

print(epsilon)



