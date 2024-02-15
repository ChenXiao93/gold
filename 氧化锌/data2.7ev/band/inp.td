CalculationMode = unocc
FromScratch = yes
PeriodicDimensions = 3
ParKPoints = auto
BoxShape = parallelepiped
ExperimentalFeatures=yes
PseudopotentialSet = pseudodojo_pbe     # 使用pbe赝势

%Spacing                 # 单位 Bohr
 0.20 | 0.20 | 0.20 
%

%LatticeVectors          
  1.0 | 0.0   | 0.0
 -0.5 | sqrt(3)/2 | 0.0
  0.0 | 0.0   | 1.0
%

a = 3.289*angstrom
b = 5.307*angstrom
%LatticeParameters
 a | a | b
%

%ReducedCoordinates
 "O"  |- 1/6 | 1/6 | -0.120238
 "O"  | 1/6 | -1/6 | 0.379762
 "Zn" | -1/6 | 1/6 | -0.500000
 "Zn" | 1/6 | -1/6 | 0.00000
%


%KPointsPath        #k点扫描路径,见“注1”
  14 |  10 | 20 | 14 | 14 | 10    #标明各个点之间的路径需要取样几个k点
 0.0 | 0.0 | 0.0  # Gamma point  #所以有7个点的坐标，只有6段路径需要标明取样点数
 0.5 | 0.0 | 0.0  # M  point
 2/3 | 1/3 | 0.0  # K point
 0.0 | 0.0 | 0.0  # Gamma point
 0.0 | 0.0 | 0.5  # A point
 0.5 | 0.0 | 0.5  # L point
 2/3 | 1/3 | 0.5  # H point
%
KPointsUseSymmetries = no

DFTULevel = dft_u_empirical              
%Species
 "Zn" | species_pseudo | hubbard_l | 2 | hubbard_u | 12.8*eV
 "O"  | species_pseudo | hubbard_l | 1 | hubbard_u | 5.29*eV
%

Output = dos
ExtraStates = 8              
ExtraStatesToConverge = 6  
