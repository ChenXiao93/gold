CalculationMode = gs
FromScratch = yes
PeriodicDimensions = 3
ParStates = no
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
 "O"  | 1/3 | 2/3 | 0.379214
 "O"  | 2/3 | 1/3 | 0.879214
 "Zn" | 1/3 | 2/3 | 0.000000
 "Zn" | 2/3 | 1/3 | 0.500000
%

%KPointsGrid
  12  |  12  |  7
%
%SymmetryBreakDir
 0 | 0 | 1
%

KPointsUseSymmetries = yes

DFTULevel = dft_u_empirical              
%Species
 "Zn" | species_pseudo | hubbard_l | 2 | hubbard_u | 12.8*eV
 "O"  | species_pseudo | hubbard_l | 1 | hubbard_u | 5.29*eV
%

