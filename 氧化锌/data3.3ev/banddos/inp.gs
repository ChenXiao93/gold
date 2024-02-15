CalculationMode = gs
FromScratch = yes
PeriodicDimensions = 3
ParKPoints = auto
ParStates = no  # n_excited_el is not implemented for parallel in states
BoxShape = parallelepiped
ExperimentalFeatures=yes
PseudopotentialSet = pseudodojo_pbe   

s =0.20 
%Spacing                 # 单位 Bohr
 s | s | s
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

KPointsUseSymmetries = yes
%SymmetryBreakDir
 0 | 0 | 1
%

DFTULevel = dft_u_empirical
%Species
 "Zn" | species_pseudo | hubbard_l | 2 | hubbard_u | 30*eV
 "O"  | species_pseudo | hubbard_l | 1 | hubbard_u | 6*eV
%
