CalculationMode = gs
FromScratch = yes
PeriodicDimensions = 3
ParKPoints = auto
BoxShape = parallelepiped
ExperimentalFeatures=yes
XCFunctional = mgga_x_tb09 + lda_c_pw
FilterPotentials = filter_none

Spacing= 0.4   # Bohr

%LatticeVectors
  0.0 | 0.5 | 0.5
  0.5 | 0.0 | 0.5
  0.5 | 0.5 | 0.0
%

a = 6.7556   # Bohr
%LatticeParameters
 a | a | a
%

%ReducedCoordinates
 "C" | 0.0 | 0.0 | 0.0
 "C" | 1/4 | 1/4 | 1/4
%

nk = 16
%KPointsGrid
  nk |  nk |  nk
 0.5 | 0.5 | 0.5
 0.5 | 0.0 | 0.0
 0.0 | 0.5 | 0.0
 0.0 | 0.0 | 0.5
%

KPointsUseSymmetries = yes
%SymmetryBreakDir
  1 | 0 | 0
%
