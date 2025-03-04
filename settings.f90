!-----------------------------------------------------------------------------------------------------------------------------------------!
!--------- USER INPUT FILE ---------------------------------------------------------------------------------------------------------------! 
!-----------------------------------------------------------------------------------------------------------------------------------------!
! - compile solver with: gfortran.exe -fopenmp -O3 permittivitySolver.f90 -o solverRun.exe
! - number of cpus/threads has to be adjusted with environment variable OMP_NUM_THREADS, e.g. on windows 11 --> $env:OMP_NUM_THREADS=8     

integer, parameter :: nx = 50                                                                                                  !Image width
integer, parameter :: ny = 50                                                                                                  !Image height
integer, parameter :: nz = 50                                                                                                  !Image depth
integer, parameter :: voxelRes = 1 !DO NOT CHANGE                                                                                 !VoxelResolution of input-raw data --> SHOULD STAY 1, VALUE DOES NOT REALLY HAVE INFLUENCE ON RESULT, MIGHT BE CHANGED TO TYPE REAL IN FUTURE VERSION
integer, parameter :: realPrecision = 2                                                                                          !floating-point precision: 1 (single) or 2 (double)
integer, parameter :: nIteration = 100000                                                                                         !max. number of iterations
real(kind=4*realPrecision), parameter :: targetResidual = 1e-30                                                                    !target residual    
integer(kind=1), parameter :: bitLength = 8                                                                                     !Bites per pixel, typically grey-values in raw images are stored as either 8-bit or 16-bit integers
integer(kind=1), parameter :: nPhases = 2 
integer, dimension(2) :: rawValues = [0,255]                                                                                      !bit-grayvalues of each phase (unsigned format)
real(kind=4*realPrecision), dimension(nPhases) :: correspConductivity = [1.0,1e-6]                                                       !conductivity of each phase (only double values, e.g. 1.0)
real(kind=4*realPrecision), dimension(nPhases) :: correspPermittivity = [79.0,4.5]                                                       !conductivity of each phase (only double values, e.g. 1.0)
real(kind=4*realPrecision), parameter :: frequency = 1.01e8                                                                     !frequency   
character(512) :: casePath = 'C:\Users\Mirko\Desktop\solver\solver_permittivity\samples\'       !END STRING WITH SLASH             !path to case folder (place raw-data in this folder, results will as be saved there as well)
character(512) :: sampleName = 'veldsteen_poreus_50cube.raw'                                                                   !input file (veldsteen_poreus_250cube, Berea_2d25um_binary_400cube, spheres_150cube)
character(1) :: evalDirection = 'Z'                                                                                             !direction of calculation (X or Y or Z)
logical :: writeTempertureField = .FALSE.                                                                                         !write calculated temperature-field (.TRUE. or .FALSE.)
logical :: writeFluxField = .FALSE.                                                                                               !write calculated flux-field (.TRUE. or .FALSE.)

!----------------------------------------------------------------------------------------------------------------------------------------!
!--------- Feel free to contact me. I will be happy to help: mirko.siegert@web.de -------------------------------!
!----------------------------------------------------------------------------------------------------------------------------------------!
