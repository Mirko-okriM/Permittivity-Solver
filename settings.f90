!----------------------------------------------------------------------------------------------------------------!
!--------- USER INPUT FILE --------------------------------------------------------------------------------------! 
!----------------------------------------------------------------------------------------------------------------!

integer, parameter :: nx = 400                                             		! Image width
integer, parameter :: ny = 400                                             		! Image height
integer, parameter :: nz = 400                                             		! Image depth
integer, parameter :: voxelRes = 1 !DO NOT CHANGE                           	! VoxelResolution of input-raw data --> SHOULD STAY 1, VALUE DOES NOT REALLY HAVE INFLUENCE ON RESULT, MIGHT BE CHANGED TO TYPE REAL IN FUTURE VERSION
integer, parameter :: realPrecision = 2                                    		! Floating-point precision: 1 (single) or 2 (double)
integer, parameter :: nIteration = 500000                                      	! Max. number of iterations
real(kind=4*realPrecision), parameter :: targetResidual = 1e-9            		! Target residual    
integer(kind=1), parameter :: bitLength = 8                                 	! Bites per pixel, typically grey-values in raw images are stored as either 8-bit or 16-bit integers
integer(kind=1), parameter :: nPhases = 5 										! Number of phases in model
integer, dimension(nPhases) :: rawValues = [1, 2, 3, 4, 5]                      ! Bit-grayvalues of each phase (unsigned format)
real(kind=4*realPrecision), dimension(nPhases) :: correspConductivity = [3.98E-15, 3.89E-05, 3.85E-05, 1.25E-05, 5.70E-12]  ! Conductivity of each phase (only double values, e.g. 1.0)
real(kind=4*realPrecision), dimension(nPhases) :: correspPermittivity = [1.00, 9.65, 3.42, 2.74, 3.49]  					! Conductivity of each phase (only double values, e.g. 1.0)
real(kind=4*realPrecision), parameter :: frequency = 1e9                  		! Frequency   
character(512) :: casePath = "/Path/to/case/folder/" ! END STRING WITH SLASH    ! Path to case folder (place raw-data in this folder)
character(512) :: sampleName = "RotondoGranite_Cube400_8bit.raw"                ! Input file
character(1) :: evalDirection = "X"                                         	! Direction of calculation (X or Y or Z)
character(len=20) :: precondType="gamg"											! Choose preconditioner (see https://petsc.org/release/manualpages/PC/PCType/)
character(len=20) :: solverType="bcgs"											! Choose solver (see https://petsc.org/main/manualpages/KSP/KSPType/)

!----------------------------------------------------------------------------------------------!
!--------- Feel free to contact us. We will be happy to help: 	-------------------------------!
!--------- mirko.siegert@web.de 								-------------------------------!
!--------- nils.kerkmann@ieg.fraunhofer.de		 				-------------------------------!
!----------------------------------------------------------------------------------------------!
        
