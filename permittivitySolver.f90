! Numerical solver for the calculation of frequency-dependent 
! electrical conductivity kEff and relative permittivity eEff in porous media
! - diffusion type problem div(k*grad(T))=0 with k=k(x,y,z) in C is solved 
!   utilizing the finite volume method
! - two opposing sides are prescribed with a fixedValue Dirichlet boundary 
!   condition (of dT=1), remaining sides are zeroGradient Neumann boundary 
!   conditions
! - calculation of k at cell faces is performed with i) harmonic mean 
!   or ii) arithmetic mean (adjust in 'function kEval')
! - linear system is solved with the biconjugate gradient stabilized method 
! - code is parallelized with openmp (number of processors has to be adjusted 
!   with environmental variable 'OMP_NUM_THREADS', e.g. powershell 
!   '$env:OMP_NUM_THREADS=8')
! - compile with                 
!              gfortran.exe -fopenmp -O3 conductivitySolver.f90 -o solverRun.exe
! - contact: mirko.siegert@web.de
! - Version 1.0 (xx.08.24)
! - Changelog: 
!   ...
!

!-------------------------------------------------------------------------------
!-----------MAIN----------------------------------------------------------------
!-------------------------------------------------------------------------------

program conductivitySolver
    use omp_lib 
    implicit none
    !USER INPUT
    !integer, parameter :: chunkSize=100000 !activate if schedule(dynamic, chunkSize) is activated during parallel-looping
    !integer, parameter :: nProcessors=8
    include './settings.f90'    
    !additional variables
    character(512) :: pathTemperature, pathFlux, pathPost, pathResult                          
    character(512) :: rawDataPath
    character(2000) :: loadDataFile, dimString0, dimString0Vector, dimString1, nxString, nyString, nzString
    integer, parameter :: dx=voxelRes
    integer, parameter :: dy=voxelRes
    integer, parameter :: dz=voxelRes
    integer(kind=bitLength/8), allocatable :: currSample(:)
    !real(kind=4), allocatable :: currSample(:)
    
    integer :: x, y, z, i, j, n, currCell
    integer(kind=16), parameter :: nCells=nx*ny*nz
    integer, parameter :: posODp1=1                 !cellPos(3,2,2)-cellPos(2,2,2)
    integer, parameter :: posODm1=-1                !cellPos(1,2,2)-cellPos(2,2,2)
    integer, parameter :: posODp2=nx                !cellPos(2,3,2)-cellPos(2,2,2)
    integer, parameter :: posODm2=-nx               !cellPos(2,1,2)-cellPos(2,2,2)
    integer, parameter :: posODp3=nx*ny             !cellPos(2,2,3)-cellPos(2,2,2)
    integer, parameter :: posODm3=-nx*ny            !cellPos(2,2,1)-cellPos(2,2,2)
    complex(kind=4*realPrecision), allocatable :: MD(:), b(:), ODp1(:), ODp2(:), ODp3(:)!, ODm1(:), ODm2(:), ODm3(:) not necessary since Matrix is symmetric
    complex(kind=4*realPrecision) :: kF, kB, kE, kW, kN, kS
    complex(kind=4*realPrecision) :: aP, aF, aB, aE, aW, aN, aS, bCell
    real(kind=4*realPrecision) :: T_P, T_F, T_B, T_E, T_W, T_N, T_S
    real(kind=4*realPrecision) :: currResidual
    complex(kind=4*realPrecision) :: totalFluxFront, totalFluxBack, totalFluxEast, totalFluxWest, totalFluxNorth, totalFluxSouth
    real(kind=4*realPrecision) :: kAveFluxX, kAveFluxY, kAveFluxZ, kEffX, kEffY, kEffZ
    real(kind=4*realPrecision) :: eAveFluxX, eAveFluxY, eAveFluxZ, eEffX, eEffY, eEffZ
    integer(kind=1) BC_F, BC_B, BC_E, BC_W, BC_N, BC_S
    character(8) :: patchTypeFront, patchTypeBack, patchTypeEast, patchTypeWest, patchTypeNorth, patchTypeSouth
    complex(kind=4*realPrecision), allocatable :: r(:), rHat(:), p(:), pHat(:), xField(:), v(:), h(:), s(:), t(:) !variables for conjugate gradient method
    complex(kind=4*realPrecision) :: alpha, beta, omegaHat, rho, rhoOld
    real(kind=4*realPrecision) :: residualFactor
    real(kind=4*realPrecision), allocatable :: qFlux(:) !variables for flux calculation (postprocessing)
    integer :: arrayPos
    
    !constants for kStar
    !real(kind=4*realPrecision), parameter :: pi = atan(1.0)*4.0                        !define pi
    real(kind=4*realPrecision), parameter :: pi = 3.1415926535897932    
    real(kind=4*realPrecision), parameter :: e0 = 8.8541878188e-12                     !permittivity of vacuum (DO NOT CHANGE)
    real(kind=4*realPrecision), parameter :: omega = 2*pi*frequency                    !angular frequency
    complex(kind=4*realPrecision), dimension(nPhases) :: kStar 
    real(kind=4*realPrecision) :: tmpReal, tmpImag   
    
    do n = 1, size(kStar)
        tmpReal = correspConductivity(n)
        tmpImag = omega*e0*correspPermittivity(n)
        kStar(n) = cmplx(tmpReal,tmpImag)
    end do
    
    !initialize path variables
    write(rawDataPath,"(2A)") TRIM(casePath), TRIM(sampleName)
    write(pathTemperature,"(2A)") TRIM(casePath), TRIM('resTemp.raw')
    write(pathFlux,"(2A)") TRIM(casePath), TRIM('resFlux.raw')
    write(pathPost,"(2A)") TRIM(casePath), TRIM('loadData_paraview.xdmf')
    write(pathResult,"(5A)") TRIM(casePath), TRIM(sampleName), TRIM('_conductivity_'), TRIM(evalDirection), TRIM('.csv')

    !openmp might not work if giant arrays are not defined as allocatable variables
    allocate(currSample(nCells))
    allocate(MD(nCells))
    allocate(b(nCells))
    allocate(ODp1(nCells-posODp1))
    allocate(ODp2(nCells-posODp2))
    allocate(ODp3(nCells-posODp3))
    allocate(r(nCells))
    allocate(rHat(nCells))
    allocate(p(nCells))
    allocate(pHat(nCells))
    allocate(xField(nCells))
    allocate(v(nCells))
    allocate(h(nCells))
    allocate(s(nCells))
    allocate(t(nCells))

    !call omp_set_num_threads(nProcessors)

    !Setup boundary logic
    patchTypeFront='wall'   !x+
    patchTypeBack='wall'    !x-
    patchTypeEast='wall'    !y+
    patchTypeWest='wall'    !y-
    patchTypeNorth='wall'   !z+
    patchTypeSouth='wall'   !z-
    if (evalDirection=='Z') then
        patchTypeNorth='patch'  !z+
        patchTypeSouth='patch'  !z-
    elseif (evalDirection=='X') then
        patchTypeFront='patch'  !x+
        patchTypeBack='patch'   !x-
    elseif (evalDirection=='Y') then
        patchTypeEast='patch'   !y+
        patchTypeWest='patch'   !y-
    endif
    BC_F=0  !fixed front bc value if face is 'patch', if value is changed, than fluxCalculation might not be correct anymore
    BC_B=1  !fixed back bc value if face is 'patch', if value is changed, than fluxCalculation might not be correct anymore
    BC_E=0  !fixed east bc value if face is 'patch', if value is changed, than fluxCalculation might not be correct anymore
    BC_W=1  !fixed west bc value if face is 'patch', if value is changed, than fluxCalculation might not be correct anymore
    BC_N=0  !fixed north bc value if face is 'patch', if value is changed, than fluxCalculation might not be correct anymore
    BC_S=1  !fixed south bc value if face is 'patch', if value is changed, than fluxCalculation might not be correct anymore

    ! Open the file for binary reading
    write(*,*) '1) Reading raw-data ...'
    open(10, file=rawDataPath, form='unformatted', access='stream', status='old')
    !open(10, file=rawDataPath, form='unformatted', access='stream', status='old', convert='big_endian')
    read(10) currSample
    close(10)
    write(*,*) '... done'
    
    write(*,*) 'Wert von Position(1,2,1):', currSample(cellPos(1,2,1))

    !Convert rawValues to signed values
    do i = 1, size(rawValues)
        rawValues(i)=unsignedToSigned(rawValues(i))
        write(*,*) 'Wert ',i,': ',rawValues(i)
    end do

    !determination of matrix elements for each cell
    write(*,*) '2) Setting up linear system ...'
    include './src/buildLinearSystem.f90'
    write(*,*) '... done' 
    
    !export linear system
!    write(*,*) '2b) Export linear system ...'
!    include './src/exportLinearSystem.f90'
!    write(*,*) '... done'    
    
    !start BiCGSTAB solver
    write(*,*) '3) Solving linear system ...'    
    include './src/BiCGSTAB.f90'    
    write(*,*) '... done'
    write(*,*) ''
    
    write(*,*) 'Final result:'
    open(1, file=pathResult, access='stream', status='replace', form='formatted')
    write(1,*) 'Iteration ', i,' residual (L1-Norm): ', currResidual
    if (evalDirection=='X') then
        write(*,*) 'totalFluxFront= ', totalFluxFront, ', totalFluxBack= ', totalFluxBack, ', kEffX= ', kEffX, ', eEffX= ', eEffX    
        write(1,*) 'totalFluxFront= ', totalFluxFront, ', totalFluxBack= ', totalFluxBack, ', kEffX= ', kEffX, ', eEffX= ', eEffX
    else if (evalDirection=='Y') then
        write(*,*) 'totalFluxEast= ', totalFluxEast, ', totalFluxWest= ', totalFluxWest, ', kEffY= ', kEffY, ', eEffY= ', eEffY
        write(1,*) 'totalFluxEast= ', totalFluxEast, ', totalFluxWest= ', totalFluxWest, ', kEffY= ', kEffY, ', eEffY= ', eEffY
    else if (evalDirection=='Z') then
        write(*,*) 'totalFluxNorth= ', totalFluxNorth, ', totalFluxSouth= ', totalFluxSouth, ', kEffZ= ', kEffZ, ', eEffZ= ', eEffZ
        write(1,*) 'totalFluxNorth= ', totalFluxNorth, ', totalFluxSouth= ', totalFluxSouth, ', kEffZ= ', kEffZ, ', eEffZ= ', eEffZ
    end if
    close(1)

    !postprocessing (if activated)
    if (writeTempertureField) then
        write(*,*) 'Writing temperature field ...'
        open(1, file=pathTemperature, access='stream', status='replace', form='unformatted')
        write(1) real(xField,4)
        close(1)
        write(*,*) '... done'
    end if
    
    if (writeFluxField) then
        write(*,*) 'Writing flux field ...'         
        include './src/fluxCalculation.f90'
        open(1, file=pathFlux, access='stream', status='replace', form='unformatted')
        write(1) real(qFlux,4)
        close(1)
        write(*,*) '... done'
    end if    
 
    include './src/writeParaviewFile.f90'    
    open(1, file=pathPost, access='stream', status='replace', form='formatted')
    write(1,'(A)') loadDataFile
    close(1)

    !write(*,*) "Press Enter to exit..."
    !read(*,*)
    
    contains
    include './src/mainFunctions.f90'
end program conductivitySolver
