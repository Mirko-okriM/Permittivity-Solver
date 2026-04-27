! Numerical solver for the calculation of frequency-dependent 
! electrical conductivity kEff and relative permittivity eEff in porous media
! - diffusion type problem div(k*grad(T))=0 with k=k(x,y,z) in C is solved 
!   utilizing the finite volume method
! - two opposing sides are prescribed with a fixedValue Dirichlet boundary 
!   condition (of dT=1), remaining sides are zeroGradient Neumann boundary 
!   conditions
! - calculation of k at cell faces is performed with i) harmonic mean 
!   or ii) arithmetic mean (adjust in 'function kEval')
! - linear system is solved with PETSc and by default with the biconjugate gradient stabilized method
! - code is parallelized with openmp (number of processors has to be adjusted 
!   with environmental variable 'OMP_NUM_THREADS' (e.g., export OMP_NUM_THREADS=8')
! - compile with:         
!               ${PETSC_DIR}/${PETSC_ARCH}/bin/mpif90 -cpp -fopenmp -O3 \
!					-I${PETSC_DIR}/include \
!					-I${PETSC_DIR}/${PETSC_ARCH}/include \
!					permittivitySolver.f90 \
!					-L${PETSC_DIR}/${PETSC_ARCH}/lib -lpetsc -lm \
!					-Wl,-rpath,${PETSC_DIR}/${PETSC_ARCH}/lib \
!					-o solverRun
! - contact: mirko.siegert@web.de / nils.kerkmann@ieg.fraunhofer.de
! - Changelog: 
! 	- Version 1.0.0 (xx.08.24) MS: 
!		- initial release
! 	- Version 1.1.0 (xx.10.25) NK:
!		- added CSR matrix format
!		- added PETSc solver
!	- Version 1.1.1 (10.12.25) NK:
!		- restructuring for cleaner case studies
!	- Version 1.1.2 (23.01.26) NK:
!		- added MPI-safe PETSc
!-------------------------------------------------------------------------------
!-----------MAIN----------------------------------------------------------------
!-------------------------------------------------------------------------------

program conductivitySolver
#include <petsc/finclude/petscksp.h>
    use petscksp
    use omp_lib
	implicit none
	
    ! USER INPUT
#include "settings.F90"

    ! Additional variables
	integer, parameter :: dx=voxelRes
    integer, parameter :: dy=voxelRes
    integer, parameter :: dz=voxelRes
	integer :: x, y, z, i, j, n, currCell, nCounter
    integer(kind=16), parameter :: nCells=nx*ny*nz
	integer :: entries_in_row  ! Number of entries per line
	integer(kind=1) BC_F, BC_B, BC_E, BC_W, BC_N, BC_S
	integer :: arrayPos
	integer(kind=bitLength/8), allocatable :: currSample(:)
	integer(kind=4), allocatable :: colInd(:), rPtr(:) !arrays for CSR-Format
	integer, allocatable :: temp_colInd(:,:) !temp int for openMP
 
	real(kind=4*realPrecision) :: kEff, eEff
	real(kind=4*realPrecision), allocatable :: qFlux(:) !variables for flux calculation (postprocessing)
	
	real(kind=4*realPrecision) :: T_P, T_F, T_B, T_E, T_W, T_N, T_S
    real(kind=4*realPrecision) :: currResidual

    complex(kind=4*realPrecision) :: kF, kB, kE, kW, kN, kS
    complex(kind=4*realPrecision) :: aP, aF, aB, aE, aW, aN, aS, bCell
	complex(kind=4*realPrecision) :: totalFlux_0, totalFlux_1
	complex(kind=4*realPrecision), allocatable :: A_mat(:), b(:) !arrays for CSR-Format
	complex(kind=4*realPrecision), allocatable :: xField(:)
    complex(kind=4*realPrecision), allocatable :: temp_coeffs(:,:) !temp array for openMP
    
    character(512) :: pathTemperature, pathFlux, pathPost, pathResult                          
    character(512) :: rawDataPath
    character(2000) :: loadDataFile, dimString0, dimString0Vector, dimString1, nxString, nyString, nzString
    character(8) :: patchTypeFront, patchTypeBack, patchTypeEast, patchTypeWest, patchTypeNorth, patchTypeSouth

    ! Constants for kStar
    real(kind=4*realPrecision), parameter :: pi = 3.1415926535897932d0    
    real(kind=4*realPrecision), parameter :: e0 = 8.8541878188d-12                     !permittivity of vacuum (DO NOT CHANGE)
    real(kind=4*realPrecision), parameter :: omega = 2.0d0*pi*frequency                !angular frequency
    real(kind=4*realPrecision) :: tmpReal, tmpImag
    complex(kind=4*realPrecision), dimension(nPhases) :: kStar 
	
	! PETSc variables
    KSP :: ksp
    PetscErrorCode :: petsc_ierr
	
	! MPI variables
	PetscMPIInt :: rank, ierr_mpi
	PetscMPIInt :: rankSize, ierr_test

    do n = 1, size(kStar)
        tmpReal = correspConductivity(n)
        tmpImag = omega*e0*correspPermittivity(n)
        kStar(n) = cmplx(tmpReal,tmpImag)
    end do
	
	! Initialize PETSc
    call PetscInitialize(PETSC_NULL_CHARACTER, petsc_ierr)
	call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr_mpi)
	
	! Debug check correct MPI ranks
	!call MPI_Comm_size(PETSC_COMM_WORLD, rankSize, ierr_test)
	
	!write(*,*) 'My rank: ', rank, ' of total ', rankSize, ' processes'
	!call flush(6)
	!call MPI_Barrier(PETSC_COMM_WORLD, ierr_mpi)
       
    ! Initialize path variables
    write(rawDataPath,"(2A)") TRIM(casePath), TRIM(sampleName)
    write(pathTemperature,"(2A)") TRIM(casePath), TRIM('resTemp.raw')
    write(pathFlux,"(2A)") TRIM(casePath), TRIM('resFlux.raw')
    write(pathPost,"(2A)") TRIM(casePath), TRIM('loadData_paraview.xdmf')
    write(pathResult,"(5A)") TRIM(casePath), TRIM(sampleName), TRIM('_conductivity_'), TRIM(evalDirection), TRIM('.csv')

    ! Openmp might not work if giant arrays are not defined as allocatable variables
    allocate(currSample(nCells))
	allocate(colInd(nCells*7))
    allocate(A_mat(nCells*7))
    allocate(b(nCells))
    allocate(xField(nCells))
	allocate(rPtr(1:nCells+1))
	allocate(temp_coeffs(nx*ny*nz, 7))
	allocate(temp_colInd(nx*ny*nz, 7))
	
	! Initialize arrays with 0 for non-zero ranks 
    if (rank /= 0) then
        A_mat = cmplx(0.0d0, 0.0d0)
        b = cmplx(0.0d0, 0.0d0)
        colInd = 0
        rPtr = 1  ! Important for CSR-Format
    end if

    ! Setup boundary logic
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
    BC_F=0  !fixed front bc value if face is 'patch', if value is changed, than flux calculation might not be correct anymore
    BC_B=1  !fixed back bc value if face is 'patch', if value is changed, than flux calculation might not be correct anymore
    BC_E=0  !fixed east bc value if face is 'patch', if value is changed, than flux calculation might not be correct anymore
    BC_W=1  !fixed west bc value if face is 'patch', if value is changed, than flux calculation might not be correct anymore
    BC_N=0  !fixed north bc value if face is 'patch', if value is changed, than flux calculation might not be correct anymore
    BC_S=1  !fixed south bc value if face is 'patch', if value is changed, than flux calculation might not be correct anymore

	! Only rank 0 reads data and builds matrix
	if (rank == 0) then	! MPI
		write(*,'(A,ES12.5,A)') 'Frequency: ', frequency, ' Hz'
		write(*,*) ''
		
		! Open the file for binary reading
		write(*,*) '1) Reading raw-data ...'
		open(10, file=rawDataPath, form='unformatted', access='stream', status='old')
		read(10) currSample
		close(10)
		write(*,*) '... done'
		
		write(*,*) 'Wert von Position(1,2,1):', currSample(cellPos(1,2,1))

		! Convert rawValues to signed values
		do i = 1, size(rawValues)
			rawValues(i)=unsignedToSigned(rawValues(i))
			write(*,*) 'Wert ',i,': ',rawValues(i)
		end do

		! Determination of matrix elements for each cell
		write(*,*) '2) Setting up linear system ...'
		include './src/buildLinearSystem.f90'
		write(*,*) '... done'   

		write(*,*) '3) Solving linear system with PETSc'
		write(*,*) 'Preconditioner: ', precondType, 'Solver Type: ', solverType
	end if	! MPI
    
	! Synchronize before solving
    call MPI_Barrier(PETSC_COMM_WORLD, ierr_mpi) ! MPI
	
    call solvePETScKSP(A_mat, colInd, rPtr, b, xField, nCells, targetResidual, &
					   nIteration, petsc_ierr, rank)

	if (rank == 0) then	! MPI		
		! Final results
		call compute_results(xField, kEff, eEff, totalFlux_0, totalFlux_1)
		if (evalDirection=='X') then
			write(*,'(A,2ES15.8,A,2ES15.8,A,ES15.8,A,ES15.8)') &
				'totalFluxFront= ', totalFlux_0, &
				', totalFluxBack= ', totalFlux_1, &
				', kEffX= ', kEff, &
				', eEffX= ', eEff    
		else if (evalDirection=='Y') then
			write(*,'(A,2ES15.8,A,2ES15.8,A,ES15.8,A,ES15.8)') &
				'totalFluxEast= ', totalFlux_0, &
				', totalFluxWest= ', totalFlux_1, &
				', kEffY= ', kEff, &
				', eEffY= ', eEff
		else if (evalDirection=='Z') then
			write(*,'(A,2ES15.8,A,2ES15.8,A,ES15.8,A,ES15.8)') &
				'totalFluxNorth= ', totalFlux_0, &
				', totalFluxSouth= ', totalFlux_1, &
				', kEffZ= ', kEff, &
				', eEffZ= ', eEff
		end if
	end if	! MPI
    
	call PetscFinalize(petsc_ierr)
	
    contains
    include './src/mainFunctions.f90'
	
	!#################### PETSc Functions ####################
	subroutine solvePETScKSP(A_vals, colInd, rowP, b_vals, x_vals, nCells, tolerance, &
							 maxIter, ierr, rank)
		use petscksp
		implicit none
		
		integer(kind=4), intent(in) :: colInd(:), rowP(:), maxIter
		integer(kind=16), intent(in) :: nCells
		PetscMPIInt, intent(in) :: rank
		integer, intent(out) :: ierr
		integer :: i, j, row_start, row_end
		integer(kind=4), allocatable :: indices(:)
		real(kind=4*realPrecision), intent(in) :: tolerance
		complex(kind=4*realPrecision), intent(in) :: A_vals(:), b_vals(:)
		complex(kind=4*realPrecision), intent(inout) :: x_vals(:)
		complex(kind=4*realPrecision), allocatable :: temp_values(:)
		
		! PETSc variables
		Mat :: A
		Vec :: b, x, x_seq
		VecScatter :: scatter_ctx
		KSP :: ksp
		PC :: pc
		PetscErrorCode :: local_ierr
		PetscInt :: final_iteration
		PetscReal :: final_residual
		
		! Create matrix
		call MatCreate(PETSC_COMM_WORLD, A, local_ierr)
		call MatSetSizes(A, PETSC_DECIDE_INTEGER, PETSC_DECIDE_INTEGER, int(nCells, 4), int(nCells, 4), local_ierr)
		call MatSetType(A, MATMPIAIJ, local_ierr)
		call MatSetUp(A, local_ierr)
		
		! Fill matrix (only rank 0)
		if (rank == 0) then
			do i = 1, nCells
				row_start = rowP(i)
				row_end = rowP(i+1) - 1
				do j = row_start, row_end
					call MatSetValue(A, i-1, colInd(j)-1, A_vals(j), INSERT_VALUES, local_ierr)
				end do
			end do
		end if
		
		call MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY, local_ierr)
		call MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY, local_ierr)
		
		! Create vectors
		call VecCreate(PETSC_COMM_WORLD, b, local_ierr)
		call VecSetSizes(b, PETSC_DECIDE_INTEGER, int(nCells, 4), local_ierr)
		call VecSetType(b, VECMPI, local_ierr)
		call VecDuplicate(b, x, local_ierr)
		
		! Fill RHS (only rank 0)
		if (rank == 0) then
			do i = 1, nCells
				call VecSetValue(b, i-1, b_vals(i), INSERT_VALUES, local_ierr)
			end do
		end if
		
		call VecAssemblyBegin(b, local_ierr)
		call VecAssemblyEnd(b, local_ierr)
		
		! Starting vector
		call VecSet(x, cmplx(0.0d0, 0.0d0, kind=4*realPrecision), local_ierr)

		
		! Solver setup
		call KSPCreate(PETSC_COMM_WORLD, ksp, local_ierr)
		call KSPSetOperators(ksp, A, A, local_ierr)
		call KSPSetType(ksp, solverType, local_ierr)
		call KSPSetTolerances(ksp, 0.0d0, tolerance, PETSC_DEFAULT_REAL, maxIter, local_ierr)
		
		call KSPGetPC(ksp, pc, local_ierr)
		call PCSetType(pc, precondType, local_ierr)
		call KSPSetNormType(ksp, KSP_NORM_DEFAULT, local_ierr)
		call KSPSetFromOptions(ksp, local_ierr)
		call KSPConvergedDefaultSetUIRNorm(ksp, local_ierr)
		
		call KSPMonitorSet(ksp, conductivity_monitor, PETSC_NULL_POINTER, PETSC_NULL_FUNCTION, local_ierr)
		
		if (rank == 0) then
			write(*,*) 'Solving system...'
		end if
		
		! Solve
		call KSPSolve(ksp, b, x, local_ierr)
		
		! Total iterations and final residual
		call KSPGetIterationNumber(ksp, final_iteration, local_ierr)
		call KSPGetResidualNorm(ksp, final_residual, local_ierr)
		
		if (rank == 0) then
			write(*,*) '... done'
			write(*,*) ''
			write(*,'(A,I0)') ' Total iterations: ', final_iteration
			write(*,'(A,ES15.8)') ' Final residual: ', final_residual
			write(*,*) ''
		end if
		
		! Extract solution 
		call VecScatterCreateToZero(x, scatter_ctx, x_seq, local_ierr)
		call VecScatterBegin(scatter_ctx, x, x_seq, INSERT_VALUES, SCATTER_FORWARD, local_ierr)
		call VecScatterEnd(scatter_ctx, x, x_seq, INSERT_VALUES, SCATTER_FORWARD, local_ierr)
		
		! Nur Rank 0 extrahiert Werte
		if (rank == 0) then
			allocate(indices(nCells))
			allocate(temp_values(nCells))
			
			do i = 1, nCells
				indices(i) = i - 1
			end do
			
			call VecGetValues(x_seq, int(nCells, 4), indices, temp_values, local_ierr)
			
			do i = 1, nCells
				x_vals(i) = temp_values(i)
			end do
			
			deallocate(indices, temp_values)
		end if
		
		! Cleanup
		call VecDestroy(x_seq, local_ierr)
		call VecScatterDestroy(scatter_ctx, local_ierr)
		
		call KSPDestroy(ksp, local_ierr)
		call VecDestroy(b, local_ierr)
		call VecDestroy(x, local_ierr)
		call MatDestroy(A, local_ierr)
	end subroutine solvePETScKSP
	
	subroutine conductivity_monitor(ksp, iteration, residual_norm, ctx, ierr)
		use petscksp
		implicit none

		integer, parameter :: nCells = nx*ny*nz
		integer :: i
		integer(kind=4), allocatable, save :: indices(:)	! save = will not be deleted after 
		logical, save :: first_call = .true.
		logical, save :: rank_initialized = .false.
		real(kind=4*realPrecision) :: kEff, eEff
		complex(kind=4*realPrecision) :: totalFlux_0, totalFlux_1
		complex(kind=4*realPrecision), allocatable :: temp_solution(:)
		
		! PETSc variables
		KSP :: ksp
		PetscInt :: iteration
		PetscReal :: residual_norm  
		PetscVoid :: ctx
		PetscErrorCode :: ierr
		Vec :: x_vec, x_seq
		VecScatter :: scatter_ctx
		PetscMPIInt, save :: rank
		PetscMPIInt :: ierr_mpi
		
		! Calculate interim results only X steps
		if (mod(iteration, 10) /= 0 .and. iteration /= 0) return
		
		! Get MPI ranks
		! call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr_mpi)
		if (.not. rank_initialized) then
			call MPI_Comm_rank(PETSC_COMM_WORLD, rank, ierr_mpi)
			rank_initialized = .true.
		end if
		
		! Get current solution from solver (MPI distributed)
		call KSPGetSolution(ksp, x_vec, ierr)
		
		! These calls must be made on ALL processes
		call VecScatterCreateToZero(x_vec, scatter_ctx, x_seq, ierr)
		call VecScatterBegin(scatter_ctx, x_vec, x_seq, INSERT_VALUES, SCATTER_FORWARD, ierr)
		call VecScatterEnd(scatter_ctx, x_vec, x_seq, INSERT_VALUES, SCATTER_FORWARD, ierr)
		
		! Only rank 0 works with x_seq
		if (rank == 0) then
			! Allocate arrays
			allocate(temp_solution(nCells))
			
			! Create index array (0-based for PETSc) only once
			if (first_call) then
				allocate(indices(nCells))
				do i = 1, nCells
					indices(i) = i - 1
				end do
				first_call = .false.
			end if
			
			! Get all values from the sequential vector (only valid on rank 0)
			call VecGetValues(x_seq, int(nCells, 4), indices, temp_solution, ierr)
			
			! Calculate results
			call compute_results(temp_solution, kEff, eEff, totalFlux_0, totalFlux_1)
			
			! Output
			if (evalDirection=='X') then
				write(*,'(A,I12,A,ES15.8,A,ES15.8,A,ES15.8)') &
					' Iteration: ', iteration, &
					', residual: ', residual_norm, &
					', kEffX= ', kEff, &
					', eEffX= ', eEff
			elseif (evalDirection=='Y') then
				write(*,'(A,I12,A,ES15.8,A,ES15.8,A,ES15.8)') &
					' Iteration: ', iteration, &
					', residual: ', residual_norm, &
					', kEffY= ', kEff, &
					', eEffY= ', eEff
			elseif (evalDirection=='Z') then
				write(*,'(A,I12,A,ES15.8,A,ES15.8,A,ES15.8)') &
					' Iteration: ', iteration, &
					', residual: ', residual_norm, &
					', kEffZ= ', kEff, &
					', eEffZ= ', eEff
			endif
			
			! Cleanup of local temp array
			deallocate(temp_solution)
		end if
		
		! Cleanup must be performed on all processes
		call VecScatterDestroy(scatter_ctx, ierr)
		if (rank == 0) then
			call VecDestroy(x_seq, ierr)
		end if

	end subroutine conductivity_monitor

end program conductivitySolver

