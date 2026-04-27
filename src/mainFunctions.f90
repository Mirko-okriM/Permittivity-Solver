    integer(kind=bitLength/8) function unsignedToSigned(greyValue) !should work for input integers kind<4 (i.e. 32 bit)
        integer :: greyValue
        if (0 <= greyValue .AND. greyValue <= (2**bitLength)/2-1) then
            unsignedToSigned=greyValue
        else
            unsignedToSigned=greyValue-(2**bitLength)
        end if
    end function unsignedToSigned

    integer function cellPos(x,y,z)
        integer :: x,y,z
        cellPos=nx*ny*(z-1)+nx*(y-1)+x
    end function

    logical function isInnerCell(x,y,z)
        integer :: x,y,z
        isInnerCell=.TRUE.
        if (x == 1 .OR. x == nx .OR. y == 1 .OR. y == ny .OR. z == 1 .OR. z == nz) then
            isInnerCell=.FALSE.
        end if
    end function

    complex(kind=4*realPrecision) function kLookUp(x,y,z)
        integer :: x,y,z,i
        do i = 1, size(rawValues)
            if (currSample(cellPos(x, y, z)) == rawValues(i)) then
              kLookUp = kStar(i)
              exit ! Exit the loop when the condition is met
            end if            
        end do
        !kLookUp=currSample(cellPos(x, y, z))
    end function

    complex(kind=4*realPrecision) function kEval(x1,y1,z1,x2,y2,z2)
        integer :: x1,y1,z1,x2,y2,z2
        complex(kind=4*realPrecision) :: kP, kNB
        kP=kLookUp(x1,y1,z1)
        kNB=kLookUp(x2,y2,z2)
        kEval=2*kP*kNB/(kP+kNB)  !harmonic
        !kEval=(kP+kNB)/2         !arithmetic
    end function

    !functions for postprocessing (determination of flux)
    complex(kind=4*realPrecision) function calcFluxX(xTarget,BC,xField)
        integer :: y, z, xTarget, currCell
        integer(kind=1) :: BC
        complex(kind=4*realPrecision) :: currK
        complex(kind=4*realPrecision), dimension(:) :: xField
        calcFluxX=0.0
        !$omp parallel do collapse(2) private(currK, currCell) reduction(+:calcFluxX)
        do z=1,nz
            do y=1,ny
                currK=kLookUp(xTarget,y,z)
                currCell=cellPos(xTarget,y,z)
                calcFluxX=calcFluxX-currK*(BC-xField(currCell))
            end do
        end do              
        !$omp end parallel do
        calcFluxX=-calcFluxX*(2*dy*dz/dx);
    end function

    complex(kind=4*realPrecision) function calcFluxY(yTarget,BC,xField)
        integer :: x, z, yTarget, currCell
        integer(kind=1) :: BC
        complex(kind=4*realPrecision) :: currK
        complex(kind=4*realPrecision), dimension(:) :: xField
        calcFluxY=0.0
        !$omp parallel do collapse(2) private(currK, currCell) reduction(+:calcFluxY)
        do z=1,nz
            do x=1,nx
                currK=kLookUp(x,yTarget,z)
                currCell=cellPos(x,yTarget,z)
                calcFluxY=calcFluxY-currK*(BC-xField(currCell))
            end do
        end do              
        !$omp end parallel do
        calcFluxY=-calcFluxY*(2*dx*dz/dy);
    end function
    
    complex(kind=4*realPrecision) function calcFluxZ(zTarget,BC,xField)
        integer :: x, y, zTarget, currCell
        integer(kind=1) :: BC
        complex(kind=4*realPrecision) :: currK
        complex(kind=4*realPrecision), dimension(:) :: xField
        calcFluxZ=0.0
        !$omp parallel do collapse(2) private(currK, currCell) reduction(+:calcFluxZ)
        do y=1,ny
            do x=1,nx
                currK=kLookUp(x,y,zTarget)
                currCell=cellPos(x,y,zTarget)
                calcFluxZ=calcFluxZ-currK*(BC-xField(currCell))
            end do
        end do      
        !$omp end parallel do
        calcFluxZ=-calcFluxZ*(2*dx*dy/dz);
    end function
	
	!############### Flux Calculation for PETSc ###############
	subroutine compute_results(solution, kEff, eEff, totalFlux_0, totalFlux_1)
		implicit none
		
		real(kind=4*realPrecision) :: kAveFlux, eAveFlux
		complex(kind=4*realPrecision), intent(in) :: solution(:)
		real(kind=4*realPrecision), intent(out) :: kEff, eEff
		complex(kind=4*realPrecision), intent(out) :: totalFlux_0, totalFlux_1
		
		if (evalDirection=='X') then
			totalFlux_0 = calcFluxX(nx, BC_F, solution) ! totalFluxFront
			totalFlux_1 = calcFluxX(1, BC_B, solution) ! totalFluxBack
			kAveFlux = (abs(real(totalFlux_1)) + abs(real(totalFlux_0))) / 2.0d0
			eAveFlux = (abs(imag(totalFlux_1)) + abs(imag(totalFlux_0))) / 2.0d0
			kEff=kAveFlux/abs(BC_F-BC_B)*(nx*dx)/(ny*dy*nz*dz)
			eEff=eAveFlux/abs(BC_F-BC_B)*(nx*dx)/(ny*dy*nz*dz)/(omega*e0)
		else if (evalDirection=='Y') then
			totalFlux_0 = calcFluxY(ny, BC_E, solution) ! totalFluxEast
			totalFlux_1 = calcFluxY(1, BC_W, solution) ! totalFluxWest
			kAveFlux = (abs(real(totalFlux_1)) + abs(real(totalFlux_0))) / 2.0d0
			eAveFlux = (abs(imag(totalFlux_1)) + abs(imag(totalFlux_0))) / 2.0d0
			kEff=kAveFlux/abs(BC_E-BC_W)*(ny*dy)/(nx*dx*nz*dz)
			eEff=eAveFlux/abs(BC_E-BC_W)*(ny*dy)/(nx*dx*nz*dz)/(omega*e0)
		else if (evalDirection=='Z') then
			totalFlux_0 = calcFluxZ(nz, BC_N, solution) ! totalFluxNorth
			totalFlux_1 = calcFluxZ(1, BC_S, solution) ! totalFluxSouth
			kAveFlux = (abs(real(totalFlux_1)) + abs(real(totalFlux_0))) / 2.0d0
			eAveFlux = (abs(imag(totalFlux_1)) + abs(imag(totalFlux_0))) / 2.0d0
			kEff=kAveFlux/abs(BC_N-BC_S)*(nz*dz)/(nx*dx*ny*dy)
			eEff=eAveFlux/abs(BC_N-BC_S)*(nz*dz)/(nx*dx*ny*dy)/(omega*e0)
		end if
	end subroutine compute_results
