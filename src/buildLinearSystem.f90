!determination of matrix elements for each cell
!matrix is symmetric, stored in full CSR format for PETSc

! Build rPtr in series first
rPtr(1) = 1

do currCell = 2, nx*ny*nz + 1
    if (isInnerCell(mod(currCell-2,nx)+1, mod((currCell-2)/nx,ny)+1, (currCell-2)/(nx*ny)+1)) then
        entries_in_row = 7
    else
        ! Count entries for boundary cells
        x = mod(currCell-2,nx)+1
        y = mod((currCell-2)/nx,ny)+1
        z = (currCell-2)/(nx*ny)+1
        
        entries_in_row = 1  ! diagonal
        if (z>1) entries_in_row = entries_in_row + 1
        if (y>1) entries_in_row = entries_in_row + 1
        if (x>1) entries_in_row = entries_in_row + 1
        if (x<nx) entries_in_row = entries_in_row + 1
        if (y<ny) entries_in_row = entries_in_row + 1
        if (z<nz) entries_in_row = entries_in_row + 1
    endif
    rPtr(currCell) = rPtr(currCell-1) + entries_in_row
end do

! Determine matrix elements in parallel
! Each cell uses its row start position from rPtr
!$omp parallel do collapse(3) private(x,y,z,currCell, &
!$omp   kF,kB,kE,kW,kN,kS,aP,aF,aB,aE,aW,aN,aS,bCell, &
!$omp   entries_in_row, nCounter) schedule(static)
do z = 1, nz
    do y = 1, ny
        do x = 1, nx
            currCell = cellPos(x,y,z)

            if (isInnerCell(x,y,z)) then
                ! Inner cell: always 7 entries
                kF = kEval(x,y,z,x+1,y,z)
                kB = kEval(x,y,z,x-1,y,z)
                kE = kEval(x,y,z,x,y+1,z)
                kW = kEval(x,y,z,x,y-1,z)
                kN = kEval(x,y,z,x,y,z+1)
                kS = kEval(x,y,z,x,y,z-1)
                aP = -(kB+kW+kS+kF+kE+kN)
                b(currCell) = (0,0)

                ! Start index of the line from rPtr:
                nCounter = rPtr(currCell)

                A_mat(nCounter    ) = kS
                colInd(nCounter   ) = cellPos(x,y,z-1)

                A_mat(nCounter + 1) = kW
                colInd(nCounter + 1) = cellPos(x,y-1,z)

                A_mat(nCounter + 2) = kB
                colInd(nCounter + 2) = cellPos(x-1,y,z)

                A_mat(nCounter + 3) = aP
                colInd(nCounter + 3) = currCell

                A_mat(nCounter + 4) = kF
                colInd(nCounter + 4) = cellPos(x+1,y,z)

                A_mat(nCounter + 5) = kE
                colInd(nCounter + 5) = cellPos(x,y+1,z)

                A_mat(nCounter + 6) = kN
                colInd(nCounter + 6) = cellPos(x,y,z+1)

            else
                ! Outer cell
                aF = (0,0); aB = (0,0); aE = (0,0)
                aW = (0,0); aN = (0,0); aS = (0,0)
                aP = (0,0)
                bCell = (0,0)

                ! evalSouth
                if (z > 1) then
                    kS = kEval(x,y,z,x,y,z-1)
                    aP = aP - kS
                    aS = kS
                else
                    if (patchTypeSouth == 'patch') then
                        kS  = kLookUp(x,y,z)
                        aP  = aP - 2*kS
                        bCell = -2*kS*BC_S
                    end if
                end if

                ! evalWest
                if (y > 1) then
                    kW = kEval(x,y,z,x,y-1,z)
                    aP = aP - kW
                    aW = kW
                else
                    if (patchTypeWest == 'patch') then
                        kW  = kLookUp(x,y,z)
                        aP  = aP - 2*kW
                        bCell = -2*kW*BC_W
                    end if
                end if

                ! evalFront
                if (x < nx) then
                    kF = kEval(x,y,z,x+1,y,z)
                    aP = aP - kF
                    aF = kF
                else
                    if (patchTypeFront == 'patch') then
                        kF  = kLookUp(x,y,z)
                        aP  = aP - 2*kF
                        bCell = -2*kF*BC_F
                    end if
                end if

                ! evalBack
                if (x > 1) then
                    kB = kEval(x,y,z,x-1,y,z)
                    aP = aP - kB
                    aB = kB
                else
                    if (patchTypeBack == 'patch') then
                        kB  = kLookUp(x,y,z)
                        aP  = aP - 2*kB
                        bCell = -2*kB*BC_B
                    end if
                end if

                ! evalEast
                if (y < ny) then
                    kE = kEval(x,y,z,x,y+1,z)
                    aP = aP - kE
                    aE = kE
                else
                    if (patchTypeEast == 'patch') then
                        kE  = kLookUp(x,y,z)
                        aP  = aP - 2*kE
                        bCell = -2*kE*BC_E
                    end if
                end if

                ! evalNorth
                if (z < nz) then
                    kN = kEval(x,y,z,x,y,z+1)
                    aP = aP - kN
                    aN = kN
                else
                    if (patchTypeNorth == 'patch') then
                        kN  = kLookUp(x,y,z)
                        aP  = aP - 2*kN
                        bCell = -2*kN*BC_N
                    end if
                end if

                b(currCell) = bCell

                ! Write matrix entries to CSR
                nCounter = rPtr(currCell)

                ! South
                if (z > 1) then
                    A_mat(nCounter)   = aS
                    colInd(nCounter)  = cellPos(x,y,z-1)
                    nCounter = nCounter + 1
                end if

                ! West
                if (y > 1) then
                    A_mat(nCounter)   = aW
                    colInd(nCounter)  = cellPos(x,y-1,z)
                    nCounter = nCounter + 1
                end if

                ! Back
                if (x > 1) then
                    A_mat(nCounter)   = aB
                    colInd(nCounter)  = cellPos(x-1,y,z)
                    nCounter = nCounter + 1
                end if

                ! Diagonal
                A_mat(nCounter)   = aP
                colInd(nCounter)  = currCell
                nCounter = nCounter + 1

                ! Front
                if (x < nx) then
                    A_mat(nCounter)   = aF
                    colInd(nCounter)  = cellPos(x+1,y,z)
                    nCounter = nCounter + 1
                end if

                ! East
                if (y < ny) then
                    A_mat(nCounter)   = aE
                    colInd(nCounter)  = cellPos(x,y+1,z)
                    nCounter = nCounter + 1
                end if

                ! North
                if (z < nz) then
                    A_mat(nCounter)   = aN
                    colInd(nCounter)  = cellPos(x,y,z+1)
                    nCounter = nCounter + 1
                end if

            end if
        end do
    end do
end do
!$omp end parallel do
