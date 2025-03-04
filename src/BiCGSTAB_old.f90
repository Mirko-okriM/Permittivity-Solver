!preconditioned conjugate gradient with Jacobi-Diagonal-Preconditioner
r=b                                                                     !original logic would be r=b-A*x, since x0=zeroVector is chosen r=b follows
deallocate(b)                                                           !b while not be used anymore    
rHat=r
rho=dotProduct(rHat,r) 
p=r
residualFactor=sumOfAbsComponents(r)
currResidual=sumOfAbsComponents(r)/residualFactor
write(*,*) 'debug-bla:', r(1)
write(*,*) 'debug-bla:', residualFactor
write(*,*) 'Iteration 0 residual (L1-Norm): ', currResidual
i=1
do while (i<=nIteration .AND. currResidual>targetResidual)
    call symHeptaMatrixTimesVector(MD, ODp1, ODp2, ODp3, p, posODp1, posODp2, posODp3, v)           ! calc v=A*p
    alpha=rho/dotProduct(rHat,v)                                                                    ! calc alpha=rho/dot(rHat,v)
    call vectorPlusScalarTimesVector(xField,alpha,p,h)                                              ! calc h=x+alpha*p
    call vectorPlusScalarTimesVector(r,-alpha,v,s)                                                  ! calc s=r-alpha*v
    call symHeptaMatrixTimesVector(MD, ODp1, ODp2, ODp3, s, posODp1, posODp2, posODp3, t)           ! calc t=A*s    
    omegaHat=dotProduct(t,s)/dotProduct(t,t)  
    call vectorPlusScalarTimesVector(h,omegaHat,s,xField)                                           ! calc x=h+omegaHat*s
    call vectorPlusScalarTimesVector(s,-omegaHat,t,r)                                               ! calc r=s-omegaHat*t
    rhoOld=rho
    rho=dotProduct(rHat,r)
    beta=rho/rhoOld*alpha/omegaHat
    call vectorPlusScalarTimesVector(p,-omegaHat,v,pHat)                                            ! calc pHat=p-omegaHat*v
    call vectorPlusScalarTimesVector(r,beta,pHat,p)                                                 ! calc p=r+beta*pHat    
    currResidual=sumOfAbsComponents(r)/residualFactor
        !calc current conductivity (code should be improved ... if-condition in every loop is not very smart)
    if (evalDirection=='X') then
        totalFluxBack=calcFluxX(1,BC_B,xField)
        totalFluxFront=calcFluxX(nx,BC_F,xField)
        kAveFluxX=(abs(real(totalFluxBack))+abs(real(totalFluxFront)))/2
        eAveFluxX=(abs(imag(totalFluxBack))+abs(imag(totalFluxFront)))/2
        kEffX=kAveFluxX/abs(BC_F-BC_B)*(nx*dx)/(ny*dy*nz*dz)
        eEffX=eAveFluxX/abs(BC_F-BC_B)*(nx*dx)/(ny*dy*nz*dz)/(omega*e0)
        write(*,*) 'Iteration ', i,' residual (L1-Norm): ', currResidual, ', kEffX= ', kEffX, ', eEffX= ', eEffX
    else if (evalDirection=='Y') then
        totalFluxWest=calcFluxY(1,BC_W,xField)
        totalFluxEast=calcFluxY(ny,BC_E,xField)
        kAveFluxY=(abs(real(totalFluxWest))+abs(real(totalFluxEast)))/2
        eAveFluxY=(abs(imag(totalFluxWest))+abs(imag(totalFluxEast)))/2
        kEffY=kAveFluxY/abs(BC_E-BC_W)*(ny*dy)/(nx*dx*nz*dz)
        eEffY=eAveFluxY/abs(BC_E-BC_W)*(ny*dy)/(nx*dx*nz*dz)/(omega*e0)
        write(*,*) 'Iteration ', i,' residual (L1-Norm): ', currResidual, ', kEffY= ', kEffY, ', eEffY= ', eEffY
    else if (evalDirection=='Z') then
        totalFluxSouth=calcFluxZ(1,BC_S,xField)
        totalFluxNorth=calcFluxZ(nz,BC_N,xField)
        kAveFluxZ=(abs(real(totalFluxSouth))+abs(real(totalFluxNorth)))/2
        eAveFluxZ=(abs(imag(totalFluxSouth))+abs(imag(totalFluxNorth)))/2
        kEffZ=kAveFluxZ/abs(BC_N-BC_S)*(nz*dz)/(nx*dx*ny*dy)
        eEffZ=eAveFluxZ/abs(BC_N-BC_S)*(nz*dz)/(nx*dx*ny*dy)/(omega*e0)
        write(*,*) 'Iteration ', i,' residual (L1-Norm): ', currResidual, ', kEffZ= ', kEffZ, ', eEffZ= ', eEffZ
    end if
    i=i+1
end do  
