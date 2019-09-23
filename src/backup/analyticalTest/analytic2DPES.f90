!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Purpose of Program:                                       !
!   TODO:
!                                                           !
! Reference:                                                !
!   1. Model surface 1: symmetric PES                       !
!   Theor Chem Acc 1998, 100, 285-299.                      !
! https://f1000.com/work/item/3809905/resources/2861282/pdf !
!   2. Model surface 2: both sym. and asym. PES             !
!   Theor Chem Acc 2004, 112, 40-51.                        !
! https://f1000.com/work/item/5637171/resources/4634218/pdf !
!                                                           !
! Output:                                                   !
!   1. IRC.*.dat                                            !
!   2. ModifiedIRC.*.dat                                    !
!                                                           !
! History:                                                  !
!   2018/07/20, Grace                                       !
!   2018/07/23, Grace, modify the structure for 2 surfaces  !
!   2018/08/10, Grace, add model potential 2 and 3          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program main
    implicit none
    real(8) :: E,x,y,ds,theta
    real(8),dimension(2)    :: G,modifiedG
    real(8),dimension(2)    :: EigValue,EigVector !smallest eigenvector 
    ! real(8),dimension(2)    :: unitNegEigV,unitOrthoNegEigV
    integer(4)  :: modelNum,i,j,npt,num_path
    real(8),parameter   :: zero=0.0001D0,ngrid=2
    character(len=100)  :: filenum,filename

    call printFun(modelNum)
    ds=0.01D0
    npt=1000
    num_path=50

    ! Running IRC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do i = 1,num_path 
    !     if ( i .eq. 1 ) then 
    !         call setVRI(x,y)
    !     else
    !         call setCoord(x,y)
    !     end if
    !     write(filenum,*) i
    !     filename='IRC.'//TRIM(ADJUSTL(filenum))//'.dat'
    !     open(100,file=filename,status='replace')
    !     do j=0,npt
    !         call getE(modelNum,E,x,y)
    !         call getG(modelNum,2,G,x,y) 
    !         call normVec(2,G)
    !         write(100,"(5(F20.9,1X))") x,y,E,G(1),G(2)
    !         if ( ( NORM2(G) .gt. - zero ) .and. ( NORM2(G) .lt. zero ) ) then
    !             call getEigV(modelNum,2,EigValue,EigVector,x,y)
    !             G(1) = - EigVector(1)
    !             G(2) = - EigVector(2)
    !         end if
    !         x = x - G(1) * ds
    !         y = y - G(2) * ds     
    !     end do
    !     close(100)
    ! end do
    
    ! Running modified IRC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do i=1,num_path
    ! ! do i=1,1
    !     if ( i .eq. 1 ) then 
    !         call setVRI(x,y)
    !     else
    !         call setCoord(x,y)
    !     end if
    !     write(filenum,*) i
    !     filename='ModifiedIRC.'//TRIM(ADJUSTL(filenum))//'.dat'
    !     open(100,file=filename,status='replace')
    !     do j=1,npt
    !         call getE(modelNum,E,x,y)
    !         call getG(modelNum,2,G,x,y)
    !         call normVec(2,G)
    !         call getEigV(modelNum,2,EigValue,EigVector,x,y)
    !         call normVec(modelNum,2,EigVector)

    !         modifiedG(1) = G(1) - 2*(DOT_PRODUCT(G,EigVector)/DOT_PRODUCT(EigVector,EigVector))*EigVector(1)
    !         modifiedG(2) = G(2) - 2*(DOT_PRODUCT(G,EigVector)/DOT_PRODUCT(EigVector,EigVector))*EigVector(2)
    !         call normVec(2,modifiedG)

    !         x = x - modifiedG(1) * ds
    !         y = y - modifiedG(2) * ds     
    !         write(100,"(6(F20.9,1X))") x,y,E,modifiedG(1),modifiedG(2),DOT_PRODUCT(G,EigVector)
    !     end do
    !     close(100)
    ! end do

    ! Finding VRI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    open(100,file='Path1.dat',status='replace')
    open(200,file='Path2.dat',status='replace')
    call setVRI(x,y)
    ! temp condition, because it is bad
    if ( x .eq. 0.0D0 .and. y .eq. 0.0D0 ) then 
        !write(*,'(A)') 'Starts from VRI'
        ! Path 1
        do i=1,npt
            call getE(modelNum,E,x,y)
            call getG(modelNum,2,G,x,y)
            if ( i .eq. 1) then
                call getEigV(modelNum,2,EigValue,EigVector,x,y)
                G(1) = G(1) + EigVector(1)
                G(2) = G(2) + EigVector(2) 
            end if
            write(100,"(5(F15.9,1X))") x,y,E,G(1),G(2)
            x = x - G(1) * ds
            y = y - G(2) * ds    
        end do
        call setVRI(x,y)
        ! Path 2
        do i=1,npt
            call getE(modelNum,E,x,y)
            call getG(modelNum,2,G,x,y)
            if ( i .eq. 1) then
                call getEigV(modelNum,2,EigValue,EigVector,x,y)
                G(1) = G(1) - EigVector(1)
                G(2) = G(2) - EigVector(2) 
            end if
            write(200,"(5(F15.9,1X))") x,y,E,G(1),G(2)
            x = x - G(1) * ds
            y = y - G(2) * ds    
        end do
    end if
    close(100)
    close(200)
end program main

subroutine printFun(modelNum)
    implicit none
    integer(4),intent(out)  :: modelNum
    write(*,'(A)') 'Select a model potential'
    write(*,'(A)') 'Symmetric PES:'
    write(*,'(A)') 'Model 1: E(x,y) = 2y + y^2 +(y+0.4x^2)*x^2'
    write(*,'(A)') 'Model 2: E(x,y) = 0.5(xy^2-yx^2-2x+2y)+(1/30)(x^4+y^4)' 
    write(*,'(A)') 'Asymmetric PES:'
    write(*,'(A)') 'Model 3: E(x,y) = 0.5(xy^2-yx^2-1x+2y)+(1/30)(x^4+y^4)' 
    write(*,'()')
    write(*,'(A)',advance='no') 'Please key-in 1-3: '
    read(*,*) modelNum
    write(*,'()')
    write(*,'(A)') 'Running a lot paths with different initial points,'
    write(*,'(A)') 'VRI point is the first point.'
    write(*,'()')
    return
end subroutine printFun

subroutine getE(modelNum,E,x,y)
    implicit none
    integer(4),intent(in)   :: modelNum
    real(8),intent(out) :: E
    real(8),intent(in)  :: x,y
    select case (modelNum)
    case (1)
        E = 2*y + y**2 + ( y + 0.4*x**2 ) * x**2
    case (2)
        E = 0.5*(x*y**2-y*x**2-2*x+2*y)+(1/30)*(x**4+y**4)
    case (3)
        E = 0.5*(x*y**2-y*x**2-1*x+2*y)+(1/30)*(x**4+y**4)
    end select
    return
end subroutine getE

subroutine setCoord(x,y)
    implicit none
    real(8),intent(out) :: x,y
    ! local variable
    integer(4)  :: seed
    real(8) :: random
    call random_seed(seed)  ! for generate random number
                            ! -1 <= x <= 1, -1 <= y <= 1
    call random_number(random)
    x = -0.5 + random ! -1
    call random_number(random)
    y = -0.5 + random
    ! write(*,'(2(F9.4))') x,y
    return
end subroutine setCoord

subroutine setVRI(x,y)
    implicit none
    real(8),intent(out) :: x,y
        x = 0.0D0
        y = 0.0D0
    return
end subroutine setVRI

subroutine getG(modelNum,dimen,G,x,y)
    implicit none
    integer(4),intent(in)   :: dimen,modelNum
    real(8),dimension(2),intent(inout)    :: G
    real(8),intent(in)  :: x,y
    G = 0.0D0
    select case (modelNum)
    case (1)
        G(1) = 0.8*x**3+2*x*(0.4*x**2+y)
        G(2) = 2+x**2+2*y
    case (2)
        G(1) = (2/15)*x**3-y*x+0.5*y**2-2
        G(2) = -0.5*x**2+x*y+(2/15)*y**3+2
    case (3)
        G(1) = (2/15)*x**3-y*x+0.5*y**2-1
        G(2) = -0.5*x**2+x*y+(2/15)*y**3+2
    end select
    return
end subroutine getG

subroutine normVec(dimen,G)
    implicit none
    integer(4),intent(in)   :: dimen
    real(8),dimension(2),intent(inout)   :: G
    ! local variable
    real(8),dimension(dimen)    :: notNormG
    real(8),parameter   ::  zero=0.00001D0
    integer(4)  :: i
    notNormG = G
    ! if it is not null vector
    if ( (NORM2(notNormG) .lt. -zero) .and. (NORM2(notNormG) .gt. zero) ) then
        do i=1,dimen
            G(i) = notNormG(i)/NORM2(notNormG)
        end do
    end if
    return
end subroutine normVec

subroutine getEigV(modelNum,dimen,EigValue,EigVector,x,y)
    implicit none
    integer(4),intent(in)   :: dimen,modelNum
    real(8),dimension(dimen),intent(inout)  :: EigVector
    real(8),intent(in)  :: x,y
    real(8),dimension(dimen),intent(inout)  :: EigValue
    ! local variables
    real(8),dimension(dimen,dimen)  :: H,eigV
    integer(4)  :: i,j
    call getH(modelNum,2,H,x,y)
    call dia(2,H,EigValue,eigV)
    ! write(*,'(A)') 'test'
    do j=1,dimen
        EigVector(j) = eigV(j,1)
    end do

    ! do i=1,dimen
    !     ! write(*,*) eigE(i)
    !     ! write(*,*) eigV(i,:)
    !     if ( EigValue(i) .le. 0.0D0 ) then
    !         ! select the negative eigenvalue
    !         ! write(*,'(A)',advance='no') 'The negative eigenvalue is '
    !         ! write(*,*)  eigE(i)
    !         do j=1,dimen
    !             NegEigV(j) = eigV(i,j)
    !         end do
    !     else 
    !         ! select the smallest eigenvalue
    !         do j=1,dimen
    !             NegEigV(j) = eigV(1,j)
    !         end do
    !     end if
    ! end do
    return
end subroutine getEigV

subroutine getH(modelNum,dimen,H,x,y)
    implicit none
    integer(4),intent(in)   :: dimen,modelNum
    real(8),dimension(dimen,dimen),intent(inout)    :: H
    real(8),intent(in)  :: x,y
    select case (modelNum)
    case (1)
        H(1,1) = 4*x**2 + 2*(0.4*x**2 +y)
        H(1,2) = 2*x
        H(2,1) = 2*x
        H(2,2) = 2
    case(2)
        H(1,1) = (2/5)*x**2 - y
        H(1,2) = -x + y
        H(2,1) = -x +y
        H(2,2) = x+(2/5)*y**2
    case(3)
        H(1,1) = (2/5)*x**2 - y
        H(1,2) = -x + y
        H(2,1) = -x +y
        H(2,2) = x+(2/5)*y**2
    end select
    return
end subroutine getH

subroutine dia(ngrid,H,eigE,eigV)
    implicit none
    integer(4),intent(in)   :: ngrid
    real(8),dimension(ngrid,ngrid),intent(inout)   :: H,eigV
    real(8),dimension(ngrid),intent(inout)  :: eigE
    ! local variables
    integer(4),parameter    :: LWMAX=10
    real(8),dimension(LWMAX)    :: WORK
    integer(4)  :: INFO,LWORK
    LWORK = -1
    eigE(:) = 0.0D0
    eigV(:,:) = 0.0D0
    call DSYEV('Vectors','Upper',ngrid,H,ngrid,eigE,WORK,LWORK,INFO)
    LWORK = MIN(LWMAX,INT(WORK(1)))
    ! solve eigenproblem
    call DSYEV('Vectors','Upper',ngrid,H,ngrid,eigE,WORK,LWORK,INFO)
    eigV = H
    ! check for convergence
    if ( INFO .gt. 0 ) then
        write(*,'(A)') 'The algorithm,DSYEV,failed to compute eigenvalues.'
        stop
    end if
    return
end subroutine dia

subroutine getUnit(dimen,NegEigV,unitNegEigV)
    implicit none
    integer(4),intent(in)   :: dimen
    real(8),dimension(dimen),intent(in) :: NegEigV
    real(8),dimension(dimen),intent(inout)    :: unitNegEigV!,unitOrthoNegEigV
    ! local variables
    integer(4)   :: i
    real(8),parameter :: angle=1.57 ! temp angle in radian
    real(8),dimension(dimen)    :: tmpUnit
    do i=1,dimen
        unitNegEigV(i) = NegEigV(i)/NORM2(NegEigV)
        ! write(*,*) unitNegEigV(i),NegEigV(i)
    end do
    ! Rotate vector for 90 degree
    !tmpUnit(1) = NegEigV(1)*COS(angle) - NegEigV(2)*SIN(angle)
    !tmpUnit(2) = NegEigV(1)*SIN(angle) + NegEigV(2)*COS(angle)
    ! write(*,'(A)') 'test'
    !do i=1,dimen
    !    unitOrthoNegEigV(i) = tmpUnit(i)/Norm2(tmpUnit)
    !    ! write(*,*) G(i),unitNegEigV(i),unitOrthoNegEigV(i)
    !end do
    ! write(*,*)  DOT_PRODUCT(unitNegEigV,unitOrthoNegEigV)
    return
end subroutine getUnit