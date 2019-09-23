!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program:                                                                          !
!   Generate new strcutures by following specific gradient direction                !                                                       
!                                                                                   !
! Reference:                                                                        !
!                                                                                   !  
! Input:                                                                            !
!   $1 = # of grid points in 1D                                                     !
!   $2 = NAtoms                                                                     !
!                                                                                   !
! Output:                                                                           !
!   newCoord.xyz                                                                    !
!                                                                                   !
! History:                                                                          !
!   2019/06/23, Grace,                                                              !
!   2019/07/12, Grace, modify the input gradients form extracting the central       ! 
!       gradient to nearby gradient.                                                !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program main
    implicit none
    character(len=100)   :: file
    integer(4)  :: num1DPts, NAtoms, i, j, k
    real(8) :: ds
    real(8),allocatable,dimension(:,:)  :: unitX, unitY, G, coord, newCoord
    integer(4),allocatable,dimension(:) :: AN
    real(8),allocatable,dimension(:)    :: AW

    ! Step 1. Use command-line argument to allocate the size of array
    ! (1) amount of grid along 1D (2) amount of atoms in one molecule
        call GETARG(1,file)
        read(file,'(I)') num1DPts
        call GETARG(2,file)
        read(file,'(I)') NAtoms
        allocate(unitX(NAtoms,3))
        allocate(unitY(NAtoms,3))
        allocate(G(NAtoms,3))
        allocate(coord(NAtoms,3))
        allocate(newCoord(NAtoms,3))
        allocate(AN(NAtoms))
        allocate(AW(NAtoms))
        ! Change atomic number to atomic weight
        open(100, file='AN.dat', status='old', action='read')
        do i = 1, NAtoms
            read(100,*) AN(i)
            call AN2AW(AN(i),AW(i))
        end do
        close(100)

    ! Main part: Use two loops to assign the filename of selected points
    do i = 1, num1DPts - 1 ! x coord; exclude the TS2 IRC
        do j = 1, (num1DPts - 1 ) / 2 ! y coord; exclude the TS1 IRC 
        ! Step 2. Define gradient unit vectors of xhat and yhat
            call getXandY(NAtoms, num1DPts, i, j,  AW, unitX, unitY)
        ! Step 3. Extract the gradient and structure of selected point
            call getSelectPts(NAtoms, i, j, AW, G, coord)
        ! Step 4. Generate a new gradient from substracting the 
            !   vector projection of g'
            call genNewGrad(NAtoms, unitX, unitY, G)
        ! Step 5. Shift the selected structure along the new gradient
            ! 5.1 Test range of ds
            ! do k = -14,14,4  
            do k = -30,-18,4
                ds = k * 0.1D0
                call genNewStruc(NAtoms, coord, G, ds, newCoord)
                call writeNewStruc(NAtoms, AN, newCoord, i, j, ds)
            end do
        end do
    end do
    deallocate(unitX)
    deallocate(unitY)
    deallocate(G)
    deallocate(coord)
    deallocate(newCoord)
    deallocate(AN)
    deallocate(AW)
stop
end program main

subroutine getXandY(NAtoms, num1DPts, indexX, indexY, AW, unitX, unitY)
    implicit none
    integer(4),intent(in)  :: NAtoms, num1DPts, indexX, indexY
    real(8),dimension(NAtoms,3),intent(inout)  :: unitX, unitY
    real(8),dimension(NAtoms),intent(in)    :: AW
    ! local variables
    integer(4)  :: i
    character(len=100)  :: buffer1, buffer2, file1, file2
    ! Step 2. Define gradient unit vectors of xhat and yhat
    ! 2.1 Extract the correct files; assign proper filenames
        ! 2019/07/12, extract the gradient from its neighboor
        write(buffer1,*) indexX
        write(buffer2,*) indexY - 1 
        file1 = TRIM(ADJUSTL(buffer1))//'_'//TRIM(ADJUSTL(buffer2))//'.grad'
        write(buffer1,*) indexX + 1
        write(buffer2,*) indexY
        file2 = TRIM(ADJUSTL(buffer1))//'_'//TRIM(ADJUSTL(buffer2))//'.grad'
        open(10, file=file1, status='old', action='read')
        open(20, file=file2, status='old', action='read')

        do i = 1, NAtoms
            read(10,*) unitX(i,1:3)
            read(20,*) unitY(i,1:3)
        end do
        close(10)
        close(20)
    ! 2.2 Change unit from Eh/Bohr to Eh/Mass-weighted; B2MW()
        call B2MW(NAtoms, 3, AW, unitX)
        call B2MW(NAtoms, 3, AW, unitY)
    ! 2.3 Normalize; NormV()
        call NormV(NAtoms, 3, unitX)
        call NormV(NAtoms, 3, unitY)
end subroutine getXandY

subroutine B2MW(dimen1, dimen2, AW, unitVec)
    implicit none
    integer(4), intent(in)  :: dimen1, dimen2
    real(8), dimension(dimen1)   :: AW
    real(8), dimension(dimen1, dimen2), intent(inout)    :: unitVec
    ! local variables
    real(8), parameter  :: Bohr2Ang=0.52977249D0
    integer(4)  :: i, j

    do j = 1, dimen2
        do  i = 1, dimen1
            unitVec(i,j) = unitVec(i,j) / (Bohr2Ang * SQRT(AW(i)) )
        end do
    end do
    return 
end subroutine B2MW

subroutine AN2AW(AN,AW)
    implicit none
    integer(4),intent(in)   :: AN
    real(8),intent(out) :: AW
    ! change atomic number to standard atomic weight
    !ref: https://en.wikipedia.org/wiki/Periodic_table
    Select case (AN)
    case(1) ! hydrogen
        AW=1.0008
    case(6) ! carbon
        AW=12.011
    case (7) ! nitrogen
        AW=14.007
    case (8) ! oxygen 
        AW=15.999
    case (9) ! florine
        AW=18.998
    case (15) ! phosphorus
        AW=30.974
    case default
        write(*,'(A)') 'Modified the subroutine "ANtoAW"'
        write(*,'(A)') 'Stop the program, GenGradStruc.f90.'
        Stop
    end select 
end subroutine AN2AW

subroutine NormV(dimen1, dimen2, Vec)
    implicit none
    integer(4),intent(in)   :: dimen1,dimen2
        real(8),dimension(dimen1,dimen2),intent(inout)   :: Vec
        ! local variable
        real(8),dimension(dimen1*dimen2)    :: arrayVec
        real(8),dimension(dimen1,dimen2)    :: originalVec
        real(8),parameter   ::  zero=0.0000001D0
        integer(4)  :: i,j,index
        index=1
        do j=1,dimen2
            do i=1,dimen1
                arrayVec(index) = originalVec(i,j)
                index = index + 1
            end do
        end do
        originalVec = Vec
        ! if it is not null vector
        if ( (NORM2(arrayVec) .lt. -zero) .or. (NORM2(arrayVec) .gt. zero) ) then
            do j=1,dimen2
                do i=1,dimen1
                    Vec(i,j) = originalVec(i,j)/NORM2(arrayVec)
                end do
            end do
        end if
    return
end subroutine NormV

subroutine getSelectPts(NAtoms,indexX,indexY,AW,Grad,coord)
    implicit none
    integer(4), intent(in)  :: NAtoms, indexX, indexY
    real(8), dimension(NAtoms), intent(in)  :: AW
    real(8), dimension(NAtoms, 3), intent(inout)    :: Grad, coord
    ! local variables
    integer(4)  :: i
    character(len=100)  :: buffer1, buffer2, filename1, filename2
    
    ! Step 3. Extract the gradient and structure of selected point
    write(buffer1,*) indexX
    write(buffer2,*) indexY
    filename1 = TRIM(ADJUSTL(buffer1))//'_'//TRIM(ADJUSTL(buffer2))//'.grad'
    filename2 = TRIM(ADJUSTL(buffer1))//'_'//TRIM(ADJUSTL(buffer2))//'.xyz'
    open(10, file=filename1, status='old', action='read')
    open(20, file=filename2, status='old', action='read')
    read(20,*) buffer1
    read(20,*) buffer1
    do i = 1, NAtoms
        read(10,*) Grad(i,1:3)
        read(20,*) buffer1, coord(i,1:3)
    end do
    close(10)
    close(20)
    ! 3.1 Change unit; B2MW()
    call B2MW(NAtoms, 3, AW, Grad)
    ! 3.2 Normalize; NormV()
    call NormV(NAtoms, 3, Grad)
    return
end subroutine getSelectPts

subroutine genNewGrad(NAtoms,unitX,unitY,G)
    implicit none
    integer(4), intent(in)  :: NAtoms
    real(8), dimension(NAtoms,3), intent(in)    :: unitX, unitY
    real(8), dimension(NAtoms,3), intent(inout) :: G
    ! local variables
    integer(4)  :: i, j
    real(8), dimension(NAtoms,3)    :: progX, progY
    real(8) :: XdotG, YdotG
    ! Step 4. Generate a new gradient from substracting the 
        !   vector projection of g'
        ! 4.1 vector projection of g' on xhat and yhat
        ! x' = (g' * xhat) * xhat 
        ! y' = (g' * yhat) * yhat
            XdotG = 0.0D0
            YdotG = 0.0D0
            do j = 1, 3
                do  i = 1, NAtoms
                    XdotG = XdotG + unitX(i,j) * G(i,j)
                    YdotG = XdotG + unitY(i,j) * G(i,j)
                end do
            end do
            do j = 1, 3
                do i = 1, NAtoms
                    progX(i,j) = XdotG * G(i,j)
                    progY(i,j) = YdotG * G(i,j)
                end do
            end do
        ! 4.2 Generate a new gradient; g'' = g' - x' - y'
            do j =  1, 3
                do i = 1, NAtoms 
                    G(i,j) = G(i,j) - progX(i,j) - progY(i,j)
                end do
            end do
        ! 4.3 Normalize; NormV()
            call NormV(NAtoms, 3, G)
    return
end subroutine genNewGrad

subroutine genNewStruc(NAtoms, coord, grad, ds, newCoord)
    implicit none
    integer(4), intent(in)  :: NAtoms
    real(8), dimension(NAtoms, 3), intent(in)   :: coord, grad
    real(8), dimension(NAtoms, 3), intent(out)   :: newCoord
    real(8), intent(in) :: ds
    ! local variables
    integer(4)  :: i, j

    newCoord = 0.0D0
    do j = 1, 3
        do i = 1, NAtoms 
             newCoord(i, j) = coord(i,j) + ds * grad(i,j)
        end do
    end do
    return
end subroutine genNewStruc

subroutine writeNewStruc(NAtoms, AN, coord, indexX, indexY, ds)
    implicit none
    integer(4), intent(in)  :: NAtoms, indexX, indexY
    integer(4), dimension(NAtoms), intent(in)   :: AN
    real(8), dimension(NAtoms, 3), intent(in)   :: coord
    real(8), intent(in) :: ds
    ! local variables
    integer(4)  :: i
    character(len=100)  :: buffer1, buffer2, buffer3, filename

    write(buffer1, *) indexX
    write(buffer2, *) indexY
    if ( ds .gt. 0) then 
        write(buffer3, '(F3.1)') ds
    else 
        write(buffer3, '(F4.1)') ds
    end if
    filename = TRIM(ADJUSTL(buffer1))//'_' &
        //TRIM(ADJUSTL(buffer2))//'_'//TRIM(ADJUSTL(buffer3)) &
        //'.xyz'
    open(10, file=filename, status='replace')
    write(10,'(I2)') NAtoms
    write(10,'(A)') TRIM(ADJUSTL(filename))
    do i = 1, NAtoms
        write(10,'(I2,3(1X,F9.6))') AN(i),coord(i,1:3)
    end do
    close(10)
    return 
end subroutine writeNewStruc