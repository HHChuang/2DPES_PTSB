!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program:                                                                          !
!   Searching the path from near the shoulder (i.e. shift away from TS1) of a       !
!   bifurcation reaction to TS2                                                     !
!                                                                                   !
! Reference:                                                                        !
!                                                                                   !  
! Input:                                                                            !
!   $1 = *.xyz                                                                      !
!   $2 = F.dat                                                                      !
!   $3 = Eigenvector of TSS2; TS2EigV.dat                                           !
!                                                                                   !
! Output:                                                                           !
!   newCoord.xyz                                                                    !
!                                                                                   !
! History:                                                                          !
!   2018/08/11, Grace                                                               !
!   2019/04/20, Grace,                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program main
    implicit none
        character(len=100)  :: filename
        integer(4)  :: NAtoms
        integer(4),allocatable,dimension(:)    :: AtomName
        real(8),allocatable,dimension(:,:)  :: iniCoord,finCoord
        real(8),allocatable,dimension(:,:)  :: F,TS2EigV,ModifiedG
        ! from benchmark information ds = 0.001. All in MW unit.
        real(8),parameter   :: ds = 0.005D0 

        ! Step 1. Using variable 'NAtoms' to allocate memory
            call GETARG(1,filename)
            call getNatoms(filename,NAtoms)
            allocate(AtomName(NAtoms))
            allocate(iniCoord(NAtoms,3))
            allocate(finCoord(NAtoms,3))
            allocate(F(NAtoms,3))
            allocate(TS2EigV(NAtoms,3))
            allocate(ModifiedG(NAtoms,3))
        ! Step 2. Extract gradient and eigenvector and all other necessary info. 
        !         for writing input 
            call getCoord(filename,NAtoms,AtomName,iniCoord)
            call GETARG(2,filename) ! $2 = F.dat 
            call getVec(filename,NAtoms,F) 
            call GETARG(3,filename) ! $3 = TS2EigV.dat
            call getVec(filename,NAtoms,TS2EigV)
            call changeCoord(NAtoms,AtomName,TS2EigV)
            call normVec(NAtoms,3,TS2EigV)
        ! Step 3. Main part: generate the next geometry
            ! the region between TSS1 and TSS2
            call getModifiedG(NAtoms,F,TS2EigV,ModifiedG)
            call genGeom(NAtoms,ds,ModifiedG,AtomName,iniCoord)
end program main

subroutine getNAtoms(filename,NAtoms)
    implicit none
    character(len=100),intent(in)   :: filename
    integer(4),intent(inout)  :: NAtoms
    open(10,file=TRIM(ADJUSTL(filename)),status='old')
    read(10,*) NAtoms
    close(10)
    return
end subroutine getNAtoms

subroutine getCoord(filename,NAtoms,AtomName,Coord)
    implicit none
    character(len=100),intent(in)   :: filename
    integer(4),intent(in)   :: NAtoms
    integer(4),dimension(NAtoms),intent(out)    :: AtomName
    real(8),dimension(NAtoms,3),intent(out) :: Coord
    ! local variables
    integer(4)  :: i
    character(len=100)  :: buffer
    open(10,file=TRIM(ADJUSTL(filename)),status='old')
    read(10,*) buffer
    read(10,*) buffer
    do i=1,NAtoms
        read(10,*) AtomName(i),Coord(i,1:3)
    end do
    close(10)
    return
end subroutine getCoord

subroutine getVec(filename,NAtoms,Vec)
    implicit none
    character(len=100),intent(in)   :: filename
    integer(4),intent(in)   :: NAtoms
    real(8),dimension(NAtoms,3),intent(out)   :: Vec
    ! local variable
    integer(4)  :: i,j
    open(10,file=TRIM(ADJUSTL(filename)),status='old')
    do i=1,NAtoms
        read(10,*) Vec(i,1:3)
    end do
    close(10)
    return
end subroutine getVec

subroutine normVec(dimen1,dimen2,Vec)
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
end subroutine normVec

subroutine changeCoord(NAtoms,AtomName,TS2EigV)
    implicit none
    integer(4),intent(in)   :: NAtoms
    integer(4),dimension(NAtoms),intent(in) :: AtomName
    real(8),dimension(NAtoms,3),intent(inout)   :: TS2EigV
    !local variable 
    integer(4)  :: i,j
    real(8),dimension(NAtoms)   :: AW
    do i=1,NAtoms
        call ANtoAW(AtomName(i),AW(i))
    end do
    do j=1,3 ! column-major
        do i=1,NAtoms
            TS2EigV(i,j)=TS2EigV(i,j)*SQRT(AW(i))
        end do 
    end do
end subroutine changeCoord

subroutine ANtoAW(AN,AW)
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
        write(*,'(A)') 'Stop the program, GenArticStruc.f90.'
        Stop
    end select 
end subroutine ANtoAW

subroutine getModifiedG(NAtoms,F,eigenVec,ModifiedG)
    implicit none
    integer(4),intent(in)   :: NAtoms
    real(8),dimension(NAtoms,3),intent(in)   :: F,eigenVec
    real(8),dimension(NAtoms,3),intent(out) :: ModifiedG
    ! local variable
    integer(4)  :: i,j,index
    real(8),dimension(3*NAtoms) :: arrayG,arrayEigV,arrayModifiedG
    arrayG = 0.0D0 
    arrayEigV = 0.0D0 
    arrayModifiedG = 0.0D0
    index = 1
    do i = 1,NAtoms
        do j = 1,3
            arrayG(index) = -F(i,j)
            arrayEigV(index) = eigenVec(i,j)
            index = index + 1
        end do
    end do
    ! Main part 
    ! wiki: https://en.wikipedia.org/wiki/Vector_projection
    do i = 1,3*NAtoms
        arrayModifiedG(i) = arrayG(i) &
            - 2*(DOT_PRODUCT(arrayG,arrayEigV)/DOT_PRODUCT(arrayEigV,arrayEigV)) &
            *arrayEigV(i)
    end do
    ! fold 1D array to 2D matrix
    index = 1
    do i=1,NAtoms
        do j = 1,3
            ModifiedG(i,j) = arrayModifiedG(index)
            index = index + 1
        end do
    end do
    return
end subroutine getModifiedG

subroutine genGeom(NAtoms,ds,ModifiedG,AtomName,iniCoord)
    implicit none
    integer(4),intent(in)   :: NAtoms
    integer(4),dimension(NAtoms),intent(in)   :: AtomName
    real(8),intent(in)  :: ds
    real(8),dimension(NAtoms,3),intent(in)  :: ModifiedG,iniCoord
    ! local variable
    integer(4)  :: i,j
    real(8),dimension(NAtoms,3) :: finCoord
    finCoord = 0.0D0
    do i = 1,NAtoms
        do j = 1,3
            finCoord(i,j) = iniCoord(i,j) - ModifiedG(i,j) * ds
        end do
    end do
    open(10,file='newCoord.xyz',status='replace')
    write(10,'(I2)') NAtoms
    write(10,'(A)') 'newCoord generate by artifical rxn. coord'
    do i=1,NAtoms
        write(10,'(I2,3(1X,F10.6))') AtomName(i),finCoord(i,1:3)
    end do
    close(10)
    return
end subroutine genGeom