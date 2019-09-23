! Debug for MoveStruc.f90 
module constants  
    implicit none 
    real(8),parameter   :: pi = 3.1415926536D0
    real(8),parameter   :: zero=0.0001D0
end module constants 

Program main
    use constants
    implicit none
    integer(4)      :: i,NAtoms
    integer(4),allocatable,dimension(:) :: AtomName
    real(8),allocatable,dimension(:,:)  :: Coord
    real(8) :: angle=-4.0D0*pi/7.0D0 !rad
    real(8),dimension(3)    :: n=(/ 1.0D0,0.0D0,0.0D0 /)

    ! Step 1. Check the input argument and read the files
        call get_NAtoms(NAtoms)
        allocate(AtomName(NAtoms))
        allocate(Coord(NAtoms,3))
        call get_coord(1,NAtoms,AtomName,Coord)
    ! Step 2. Translation: move two structures to the center of mass
        call get_move_to_center(NAtoms,AtomName,Coord)
    ! Step 3. Rotation
        call do_rotation(NAtoms,AtomName,angle,n,Coord) 
        call print_moveCoord(NAtoms,AtomNAme,Coord)
        call print_coord(NAtoms,AtomName,Coord)
    stop
end program main
 
!FIXME: written for debug
subroutine print_coord(NAtoms,AtomName,Coord)
    implicit none
    integer(4),intent(in)   :: NAtoms
    integer(4),dimension(NAtoms)    :: AtomName
    real(8),dimension(NAtoms,3) :: Coord
    integer(4)  :: i
    write(*,'()')
    do i=1,NAtoms
        write(*,'(I2,3(1X,F14.10))') AtomName(i),Coord(i,1:3)
    end do
    write(*,'()')
    return
end subroutine print_coord

subroutine print_moveCoord(NAtoms,AtomNAme,Coord)
    implicit none
    integer(4),intent(in)   :: NAtoms
    integer(4),dimension(NAtoms)    :: AtomName
    real(8),dimension(NAtoms,3),intent(in)  :: Coord
    ! local variable
    integer(4)  :: i
    character(len=100)  :: struc
    call GETARG(1,struc)
    open(10,file='move.'//TRIM(struc),status='replace')
    write(10,'(I2)') NAtoms
    write(10,'(A)') 'adjustable coordinate'
    do i=1,NAtoms
        write(10,'(I2,3(1X,F14.10))') AtomName(i),Coord(i,1:3)
    end do
    return
end subroutine print_moveCoord

subroutine get_NAtoms(NAtoms)
    implicit none
    integer(4),intent(out)  :: NAtoms
    character(len=100)  :: struc
    logical :: filestat
    call GETARG(1,struc)
    INQUIRE(file=struc,exist=filestat)
    if (filestat) then
        open(10,file=struc,status='old')
    else
        write(*,'(A)') TRIM(struc)//" doesn't exist"
        stop
    end if
    read(10,*) NAtoms
    close(10)
    return
end subroutine get_NAtoms

subroutine get_coord(fid,NAtoms,AtomName,Coord)
    implicit none
    integer(4),intent(in)   :: fid,NAtoms
    integer(4),dimension(NAtoms),intent(out)    :: AtomName
    real(8),dimension(NAtoms,3),intent(out) :: Coord
    character(len=100)  :: struc
    logical :: filestat
    integer(4)  :: i
    call GETARG(fid,struc)
    INQUIRE(file=struc,exist=filestat)
    if (filestat)   then
        open(10,file=struc,status='old',action='read')
    else
        write(*,'(A)') TRIM(struc)//" doesn't exist"
        stop
    end if
    read(10,*) struc ! buffer
    read(10,*) struc ! buffer
    do i=1,NAtoms
        read(10,*) AtomName(i),Coord(i,1:3)
    end do
    close(10)
    return
end subroutine get_coord

subroutine get_move_to_center(NAtoms,AtomName,Coord)
    implicit none
    integer(4),intent(in)   :: NAtoms
    integer(4),dimension(NAtoms),intent(in) :: AtomName
    real(8),dimension(NAtoms,3),intent(inout)  :: Coord
    real(8),dimension(3)    :: com
    integer(4)  :: i,j,totalmass
    ! Calculate center of mass (only atomic number; # of protons, are extracted)
    ! Since we only have the information of proton, there is no need to multiply
    ! a constant to transform the value to a.u.
    totalmass=0
    do i=1,NAtoms
        totalmass=totalmass+AtomName(i)
    end do
    com(:)=0.0D0
    do i=1,3
        do j=1,NAtoms
            com(i)=com(i)+AtomName(j)*Coord(j,i)
        end do
        com(i)=com(i)/totalmass
    end do
    do i=1,NAtoms
        do j=1,3
            Coord(i,j)=Coord(i,j)-com(j)
        end do
    end do
    return
end subroutine get_move_to_center

subroutine do_rotation(NAtoms,AtomName,angle,n,Coord)
    use constants
    implicit none
    integer(4),intent(in)   :: NAtoms
    integer(4),dimension(NAtoms),intent(in) :: AtomName
    real(8),dimension(NAtoms,3),intent(inout)  :: Coord
    real(8),intent(in)  :: angle
    real(8),dimension(3),intent(in)  :: n
    ! local variable
    real(8),dimension(3,3)  ::rotM
        call get_rotM(n,angle,rotM)
        call operateR(NAtoms,rotM,Coord) 
    return
end subroutine do_rotation

subroutine operateR(NAtoms,rotM,coord)
    implicit none
    integer(4),intent(in)   :: NAtoms
    real(8),dimension(3,3)  :: rotM
    real(8),dimension(NAtoms,3),intent(inout)    :: coord
    ! local variable
    real(8),dimension(3,NAtoms) :: coord_t,tmp
    integer(4)  :: i,j,k
    ! transpose 
    do i=1,NAtoms
        do j=1,3
            coord_t(j,i)=coord(i,j)
        end do
    end do
    tmp(:,:)=0.0D0
    do i=1,3
        do j=1,3
            do k=1,NAtoms
                tmp(i,k)=tmp(i,k)+rotM(i,j)*coord_t(j,k)
            end do
        end do
    end do
    ! transpose
    do i=1,NAtoms
        do j=1,3
            coord(i,j)=tmp(j,i)
        end do
    end do
    return
end subroutine operateR      

subroutine get_rotM(rotAxis,angle,RotMatrix)
    use constants
    implicit none
    real(8),dimension(3),intent(in) :: rotAxis
    real(8),intent(in) :: angle
    real(8),dimension(3,3),intent(out)  :: RotMatrix
    ! local variable
    integer(4)  :: i
    real(8) :: length
    real(8),dimension(4)    :: q
    ! calculate quaternion
    q(1)=COS(angle/2)
    do i=1,3
        q(i+1)=SIN(angle/2)*rotAxis(i)
    end do
    ! normalize q to get unit quaternion
    length =0.0D0
    do i=1,4
        length = length + q(i)**2 
    end do
    if ( length .gt. zero ) then
        do i=1,4
            q(i) = q(i)/length
        end do
    end if
    write(*,*) q(:)
    ! create rotation matrix
                    ! (1,1)
    RotMatrix=RESHAPE(SOURCE=(/ 1-2*( q(3)**2 + q(4)**2 ),&
                    ! (1,2)
                    2*( q(2)*q(3) - q(4)*q(1) ), &
                    ! (1,3)
                    2*( q(2)*q(4) + q(3)*q(1) ), &
                    1-2*( q(2)**2 + q(4)**2 ),    &
                    ! (2,3)
                    2*( q(3)*q(4) - q(2)*q(1) ), &
                    ! (2,1)
                    2*( q(2)*q(3) + q(4)*q(1) ), &
                    ! (2,2)
                    ! (3,1)
                    2*( q(2)*q(4) - q(3)*q(1) ), &
                    ! (3,2)
                    2*( q(3)*q(4) + q(2)*q(1) ), &
                    ! (3,3)
                    1-2*( q(2)**2 + q(3)**2 ) /), &
                SHAPE=(/ 3,3 /) ) 
    ! plot the rotational matrix
    write(*,'()')
    do i=1,3
        write(*,'(3(F5.2,1X))') RotMAtrix(1,i),RotMAtrix(2,i),RotMAtrix(3,i)
    end do
    write(*,'()')
    return
end subroutine get_rotM
