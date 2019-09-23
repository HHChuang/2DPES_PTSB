! $1 = traj
! $2 = rotM.ts.$1
! 2018/12/10, Grace

Program main
    implicit none
    integer(4)      :: NAtoms,npts,i,j
    integer(4),allocatable,dimension(:) :: AtomName
    real(8),allocatable,dimension(:,:)  :: Coord
    real(8),dimension(3,3)  :: rotM
    character(len=100)  :: traj

    ! Step 1. Check the input argument and read the files
        call get_NAtoms(NAtoms)
        call get_npts(NAtoms,npts)
        call get_rotM(rotM)
        ! call get_rotM(3,rotM3)
        allocate(AtomName(NAtoms))
        allocate(Coord(NAtoms,3))
    ! Step 2. Rotation
        call GETARG(1,traj)
        open(100,file='rot.'//TRIM(traj),status='replace')
        do i=1,npts
            call get_Coord(i,NAtoms,AtomName,Coord)
            call get_move_to_center(NAtoms,AtomName,Coord)
            call operateR(NAtoms,rotM,Coord)
            write(100,'(I2)') NAtoms
            write(100,'(I5)') i
            do j=1,NAtoms
                write(100,'(I2,3(1X,F14.10))') AtomName(j),Coord(j,1:3)
            end do
        end do
    stop
end program main
 
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

subroutine get_npts(NAtoms,npts)
    implicit none
    integer(4),intent(in)  :: NAtoms
    integer(4),intent(out)  :: npts
    integer(4)  :: totline
    character(len=100)  :: struc

    call GETARG(1,struc)
    call system('wc -l '//TRIM(adjustl(struc))//" | awk '{print $1}' > totline" )
    open(10,file='totline',status='old')
    read(10,*) totline
    npts = totline / (NAtoms +2)
    close(10,status='delete')
    return
end subroutine get_npts

subroutine get_rotM(rotM)
    implicit none
    real(8),dimension(3,3),intent(out)  :: rotM
    integer(4)  :: i
    character(len=100)  :: struc

    call GETARG(2,struc)
    open(10,file=TRIM(struc),status='old')
    do i=1,3
        read(10,*) rotM(i,1:3)
    end do
    close(10)
    return
end subroutine get_rotM

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

subroutine get_Coord(order,NAtoms,AtomName,Coord)
    implicit none
    integer(4),intent(in)   :: order,NAtoms
    integer(4),dimension(NAtoms),intent(out)    :: AtomName
    real(8),dimension(NAtoms,3),intent(out) :: Coord
    ! local variable
    character(len=100)  :: struc,ini_char,fin_char
    integer(4)  :: a,ini,fin

    call GETARG(1,struc)
    open(10,file=TRIM(ADJUSTL(struc)),status='old')
    ini = (NAtoms+2) * (order-1) + 1
    fin = (NAtoms+2) * order
    write(ini_char,*) ini
    write(fin_char,*) fin
    call system("sed -n '"//TRIM(ADJUSTL(ini_char)) &
        //','//TRIM(ADJUSTL(fin_char))//" p' " &
        //TRIM(ADJUSTL(struc))//' > coord' )
    open(10,file='coord',status='old')
    read(10,*) struc ! buffer
    read(10,*) struc ! buffer
    do a=1,NAtoms
        read(10,*) AtomName(a),Coord(a,1:3)
    end do
    close(10,status='delete')
    return
end subroutine get_Coord