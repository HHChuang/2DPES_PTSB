!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program: Assume the energy on this PES is the fingure-print of the corresponding structure; i.e., each point on this PES has unique energy. And compare the RMSD of the whole trajectories between itself and this PES.
!
! Input: 
!   1. rot.Traj#
!   2. coord.rot.Traj#
!   3. #_Struc.*.xyz
!
! Output: 
!   RMSD.Traj# TODO:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Program main
    implicit none
    integer(4)  :: NAtoms, nTrajPts
    real(8),allocatable,dimension(:,:)  :: coord_Traj, coord_PES
    integer(4)  :: i
    real(8) :: rmsd

    ! 1. 
    call get_NAtoms(NAtoms)
    allocate(coord_Traj(NAtoms,3))
    allocate(coord_PES(NAtoms,3))
    call get_nTrajPts(NAtoms,nTrajPts)

    open(100,file='RMSD.dat',status='replace')
    do i = 1, nTrajPts
    ! 2. Get the structure of trajectories 
        call get_coordTraj(i,NAtoms,coord_Traj)
    ! 3. Get the structure of 2D-PES
        call get_coordPES(i,NAtoms,coord_PES)
    ! 4. Calculate RMSD
        call get_RMSD(NAtoms,coord_Traj,coord_PES,rmsd)
        write(*,*) rmsd
        write(100,*) rmsd
    end do
    close(100)
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

subroutine get_nTrajPts(NAtoms,nTrajPts)
    implicit none
    integer(4),intent(in)   :: NAtoms
    integer(4),intent(out)   :: nTrajPts
    ! local variable
    integer(4)  :: tmp
    character(len=100)  :: filename

    call GETARG(1,filename)
    call system('wc -l '//TRIM(filename)//" | awk '{print $1}' > totalLine.dat")
    open(10,file='totalLine.dat',status='old')
    read(10,*)  tmp
    nTrajPts = tmp/(NAtoms+2)
    close(10,status='delete')
    return
end subroutine get_nTrajPts

subroutine get_coordTraj(order,NAtoms,coord)
    implicit none
    integer(4),intent(in)   :: order, NAtoms
    real(8),dimension(NAtoms,3),intent(out) :: coord
    ! local variable 
    integer(4)  :: fileL,iniL,finL
    character(len=100)   :: filename,ini_char, fin_char

    call GETARG(1, filename)
    fileL = NAtoms + 2
    iniL = 1 + (order - 1) * fileL 
    finL = order * fileL
    write(ini_char,*) iniL
    write(fin_char,*) finL
    call system('sed -n "'//TRIM(ADJUSTL(ini_char))//','//TRIM(ADJUSTL(fin_char))//' p" '//TRIM(ADJUSTL(filename))//' > tmp.xyz')
    filename='tmp.xyz'
    call readStruc(filename,NAtoms,coord)
    return 
end subroutine get_coordTraj

subroutine readStruc(filename,NAtoms,coord)
    implicit none
    character(len=100),intent(in)   :: filename
    integer(4),intent(in)   :: NAtoms
    real(8),dimension(NAtoms,3),intent(out) :: coord
    ! locat variable 
    integer(4)  :: i
    character(len=100)  :: buffer

    open(10,file=TRIM(ADJUSTL(filename)),status='old',action='read')
    read(10,*) buffer
    read(10,*) buffer
    do i=1,NAtoms
        read(10,*) buffer,coord(i,1:3)
    end do
    close(10,status='delete')
    return
end subroutine readStruc

subroutine get_coordPES(order,NAtoms,coord)
    implicit none 
    integer(4),intent(in)   :: order, NAtoms
    real(8),dimension(NAtoms,3),intent(out) :: coord
    ! local variable 
    character(len=100)  :: tmpL,energy,filename

    ! pick the position of corresponding structure 
    call GETARG(2, filename)
    write(tmpL,*) order
    call system('sed -n "'//TRIM(ADJUSTL(tmpL))//','//TRIM(ADJUSTL(tmpL))//' p" '//TRIM(ADJUSTL(filename))//" | awk '{print $3}' | sed 's/-//g' > tmp")
    open(10,file='tmp',status='old')
    read(10,*) filename 
    energy=filename(1:10)
    ! write(*,*) tmpL,energy
    close(10,status='delete')

    ! extract the structure on this PES
    call GETARG(3, filename)
    write(tmpL,*) NAtoms
    call system('grep -B1 -A'//TRIM(ADJUSTL(tmpL))//' "'//TRIM(ADJUSTL(energy))//'" '//TRIM(ADJUSTL(filename))//' > tmp.xyz ')
    filename='tmp.xyz'
    call readStruc(filename,NAtoms,coord)
    return 
end subroutine get_coordPES

subroutine get_RMSD(NAtoms,Coord1,Coord2,rmsd)
    implicit none
    integer(4), intent(in)   :: NAtoms
    real(8), dimension(NAtoms,3), intent(in)  :: Coord1,Coord2
    real(8), intent(out) :: rmsd
    ! local variable
    real(8) :: sumDiff
    integer(4)  :: i,j
    sumDiff = 0.0D0
    do j = 1, 3 
        do i = 1, Natoms 
            sumDiff = sumDiff + ( coord1(i,j) - coord2(i,j) )**2
        end do
    end do
    rmsd = SQRT(sumDiff/(NAtoms*1.0D0)) ! NAtoms must larger then 0
    return
end subroutine get_RMSD