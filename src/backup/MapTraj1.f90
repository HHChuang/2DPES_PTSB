!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Program : Mapping the trajectories to numerical PES by RMSD         !
!                                                                       !
!   Pre-requisted shell script  :                                       !
!       1. getPESwStruc.sh                                              !
!               format:                                                 !
!                       # of atom                                       !
!                       coord_x coord_y E                               !
!                       #atom   x   y   z                               !
!                                                                       !
!   Input :                                                             !
!           $1 = trajectory                                             !
!               format:                                                 !
!                       # of atom                                       !
!                       comment                                         !
!                       #atom       x   y   z                           !
!           $2 = orderList.dat                                          !
!                and a losts *.xyz; all the structures on this 2D-PES   !
!                                                                       !
!                                                                       !
!   Output :                                                            !
!           coord.$1                                                    !
!                                                                       !
!   History:                                                            !
! 2018/10/12, Grace                                                     !
! 2019/06/04, Grace, modify the structure of code. Map trajectory from  !
!   the transition-state structure.                                     !
! 2019/06/10, Grace, Turn the coordinate of trajectory; add function    !
!   tuneTraj(). Change the coordinate from the amount of grid points to !
!   mass-weighted coordinate; add function gp2mw().                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Program main
    implicit none
    integer(4)  :: NAtoms,order_pes_ts,num_pes
    integer(4),allocatable,dimension(:,:)   :: orderList
    real(8),allocatable,dimension(:)    :: energyList
    character(len=100)  :: filename

    ! Step 1. I/O import trajectory and numerical PES
        call print_purpose()
        call get_NAtoms(NAtoms)
        call get_numPts(2,NAtoms,num_pes)
        allocate(orderList(num_pes,2))
        allocate(energyList(num_pes))
        call get_orderList(num_pes,orderList,energyList)
        

    ! Step 2. Find the initial coordinate; compare the rmsd by bubble sorting concept
    !   most time consuming step 
        call GETARG(1,filename)
        call find1stPts(NAtoms,num_pes,orderList,order_pes_ts) ! output: order_pes_ts
    ! Step 3. Compare the RMSD to its neighborhood, and find the next point
        open(100,file='coord.'//TRIM(ADJUSTL(filename)),status='replace')
        ! Reverse path of trajectory; second argument = .true.
        call halfPts(100,.true.,NAtoms,order_pes_ts,num_pes,orderList,energyList)
        ! tss of trajectory
        write(100,'(2(I4,1X),F12.7)') orderList(order_pes_ts,1:2),energyList(order_pes_ts)
        ! Forward path of trajectory; second argument = .false.
        call halfPts(100,.false.,NAtoms,order_pes_ts,num_pes,orderList,energyList)
        close(100)
    stop
end program main

subroutine print_purpose()
    implicit none
    character(len=100)   :: filename

    call GETARG(1,filename)

    write(*,'()') 
    write(*,'(A)') '---------------------------------------------------------------------'
    write(*,'(A)') '|  Purpose:                                                         |'
    write(*,'(A)') '|   Mapping the trajectories to numerical PES by RMSD               |'
    write(*,'(A)') '|                                                                   |'
    write(*,'(A)') '|  Limitation:                                                      |'
    write(*,'(A)') '|   1. the potential should be square; amount of point from x = y   |'
    write(*,'(A)') '|   2. atomic number in trajectory must be integer                  |'
    write(*,'(A)') '|                                                                   |'
    write(*,'(A)') '|  Output file:                                                     |'
    write(*,'(A)') '|   coord.'//TRIM(ADJUSTL(filename))
    write(*,'(A)') '---------------------------------------------------------------------'
    write(*,'()') 
    return
end subroutine print_purpose

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

subroutine get_numPts(fid,NAtoms,numPts)
    implicit none
    integer(4),intent(in)   :: fid,NAtoms
    integer(4),intent(out)   :: numPts
    ! local variable
    integer(4)  :: totalLine
    character(len=100)  :: struc
    call GETARG(fid,struc)
    call system('wc -l '//TRIM(struc)//" | awk '{print $1}' > totalLine.dat")
    open(10,file='totalLine.dat',status='old')
    read(10,*)  totalLine
    close(10,status='delete')
    if (fid == 1) then 
        numPts=totalLine/(NAtoms+2)
    else 
        numPts=totalLine
    end if 
    
    return
end subroutine get_numPts

subroutine get_orderList(num_pes,orderList,energyList)
    implicit none
    integer(4),intent(in)   :: num_pes
    integer(4),dimension(num_pes,2),intent(out) :: orderList
    real(8),dimension(num_pes),intent(inout)    :: energyList
    ! local variable
    character(len=100)   :: filename
    logical :: filestat
    integer(4)  :: i

    call GETARG(2,filename)
    INQUIRE(file=filename,exist=filestat)
    open(10,file=filename,status='old',action='read')
    do i=1,num_pes
        read(10,*) orderList(i,1:2),energyList(i)
    end do
    close(10)
    return
end subroutine

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
    close(10)
return
end subroutine readStruc

subroutine get_PESfileName(order_pes,num_pes,orderList,filename)
    implicit none
    integer(4),intent(in)   :: order_pes,num_pes
    integer(4),dimension(num_pes,2),intent(in)  :: orderList
    character(len=100),intent(inout)    :: filename
    ! local variable
    character(len=50)   :: x_coord,y_coord

    write(x_coord,*) orderList(order_pes,1)
    write(y_coord,*) orderList(order_pes,2)
        filename=TRIM(ADJUSTL(x_coord))//'_' &
            //TRIM(ADJUSTL(y_coord))//'.xyz'
    return
end subroutine get_PESfileName

subroutine get_RMSD(NAtoms,Coord1,Coord2,rmsd)
    implicit none
    integer(4),intent(in)   :: NAtoms
    real(8),dimension(NAtoms,3),intent(in)  :: Coord1,Coord2
    real(8),intent(out) :: rmsd
    ! local variable
    real(8) ::sumDiff
    integer(4)  :: i,j
    sumDiff=0.0D0
    do i=1,NAtoms
        do j=1,3
            sumDiff = sumDiff + (coord1(i,j)-coord2(i,j))**2
        end do
    end do
    rmsd = SQRT(sumDiff/(NAtoms*1.0D0)) !NAtoms must larger then 0
    return
end subroutine get_RMSD

subroutine find1stPts(NAtoms,num_pes,orderList,order_pes)
    implicit none
    integer(4),intent(in)   :: NAtoms,num_pes
    integer(4),dimension(num_pes,3),intent(in)  :: orderList
    integer(4),intent(inout)    :: order_pes
    ! local variable
    character(len=100)   :: filename,buffer
    integer(4)  :: i
    real(8),dimension(NAtoms,3) :: coord_traj,coord_pes
    character(len=100)  :: name_Coord_pes
    real(8) :: rmsd1,rmsd2

    call GETARG(1,filename)
    buffer='tss.'//TRIM(ADJUSTL(filename))//'.xyz'
    call readStruc(buffer,NAtoms,coord_traj)
    call get_PESfileName(1,num_pes,orderList,name_Coord_pes)
    call readStruc(name_Coord_pes,NAtoms,Coord_pes)
    call get_RMSD(NAtoms,coord_traj,coord_pes,rmsd1)
    do i=2,num_pes
        call get_PESfileName(i,num_pes,orderList,name_Coord_pes)
        call readStruc(name_Coord_pes,NAtoms,coord_pes)
        call get_RMSD(NAtoms,coord_traj,coord_pes,rmsd2)
        if ( rmsd1 .gt. rmsd2 ) then
            rmsd1 = rmsd2
            order_pes=i
        end if
    end do
    return
end subroutine find1stPts

subroutine halfPts(fileID,direction,NAtoms,order_ts,num_pes,orderList,energyList)
    implicit none
    logical,intent(in)  :: direction
    integer(4),intent(in)   :: fileID,NAtoms,num_pes
    integer(4),intent(in)    :: order_ts
    integer(4),dimension(num_pes,2),intent(in)  :: orderList
    real(8),dimension(num_pes),intent(in)   :: energyList
    ! local variable
    integer(4)  :: order_center
    integer(4)  :: njobs,i,order_neighbor,write_ini,write_fin,write_int
    character(len=100)  :: jobpre,jobnum,jobname
    real(8),dimension(NAtoms,3) :: coord_traj
    integer(4),allocatable,dimension(:)   :: coordList
    real(8), allocatable,dimension(:,:)   :: tuneCoord

    if (direction) then
        ! reverse direction
        call system('ls | grep reverse | grep -c xyz > njobs.tmp')
        open(10,file='njobs.tmp',status='old',action='read')
        read(10,*) njobs
        close(10,status='delete')
        jobpre='reverse.'
        write_ini=njobs
        write_fin=1
        write_int=-1
    else
        ! forward direction
        call system('ls | grep forward | grep -c xyz > njobs.tmp')
        open(10,file='njobs.tmp',status='old',action='read')
        read(10,*) njobs
        close(10,status='delete')
        jobpre='forward.'
        write_ini=1
        write_fin=njobs
        write_int=1
    end if

    allocate(coordList(njobs)) 
    allocate(tuneCoord(njobs,2))

    ! Searching the coordinate from its neighbors
    order_center = order_ts
    do i = 1,njobs
        write(jobnum,*) i
        jobname=TRIM(ADJUSTL(jobpre))//TRIM(ADJUSTL(jobnum))//'.xyz'
        call readStruc(jobname, NAtoms, coord_traj)
        call findNextPts(NAtoms, coord_traj, order_center, num_pes, orderList, energyList, order_neighbor)
        call tuneTraj(NAtoms, num_pes, orderList, order_neighbor, coord_traj, tuneCoord(i,1:2))
        order_center = order_neighbor
        coordList(i) = order_neighbor
    end do

    ! test the boundary conditions
    ! do i=1,num_pes
    !     call findNextPts(NAtoms, coord_traj, i, num_pes, orderList, energyList, order_neighbor)
    ! end do
    ! stop

    ! write(*,*) write_ini,write_fin,write_int
    ! import those coordinates 
    do i = write_ini,write_fin,write_int
        ! original coordinates
        ! write(fileID,'(2(I4,1X),F12.7)') orderList(coordList(i),1:2),energyList(coordList(i))

        ! shift coordinate
        write(fileID,'(3(F12.7,1X))') tuneCoord(i,1:2), energyList(coordList(i))
        ! stop
    end do
    return
end subroutine halfPts


subroutine findNextPts(NAtoms,coord_traj,order_center,num_pes,orderList,energyList,order_neighbor)
    implicit none
    integer(4),intent(in)   :: NAtoms,num_pes,order_center
    real(8),dimension(NAtoms,3),intent(in)  :: coord_traj
    integer(4),dimension(num_pes,2),intent(in) :: orderList
    real(8),dimension(num_pes),intent(in)   :: energyList
    integer(4),intent(out)   :: order_neighbor
    ! local variable 
    integer(4)  :: i,j
    integer(4)  :: num_1D,num_column
    real(8) :: rmsd_small,rmsd
    character(len=100)  :: name_coord_pes
    real(8),dimension(NAtoms,3) :: coord_pes

    num_1D = INT(SQRT(num_pes*1.0D0))
    order_neighbor=0
    
    rmsd_small = 100000000.D0 ! vary big rmsd to be replaced
    do i = order_center - num_1D, order_center + num_1D, num_1D
        if ( ( i .lt. 1 ) .or. ( i .gt. num_pes ) ) cycle
        if ( MODULO(i,num_1D) .eq. 0 ) then
            num_column = i/num_1D
        else
            num_column = i/num_1D + 1
        end if
        do j = i - 1, i +1
            ! avoid the original point
                if ( j .eq. order_center ) cycle 
            ! set the boudary of one column in the OrderList
                if (     ( j .lt. (num_column-1) * num_1D+1 ) &
                    .or. ( j .gt. num_column*num_1D       )     ) cycle 
            ! call get_coord_pes(2,j,NAtoms,x,y,E,AtomName,Coord_pes)  
                call get_PESfileName(j, num_pes, orderList, name_coord_pes)
                call readStruc(name_coord_pes, NAtoms, coord_pes)

            ! write(*,*) j, orderList(j,1:2)
            ! write(*,'(2(I3,1X))') orderList(j,1:2)

            call get_RMSD(NAtoms, coord_traj, coord_pes, rmsd)
            if ( rmsd_small .gt. rmsd ) then
                rmsd_small = rmsd
                order_neighbor = j
            end if

        end do
     end do
    return
end subroutine findNextPts

! 2019/06/10, Grace, shift the original coordinate of trajectory
subroutine get_coord_pes(NAtoms, x, y, coord_pes) 
    implicit none
    integer(4), intent(in)   :: NAtoms, x,y
    real(8), dimension(NAtoms,3), intent(out)   :: coord_pes
    ! local variables
    character(len=100)  :: coord_x, coord_y, filename

    write(coord_x,*) x
    write(coord_y,*) y
    filename=TRIM(ADJUSTL(coord_x))// &
        '_'//TRIM(ADJUSTL(coord_y))//'.xyz'
    call readStruc(filename,NAtoms,coord_pes)
end subroutine get_coord_pes

subroutine get_vec(NAtoms, num_pes, orderList, order_nearby, coord_original, vec)
    implicit none
    integer(4), intent(in)  :: NAtoms, num_pes, order_nearby
    integer(4),dimension(num_pes,2),intent(in)  :: orderList
    real(8), dimension(NAtoms, 3), intent(in)   :: coord_original
    real(8), dimension(NAtoms, 3), intent(inout)  :: vec
    ! local variable
    integer(4)  :: i, j
    real(8), dimension(NAtoms, 3)   :: coord_pes
    real(8) :: norm

    call get_coord_pes(NAtoms, orderList(order_nearby, 1), orderList(order_nearby, 2), coord_pes)
    ! Fortran is row-major
     norm = 0.0D0
    do i = 1, NAtoms
        do j = 1, 3
            vec(i, j) = coord_pes(i,j) - coord_original(i,j)
            ! Normalize 
            norm = norm + vec(i,j) * vec(i,j)
        end do
    end do
    norm = SQRT(norm) 
    ! Helf length FIXME: may have bug
    do i = 1, NAtoms
        do j = 1, 3 
            ! vec(i, j) = vec(i, j)*0.5D0 / norm
            vec(i, j) = vec(i, j) / norm
        end do
    end do
    return
end subroutine get_vec

subroutine get_nTimes(NAtoms, vec_pes, vec_traj, nTimes)
    implicit none
    integer(4), intent(in)  :: NAtoms
    real(8), dimension(NAtoms, 3), intent(in)   :: vec_pes, vec_traj
    real(8), intent(inout)    :: nTimes
    ! local variables
    integer(4)  :: i, j

    nTimes = 0
    do i = 1, NAtoms
        do j = 1, 3
            nTimes = vec_pes(i, j) * vec_traj(i, j)
        end do
    end do

    return
end subroutine get_nTimes

subroutine tuneTraj(NAtoms,num_pes,orderList,order_2DPES,coord_traj, tuneCoord)
    implicit none
    integer(4),intent(in)   :: NAtoms, num_pes
    integer(4),dimension(num_pes, 2),intent(in)  :: orderList
    integer(4),intent(in)   :: order_2DPES
    real(8), dimension(NAtoms, 3), intent(in)    :: coord_traj
    real(8), dimension(2), intent(inout)  :: tuneCoord
    ! local variable
    integer(4)  :: i, j 
    integer(4)  :: num_1D, order_up, order_below, order_right, order_left
    real(8), dimension(4)   :: nTimes ! nTimes=(n_up, n_below, n_right, n_left)
    real(8),dimension(NAtoms,3) :: coord_original, coord_pes
    real(8),dimension(NAtoms,3) :: vec_up, vec_below, vec_right, vec_left, vec_traj

    num_1D = INT(SQRT(num_pes*1.0D0))

    ! Step 1. Extract the indeice of 4 structures from 
    !   this 2D-PES which are arounds the original selected grid point. 
        order_up = order_2DPES + 1
        order_below = order_2DPES - 1
        order_right = order_2DPES + num_1D
        order_left = order_2DPES - num_1D
    ! Step 2. Calculate these vectors; the structure difference between
    !   one of the four structures and the original selected structure.
        call get_coord_pes(NAtoms, orderList(order_2DPES, 1), orderList(order_2DPES, 2), coord_original)
        call get_vec(NAtoms, num_pes, orderList, order_up, coord_original, vec_up)
        call get_vec(NAtoms, num_pes, orderList, order_below, coord_original, vec_below)
        call get_vec(NAtoms, num_pes, orderList, order_right, coord_original, vec_right)
        call get_vec(NAtoms, num_pes, orderList, order_left, coord_original, vec_left)
        
    ! Step 3. Project the coordinate of trajectory on this tiny grid.
        ! 1. Structure difference between the trajectory and the original
        !       grid point on this 2D-PES.
        do i = 1, NAtoms
            do j = 1, 3
                vec_traj(i, j) = coord_original(i, j) - coord_traj(i, j)
            end do
        end do
        ! 2. vector analysis
        call get_nTimes(NAtoms, vec_up, vec_traj, nTimes(1))
        call get_nTimes(NAtoms, vec_below, vec_traj, nTimes(2))
        call get_nTimes(NAtoms, vec_right, vec_traj, nTimes(3))
        call get_nTimes(NAtoms, vec_left, vec_traj, nTimes(4))

        tuneCoord(1) = orderList(order_2DPES,1) * 1.0D0  &
            + nTimes(3) * 0.5D0  &
            - nTimes(4) * 0.5D0
        tuneCoord(2) = orderList(order_2DPES,2) * 1.0D0  &
            + nTimes(1) * 0.5D0  &
            - nTimes(2) * 0.5D0
        ! write(*,*) nTimes(1:4)
    ! Step 4. Check the boundary condition. TODO: 
    
    return
end subroutine tuneTraj

! 2019/06/10, Grace, change the coordinate from the amount of grid point (integer)
!   to mass-weighted coordinate (real)
! subroutine gp2mw()
!     implicit none

!     return
! end subroutine gp2mw