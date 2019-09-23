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
!                       #atom       x   y   z                               !
!           $2 = numerical PES                                          !
!               format:                                                 !
!                       # of atom                                       !
!                       coord_x coord_y E                               !
!                       #atom   x   y   z                               !
!                                                                       !
!   Output :                                                            !
!           coord.$1                                                    !
!                                                                       !
!   History:                                                            !
! 2018/10/12, Grace                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Program main
    implicit none
    integer(4)  :: i,NAtoms
    integer(4)  :: order_pes,tra,num_tra,num_pes,num_1D
    character(len=100)  :: struc
    integer(4),allocatable,dimension(:,:)   :: orderList
    real(8),allocatable,dimension(:)    :: energyList

    ! Step 1. I/O import trajectory and numerical PES
        call get_NAtoms(NAtoms)
        call get_numPts(1,NAtoms,num_tra)
        call get_numPts(2,NAtoms,num_pes)
        allocate(orderList(num_pes,3))
        allocate(energyList(num_pes))
        num_1D = INT(SQRT(num_pes*1.0D0))
        call GETARG(1,struc)
        call print_purpose(struc)
        open(100,file='coord.'//TRIM(struc),status='replace')
    ! Step 2. Find the initial coordinate; compare the rmsd by bubble sorting concept
    !   most time consuming step 
        call get_orderList(NAtoms,num_pes,orderList,energyList)
        stop
        call get_iniPts(NAtoms,num_pes,order_pes) !only searching the 10% region of 2D-PES
        write(100,'(2(I4,1X),F12.7)') orderList(order_pes,2:3),energyList(order_pes)
    ! Step 3. Compare the RMSD to its neighborhood, and find the next point
        ! do i=2,num_tra 
        !     call get_nextPts(NAtoms,num_pes,i,num_1D,order_pes)
        !     write(100,'(2(I4,1X),F12.7)') orderList(order_pes,2:3),energyList(order_pes)
        ! end do
        close(100)
    stop
end program main

subroutine print_purpose(file)
    implicit none
    character(len=100),intent(in)   :: file
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
    write(*,'(A)') '|   coord.'//TRIM(ADJUSTL(file))
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
    call system('nl '//TRIM(struc)//" | tail -n -1 | awk '{print $1}' > totalLine.dat")
    open(10,file='totalLine.dat',status='old')
    read(10,*)  totalLine
    close(10,status='delete')
    numPts=totalLine/(NAtoms+2)
    return
end subroutine get_numPts

subroutine get_coord_tra(fid,order,NAtoms,AtomName,Coord)
    implicit none
    integer(4),intent(in)   :: fid,order,NAtoms
    integer(4),dimension(NAtoms),intent(inout)    :: AtomName
    real(8),dimension(NAtoms,3),intent(inout) :: Coord
    ! local variable
    integer(4)  :: i,line_i,line_f
    character(len=100)  :: struc,char_line_i,char_line_f
    logical :: filestat
    
    call GETARG(fid,struc)
    INQUIRE(file=struc,exist=filestat)
    if (filestat)   then
        line_i = (NAtoms+2)*order
        line_f = NAtoms+2
        write(char_line_i,*) line_i
        write(char_line_f,*) line_f
        call system('head -n '//TRIM(ADJUSTL(char_line_i))//' '//TRIM(struc)//' | tail -n '//TRIM(ADJUSTL(char_line_f))//' > tmp.dat')
        open(10,file='tmp.dat',status='old',action='read')
    else
        write(*,'(A)') TRIM(struc)//" doesn't exist"
        stop
    end if
    read(10,*) struc ! buffer
    read(10,*) struc ! buffer
    do i=1,NAtoms
        read(10,*) AtomName(i),Coord(i,1:3)
    end do
    close(10,status='delete')
    return
end subroutine get_coord_tra

subroutine get_coord_pes(fid,order,NAtoms,x,y,E,AtomName,Coord)
    implicit none
    integer(4),intent(in)   :: fid,order,NAtoms
    integer(4),intent(out)  :: x,y
    real(8),intent(out) :: E
    integer(4),dimension(NAtoms),intent(out)    :: AtomName
    real(8),dimension(NAtoms,3),intent(out) :: Coord
    ! local variable
    integer(4)  :: i,line_i,line_f
    character(len=100)  :: struc,char_line_i,char_line_f
    logical :: filestat
    
    call GETARG(fid,struc)
    INQUIRE(file=struc,exist=filestat)
    if (filestat)   then
        line_i = (NAtoms+2)*order
        line_f = NAtoms+2
        write(char_line_i,*) line_i
        write(char_line_f,*) line_f
        call system('head -n '//TRIM(char_line_i)//' '//TRIM(struc)//' | tail -n '//TRIM(char_line_f)//"  > tmp.dat")
        open(10,file='tmp.dat',status='old',action='read')
    else
        write(*,'(A)') TRIM(struc)//" doesn't exist"
        stop
    end if
    read(10,*) struc ! buffer
    read(10,*) x,y,E
    do i=1,NAtoms
    ! do i=1,20
        read(10,*) AtomName(i),Coord(i,1:3)
    end do
    close(10,status='delete')
    return
end subroutine get_coord_pes

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

subroutine get_orderList(NAtoms,num_pes,orderList,energyList)
    implicit none
    integer(4),intent(in)   :: NAtoms,num_pes
    integer(4),dimension(num_pes,3),intent(out) :: orderList
    real(8),dimension(num_pes),intent(out)  :: energyList
    ! local variable
    integer(4)  :: i,x,y
    real(8) :: E
    integer(4),dimension(NAtoms)    :: AtomName ! redundant
    real(8),dimension(NAtoms,3) :: Coord ! redundant
    do i=1,num_pes
        call get_coord_pes(2,i,NAtoms,x,y,E,AtomName,Coord)
        orderList(i,1)=i
        orderList(i,2)=x
        orderList(i,3)=y
        energyList(i)=E
    end do
    return
end subroutine get_orderList

subroutine get_iniPts(NAtoms,num_pes,order_pes)
    implicit none
    integer(4),intent(in)   :: NAtoms,num_pes
    integer(4),intent(inout)   :: order_pes
    ! local variable
    integer(4)  :: i,x1,x2,y1,y2
    real(8) :: E1,E2,rmsd1,rmsd2
    integer(4),dimension(NAtoms)    :: AtomName
    real(8),dimension(NAtoms,3) :: Coord_tra,Coord_pes1,Coord_pes2

    order_pes=1 ! compare with the first point
    call get_coord_tra(1,1,NAtoms,AtomName,Coord_tra)
    call get_coord_pes(2,order_pes,NAtoms,x1,y1,E1,AtomName,Coord_pes1)  
    call get_RMSD(NAtoms,Coord_tra,Coord_pes1,rmsd1)
    ! do i=2,num_pes
    do i=2,num_pes/10  !only searching the 10% region of 2D-PES
        call get_coord_pes(2,i,NAtoms,x2,y2,E2,AtomName,Coord_pes2)  
        call get_RMSD(NAtoms,Coord_tra,Coord_pes2,rmsd2)
        write(*,*) rmsd1
        if ( rmsd1 .gt. rmsd2 ) then
            rmsd1 = rmsd2
            order_pes=i
        end if
    end do
    return
end subroutine get_iniPts

subroutine get_nextPts(NAtoms,num_pes,tra,num_1D,order_pes)
    implicit none
    integer(4),intent(in)   :: NAtoms,num_pes,tra,num_1D
    integer(4),intent(inout)   :: order_pes
    ! local variable
    integer(4)  :: i,j,x,y,order1,num_column
    real(8) :: E,rmsd_small,rmsd
    integer(4),dimension(NAtoms)    :: AtomName
    real(8),dimension(NAtoms,3) :: Coord_tra,Coord_pes

    order1 = order_pes
    call get_coord_tra(1,tra,NAtoms,AtomName,Coord_tra)
    rmsd_small = 100000000.D0 ! vary big rmsd to be replaced
    do i = order1 - num_1D, order1 + num_1D, num_1D
        if ( ( i .lt. 1 ) .or. ( i .gt. num_pes ) ) cycle
        if ( MODULO(i,num_1D) .eq. 0 ) then
            num_column = i/num_1D
        else
            num_column = i/num_1D + 1
        end if
        do j = i - 1, i +1
            ! avoid the original point
                if ( j .eq. order1 ) cycle 
            ! set the boudary of one column in the OrderList
                if (     ( j .lt. (num_column-1)*num_1D+1 ) &
                    .or. ( j .gt. num_column*num_1D       )     ) cycle 
            call get_coord_pes(2,j,NAtoms,x,y,E,AtomName,Coord_pes)  
            ! write(*,'(3(I3,1X))') x,y
            call get_RMSD(NAtoms,Coord_tra,Coord_pes,rmsd)
            if ( rmsd_small .gt. rmsd ) then
                rmsd_small = rmsd
                order_pes = j
            end if
        end do
     end do
    return
end subroutine get_nextPts
