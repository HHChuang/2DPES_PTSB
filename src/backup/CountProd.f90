!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program :                                                 !
!   Analyse the last point from trajectory, and then        !
!   assign the product name by user define.                 !
!                                                           !
! Input :                                                   !
!   $1 = name of trajectory                                 !
!   $2 = define the reactant                                !
!   $3 = define the product 1                               !
!   $4 = define the product 2                               !
!   format:                                                 !
!       # Distance criteria # do not remove this line       !
!       amount of distance                                  !
!       1st atom, 2nd atom                                  !
!       # Angle criteria # do not remove this line          !
!       amount of angle                                     !
!       1st atom, 2nd atom, 3rd atom                        !
!       # Dihedral angle criteria # do not remove this line !
!       amount of dihedral angle                            !
!       1st atom, 2nd atom, 3rd atom, 4th atom              !
!                                                           !
! Output :                                                  !
!   std-out                                                 !
!                                                           !
! Reference code:                                           !
!   proganal                                                !
!                                                           !
! 2018/11/30, Grace                                         !
! 2018/12/04, Grace, separate the trajectory to 2 direction !
! and then add the 4th input argument to include the        !
! user-defined reactant.                                    !
! 2018/12/11, Grace, remove the redundant part, and classify!
! all trajectories into 10 categories.                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Program main 
    implicit none
    integer(4)  :: NAtoms
    real(8),allocatable,dimension(:,:) :: coord1,coord2
    character(len=100)  :: point1,point2

    ! Step 1. Std-out purpose and extract the coordinate 
    !   of the last point from the trajectory
        ! call print_purpose()
        call getNAtoms(NAtoms)
        allocate(coord1(NAtoms,3))
        allocate(coord2(NAtoms,3))
        call getLastPoint(1,NAtoms,coord1)
        call getLastPoint(2,NAtoms,coord2)
    ! Step 2. Assign the product as P1 or P2
        call getPoint(NAtoms,coord1,point1) ! call getNumProd
        call getPoint(NAtoms,coord2,point2) ! call getNumProd
    ! Step 3. - Reorder trajectory
    !         - Only 10 categories exist now. and swap point1 
    !           for some trajectories. 
    !         - Remove redundant part
        call reorder_cut_Traj(NAtoms,point1,point2) ! output: reorder.$1
    ! Step4. Std-out
        write(*,'(A)') TRIM(ADJUSTL(point1))//' '//TRIM(ADJUSTL(point2))
stop
end program main

subroutine print_purpose()
    implicit none
    write(*,'(A)') '------------------------------------------'
    write(*,'(A)') 'Analyse the last point from trajectory,' 
    write(*,'(A)') 'and then assign the product name by' 
    write(*,'(A)') 'user-defined criteria.'
    write(*,'()')
    write(*,'(A)') 'Reorder all trajectories, and then remove'
    write(*,'(A)') 'redundant part of trajectory in RP1 and '
    write(*,'(A)') 'RP2 categories.'
    write(*,'(A)') '------------------------------------------'
    return
end subroutine print_purpose

subroutine getPoint(NAtoms,coord,point)
    implicit none
    integer(4),intent(in)   :: NAtoms
    real(8),dimension(NAtoms,3),intent(in)  :: coord
    character(len=100),intent(out)  :: point
    ! local variables
    logical :: NumProd0,NumProd1, NumProd2

    call getNumProd(1,NAtoms,coord,NumProd0) ! reactant
    call getNumProd(2,NAtoms,coord,NumProd1) ! product 1
    call getNumProd(3,NAtoms,coord,NumProd2) ! product 2
    if ( NumProd0 ) then
        point='R'
    else if ( NumProd1 ) then  
        point='P1'
    else if ( NumProd2 ) then
        point='P2'
    else
        point='none'
    end if
    return
end subroutine getPoint

subroutine getNAtoms(NAtoms)
    implicit none
    integer(4),intent(out)  :: NAtoms
    ! local variable
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
end subroutine getNAtoms

subroutine getLastPoint(index,NAtoms,coord)
    implicit none
    integer(4),intent(in)   :: index,NAtoms
    real(8),dimension(NAtoms,3),intent(out) :: coord
    ! local variable
    character(len=100)  :: struc,NAtoms_char,buffer
    character(len=100)  :: ini_char,fin_char
    integer(4)  :: i,ini,fin

    write(NAtoms_char,*) NAtoms
    call GETARG(1,struc)
    select case(index)
    case(1) 
        call system("grep -n 'runpoint 1' "//TRIM(ADJUSTL(struc))&
            //" | tail -n 1 | cut -d ':' -f 1 > buffer.txt" )
        open(10,file='buffer.txt',status='old')
        read(10,*) i
        close(10,status='delete')
        ini=i-1-NAtoms
        fin=i-2
        write(ini_char,*) ini
        write(fin_char,*) fin
        call system("sed -n '"//TRIM(ADJUSTL(ini_char))&
            //','//TRIM(ADJUSTL(fin_char))//&
            " p' "//TRIM(ADJUSTL(struc))//' > lastpoint.txt')
    case(2)
        call system('tail -n '//TRIM(ADJUSTL(NAtoms_char))&
            //' '//TRIM(ADJUSTL(struc))//' > lastpoint.txt')
    end select

    open(10,file='lastpoint.txt',status='old')
    do i=1,NAtoms
        read(10,*) buffer,coord(i,1:3)
    end do
    close(10,status='delete')
    return
end subroutine getLastPoint

subroutine getNumProd(order,NAtoms,coord,NumProd)
    implicit none
    integer(4),intent(in)   :: order ! order = 1, 2 or 3
    integer(4),intent(in)   :: NAtoms
    real(8),dimension(NAtoms,3),intent(in) :: coord
    logical,intent(out)  :: NumProd
    ! local variable
    integer(4)  :: i
    integer(4)  :: NumDis,NumAng,NumDih
    integer(4),allocatable,dimension(:,:) :: IndexDis,IndexAng,IndexDih
    real(8),allocatable,dimension(:) :: CriDis,CriAng,CriDih
    character(len=10),allocatable,dimension(:)   :: CondDis,CondAng,CondDih
    character(len=10),allocatable,dimension(:)   :: LogDis,LogAng,LogDih
    real(8),allocatable,dimension(:) :: LPDis,LPAng,LPDih ! LP=last point
    ! temporary indicator
    integer(4) :: tmpInd, sumInd, Ind, totInd

    ! 1. get the user-defined criteria
        call getNumDefProd(order,NumDis,NumAng,NumDih)
        allocate(IndexDis(NumDis,2))
        allocate(IndexAng(NumAng,3))
        allocate(IndexDih(NumDih,4))
        allocate(CriDis(NumDis))
        allocate(CriAng(NumAng))
        allocate(CriDih(NumDih))
        allocate(CondDis(NumDis))
        allocate(CondAng(NumAng))
        allocate(CondDih(NumDih))
        allocate(LogDis(NumDis))
        allocate(LogAng(NumAng))
        allocate(LogDih(NumDih))
        allocate(LPDis(NumDis))
        allocate(LPAng(NumAng))
        allocate(LPDih(NumDih))
        call getDefProd(order,NumDis,NumAng,NumDih,IndexDis,IndexAng,IndexDih,CriDis,CriAng,CriDih,CondDis,CondAng,CondDih,LogDis,LogAng,LogDih)
    ! 2. calculate the geometric parameter of the last point of the trajectory
        do i=1,NumDis
            call getDistance(NAtoms,coord,IndexDis(i,1),IndexDis(i,2),LPDis(i))
        end do
        do i=1,NumAng
            call getAngle(NAtoms,coord,IndexAng(i,1),IndexAng(i,2),IndexAng(i,3),LPAng(i))
        end do
        do i=1,NumDih
            call getDihedralAngle(NAtoms,coord,IndexDih(i,1),IndexDih(i,2),IndexDih(i,3),IndexDih(i,4),LPDih(i))
        end do
    ! 3. compare the geometric parameter with user-defined criteria
    !       to assign product
        Ind = 0
        totInd = 0
        if ( NumDis .ne. 0 ) then
            do i=1,NumDis
                call getTmpInd(LPDis(i),CriDis(i),CondDis(i),LogDis(i),tmpInd,sumInd)
                Ind = Ind + tmpInd
                totInd = totInd + sumInd
            end do
        end if 
        if ( NumAng .ne. 0 ) then
            do i=1,NumAng
                call getTmpInd(LPAng(i),CriAng(i),CondAng(i),LogAng(i),tmpInd,sumInd)
                Ind = Ind + tmpInd
                totInd = totInd + sumInd
            end do
        end if
        if ( NumDih .ne. 0 ) then
            do i=1,NumDih
                call getTmpInd(LPDih(i),CriDih(i),CondDih(i),LogDih(i),tmpInd,sumInd)
                Ind = Ind + tmpInd
                totInd = totInd + sumInd
            end do
        end if
        if ( Ind .eq. totInd ) then
            NumProd = .true.
        else
            NumProd = .false.
        end if
    return
end subroutine getNumProd

subroutine getNumDefProd(order,NumDis,NumAng,NumDih)
    implicit none
    integer(4),intent(in)   :: order
    integer(4),intent(out)  :: NumDis,NumAng,NumDih
    ! local variable
    character(len=100)  :: struc
    logical :: filestat

    call GETARG(order+1,struc) ! get the order+1 argument
    INQUIRE(file=struc,exist=filestat)
    if (filestat) then
        open(10,file=struc,status='old')
        call system('grep -A 1 \# '//TRIM(struc) &
        //'| sed "/--/ d" | sed -n "2,2 p" > NumDis.tmp')
        call system('grep -A 1 \# '//TRIM(struc) &
        //'| sed "/--/ d" | sed -n "4,4 p" > NumAng.tmp')
        call system('grep -A 1 \# '//TRIM(struc) &
        //'| sed "/--/ d" | sed -n "6,6 p" > NumDih.tmp')
        open(10,file='NumDis.tmp',status='old')
        read(10,*) NumDis
        close(10,status='delete')
        open(10,file='NumAng.tmp',status='old')
        read(10,*) NumAng
        close(10,status='delete')
        open(10,file='NumDih.tmp',status='old')
        read(10,*) NumDih
        close(10,status='delete')
    else
        write(*,'(A)') TRIM(struc)//" doesn't exist"
        stop
    end if
    close(10)
    return
end subroutine getNumDefProd

subroutine getDefProd(order,NumDis,NumAng,NumDih,IndexDis,IndexAng,IndexDih,CriDis,CriAng,CriDih,CondDis,CondAng,CondDih,LogDis,LogAng,LogDih)
    implicit none
    integer(4),intent(in)   :: order,NumDis,NumAng,NumDih
    integer(4),dimension(NumDis,2),intent(out)  :: IndexDis
    integer(4),dimension(NumAng,3),intent(inout)  :: IndexAng
    integer(4),dimension(NumDih,4),intent(out)  :: IndexDih
    real(8),dimension(NumDis),intent(out)   :: CriDis
    real(8),dimension(NumAng),intent(inout)   :: CriAng
    real(8),dimension(NumDih),intent(inout)   :: CriDih
    character(len=10),dimension(NumDis),intent(out)   :: CondDis,LogDis
    character(len=10),dimension(NumAng),intent(inout)   :: CondAng,LogAng
    character(len=10),dimension(NumDih),intent(inout)   :: CondDih,LogDih
    !local variable
    character(len=100)  :: struc
    character(len=100)  :: NumDis_char,NumAng_char,NumDih_char
    character(len=100)  :: NumDisPlus1_char,NumAngPlus1_char,NumDihPlus1_char
    integer(4)  :: i

    call GETARG(order+1,struc)
    write(NumDis_char,*) NumDis
    write(NumAng_char,*) NumAng
    write(NumDih_char,*) NumDih
    write(NumDisPlus1_char,*) NumDis+1
    write(NumAngPlus1_char,*) NumAng+1
    write(NumDihPlus1_char,*) NumDih+1

    if ( NumDis .ne. 0 ) then
        call system('grep -A '//TRIM(ADJUSTL(NumDisPlus1_char))&
            //" '# Distance criteria' "&
            //TRIM(struc)//' | tail -n '//TRIM(ADJUSTL(NumDis_char))&
            //' > Dis.tmp')
        open(10,file='Dis.tmp',status='old')
        do i=1,NumDis
            read(10,*) IndexDis(i,1:2),CriDis(i),CondDis(i),LogDis(i)
        end do
        close(10,status='delete')
    end if
    if ( NumAng .ne. 0 ) then
        call system('grep -A '//TRIM(ADJUSTL(NumAngPlus1_char))&
            //" '# Angle criteria ' "&
            //TRIM(struc)//' | tail -n '//TRIM(ADJUSTL(NumAng_char))&
            //' > Ang.tmp')
        open(10,file='Ang.tmp',status='old')
        do i=1,NumAng
            read(10,*) IndexAng(i,1:3),CriAng(i),CondAng(i),LogAng(i)
            ! unit convert, from degree to redian
            CriAng(i) = CriAng(i) * 3.14159D0 / 180.0D0
        end do
        close(10,status='delete')
    end if
    if ( NumDih .ne. 0 ) then
        call system('grep -A '//TRIM(ADJUSTL(NumDihPlus1_char))&
            //" '# Dihedral angle criteria' "&
            //TRIM(struc)//' | tail -n '//TRIM(ADJUSTL(NumDih_char))&
            //' > Dih.tmp')
        open(10,file='Dih.tmp',status='old')
        do i=1,NumDih
            read(10,*) IndexDih(i,1:4),CriDih(i),CondDih(i),LogDih(i)
        end do
        close(10,status='delete')
    end if
    return
end subroutine getDefProd

subroutine getDistance(NAtoms,coord,Index1,Index2,length)
    implicit none
    integer(4),intent(in)   :: NAtoms,Index1,Index2
    real(8),dimension(NAtoms,3),intent(in)  :: coord
    real(8),intent(inout) :: length
    ! local variable
    integer(4)  :: i
    real(8) :: count
    length=0.0D0
    do i=1,3
        count=( coord(Index1,i) - coord(Index2,i) )**2
        length=length+count
    end do
    length=sqrt(length)
    return
end subroutine getDistance

subroutine getAngle(NAtoms,coord,Index1,Index2,Index3,angle)
    implicit none
    integer(4),intent(in)   :: NAtoms,Index1,Index2,Index3
    real(8),dimension(NAtoms,3),intent(in)  :: coord
    real(8),intent(inout) :: angle ! unit: radian
    ! local variable
    integer(4)  :: i
    real(8) :: dist12,dist13,dist23,value

    call getDistance(NAtoms,coord,Index1,Index2,dist12)
    call getDistance(NAtoms,coord,Index1,Index3,dist13)
    call getDistance(NAtoms,coord,Index2,Index3,dist23)
    value=( (-dist13**2 + dist12**2 + dist23**2) &
            / (2*dist12*dist23) &
          )
    angle=ACOS(value) 
    return
end subroutine getAngle

! reference for formula: https://en.wikipedia.org/wiki/Dihedral_angle
subroutine getDihedralAngle(NAtoms,coord,Index1,Index2,Index3,Index4,dihedral)
    implicit none
    integer(4),intent(in)   :: NAtoms,Index1,Index2,Index3,Index4
    real(8),dimension(NAtoms,3),intent(in)  :: coord
    real(8),intent(out) :: dihedral
    ! local variable
    integer(4)  :: i
    real(8),dimension(3)    :: b1,b2,b3,yA
    real(8) :: modB2,termY,termX
    real(8),dimension(3)    :: CP23,CP12 !cross product 

    modB2=0.0D0
    do i=1,3 !x: i=1, y: i=2, z: i=3
        b1(i) = coord(Index2,i) - coord(Index1,i)
        b2(i) = coord(Index3,i) - coord(Index2,i)
        b3(i) = coord(Index4,i) - coord(Index3,i)
        modB2 = b2(i)**2 + modB2
    end do
    modB2 = sqrt(modB2)
    call getcp(b2,b3,CP23)
    call getcp(b1,b2,CP12)
    termY = 0.0D0
    termX = 0.0D0
    ! yA is x-coord. etc of modulus of B2 times B1
    do i=1,3
        yA(i) = modB2*b1(i)
        termY = yA(i)*CP23(i) + termY
        termX = CP12(i)*CP23(i) + termX
    end do
    dihedral = (180.0D0/3.141592) * atan2(termY,termX)
    return
end subroutine getDihedralAngle

subroutine getcp(a, b, cp)
  real(8), dimension(3),intent(out) :: cp
  real(8), dimension(3), intent(in) :: a, b

  cp(1) = a(2) * b(3) - a(3) * b(2)
  cp(2) = a(3) * b(1) - a(1) * b(3)
  cp(3) = a(1) * b(2) - a(2) * b(1)
end subroutine getcp

subroutine getTmpInd(LP,Cri,Cond,Log,tmpInd,sumInd)
    implicit none
    real(8),intent(in)  :: LP,Cri
    integer(4),intent(out)  :: tmpInd,sumInd
    character(len=4),intent(in) :: Cond, Log
    select case(Cond)
    case ('lt')
        if ( LP .lt. Cri ) then
            tmpInd = 1 ! true
        else
            tmpInd = 0 ! false
        end if
    case ('gt')
        if ( LP .gt. Cri ) then 
            tmpInd = 1 ! true
        else
            tmpInd = 0 ! false
        end if
    case default
        write(*,'(A)') 'Wrong input argument! Stop program'
        stop
    end select

    select case(Log)
    case ('and')
        sumInd=1
    case ('or')
        sumInd=0
    case('end')
        sumInd=1
    case default
        write(*,'(A)') 'Wrong input argument! Stop program'
        stop
    end select
        ! write(*,*) LP,Cri,tmpInd
    return
end subroutine getTmpInd

subroutine reorder_cut_Traj(NAtoms,point1,point2)
    implicit none
    integer(4),intent(in)   :: NAtoms
    character(len=100),intent(inout)  :: point1,point2
    ! local variable
    character(len=100)  :: arg1,point_tmp
    integer(4)  :: onefile,keyline,ini_traj1,fin_traj1,ini_traj2,fin_traj2
    character(len=100)   :: ini_traj1_char,fin_traj1_char,ini_traj2_char,fin_traj2_char
    
    ! Step1. input necessary variable, and then cut total 
    !       trajectory into two part, traj1 and traj2.
        call GETARG(1,arg1)
        onefile = NAtoms + 2
        call system("grep -n 'runpoint 1' "//TRIM(ADJUSTL(arg1)) &
                //"| tail -n 1 | cut -d ':' -f 1 > keyline")
        open(10,file='keyline',status='old')
        read(10,*) keyline
        close(10,status='delete')
        ini_traj1 = 1
        fin_traj1 = keyline - 2
        ini_traj2 = keyline + NAtoms +1 ! remove the redundant TSS
        call system("wc -l "//TRIM(ADJUSTL(arg1)) &
                //"| awk '{print $1}' > fin_traj2 ")
        open(10,file='fin_traj2',status='old')
        read(10,*) fin_traj2
        close(10,status='delete')
        write(ini_traj1_char,*) ini_traj1
        write(fin_traj1_char,*) fin_traj1
        write(ini_traj2_char,*) ini_traj2
        write(fin_traj2_char,*) fin_traj2
        call system("sed -n '"//TRIM(ADJUSTL(ini_traj1_char))//',' &
                        //TRIM(ADJUSTL(fin_traj1_char))//" p' " &
                        //TRIM(ADJUSTL(arg1))//' > traj1')
        call system("sed -n '"//TRIM(ADJUSTL(ini_traj2_char))//',' &
                        //TRIM(ADJUSTL(fin_traj2_char))//" p' " &
                        //TRIM(ADJUSTL(arg1))//' > traj2')
    ! Step2. remove the redundant part
        call cutTraj(NAtoms,point1,'traj1')
        call cutTraj(NAtoms,point2,'traj2')
    ! Step3. Reorder the total trajectory, output: reorder.$arg1
    ! 10 categories: 1. RR, 2. P1P1, 3. P2P2, 4. nonenone, 5. noneR
    !   6. noneP1, 7. noneP2, 8. P1P2, 9. RP1, 10. RP2.
        if ( &
            ( (point1 .eq. 'R') .and. (point2 .eq. 'R') ) .or. &
            ( (point1 .eq. 'P1') .and. (point2 .eq. 'P1') ) .or. &
            ( (point1 .eq. 'P2') .and. (point2 .eq. 'P2') ) .or. &
            ( (point1 .eq. 'none') .and. (point2 .eq. 'none') ) &
        ) then
            ! First group: RR, P1P1, P2P2, nonenone
            call traj12(onefile,arg1)
        else
            ! Second group: noneR, noneP1,noneP2,P1P2,RP1,RP2
            if ( (point1 .eq. 'none') .or. (point1 .eq. 'R') ) then
                call traj12(onefile,arg1)
            else if ( (point1 .eq. 'P1') .and. (point2 .eq. 'P2') ) then 
                call traj12(onefile,arg1)
            else 
                call traj21(onefile,arg1)
                point_tmp = point1
                point1 = point2
                point2 = point_tmp
            end if  
        end if
        call system('rm -f traj1 traj2')
    return
end subroutine reorder_cut_Traj

subroutine reverseTraj(onefile,traj)
    implicit none
    integer(4),intent(in)   :: onefile
    character(len=*),intent(in)   :: traj
    ! local variable
    integer(4)  :: totline,npts,i,tmp
    character(len=100)   :: onefile_char,tmp_char

    call system('wc -l '//TRIM(ADJUSTL(traj)) &
        //"| awk '{print $1}'  > totline")
    open(10,file='totline',status='old')
    read(10,*)  totline
    close(10,status='delete')
    npts = totline / onefile
    write(onefile_char,*) onefile
    call system('rm -f revtraj')
    do i=1,npts
        tmp = onefile * i
        write(tmp_char,*) tmp
        call system('tail -n '//TRIM(ADJUSTL(tmp_char))//' '//&
                TRIM(ADJUSTL(traj))//' | head -n '//&
                TRIM(ADJUSTL(onefile_char))//' >> revtraj ')
    end do
    call system('mv -f revtraj '//TRIM(ADJUSTL(traj)))
    return
end subroutine reverseTraj

subroutine traj12(onefile,tottraj)
    implicit none
    integer(4),intent(in)   :: onefile
    character(len=100),intent(in)   :: tottraj
    call reverseTraj(onefile,"traj1")
    call system('cat traj1 > reorder.'//TRIM(ADJUSTL(tottraj)) )
    call system('cat traj2 >> reorder.'//TRIM(ADJUSTL(tottraj)) )
    return
end subroutine traj12

subroutine traj21(onefile,tottraj)
    implicit none
    integer(4),intent(in)   :: onefile
    character(len=100),intent(in)   :: tottraj
    call reverseTraj(onefile,"traj2")
    call system('cat traj2 > reorder.'//TRIM(ADJUSTL(tottraj)) )
    call system('cat traj1 >> reorder.'//TRIM(ADJUSTL(tottraj)) )
    return
end subroutine traj21

subroutine cutTraj(NAtoms,point,traj)
    implicit none
    integer(4),intent(in)   :: NAtoms
    character(len=100),intent(in)  :: point
    character(len=*),intent(in) :: traj
    !local variable
    logical :: NumProd0,NumProd1, NumProd2
    real(8),dimension(NAtoms,3) :: coord
    integer(4)  :: i,j,totline,npts,tailline,headline,cutpts,cutline
    character(len=100)  :: tailline_char,headline_char,buffer,cutline_char

    if ( point .eq. 'none' ) then
        call system('wc -l '//TRIM(ADJUSTL(traj))&
                //"| awk '{print $1}' > totline")
        open(10,file='totline',status='old')
        read(10,*) totline
        close(10,status='delete')
        headline = NAtoms +2
        npts = totline / headline
        write(headline_char,*) headline
        call system('cp '//TRIM(ADJUSTL(traj))//' tmptraj' )
        do i=1,npts
            tailline = (NAtoms+2)*i
            write(tailline_char,*) tailline
            ! check the structure form the last one to the 
            ! first one (i.e. TSS)
            call system('tail -n '//TRIM(ADJUSTL(tailline_char)) &
                    //' '//TRIM(ADJUSTL(traj))//' | head -n '&
                    //TRIM(ADJUSTL(headline_char))//' > coord')
            open(10,file='coord',status='old')
            read(10,*) buffer
            read(10,*) buffer
            do j=1,NAtoms
                read(10,*) buffer,coord(j,1:3)
            end do
            close(10,status='delete')
            call getPoint(NAtoms,coord,point)
            cutpts = i 
            if ( point .ne. 'none' ) exit
        end do
        ! only reserve useful part of trajectory
        cutpts = npts - cutpts + 1
        cutline = headline * cutpts 
        write(cutline_char,*) cutline
        call system("head -n "//TRIM(ADJUSTL(cutline_char)) &
            //" tmptraj > "//TRIM(ADJUSTL(traj)) )
        call system('rm -f tmptraj')
    end if
    return
end subroutine cutTraj