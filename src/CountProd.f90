!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program :                                                         !
!   Analyse the last point from trajectory, and then assign the     !
!   product name by user define.                                    !
!                                                                   !
! Input :                                                           !
!   $1 = single molecular structure; can be used via jmol.          !
!   $2 = define the reactant                                        !
!   $3 = define the product 1                                       !
!   $4 = define the product 2                                       !
!   format:                                                         !
!       # Distance criteria # do not remove this line               !
!       amount of distance                                          !
!       1st atom, 2nd atom                                          !
!       # Angle criteria # do not remove this line                  !
!       amount of angle                                             !
!       1st atom, 2nd atom, 3rd atom                                !
!       # Dihedral angle criteria # do not remove this line         !
!       amount of dihedral angle                                    !
!       1st atom, 2nd atom, 3rd atom, 4th atom                      !
!                                                                   !
! Output :                                                          !
!   std-out                                                         !
!                                                                   !
! Reference code:                                                   !
!   proganal                                                        !
!                                                                   !
! 2018/11/30, Grace                                                 !
! 2018/12/04, Grace, separate the trajectory to 2 direction and then!
!   add the 4th input argument to include the user-defined reactant.!
! 2018/12/11, Grace, remove the redundant part, and classify all    !
!    trajectories into 10 categories.                               !
! 2019/07/25, Grace, based on the new process, remove all the       !
!   redundent parts. Back up the old version in /src/backup         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Program main 
    implicit none
    integer(4)  :: NAtoms
    real(8),allocatable,dimension(:,:) :: coord
    character(len=100)  :: defPts

    ! Step 1. Allocate memory 
        call getNAtoms(NAtoms)
        allocate(coord(NAtoms,3))
    ! Step 2. Assign the name of input geometry based on 
    !   the user-defined criteria; the 2nd, 3rd and the 4th input 
    !   arguments.
        call getCoord(NAtoms,coord)
        call getPtsName(NAtoms,coord,defPts) 
    ! Step 3. Std-out
        write(*,'(A)') TRIM(ADJUSTL(defPts))
        deallocate(coord)
stop
end program main

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

subroutine getCoord(NAtoms,coord)
    implicit none
    integer(4),intent(in)   :: NAtoms
    real(8),dimension(NAtoms,3),intent(out) :: coord
    ! local variable
    integer(4)  :: i
    character(len=100)  :: buff, struc
    logical :: filestat

    call GETARG(1,struc)
    open(10,file=struc,status='old')
    read(10,*) buff
    read(10,*) buff
    do i = 1, NAtoms
        read(10,*) buff, coord(i,1:3)
    end do
    close(10)
    return
end subroutine getCoord

subroutine getPtsName(NAtoms,coord,defPts)
    implicit none
    integer(4),intent(in)   :: NAtoms
    real(8),dimension(NAtoms,3),intent(in)  :: coord
    character(len=100),intent(out)  :: defPts
    ! local variables
    logical :: NumProd0,NumProd1, NumProd2

    call getNumProd(1,NAtoms,coord,NumProd0) ! reactant
    call getNumProd(2,NAtoms,coord,NumProd1) ! product 1
    call getNumProd(3,NAtoms,coord,NumProd2) ! product 2
    if ( NumProd0 ) then
        defPts='R'
    else if ( NumProd1 ) then  
        defPts='P1'
    else if ( NumProd2 ) then
        defPts='P2'
    else
        defPts='none'
    end if
    return
end subroutine getPtsName

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