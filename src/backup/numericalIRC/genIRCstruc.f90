!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program:                                                                          !
!   Reproduce IRC path and then compare with gaussian result                        !
!                                                                                   !
! Reference:                                                                        !
!   The Journal of Physical Chemistry 1970, 74, 4161.                               !
!       https://f1000.com/work/item/1366252/resources/2855429/pdf                   !
!                                                                                   !  
! Input:                                                                            !
!   $1 = Q0.out   (gaussian output file w/ force and freq=hpmodes)                  !
!   $2 = 0 or 1 (0 = TS, 1 = other point)                                           !
!   $3 = *.com (gaussian input file; assign output name)                            !
!                                                                                   !
! Output:                                                                           !
!   *.com (user defined name = $4)                                                  !
!                                                                                   !
! History:                                                                          !
!   2018/08/16, Grace                                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program main
implicit none
    character(len=100)  :: filename
    character(len=100)  :: theory = 'RHF',basis = 'STO-3G',charge = '-1',multi = '1'
    integer(4)  :: NAtoms,i
    real(8) :: freq
    real(8),parameter   :: ds = 0.001D0
    integer(4),allocatable,dimension(:)    :: AtomName
    real(8),allocatable,dimension(:,:)  :: iniCoord,finCoord
    real(8),allocatable,dimension(:,:)  :: G,eigenVec

    ! Step 1. Using variable 'NAtoms' to allocate memory
        call GETARG(1,filename)
            call getNatoms(filename,NAtoms)
            allocate(AtomName(NAtoms))
            allocate(iniCoord(NAtoms,3))
            allocate(finCoord(NAtoms,3))
            allocate(G(NAtoms,3))
            allocate(eigenVec(NAtoms,3))
    ! Step 2. Extract gradient of point Q0
        call GETARG(2,filename)
            read(filename,'(I1)') i
        call GETARG(1,filename)
            if ( i .eq. 0 ) then 
                call getH(filename,NAtoms,freq,eigenVec)
                G = eigenVec
            else if ( i .eq. 1 ) then
                call getG(filename,NAtoms,G)
            end if 
            ! call changeGvec2MW(NAtoms,3,AtomName,G)
            call normVec(NAtoms,3,G)
    ! Step 3. Generate the next structure
            call getCoord(filename,NAtoms,AtomName,iniCoord)
            ! call changeCar2MW(NAtoms,3,AtomName,iniCoord)
        call GETARG(3,filename)
            call genGeom(NAtoms,ds,G,iniCoord,finCoord)
            ! call changeMW2Car(NAtoms,3,AtomName,finCoord)
            call writeG09input(NAtoms,filename,theory,basis,charge,multi,AtomName,finCoord)
end program main

subroutine getNAtoms(filename,NAtoms)
    implicit none
    character(len=100),intent(in)   :: filename
    integer(4),intent(out)  :: NAtoms
    call system('grep NAtoms '//TRIM(ADJUSTL(filename))//" | head -n 1 | awk '{print $2}' > NAtoms.tmp ")
    open(10,file='NAtoms.tmp',status='old')
    read(10,*) NAtoms
    close(10,status='delete')
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
    character(len=100)  :: buffer1,buffer2
    i = NAtoms + 2
    write(buffer1,*) i
    write(buffer2,*) NAtoms
    call system('grep -A '//TRIM(ADJUSTL(buffer1))//' Coordinates '//&
        TRIM(ADJUSTL(filename))//' | tail -n '//TRIM(ADJUSTL(buffer2))//&
        " | awk '{print $2,$4,$5,$6}' > Coord.tmp")
    open(10,file='Coord.tmp',status='old')
    do i=1,NAtoms
        read(10,*) AtomName(i),Coord(i,1:3)
    end do
    close(10,status='delete')
    return
end subroutine getCoord

subroutine getG(filename,NAtoms,G)
    implicit none
    character(len=100),intent(in)   :: filename
    integer(4),intent(in)   :: NAtoms
    real(8),dimension(NAtoms,3),intent(inout)   :: G
    ! local variable
    real(8),dimension(NAtoms,3) :: F
    integer(4)  :: i,j
    character(len=100)  :: buffer1,buffer2
    i = NAtoms + 2
    write(buffer1,*) i
    write(buffer2,*) NAtoms
    call system('grep -A '//TRIM(ADJUSTL(buffer1))//" 'Forces (Hartrees/Bohr)' "//&
        TRIM(ADJUSTL(filename))//' | tail -n '//TRIM(ADJUSTL(buffer2))//&
        " | awk '{print $3,$4,$5}' > G.tmp")
    open(10,file='G.tmp',status='old')
    do i=1,NAtoms
        read(10,*) F(i,1:3)
        do j=1,3
            G(i,j) = -F(i,j)
        end do
    end do
    close(10,status='delete')
    ! close(10)
    return
end subroutine getG

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

!FIXME: amu to au !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   amu = atomic mass unit; daltons, Da or u. 1 amu = 1.6605*10^-21 kg      !
!       wiki: https://en.wikipedia.org/wiki/Unified_atomic_mass_unit        !
!   a.u. : me = e = hbar = k_e = 1, 1 a.u. = 9.1*10^-31 kg                  !
!       wiki: https://en.wikipedia.org/wiki/Atomic_units                    !
!   1 amu = 1825 a.u. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! unit of gradient: Hartrees/Bohr
subroutine changeGvec2MW(dimen1,dimen2,AtomName,Vec)
    implicit none
    integer(4),intent(in)   :: dimen1,dimen2
    integer(4),dimension(dimen1),intent(in) :: AtomName
    real(8),dimension(dimen1,dimen2),intent(inout)   :: Vec
    ! local variable
    real(8),parameter:: amu2au = 1825.0D0
    integer(4)  :: i,j
    ! do i=1,dimen1
        ! write(*,*) AtomName(i)
    ! end do
    do i = 1,dimen1
        do j = 1,dimen2
            Vec(i,j) = Vec(i,j)/(SQRT(AtomName(i)*amu2au))
        end do
    end do
    return
end subroutine changeGvec2MW

subroutine changeCar2MW(dimen1,dimen2,AtomName,Coord)
    implicit none
    integer(4),intent(in)   :: dimen1,dimen2
    integer(4),dimension(dimen1),intent(in) :: AtomName
    real(8),dimension(dimen1,dimen2),intent(inout)   :: Coord
    ! local variable
    real(8),parameter:: amu2au = 1825.0D0, ang2bohr=1.889725989D0
    integer(4)  :: i,j
    do i = 1,dimen1
        do j = 1,dimen2
            Coord(i,j) = Coord(i,j)*ang2bohr*(SQRT(AtomName(i)*amu2au))
        end do
    end do
    return
end subroutine changeCar2MW

subroutine changeMW2Car(dimen1,dimen2,AtomName,Coord)
    implicit none
    integer(4),intent(in)   :: dimen1,dimen2
    integer(4),dimension(dimen1),intent(in) :: AtomName
    real(8),dimension(dimen1,dimen2),intent(inout)   :: Coord
    ! local variable
    real(8),parameter:: amu2au = 1825.0D0, ang2bohr=1.889725989D0
    integer(4)  :: i,j
    do i = 1,dimen1
        do j = 1,dimen2
            Coord(i,j) = Coord(i,j)/(ang2bohr*(SQRT(AtomName(i)*amu2au)))
        end do
    end do
    return
end subroutine changeMW2Car

subroutine getH(filename,NAtoms,freq,eigenVec)
    implicit none
    character(len=100),intent(in)   :: filename
    integer(4),intent(in)   :: NAtoms
    real(8),intent(out) :: freq
    real(8),dimension(NAtoms,3),intent(out) :: eigenVec
    ! local variable
    ! real(8),dimension(NAtoms*3) :: tmpV
    integer(4)  :: i,j
    character(len=100)  :: buffer1,buffer2,buffer3
    freq = 0.0D0
    eigenVec = 0.0D0

    i = 7 + 3*NAtoms
    write(buffer1,*) i
    i = 7 + 3*NAtoms + 1
    write(buffer2,*) i
    i = 3*NAtoms
    write(buffer3,*) i

    call system("grep 'Frequencies ---' "//TRIM(ADJUSTL(filename))//&
        " | head -n 1 | awk '{print $3}' > Freq.tmp")
    open(10,file='Freq.tmp',status='old')
    read(10,*) freq
    close(10,status='delete')
    call system("grep -A "//TRIM(ADJUSTL(buffer1))//&
        " 'Frequencies ---' "//TRIM(ADJUSTL(filename))//&
        " | head -n "//TRIM(ADJUSTL(buffer2))//" | tail -n "&
        //TRIM(ADJUSTL(buffer3))//" | awk '{print $4}' > eigenV.tmp")
    open(10,file='eigenV.tmp',status='old')
    
    do i=1,NAtoms
        do j=1,3
            read(10,*) eigenVec(i,j)
        end do
    end do
    close(10,status='delete')
    return
end subroutine getH

subroutine genGeom(NAtoms,ds,G,iniCoord,finCoord)
    implicit none
    integer(4),intent(in)   :: NAtoms
    real(8),intent(in)  :: ds
    real(8),dimension(NAtoms,3),intent(in)  :: G,iniCoord
    real(8),dimension(NAtoms,3),intent(out) :: finCoord
    ! local variable
    integer(4)  :: i,j
    finCoord = 0.0D0
    do i = 1,NAtoms
        do j = 1,3
            finCoord(i,j) = iniCoord(i,j) - G(i,j) * ds
        end do
    end do
    return
end subroutine genGeom

subroutine writeG09input(NAtoms,filename,theory,basis,charge,multi,AtomName,finCoord)
    implicit none
    integer(4),intent(in)   :: NAtoms
    integer(4),dimension(NAtoms)    :: AtomName
    real(8),dimension(NAtoms,3),intent(in) :: finCoord
    character(len=100),intent(in)   :: filename,theory,basis,charge,multi
    ! local variable
    integer(4)  :: i
    open(10,file=filename,status='replace')
    write(10,'(A)') '# '//TRIM(ADJUSTL(theory))//'/'//TRIM(ADJUSTL(basis))//' force'
    write(10,'()')
    write(10,'(A)') 'Generate the IRC structure by the gradient from previous structure'
    write(10,'()')
    write(10,'(A)') TRIM(ADJUSTL(charge))//' '//TRIM(ADJUSTL(multi))
    do i = 1,NAtoms
        write(10,'(I2,3(1X,F12.8))') AtomName(i),finCoord(i,1:3)
    end do
    write(10,'()')
    write(10,'(A)') '--Link1--'
    write(10,'(A)') '# '//TRIM(ADJUSTL(theory))//'/'//TRIM(ADJUSTL(basis))//' Freq=hpmodes'
    write(10,'()')
    write(10,'(A)') 'Generate the IRC structure by the gradient from previous structure'
    write(10,'()')
    write(10,'(A)') TRIM(ADJUSTL(charge))//' '//TRIM(ADJUSTL(multi))
    do i = 1,NAtoms
        write(10,'(I2,3(1X,F12.8))') AtomName(i),finCoord(i,1:3)
    end do
    write(10,'()')
    return
end subroutine writeG09input