! 2018/11/08, Grace

program main
    implicit none
    integer(4)  :: i,j
    character(len=100),dimension(22)  :: filename
    real(8),dimension(3*5)  :: TS1,TS2
    real(8),dimension(5) :: dot1,dot2
    real(8),dimension(2)    :: eigVal

    write(*,'(A)') 'Output: Dot.dat'
    write(*,'()')
    open(10,file='list',status='old')
    do i=1,22
        read(10,*) filename(i)
    end do
    close(10)
    open(10,file='TS1.vector.txt',status='old')
    open(20,file='TS2.vector.txt',status='old')
    do i=1,3*5
        read(10,*) TS1(i)
        read(20,*) TS2(i)
    end do
    close(10)
    close(20)
    
    open(10,file='Dot.dat',status='replace')
    do i=1,22
        do j=1,5
            call getdot(15,TS1,filename(i),j,dot1(j))
            call getdot(15,TS2,filename(i),j,dot2(j))
        end do
        call getEigVal(15,filename(i),MAXLOC(dot1),eigVal(1))
        call getEigVal(15,filename(i),MAXLOC(dot2),eigVal(2))
        write(10,'(A7,2(1X,I1,1X,F9.5,1X,ES14.7))') TRIM(filename(i)),MAXLOC(dot1),MAXVAL(dot1),eigVal(1),MAXLOC(dot2),MAXVAL(dot2),eigVal(2)
    end do
    close(10)
stop
end program main

subroutine getdot(dimen,TS,filename,order,dot)
    implicit none
    integer(4),intent(in)   :: dimen,order
    real(8),dimension(dimen),intent(in) :: TS
    character(len=100),intent(in)   :: filename
    real(8),intent(inout) :: dot
    ! local variable !
    character(len=100)   :: ch_order
    integer(4)  :: i
    real(8),dimension(dimen)    :: fileVec
    real(8) :: cum

    write(ch_order,*)   order
    open(100,file=TRIM(filename)//'.vector'//TRIM(ADJUSTL(ch_order))//'.dat',status='old')
    do i=1,dimen
       read(100,*) fileVec(i)
    end do
    close(100)
    dot=0.0D0
    do i=1,dimen
       cum=TS(i)*fileVec(i)
       dot=dot+cum
    end do
    dot=ABS(dot)
return
end subroutine getdot
    
subroutine getEigVal(dimen,filename,order,eigVal)
    implicit none
    integer(4),intent(in)   :: dimen,order
    character(len=100),intent(in)   :: filename
    real(8),intent(out) :: eigVal
    ! local variable
    real(8) :: mass,freq
    character(len=100)  :: ch_order
    real(8),parameter   :: zero=0.00001D0

    write(ch_order,*)   order
    open(100,file=TRIM(filename)//'.mass'//TRIM(ADJUSTL(ch_order))//'.dat',status='old')
    read(100,*)  mass
    close(100)
    open(100,file=TRIM(filename)//'.freq'//TRIM(ADJUSTL(ch_order))//'.dat',status='old')
    read(100,*)  freq
    close(100)
    ! 1amu=1840au, 1Bohr=5.2819*10^-9 cm
    eigVal=(mass*freq**2)*(1840.0D0)*(5.2819D0*10**(-9.0D0))**2
    if ( freq .lt. zero ) then 
        eigVal = -eigVal ! complex number
    end if
    return
end subroutine getEigVal