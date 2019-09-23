!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program                                                   !
!   Calculate the dot product between a given gradient and  !
!   and several eigenvectors from Hessian                   !
!                                                           !
! I/O part                                                  !
!   input:                                                  !
!       $1 = fileList.dat; name of files                    !
!                                                           !
!   output:                                                 !
!       dot.dat                                             !
!       (format)                                            !
!       $coord Eigenvalue dotProduct                        !
!                                                           !
! Flowchart                                                 !
!       1. Extract the # of Atoms in order to allocate      !
!           all vectors                                     !
!       2. calculate the dot product                        !
!       3. write the result into file, dot.dat.             !
!                                                           !
! Pre-requested script                                      !
!   getVec.sh                                               !
! History                                                   !
!   2018/07/09, Grace                                       !
!   2018/07/17, Grace, follow the definition of VRI         !
!               J. Mole. STRUC. 695-696 (2004) 95-101       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Program main
    implicit none
    character(len=100)  :: fileList,fileName
    integer(4)  :: NAtoms,NFile,NumEig
    real(8) :: Freq,dotP
    real(8),allocatable,dimension(:)    :: GVec1,EigVec2
    integer(4)  :: i,j

! 1. Extract the # of Atoms in order to allocate all vectors
!    var.: NAtoms  
    call getNAtoms(NAtoms)
    
! ! 2. calculate the dot product
    call getNFile(NFile)
    call GETARG(1,fileList)
    open(100,file=TRIM(fileList),status='old',action='read')
    open(200,file='dot.dat',status='replace')
    ! create file dot.dat
    do i=1,NFile
        read(100,*) fileName
        write(200,'(A)',advance='no') TRIM(fileName)
        call getNumEig(fileName,NumEig)
        allocate(GVec1(3*NAtoms))
        call getGVec(fileNAme,NAtoms,GVec1)
        do j=1,NumEig
            allocate(EigVec2(3*NAtoms))
            call getEig2(fileName,j,NAtoms,Freq,EigVec2)
            dotP = DOT_PRODUCT(GVec1,EigVec2)
            deallocate(EigVec2)
        ! 3. write the result into file, dot.dat. 
            if ( j .ne. NumEig ) then
                write(200,'(F12.4,1X,F13.9)',advance='no') Freq, dotP
            else
                write(200,'(F12.4,1X,F13.9)') Freq,dotP
            end if
        end do
        deallocate(GVec1)
    end do
    close(100)
    close(200)
    ! call system("sed -i 's/n//g' dot.dat")
    stop
end program main

subroutine getNAtoms(NAtoms)
    implicit none
    integer(4),intent(out)  :: NAtoms
    ! local variable
    character(len=100)  :: input,testfile
    call GETARG(1,input)
    call system("head -n 1 "//TRIM(input)//" > testfile.tmp")
    open(10,file='testfile.tmp',status='old',action='read')
    read(10,*) testfile
    close(10,status='delete')
    call system("head -n 1 "//TRIM(testfile)//".G.dat > NAtoms.tmp")
    open(10,file='NAtoms.tmp',status='old',action='read')
    read(10,*) NAtoms
    close(10,status='delete')
    return
end subroutine getNAtoms

subroutine getNFile(NFile)
    implicit none
    integer(4),intent(out)  :: NFile
    ! local variables
    character(len=100)  :: fileName
    call GETARG(1,fileName)
    open(10,file=TRIM(fileName),status='old',action='read')
    NFile=0
    do while (.not. eof(10))
        read(10,*) fileName ! buffer
        NFile=NFile+1
    end do
    close(10)
    return
end subroutine getNFile

subroutine getNumEig(fileName,NumEig)
    implicit none
    character(len=100),intent(in)   :: fileName
    integer(4),intent(out)  :: NumEig
    ! local variables
    integer(4)  :: NAtoms,count
    character(len=100)  :: buffer
    open(10,file=TRIM(fileName)//'.Eig.dat',status='old',action='read')
    read(10,*) NAtoms
    rewind(10)
    count=0
    do while (.not. eof(10))
        read(10,*) buffer
        count=count+1
    end do
    NumEig = count / ( 3 * NAtoms +2)
    close(10)
    return
end subroutine getNumEig

subroutine getGVec(fileName,NAtoms,GVec1)
    implicit none
    character(len=100),intent(in)   :: fileName
    integer(4),intent(in)   :: NAtoms
    real(8),dimension(3*NAtoms),intent(out) :: GVec1
    ! local variables
    character(len=100)  :: buffer
    integer(4)  :: i
    open(10,file=TRIM(fileName)//'.G.dat',status='old',action='read')
    read(10,*)  buffer
    do i=1,3*NAtoms
        read(10,*)  GVec1(i)
    end do
    return
end subroutine getGVec

subroutine getEig2(fileName,order,NAtoms,Freq,EigVec2)
    implicit none
    character(len=100),intent(in)   :: fileName
    integer(4),intent(in)   :: order,NAtoms
    real(8),intent(out) :: Freq
    real(8),dimension(3*NAtoms),intent(out) :: EigVec2
    ! local variables
    character(len=100)  :: buffer
    integer(4)  :: buffer_number,i
    open(10,file=TRIM(fileName)//'.Eig.dat',status='old',action='read')
    buffer_number = order - 2
    if ( buffer_number .eq. -1) then
        read(10,*) buffer ! NAtoms
        read(10,*) Freq
        do i=1,3*NAtoms
            read(10,*)  EigVec2(i)
        end do
    else
        ! read buffer
        do i=(3*NAtoms+2)*buffer_number+1,(3*NAtoms+2)*(buffer_number+1)
            read(10,*) buffer
        end do
        ! read wanted information
        read(10,*) buffer ! NAtoms
        read(10,*) Freq
        do i=1,3*NAtoms
            read(10,*)  EigVec2(i)
        end do
    end if
    close(10)
    return
end subroutine getEig2