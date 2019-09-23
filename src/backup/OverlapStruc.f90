!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Program : Take structure 1 as the reference and then move the       !
!             origin of structure 2 to structure 1, after that, rotate  !
!             structure 2 as close as possible to structure 1.          !
!                                                                       !
!   Pre-requisted shell script  :                                       !
!       1. getCoord                                                     !
!           format:                                                     !
!                   num. of atom                                        !
!                   comment                                             !
!                   #atom   x   y   z                                   !
!                                                                       !
!   Input :                                                             !
!           $1 = structure 1 (variable)                                 !
!           $2 = structure 2 (reference)                                !
!           name of the file: *.xyz                                     !
!                                                                       !
!   Output :                                                            !
!           1. overlap.$1                                               !
!           2. RMSD.$1                                                  !
!           3. rotM                                                     !
!                                                                       !
!   Visualizer :                                                        !
!           Jmol                                                        !
!                                                                       !
!   Reference:                                                          !
!       1. Video and slide: Math for Game Programmers: Understanding    !
!          Quaternions, Jim Van Verth                                   !!!!!!!!
!                                                                              !
!https://www.gdcvault.com/play/1017653/Math-for-Game-Programmers-Understanding !
!                                                                              !
!       2. Youtube: 3D Rotations in General: Rodrigues Rotation Formula !!!!!!!!
!                   and Quaternion Exponentials                         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                              !
!https://www.youtube.com/watch?annotation_id=annotation_1566057571&feature=iv&src_vid=UaK2q22mMEg&v=q-ESzg03mQc!
!                                                                                                              !
!       3. Github: 3D rotation converter                                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                       !
!           https://github.com/gaschler/rotationconverter               !
!                                                                       !
!       4. Superimposing multiple structures and exploring protein      !
!           binding sites.                                              !
!                                                                       !
!           https://is.muni.cz/th/u0qtg/thesis.pdf                      !   
!           code: http://dillgroup.stonybrook.edu/#/landscapes          !                  
!                                                                       !
!   History:                                                            !
! 2017/05/25, Grace, Cherri's rotation matrix.                          !
! 2017/08/07, Grace, Rodrigues rotation (Chun-I).                       !
! 2018/09/12, Grace, revised, use quaternion.                           !
! 2018/12/12, Grace, generate rotation matrices and RMSD files, in order!
! to conjugate with program OotTranTraj.f90                             !
! 2018/12/13, Grace, rewrite this program by referring reference 4.     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module constants  
    implicit none 
    real(8),parameter   :: pi = 3.1415926536D0
    real(8),parameter   :: zero=0.0001D0
end module constants 

Program main
    implicit none
    integer(4)      :: i,j,NAtoms
    real(8) :: rmsd
    integer(4),allocatable,dimension(:) :: AtomName
    real(8),allocatable,dimension(:,:)  :: Coord_ref,Coord_var
    real(8),dimension(3,3)  :: CovMat,RMat

    ! Step 1. Check the input argument and read the files
        ! call print_purpose() !TODO: un-comment this if this code works
        call get_NAtoms(NAtoms)
        allocate(AtomName(NAtoms))
        allocate(Coord_ref(NAtoms,3))
        allocate(Coord_var(NAtoms,3))
        call get_coord(1,NAtoms,AtomName,Coord_var)
        call get_coord(2,NAtoms,AtomName,Coord_ref)
    ! Step 2. Translation: move two structures to the center of mass
        call get_move_to_center(NAtoms,Coord_var)
        call get_move_to_center(NAtoms,Coord_ref)
    ! Step 3. 
        call get_rmsd(NAtoms,Coord_var,Coord_ref,rmsd)
        write(*,'(A)',advance='no') 'Initial RMSD: '
        write(*,*) rmsd

        ! Covariance matrix -> rotation matrix
        call get_CovMat(Natoms,Coord_var,Coord_ref,CovMat)
        ! call get_Cov_RMat(CovMat,RMat)
        call operate_Rmat_coord(NAtoms,RMat,Coord_var)
        call print_coord(NAtoms,AtomName,Coord_var)

        ! call get_rmsd(NAtoms,Coord_var,Coord_ref,rmsd)
        ! write(*,'(A)',advance='no') 'Final RMSD: '
        ! write(*,*) rmsd
        ! write(*,'(A)',advance='no') 'RMSD based on quaternion: ' 
        ! call get_QRMSD(NAtoms,coord_var,coord_ref,qEigVal,rmsd)
        ! write(*,*) rmsd
        ! call print_coord(NAtoms,AtomName,Coord_var)
    ! Step 6. Print out the adjustable coordinate to file, move.$2
        ! call print_fileCoord(NAtoms,AtomName,Coord_var)
        ! call print_filerotM(Rmat)
    stop
end program main

subroutine print_purpose()
    implicit none
    write(*,'()') 
    write(*,'(A)') '-------------------------------------------------------------'
    write(*,'(A)') '|  Purpose:                                                 |'
    write(*,'(A)') '|  This program takes the first structure as the reference, |'
    write(*,'(A)') '|  and then change the coordinate of the second structure.  |'
    write(*,'(A)') '|                                                           |'
    write(*,'(A)') '|  Limitation:                                              !'
    write(*,'(A)') '|  1. only consider the amount of protons to move c.o.m.    |'
    write(*,'(A)') '|           subroutine: get_move_to_center                  |'
    write(*,'(A)') '|  2. only consider the following atom: H,C,O,N             |'
    write(*,'(A)') '|           subroutine: get_I                               |'
    write(*,'(A)') '|                                                           |'
    write(*,'(A)') '|  Input arguments:                                         |'
    write(*,'(A)') '|       $1 = variable structure                             |'
    write(*,'(A)') '|       $2 = reference structure                            |'
    write(*,'(A)') '|                                                           |'
    write(*,'(A)') '|  Output file:                                             |'
    write(*,'(A)') '|       move.$2                                             |'
    write(*,'(A)') '-------------------------------------------------------------'
    write(*,'()') 
    return
end subroutine print_purpose

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

subroutine print_filerotM(rotM)
    implicit none
    real(8),dimension(3,3),intent(in) :: rotM
    integer(4)  :: i
    character(len=100)  :: struc

    call GETARG(1,struc)
    open(10,file='rotM.'//TRIM(struc),status='replace')
    ! write(*,'(A)') 'rotM'
    do i=1,3
        ! write(*,'(3(1X,F14.10))') rotM(i,1:3)
        write(10,'(3(1X,F14.10))') rotM(i,1:3)
    end do
    close(10)
    return
end subroutine print_filerotM

subroutine print_fileCoord(NAtoms,AtomNAme,Coord)
    implicit none
    integer(4),intent(in)   :: NAtoms
    integer(4),dimension(NAtoms)    :: AtomName
    real(8),dimension(NAtoms,3),intent(in)  :: Coord
    ! local variable
    integer(4)  :: i
    character(len=100)  :: struc
    call GETARG(1,struc)
    open(10,file='overlap.'//TRIM(struc),status='replace')
    write(10,'(I2)') NAtoms
    write(10,'(A)') 'adjustable coordinate'
    do i=1,NAtoms
        write(10,'(I2,3(1X,F14.10))') AtomName(i),Coord(i,1:3)
    end do
    return
end subroutine print_fileCoord

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
    ! local variable
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

subroutine get_move_to_center(NAtoms,Coord)
    implicit none
    integer(4),intent(in)   :: NAtoms
    real(8),dimension(NAtoms,3),intent(inout)  :: Coord
    real(8),dimension(3)    :: mean
    integer(4)  :: i,j

    mean(:)=0.0D0
    do i=1,3
        do j=1,NAtoms
            mean(i) = mean(i) + Coord(j,i)
        end do
        mean(i) = mean(i)/NAtoms
    end do
    do i=1,NAtoms
        do j=1,3
            Coord(i,j)=Coord(i,j)-mean(j)
        end do
    end do
    return
end subroutine get_move_to_center

subroutine get_CovMat(NAtoms,coord1,coord2,CovMat)
    implicit none
    integer(4),intent(in)   :: NAtoms
    real(8),dimension(NAtoms,3),intent(in)  :: coord1,coord2
    real(8),dimension(3,3),intent(out)  :: CovMat
    ! local variable
    integer(4)  :: i,j,k

    ! covariance matrix = coord1^T * coord2
    CovMat = 0.0D0
    do i = 1,3
        do j = 1,3
            do k = 1,NAtoms
                CovMat(i,j) = CovMat(i,j) + coord1(k,i)*coord2(k,j)
            end do
        end do
    end do
    ! do i=1,3
    !     write(*,*) CovMat(i,:)
    ! end do
    return
end subroutine get_CovMat

subroutine get_RMSD(NAtoms,coord1,coord2,rmsd)
    implicit none
    integer(4),intent(in)   :: NAtoms
    real(8),dimension(NAtoms,3),intent(in)  :: coord1,coord2
    real(8),intent(out) :: rmsd
    ! local variable
    real(8) ::sumDiff
    integer(4)  :: i,j
    sumDiff=0.0D0
    do i=1,NAtoms
        do j=1,3
            sumDiff = sumDiff + ( coord1(i,j) - coord2(i,j) )**2
        end do
    end do
    rmsd = SQRT(sumDiff/ DBLE(NAtoms) ) 
    return
end subroutine get_RMSD

subroutine dia(dim,Matrix,eigValue,eigVec)
    implicit none
    integer(4),intent(in)   :: dim
    real(8),intent(in),dimension(dim,dim)   :: Matrix
    real(8),intent(inout),dimension(dim,dim)    :: eigVec
    real(8),intent(inout),dimension(dim)  :: eigValue
    integer(4),parameter    :: LWMAX=10000
    real(8),dimension(LWMAX)    :: WORK
    integer(4)  :: INFO,LWORK,i
    ! DSYEV computes all eigenvalues and, optionally, eigenvectors of
    ! a real symmetric matrix A

    ! DSYEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO)
    ! Query the optimal workspace
    LWORK = -1
    eigValue(:)=0.0D0
    call DSYEV('Vectors','Upper',dim,Matrix,dim,eigValue,WORK,LWORK,INFO)
    LWORK = MIN(LWMAX,INT(WORK(1)))
    ! Solve eigenproblem
    call DSYEV('Vectors','Upper',dim,Matrix,dim,eigValue,WORK,LWORK,INFO)
    eigVec = Matrix
    ! Check for convergence
    if ( INFO .gt. 0 ) then
        write(*,'(A)') 'The algorithm,DSYEV,failed to compute eigenvalues.'
        stop
    end if
    return
end subroutine dia

subroutine operate_Rmat_coord(NAtoms,RMat,coord)
    implicit none
    integer(4),intent(in)   :: NAtoms
    real(8),dimension(3,3)  :: RMat
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
                tmp(i,k)=tmp(i,k)+RMat(i,j)*coord_t(j,k)
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
end subroutine operate_Rmat_coord 

subroutine get_Cov_RMat(CovMat,RMat)
    implicit none
    real(8),dimension(3,3),intent(in)   :: CovMat
    real(8),dimension(3,3),intent(out)  :: RMat
    ! local variable
    integer(4)  :: i,j
    real(8),dimension(3)    :: eigValue
    real(8),dimension(3,3)  :: eigVec

    call dia(3,CovMat,eigValue,eigVec)
    ! write(*,*) eigValue(:)
    RMat = 0.0D0
    ! R = v^T * v
    do i = i,3
        do j = 1,3
            RMat(i,j) = sqrt(eigVec(i,3) * eigVec(j,3))
        end do
    end do
    return
end subroutine get_Cov_RMat
