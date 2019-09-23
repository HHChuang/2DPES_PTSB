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
!           3. rotM1                                                    !
!           4. rotM2                                                    !
!           5. rotM3 (?)                                                !
!                                                                       !
!   Visualizer :                                                        !
!           Jmol                                                        !
!                                                                       !
!   Reference:                                                          !
!       1. Video and slide: Math for Game Programmers: Understanding    !
!          Quaternions, Jim Van Verth                                   !!!!!!!!
!https://www.gdcvault.com/play/1017653/Math-for-Game-Programmers-Understanding !
!       2. Youtube: 3D Rotations in General: Rodrigues Rotation Formula !!!!!!!!
!                   and Quaternion Exponentials                         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!https://www.youtube.com/watch?annotation_id=annotation_1566057571&feature=iv&src_vid=UaK2q22mMEg&v=q-ESzg03mQc!
!       3. Github: https://github.com/gaschler/rotationconverter        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       4. https://is.muni.cz/th/u0qtg/thesis.pdf                       !
!                                                                       !
!   History:                                                            !
! 2017/05/25, Grace, Cherri's rotation matrix.                          !
! 2017/08/07, Grace, Rodrigues rotation (Chun-I).                       !
! 2018/09/12, Grace, revised, use quaternion.                           !
! 2018/12/12, Grace, generate rotation matrices and RMSD files, in order!
! to conjugate with program OotTranTraj.f90                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module constants  
    implicit none 
    real(8),parameter   :: pi = 3.1415926536D0
    real(8),parameter   :: zero=0.0001D0
end module constants 

Program main
    implicit none
    integer(4)      :: i,NAtoms
    real(8) :: rmsd
    integer(4),allocatable,dimension(:) :: AtomName
    real(8),allocatable,dimension(:,:)  :: Coord_ref,Coord_var

    ! Step 1. Check the input argument and read the files
        ! call print_purpose() !TODO: un-comment this if this code works
        call get_NAtoms(NAtoms)
        allocate(AtomName(NAtoms))
        allocate(Coord_ref(NAtoms,3))
        allocate(Coord_var(NAtoms,3))
        call get_coord(1,NAtoms,AtomName,Coord_var)
        call get_coord(2,NAtoms,AtomName,Coord_ref)
    ! Step 2. Translation: move two structures to the center of mass
        call get_move_to_center(NAtoms,AtomName,Coord_var)
        ! call print_coord(NAtoms,AtomName,Coord_var)
        call get_move_to_center(NAtoms,AtomName,Coord_ref)
        ! call print_coord(NAtoms,AtomName,Coord_ref)
    ! Step 3. First rotation: using the largest moment of 
    !         inertia to rotate the variable structure
        call do_rotation(NAtoms,AtomName,Coord_ref,3,Coord_var) ! 3 = the largest one
        ! write(*,'(A)') '1st. rot.'
        ! call print_coord(NAtoms,AtomName,Coord_var)
    ! Step 4. Second rotation: using the second largest moment of 
    !         inertia to rotate the variable structure
        ! call do_rotation(NAtoms,AtomName,Coord_ref,2,Coord_var) 
        ! write(*,'(A)') '2nd. rot.'
        ! call print_coord(NAtoms,AtomName,Coord_var)
    ! Step 5. Third rotation: using the smallest moment of 
    !         inertia to rotate the variable structure
        ! call do_rotation(NAtoms,AtomName,Coord_ref,1,Coord_var) 
        ! write(*,'(A)') '3rd. rot.'
        ! call print_coord(NAtoms,AtomName,Coord_var)
    ! Step 6. Print out the adjustable coordinate to file, move.$2
        call print_moveCoord(NAtoms,AtomNAme,Coord_var)
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

!FIXME: written for RotTraj.f90
subroutine print_rotM(order,rotM)
    implicit none
    integer(4),intent(in)   :: order
    real(8),dimension(3,3),intent(in) :: rotM
    character(len=100)  :: order_char
    integer(4)  :: i

    write(order_char,*) order
    open(10,file='rotM'//TRIM(ADJUSTL(order_char)),status='replace')
    do i=1,3
        write(10,'(3(1X,F14.10))') rotM(i,1:3)
    end do
    close(10)
    return
end subroutine print_rotM

subroutine print_moveCoord(NAtoms,AtomNAme,Coord)
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
end subroutine print_moveCoord

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

subroutine do_rotation(NAtoms,AtomName,Coord_ref,order,Coord_var)
    use constants
    implicit none
    integer(4),intent(in)   :: NAtoms,order
    integer(4),dimension(NAtoms),intent(in) :: AtomName
    real(8),dimension(NAtoms,3),intent(in)  :: Coord_ref
    real(8),dimension(NAtoms,3),intent(inout)  :: Coord_var
    ! local variable
    real(8) :: angle1,angle2,rmsd1,rmsd2
    real(8),dimension(3)  :: I1,I2,n
    real(8),dimension(3,3)  ::rotM1, rotM2
    real(8),dimension(NAtoms,3) :: Coord_var1,Coord_var2

    ! Step 1. Calculate the selected moment of inertia for two structures
        call get_I(NAtoms,AtomName,Coord_var,order,I1) ! order: 3 > 2 > 1
        write(*,'(A)') 'moment of inertia'
        write(*,*) I1
        call get_I(NAtoms,AtomName,Coord_ref,order,I2)
        write(*,*) I2
    ! Step 2. Calculate rotational angle by using their dot product
        call get_angle(I1,I2,angle1)
        ! The rotational angle can be positive or negative. Or, if the
        ! angle is callculated as 0, it can also be pi.
        if ( (angle1 .gt. -zero) .and. (angle1 .lt. zero) ) then
        ! case 1.: angle = 0, pi (pi = -pi)
            angle1 = pi
            angle2 = pi
        else
            ! case 2.: angle = + or - 
            angle2 = - angle1
        end if
    ! Step 3. Calculate the rotational axis by using their cross product
        call get_cross(I1,I2,n)
        write(*,'(A)') 'rotation axis'
        write(*,*) n
    ! Step 4. Formulate two rotational matries by using quaternion
        call get_rotM(n,angle1,rotM1)
        call get_rotM(n,angle2,rotM2)
    ! Step 5. Operate the rotation matrix on the variable structure 1 and 2 twice
        Coord_var1 = Coord_var
        Coord_var2 = Coord_var
        call operateR(NAtoms,rotM1,Coord_var1) 
        call operateR(NAtoms,rotM2,Coord_var2) 
        call get_RMSD(NAtoms,Coord_var1,Coord_ref,rmsd1)
        call get_RMSD(NAtoms,Coord_var2,Coord_ref,rmsd2)
    ! Step 6. Select the smallest rmsd as the final structure
        write(*,'(A)') 'rmsd'
        write(*,*) rmsd1,rmsd2
        if ( rmsd1 .lt. rmsd2 ) then
            ! write(*,'(A)') '+'
            Coord_var = Coord_var1
            ! Print rotaiton matrix for RotTraj.f90
            call print_rotM(4-order,rotM1)
        else
            ! write(*,'(A)') '-'
            Coord_var = Coord_var2
            !  Print rotaiton matrix for RotTraj.f90
            call print_rotM(4-order,rotM2)
        end if

        call get_I(NAtoms,AtomName,Coord_var,order,I1) ! order: 3 > 2 > 1
        write(*,'(A)') 'moment of inertia '
        write(*,*) I1
    return
end subroutine do_rotation

subroutine get_I(NAtoms,AtomName,Coord,order,MI)
    implicit none
    integer(4),intent(in)   :: NAtoms,order
    integer(4),dimension(NAtoms),intent(in) :: AtomName
    real(8),dimension(NAtoms,3),intent(in)  :: Coord
    real(8),dimension(3),intent(out)  :: MI
    ! local variable
    real(8),dimension(3,3)      :: Mom,eigVec
    real(8),dimension(NAtoms)   :: mass
    real(8),dimension(3)    :: eigValue
    integer(4)  :: i,j 
    ! Step 1. Change atomic number to  mass number
        do i=1,NAtoms
            Select case ( AtomName(i) )
                case (1) ! hydrogen
                    mass(i) = 1.008D0
                case (6) ! carbon
                    mass(i) = 12.011D0
                case (7) ! nitrogen
                    mass(i) = 14.007D0
                case (8) ! oxygen
                    mass(i) = 15.999
            end select
        end do
    ! Step 2. Calculate the matrix element of the moment of inertia 
    !         (rank-2 tensor)
        Mom(:,:)=0.0D0
        do i=1,NAtoms
            ! Calculate the diagonal term
            Mom(1,1)=Mom(1,1)+(Coord(i,2)**2+Coord(i,3)**2)*mass(i) ! Ixx
            Mom(2,2)=Mom(2,2)+(Coord(i,1)**2+Coord(i,3)**2)*mass(i) ! Iyy
            Mom(3,3)=Mom(3,3)+(Coord(i,1)**2+Coord(i,2)**2)*mass(i) ! Izz
            ! Calculate the off-diagonal term
            Mom(1,2)=Mom(1,2)-(Coord(i,1)*Coord(i,2)*mass(i)) ! Ixy
            Mom(2,3)=Mom(2,3)-(Coord(i,2)*Coord(i,3)*mass(i)) ! Iyz
            Mom(1,3)=Mom(1,3)-(Coord(i,1)*Coord(i,3)*mass(i)) ! Ixz
        end do
        ! It is a symmetry matrix
        Mom(2,1)=Mom(1,2)
        Mom(3,2)=Mom(2,3)
        Mom(3,1)=Mom(1,3)
    ! Step 3. Diagonalization of moment of inertia tensor
        call dia(Mom,eigValue,eigVec) 
        ! write(*,'(A)') 'Eigenvalue and eigenvector of moment of inertia'
        ! write(*,*) eigValue(:)   
        ! write(*,'()')
        ! write(*,*) eigVec(:,1)
        ! write(*,*) eigVec(:,2)
        ! write(*,*) eigVec(:,3)
        ! write(*,'()')
        MI=eigVec(:,order)
    return
end subroutine get_I

subroutine dia(Mom,eigValue,eigVec)
    implicit none
    real(8),intent(in),dimension(3,3)   :: Mom
    real(8),intent(inout),dimension(3,3)    :: eigVec
    real(8),intent(inout),dimension(3)  :: eigValue
    integer(4),parameter    :: LWMAX=1000
    real(8),dimension(LWMAX)    :: WORK
    integer(4)  :: INFO,LWORK,i
    ! DSYEV computes all eigenvalues and, optionally, eigenvectors of
    ! a real symmetric matrix A

    ! DSYEV(JOBZ,UPLO,N,A,LDA,W,WORK,LWORK,INFO)
    ! Query the optimal workspace
    LWORK=-1
    eigValue(:)=0.0D0
    call DSYEV('Vectors','Upper',3,Mom,3,eigValue,WORK,LWORK,INFO)
    LWORK=MIN(LWMAX,INT(WORK(1)))
    ! Solve eigenproblem
    call DSYEV('Vectors','Upper',3,Mom,3,eigValue,WORK,LWORK,INFO)
    eigVec=Mom
    ! Check for convergence
    if ( INFO .gt. 0 ) then
        write(*,'(A)') 'The algorithm,DSYEV,failed to compute eigenvalues.'
        stop
    end if
    return
end subroutine dia

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

subroutine get_rotM(rotAxis,angle,RotMatrix)
    ! online convertor: https://www.andre-gaschler.com/rotationconverter/
    use constants
    implicit none
    real(8),dimension(3),intent(in) :: rotAxis
    real(8),intent(in) :: angle
    real(8),dimension(3,3),intent(out)  :: RotMatrix
    ! local variable
    integer(4)  :: i
    real(8) :: length
    real(8),dimension(4)    :: q
    ! calculate quaternion
    q(1)=COS(angle/2)
    do i=1,3
        q(i+1)=SIN(angle/2)*rotAxis(i)
    end do
    ! normalize q to get unit quaternion
    length =0.0D0
    do i=1,4
        length = length + q(i)**2 
    end do
    if ( length .gt. zero ) then
        do i=1,4
            q(i) = q(i)/length
        end do
    end if
    write(*,'(A)') 'angle'
    write(*,*) angle
    ! write(*,'(A)') 'quaternion'
    ! write(*,*) q(:)
    ! create rotation matrix
                ! (1,1)
    RotMatrix=RESHAPE(SOURCE=(/ 1-2*( q(3)**2 + q(4)**2 ),&
                    ! (1,2)
                    2*( q(2)*q(3) - q(4)*q(1) ), &
                    ! (1,3)
                    2*( q(2)*q(4) + q(3)*q(1) ), &
                    ! (2,1)
                    2*( q(2)*q(3) + q(4)*q(1) ), &
                    ! (2,2)
                    1-2*( q(2)**2 + q(4)**2 ),    &
                    ! (2,3)
                    2*( q(3)*q(4) - q(2)*q(1) ), &
                    ! (3,1)
                    2*( q(2)*q(4) - q(3)*q(1) ), &
                    ! (3,2)
                    2*( q(3)*q(4) + q(2)*q(1) ), &
                    ! (3,3)
                    1-2*( q(2)**2 + q(3)**2 ) /), &
                SHAPE=(/ 3,3 /) ) 
                write(*,'()')
                do i=1,3
                    write(*,'(3(F7.4,1X))') RotMAtrix(1,i),RotMAtrix(2,i),RotMAtrix(3,i)
                end do
                write(*,'()')
    return
end subroutine get_rotM
 
subroutine get_cross(vec1,vec2,cross)
    use constants
    implicit none
    real(8),dimension(3),intent(in) :: vec1,vec2
    real(8),dimension(3),intent(inout)  :: cross
    ! local variable
    integer(4)  :: i
    real(8) :: length
    ! cross product
    cross(1)=vec1(2)*vec2(3)-vec1(3)*vec2(2)
    cross(2)=-vec1(1)*vec2(3)+vec1(3)*vec2(1)
    cross(3)=vec1(1)*vec2(2)-vec1(2)*vec2(1)
    ! normalized, length must be positive or zero
    length=SQRT(cross(1)**2+cross(2)**2+cross(3)**2)
    if ( length .gt. zero ) then 
        do i=1,3
            cross(i)=cross(i)/length 
        end do
    end if
    return
end subroutine get_cross

subroutine get_angle(vec1,vec2,angle)
    use constants
    implicit none   
    real(8),dimension(3),intent(in) :: vec1,vec2
    real(8),intent(out) :: angle
    integer(4)  :: i
    real(8) :: dot,length_v1,length_v2
    ! dot product
    dot=0.0D0
    length_v1=0.0D0
    length_v2=0.0D0
    do i=1,3
        dot=dot+vec1(i)*vec2(i)
        length_v1=length_v1+vec1(i)**2
        length_v2=length_v2+vec2(i)**2
    end do
    if ( ( (length_v1 .gt. -zero) .and. (length_v1 .lt. zero) ) .or. &
         ( (length_v2 .gt. -zero) .and. (length_v2 .lt. zero) ) ) then
            write(*,'(A)')'angle has problem'
            stop
        else
            angle=ACOS(dot/(length_v1*length_v2))
    end if
    return
end subroutine get_angle

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
            sumDiff = sumDiff + (coord1(i,j)-coord2(i,j))**2
        end do
    end do
    rmsd = SQRT(sumDiff/(NAtoms*1.0D0)) !NAtoms must larger then 0
    return
end subroutine get_RMSD
